library(anndata)
library(stringr)
library(dendextend)
library(parallel)
library(MASS)
library(reticulate)
library(data.tree)
library(splitstackshape)

use_python("~/miniconda3/envs/renv/bin/python")

# Compute Poisson deviances
poisson_dev_batch <- function(y,x) {
  if (is.null(x)) {
    n <- Matrix::colSums(y)
    pis <- Matrix::rowSums(y)/sum(y)
    mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
    d <- 2 * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
    d[d<0] <- 0
    
    return(sqrt(d)*ifelse(y>mu,1,-1))
  } else {
    y1 <- lapply(unique(x),function(i) y[,x==i,drop=F])
    n <- lapply(y1,Matrix::colSums)
    pis <- lapply(y1,function(data) Matrix::rowSums(data)/sum(data))
    mu <- lapply(1:length(y1),function(ind)
      crossprod(array(pis[[ind]],dim=c(1,length(pis[[ind]]))),
                array(n[[ind]],dim=c(1,length(n[[ind]])))))
    d <- lapply(1:length(y1),function(ind)
      2 * (y1[[ind]] * log(ifelse(y1[[ind]] == 0, 1, y1[[ind]]/mu[[ind]])) -
             (y1[[ind]] - mu[[ind]])))
    
    res <- array(0,dim=dim(y))
    rownames(res) <- rownames(y)
    colnames(res) <- colnames(y)
    for (ind in 1:length(y1)) {
      d[[ind]][d[[ind]]<0] <- 0
      res[,x==unique(x)[ind]] <- as.matrix(sqrt(d[[ind]])*
                                             ifelse(y1[[ind]]>mu[[ind]],1,-1))
    }
    
    return(res)
  }
}

# Compute Poisson dispersion statistics
poisson_dispersion_stats <- function(y){
  n <- Matrix::colSums(y)
  pis <- Matrix::rowSums(y)/sum(y)
  mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
  y2 <- (y - mu)^2 / mu
  
  disp <- Matrix::rowSums(y2)/ncol(y2)
  
  if (!'matrix'%in%class(y2)) {
    y2 <- as.matrix(y2)
  }
  
  return(sqrt(ncol(y))*(disp-1)/sqrt(matrixStats::rowVars(y2)))
}

# Perform dimension reduction
reduce_dimension <- function(y,x,num_PCs) {
  pdev <- poisson_dev_batch(y,x)
  pdev <- t(scale(Matrix::t(pdev),scale=F))
  PCs <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev)),k=num_PCs)
  projection <- t(crossprod(PCs$vectors,pdev))
  
  return(list(PCs, projection))
}

# Compute expected sum of squares from dimension reduction scores
compute_ess <- function(redduc) {
  sum((rowSums(sweep(redduc,2,colMeans(redduc),'-')^2)))
}

# Compute test statistic
ward_linkage <- function(redduc,labels) {
  ess1 <- compute_ess(redduc[labels==1,])
  ess2 <- compute_ess(redduc[labels==2,])
  ess <- compute_ess(redduc)
  return((ess-(ess1+ess2))/length(labels))
}


# Fit model for one batch
fit_model_batch <- function(y,on_genes,num_PCs) {
  # Compute sample moments of the on genes
  on_counts <- Matrix::t(y[on_genes,])
  cov <- cov(as.matrix(on_counts))
  means <- Matrix::colMeans(on_counts)
  
  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  mus.sum <- tcrossprod(array(mus,dim=c(length(mus),1)),
                        array(1,dim=c(length(mus),1)))+
    tcrossprod(array(1,dim=c(length(mus),1)),array(mus,dim=c(length(mus),1)))
  sigmas.sum <- tcrossprod(array(sigmas,dim=c(length(sigmas),1)),
                           array(1,dim=c(length(sigmas),1)))+
    tcrossprod(array(1,dim=c(length(sigmas),1)),
               array(sigmas,dim=c(length(sigmas),1)))
  rhos <- suppressWarnings(log(cov/(exp(mus.sum+0.5*sigmas.sum))+1))
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas
  
  # Make the covariance matrix positive-definite
  on_cov_eigs <- RSpectra::eigs_sym(as.matrix(rhos),k=min(c(nrow(rhos)-1,num_PCs)))
  num_pos <- sum(on_cov_eigs$values>0)
  on_cov.sub <- on_cov_eigs$vectors[,1:num_pos]%*%
    sqrt(diag(on_cov_eigs$values[1:num_pos]))
  on_cov <- tcrossprod(on_cov.sub)
  diag(on_cov) <- diag(rhos)
  on_cov <- sfsmisc::posdefify(on_cov)
  on_cov.sqrt <- t(chol(on_cov))
  
  return(list(Matrix::rowMeans(y),mus,on_cov.sqrt))
}

# Fit model
fit_model <- function(y,on_genes,x,num_PCs) {
  on_means <- list()
  on_cov.sqrt <- list()
  lambdas <- list()
  
  for (b in unique(x)) {
    params <- fit_model_batch(y[,x==b],on_genes=on_genes,num_PCs=num_PCs)
    lambdas[[as.character(b)]] <- params[[1]]
    on_means[[as.character(b)]] <- params[[2]]
    on_cov.sqrt[[as.character(b)]] <- params[[3]]
  }
  
  return(list(lambdas,on_means,on_cov.sqrt))
}

# Generate one null sample
generate_null <- function(y,params,on_genes,x,null_gene) {
  lambdas <- params[[1]]
  on_means <- params[[2]]
  on_cov.sqrt <- params[[3]]
  
  null <- array(0,dim=dim(y))
  rownames(null) <- rownames(y)
  
  for (b in unique(x)) {
    num_gen <- min(sum(x==b),null_gene)
    names(lambdas[[as.character(b)]]) <- rownames(y)
    null[-on_genes,which(x==b)[1:num_gen]] <-
      array(rpois(num_gen*(nrow(null)-length(on_genes)),
                  lambdas[[as.character(b)]][-on_genes]),
            dim=c(nrow(null)-length(on_genes),num_gen))
    Y <- exp(sweep(on_cov.sqrt[[as.character(b)]]%*%
                     array(rnorm(num_gen*length(on_genes)),
                           dim=c(length(on_genes),num_gen)),1,
                   on_means[[as.character(b)]],'+'))
    null[on_genes,which(x==b)[1:num_gen]] <-
      array(rpois(length(Y),Y),dim=dim(Y))
  }
  
  return(list(null[,colSums(null)>0],x[colSums(null)>0]))
}

# Generate one null sample and compute test statistic
generate_null_statistic <- function(y,params,on_genes,x,num_PCs,
                                    gm,labs,posthoc,null_gene) {
  null_set <- generate_null(y,params,on_genes,x,null_gene)
  null <- null_set[[1]]
  batch_null <- null_set[[2]]
  
  if (!posthoc) {
    null_gm <- reduce_dimension(null,batch_null,num_PCs)[[2]]
    null_gm.d <- dist(null_gm)
    hc2 <- cutree(hclust(null_gm.d,method='ward.D'),2)
  } else {
    null_gm <- reduce_dimension(null,batch_null,num_PCs)[[2]]
    pdev <- poisson_dev_batch(null,batch_null)
    pdev <- t(scale(Matrix::t(pdev),scale=F))
    gm2 <- t(crossprod(gm[[1]]$vectors,pdev))
    knns <- queryKNN(gm[[2]], gm2,
                     k=15, BNPARAM=AnnoyParam())$index
    neighbor.labels <- apply(knns,1,function(x) labs[x])
    hc2 <- unlist(apply(neighbor.labels,2,function(x)
      sample(names(table(x))[as.vector(table(x))==max(table(x))],1)))
    
    if (length(unique(hc2))==1) {
      null_gm.d <- dist(null_gm)
      hc2 <- cutree(hclust(null_gm.d,method='ward.D'),2)
    }
  }
  
  Qclust2 <- sapply(unique(batch_null),function(b)
    if (length(unique(hc2[batch_null==b]))==2 &
        min(table(hc2[batch_null==b]))>=2) {
      ward_linkage(null_gm[batch_null==b,],hc2[batch_null==b])
    } else {
      0
    })
  names(Qclust2) <- unique(batch_null)
  return(median(Qclust2))
}

# Test one split
test_split <- function(data,ids1,ids2,var.genes,num_PCs,batch,
                       alpha_level,cores,posthoc,null_gene) {
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  batch <- batch[c(ids1,ids2)]
  
  # Re-run dimension reduction and calculate test statistic
  gm <- reduce_dimension(true[var.genes,],batch,num_PCs)
  gm_sub.x <- gm[[2]]
  labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
  Qclust <- sapply(unique(batch),function(b)
    ward_linkage(gm_sub.x[batch==b,],labs[batch==b]))
  names(Qclust) <- unique(batch)
  stat <- median(Qclust)
  
  # Determine the set of "on" genes
  phi_stat <- poisson_dispersion_stats(true[var.genes,])
  check_means <- matrixStats::rowMins(sapply(unique(batch),function(b)
    Matrix::rowSums(true[var.genes,batch==b])))
  on_genes <- which(pnorm(phi_stat,lower.tail=F)<0.05&check_means!=0)
  
  # Fit model
  params <- fit_model(true[var.genes,],on_genes,batch,num_PCs)
  
  # Generate null distribution of test statistics
  Qclusts2_1 <- mclapply(1:10,function(i) {
    generate_null_statistic(true[var.genes,],params,on_genes,batch,
                            num_PCs,gm,labs,posthoc,null_gene)
  },mc.cores=cores)
  
  # Quit early if p-value is much smaller or much larger than alpha level
  Qclusts2 <- unlist(Qclusts2_1)
  fit <- fitdistr(Qclusts2,'normal')
  pval <- 1-pnorm(stat,mean=fit$estimate[1],sd=fit$estimate[2])
  if (pval < 0.1*alpha_level | pval > 10*alpha_level) {
    return(pval)
  }
  
  # Otherwise, keep going
  Qclusts2_2 <- mclapply(11:50,function(i) {
    generate_null_statistic(true[var.genes,],params,on_genes,
                            batch,num_PCs,gm,labs,posthoc,null_gene)
  },mc.cores=cores)
  
  # Compute smoothed p-value
  Qclusts2 <- c(Qclusts2,unlist(Qclusts2_2))
  fit <- fitdistr(Qclusts2,'normal')
  return(1-pnorm(stat,mean=fit$estimate[1],sd=fit$estimate[2]))
}

# Full clustering pipeline with built-in hypothesis testing
scSHCmain <- function(data,batch=NULL,alpha=0.05,num_features=2500,
                  num_PCs=30,parallel=T,cores=2,null_gene=5000, stop_prop=0.05) {
  if (!parallel) {
    cores <- 1
  }
  if (is.null(batch)) {
    batch <- rep("1",ncol(data))
  }
  if (is.factor(batch)|is.numeric(batch)) {
    warning("Converting batch vector to character")
    batch <- as.character(batch)
  }
  names(batch) <- colnames(data)
  if (is.null(colnames(data))) {
    warning("Assigning cell names automatically to the columns")
    colnames(data) <- paste0('cell',1:ncol(data))
  }
  
  # Get variable features
  dev <- scry::devianceFeatureSelection(data)
  var.genes <- rownames(data)[order(dev,decreasing=T)[1:num_features]]
  
  # Dimension reduction and clustering
  gm.x <- reduce_dimension(data[var.genes,],batch,num_PCs)[[2]]
  gm.d <- dist(gm.x)
  hcl <- fastcluster::hclust(gm.d,method='ward.D')
  dend <- as.dendrogram(hcl)
  
  # Traverse the tree and store clusters and q-FWERs
  dends_to_test <- list(dend)
  clusters <- list()
  node0 <- NULL
  counter <- 0
  parents <- list('root')
  while (length(dends_to_test)> 0 ) {
    
    if (length(get_leaves_attr(dends_to_test[[1]], "members")) > ncol(data)*stop_prop) {
      # Identify the two-way split
      cuts <- dendextend::cutree(dends_to_test[[1]],k=2,order_clusters_as_data=F)
      leaves <- get_leaves_attr(dends_to_test[[1]],'label')
      ids1 <- leaves[cuts==1]
      ids2 <- leaves[cuts==2]
      alpha_level <- alpha*((length(leaves)-1)/(ncol(data)-1))
      
      # Remove any cells from batches with poor representation
      tab <- table(batch[c(ids1,ids2)],cuts)
      to.keep <- rownames(tab)[which(matrixStats::rowMins(tab)>20)]
      ids1 <- ids1[batch[ids1]%in%to.keep]
      ids2 <- ids2[batch[ids2]%in%to.keep]
      
      # Get p-value of the split
      if (length(to.keep)>0) {
        test <- test_split(data,ids1,ids2,var.genes,num_PCs,
                           batch,alpha_level,cores,posthoc=F,null_gene)
      } else {
        test <- 1
      }
      
      # Compare to significance threshold
      if (test < alpha_level) {
        # If significant: add left and right branches to testing stack
        left <- dends_to_test[[1]][[1]]
        right <- dends_to_test[[1]][[2]]
        dends_to_test[[length(dends_to_test)+1]] <- left
        dends_to_test[[length(dends_to_test)+1]] <- right
        
        # Compute q-FWER
        if (is.null(node0)) {
          node0 <- Node$new(paste0('Node 0: ',
                                   min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))
        } else {
          do.call("<-",list(paste0('node',counter),
                            eval(parse(text=parents[[1]]))$AddChild(paste0('Node ',counter,': ',
                                                                           min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))))
        }
        parents[[length(parents)+1]] <- paste0('node',counter)
        parents[[length(parents)+1]] <- paste0('node',counter)
        counter <- counter + 1
      } else {
        # If not significant: create cluster
        clusters[[length(clusters)+1]] <- leaves
        
        # Compute q-FWER
        if (is.null(node0)) {
          node0 <- Node$new(paste0('Cluster 0: ',
                                   min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))
        } else {
          do.call("<-",list(paste0('Cluster',length(clusters)),
                            eval(parse(text=parents[[1]]))$AddChild(paste0('Cluster ',length(clusters),': ',
                                                                           min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))))
        }
      }
    } else {
      leaves <- get_leaves_attr(dends_to_test[[1]],'label')
      clusters[[length(clusters)+1]] <- leaves
      if (is.null(node0)) {
        node0 <- Node$new(paste0('Cluster 0: ',
                                 min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))
      } else {
        do.call("<-",list(paste0('Cluster',length(clusters)),
                          eval(parse(text=parents[[1]]))$AddChild(paste0('Cluster ',length(clusters),': min cell' ) )))
      }
    }
    dends_to_test[[1]] <- NULL
    parents[[1]] <- NULL
  }
  
  # Produce vector of cluster labels
  cluster_labels <- rep(0,ncol(data))
  names(cluster_labels) <- colnames(data)
  for (i in 1:length(clusters)) {
    cluster_labels[clusters[[i]]] <- i
  }
  
  return(list(cluster_labels,node0))
}


gw_lst <- c("gw15", "gw20", "gw22", "gw34")
total_num <- 500000
sig_level_lst <- c(0.0005)

for (sig_level in sig_level_lst) {
  for (gw in gw_lst) {
    adata_tot <- read_h5ad(paste0(gw, ".h5ad"))
    adata_tot <- AnnData(adata_tot$raw$X, obs = adata_tot$obs, var = adata_tot$var, obsm = list(spatial = adata_tot$obsm$spatial))
    cluster_lst <- table(adata_tot$obs$H2_annotation)
    cluster_num <- as.numeric(round(cluster_lst/nrow(adata_tot)*total_num))
    names(cluster_num) <- names(cluster_lst)
    set.seed(123)
    gw_sample <- stratified(adata_tot$obs, group = "H2_annotation", size = cluster_num, keep.rownames = T)
    gw_rn <- gw_sample$rn
    adata <- adata_tot[rownames(adata_tot$obs) %in% gw_rn]
    cluster_lst <- sort(table(adata$obs$H2_annotation))
    
    setwd("scSHC_result")
    if (!dir.exists("H2")) {
      dir.create(paste0("H2"))
    }
    setwd("H2")
    if (!dir.exists(paste(as.character(sig_level), "result", sep = "_"))) {
      dir.create(paste(as.character(sig_level), "result", sep = "_"))
    }
    setwd(paste(as.character(sig_level), "result", sep = "_"))
    if (!dir.exists(gw)) {
      dir.create(gw)
    }
    setwd(gw)
    for (i in 1:length(cluster_lst)) {
      if (is.na(names(cluster_lst[i]))) {
        cluster_name <- names(cluster_lst)
      } else {
        cluster_name <- names(cluster_lst[i])
      }
      adata_h2 <- adata[adata$obs$H2_annotation == cluster_name]
      adata_count <- as(t(as.matrix(adata_h2$X)), "dgCMatrix")
      null_genenum <- min(8000, as.integer(cluster_lst[i]/3))
      try({
        adata_cluster <- scSHCmain(adata_count, alpha = sig_level, num_features = 300, parallel = FALSE, 
                                   null_gene = null_genenum, stop_prop = 0.1)
        adata_cluster_df <- data.frame(adata_cluster[[1]])
        if (str_detect(cluster_name, "/") == TRUE) {
          filename <- gsub("/", "|", cluster_name)
        } else {
          filename <- cluster_name
        }
        write.csv(adata_cluster_df, paste(filename, ".csv",sep = ""))
        write.table(adata_cluster[[2]], paste(filename, "_node.txt",sep = ""), row.names = FALSE)
      }, silent = TRUE)
      
    }
    setwd("../../../../")
  }
}



