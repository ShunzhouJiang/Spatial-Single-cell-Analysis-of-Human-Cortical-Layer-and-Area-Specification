library(reticulate)
library(dplyr)
library(scrattch.hicat)
library(stringr)
library(tibble)
library(purrr)
library(gridExtra)
library(sna)
library(Hmisc)
library(reshape2)
#library(ggalt)
library(ggforce)
library(dplyr)
library(ggrepel)
library(ggnewscale)


for (cutoff in c(0.02)){
  col_classes <- c("character", "character", "character")  
  cells.cl.df = read.csv('constellation_group_obs.csv', header = TRUE, colClasses = col_classes, row.names = 1)
  cluster_size = table(cells.cl.df$group)
  cells.cl.df$cell_name = rownames(cells.cl.df)
  rd.dat <- py_load_object('constellation_red.pickle')
  rownames(rd.dat$pca) = rownames(cells.cl.df)
  rownames(rd.dat$umap) = rownames(cells.cl.df)
  
  
  np <- import("numpy")
  knn.result = list()
  knn.result[[1]] <- np$load('constellation_indices.npy')
  knn.result[[1]] = knn.result[[1]] + 1
  knn.result[[2]] <- np$load('constellation_distances.npy')
  groups.col = 'group'
  
  
  cell_names <- cells.cl.df$cell_name
  # Cluster, or area-celltype combination.
  cl <- cells.cl.df[[groups.col]] %>% as.factor %>% magrittr::set_names(value = cell_names)
  
  cl.numeric <- as.numeric(cl) %>% magrittr::set_names(value = cell_names)
  
  cl.df <- get_cl_df(cl)
  
  cl.df$clade <- str_split_fixed(cl.df$cluster_label, "_", 2)[ ,2] %>% tolower 
  cl.df$area <- str_split_fixed(cl.df$cluster_label, "_", 2)[ ,1] %>% tolower 
  
  cl.df[names(cluster_size),'cluster_size'] = cluster_size
  cells.cl.df$cluster_label = cells.cl.df$group
  
  cells.cl.df <- cells.cl.df %>% dplyr::select(-cluster_label) %>% dplyr::rename(cluster_label = groups.col) %>%
    left_join(
      # %>% select(cell.name, groups.col, combined.cluster.2),
      cl.df, by = "cluster_label") %>%
    # Requires cells.cl.df (metadata) to have column being used for groups
    # named 'cluster_label' to match with cl_df during join.
    mutate(cluster_id = as.factor(cluster_id))
  
  
  rd.cl.center <- get_RD_cl_center(rd.dat = rd.dat$umap, cl)
  
  rd.cl.center %<>% 
    as.data.frame %>% 
    set_names(c("x", "y")) %>%
    add_column(cl = cl.df$cluster_id, .before = "x") %>%
    # add_column preserves rownames.
    # but moving rownames to column cluster_label anyway bc of left_join below.
    # Needs to be cl (not cl_id) or else you get error:
    # Error in `$<-.data.frame`(`*tmp*`, "edge.frac.within", value = numeric(0)) : 
    # replacement has 0 rows, data has 26 
    rownames_to_column("cluster_label")
  
  cl.center.df <- left_join(rd.cl.center, cl.df,
                            by = c("cluster_label"))
  
  knn.outlier.th=2
  outlier.frac.th=0.5
  row.names(knn.result[[1]]) = row.names(knn.result[[2]])=row.names(rd.dat$umap) 
  knn  = knn.result[[1]]
  knn.dist = knn.result[[2]]
  cl.knn.dist.mean = tapply(names(cl),cl, function(x) mean(knn.dist[x,-1]))
  cl.knn.dist.sd = tapply(names(cl),cl, function(x) sd(knn.dist[x,-1]))
  cl.knn.dist.th = (cl.knn.dist.mean + knn.outlier.th * cl.knn.dist.sd)
  
  knn.dist.th=cl.knn.dist.th[as.character(cl[row.names(knn)])]
  outlier = apply(knn.dist, 2, function(x) x>  knn.dist.th)
  row.names(outlier)  = row.names(knn.dist)
  knn[outlier] = NA
  select.cells = row.names(outlier)[rowMeans(outlier) < outlier.frac.th]
  
  pred.result = predict_knn(knn[select.cells,], row.names(rd.dat$umap), cl.numeric)
  
  pred.prob = pred.result$pred.prob
  knn.cell.cl.counts = round(pred.prob * ncol(knn))
  knn.cl.cl.counts = do.call("rbind",tapply(row.names(pred.prob), cl[row.names(pred.prob)], function(x)colSums(knn.cell.cl.counts[x,])))
  colnames(knn.cl.cl.counts) = rownames(knn.cl.cl.counts)
  knn.cl.df = as.data.frame(as.table(knn.cl.cl.counts))
  colnames(knn.cl.df)[1:2] = c("cl.from","cl.to")
  from.size = rowSums(knn.cl.cl.counts)
  to.size = colSums(knn.cl.cl.counts)
  total = sum(knn.cl.cl.counts)
  knn.cl.df$cl.from.total= from.size[as.character(knn.cl.df$cl.from)]
  knn.cl.df$cl.to.total = to.size[as.character(knn.cl.df$cl.to)]
  knn.cl.df = knn.cl.df[knn.cl.df$Freq > 0,]
  knn.cl.df$pval.log = knn.cl.df$odds  = 0
  for(i in 1:nrow(knn.cl.df)){
    q = knn.cl.df$Freq[i] - 1
    k = knn.cl.df$cl.from.total[i]
    m = knn.cl.df$cl.to.total[i]
    n = total - m
    knn.cl.df$pval.log[i]=phyper(q, m=m, n=n, k=k, lower.tail = FALSE, log.p=TRUE)
    knn.cl.df$odds[i] = (q + 1) / (k * m /total)
  }
  knn.cl.df$frac = knn.cl.df$Freq/knn.cl.df$cl.from.total
  knn.cl.df$cl.from.label = cl.df[as.character(knn.cl.df$cl.from),"cluster_label"]
  knn.cl.df$cl.to.label = cl.df[as.character(knn.cl.df$cl.to),"cluster_label"]
  knn.cl = list(knn.result=knn.result, pred.result=pred.result, knn.cl.df=knn.cl.df)
  filterKNN <- function(knn.cl, frac.th = 0.1) {
    knn.cl.df.filter <- knn.cl$knn.cl.df %>% dplyr::filter(frac >= frac.th) %>% 
      mutate(cl.from = as.numeric(cl.from), cl.to = as.numeric(cl.to))
  }
  
  knn.cl.df.filter <- filterKNN(knn.cl = knn.cl, frac.th = cutoff)
  st <- format(Sys.time(), "%Y%m%d.%H%M_")
  
  node.label = "cluster_id"
  exxageration = 1
  curved = TRUE
  plot.parts=FALSE
  plot.hull = NULL
  plot.height=25
  plot.width=25
  node.dodge=TRUE
  label.size=2
  max_size=10
  
  out.dir = dir
  st = format(Sys.time(), "%Y%m%d_%H%M%S_")
  
  if(!file.exists(out.dir)){
    dir.create(out.dir)
  }
  ###==== Cluster nodes will represent both cluster.size (width of point) and edges within cluster (stroke of point)
  
  # select rows that have edges within cluster
  knn.cl.same <- knn.cl.df.filter[knn.cl.df.filter$cl.from == knn.cl.df.filter$cl.to,]
  
  cl.center.df = cl.center.df[, ! names(cl.center.df) %in% c('size')]
  #append fraction of edges within to cl.center.umap for plotting of fraction as node linewidth
  cl.center.df$edge.frac.within <- knn.cl.same$frac[match(cl.center.df$cl, knn.cl.same$cl.from)]
  
  
  ###==== plot nodes
  #labels <- cl.center.df[[node.label]]
  labels <- cl.center.df$cluster_label
  
  p.nodes <-   ggplot() +
    geom_point(data=cl.center.df,
               shape=19,
               aes(x=x,
                   y=y,
                   size=cluster_size,
                   color=alpha(cluster_color, 0.8))) +
    scale_size_area(trans="sqrt", max_size=max_size, breaks = c(5000,10000,15000,20000)) +
    scale_color_identity() +
    geom_text(data=cl.center.df,
              aes(x=x,
                  y=y,
                  label=labels),
              size = label.size)
  
  g <- ggplot_build(p.nodes)
  dots <<- g[["data"]][[1]] #dataframe with geom_point size, color, coords
  
  nodes <<- left_join(cl.center.df, dots, by=c("x","y"))
  
  if (node.dodge==TRUE){
    
    #<><><># make update here to convert units by scale. check geom_mark_hull code for oneliner
    
    # dodge nodes starting at center of plot moving outward
    
    nodes$r<- (nodes$size/10)/2
    
    
    x.list <- c(mean(nodes$x), nodes$x )
    y.list <- c(mean(nodes$y), nodes$y)
    dist.test <- as.matrix(dist(cbind(x.list, y.list)))
    nodes$distance <- dist.test[2:nrow(dist.test), 1]
    nodes <- nodes[order(nodes$distance),]
    
    
    for (d1 in 1:(nrow(nodes)-1)) {
      j <- d1+1
      for (d2 in j:nrow(nodes)) {
        print(paste(d1,d2))
        
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
        
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
        
        if (distSq < radSumSq) {
          print(paste(d1,d2))
          
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
    for (d1 in 1:(nrow(nodes)-1)) {
      j <- d1+1
      for (d2 in j:nrow(nodes)) {
        print(paste(d1,d2))
        
        distSq <- sqrt(((nodes$x[d1]-nodes$x[d2])*(nodes$x[d1]-nodes$x[d2]))+((nodes$y[d1]-nodes$y[d2])*(nodes$y[d1]-nodes$y[d2])))
        
        radSumSq <- (nodes$r[d1] *1.25)+ (nodes$r[d2]*1.25) # overlapping radius + a little bit extra
        
        if (distSq < radSumSq) {
          print(paste(d1,d2))
          
          subdfk <- nodes[c(d1,d2),]
          subdfk.mod <- subdfk
          subdfd1 <- subdfk[1,]
          subdfd2  <- subdfk[2,]
          angsk <- seq(0,2*pi,length.out=nrow(subdfd2)+1)
          subdfd2$x <- subdfd2$x+cos(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfd2$y <- subdfd2$y+sin(angsk[-length(angsk)])*(subdfd1$r+subdfd2$r+0.5)#/2
          subdfk.mod[2,] <- subdfd2
          nodes[c(d1,d2),] <- subdfk.mod
        }
      }
    }
    
  }
  
  
  nodes <- nodes[order(nodes$cluster_id),]
  
  
  
  ## when printing lines to pdf the line width increases slightly. This causes the edge to extend beyond the node. Prevent this by converting from R pixels to points.
  conv.factor <- ggplot2::.pt*72.27/96
  
  
  ## line width of edge can be scaled to node point size
  nodes$node.width <- nodes$size
  
  knn.cl <- knn.cl.df.filter
  ##from knn.cl data frame remove all entries within cluster edges.
  knn.cl.d <- knn.cl[!(knn.cl$cl.from == knn.cl$cl.to),]
  nodes$cl=as.numeric(as.character(nodes$cl))
  knn.cl.d$cl.from <- as.numeric(as.character(knn.cl.d$cl.from))
  knn.cl.d$cl.to <- as.numeric(as.character(knn.cl.d$cl.to))
  
  knn.cl.d <- left_join(knn.cl.d, select(nodes, cl, node.width), by=c("cl.from"="cl"))
  colnames(knn.cl.d)[colnames(knn.cl.d)=="node.width"]<- "node.pt.from"
  knn.cl.d$node.pt.to <- ""
  knn.cl.d$Freq.to <- ""
  knn.cl.d$frac.to <- ""
  
  knn.cl.bid <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i,])
    r <- subset(knn.cl.d[i:nrow(knn.cl.d),])
    r <- r[(line$cl.from == r$cl.to & line$cl.to == r$cl.from ),]
    
    if (dim(r)[1] != 0) {
      line$Freq.to <- r$Freq
      line$node.pt.to <- r$node.pt.from
      line$frac.to <- r$frac
      knn.cl.bid <- rbind(knn.cl.bid, line)
    }
    #print(i)
  }
  
  #unidirectional
  knn.cl.uni <- NULL
  for (i in 1:nrow(knn.cl.d)) {
    
    line <- subset(knn.cl.d[i, ])
    r <- knn.cl.d[(line$cl.from == knn.cl.d$cl.to & line$cl.to == knn.cl.d$cl.from ),]
    
    if (dim(r)[1] == 0) {
      knn.cl.uni <- rbind(knn.cl.uni, line)
    }
    #print(i)
  }
  
  #min frac value = 0.01
  knn.cl.uni$node.pt.to <- nodes$node.width[match(knn.cl.uni$cl.to, nodes$cl)]
  knn.cl.uni$Freq.to <- 1
  knn.cl.uni$frac.to <- 0.01
  knn.cl.lines <- rbind(knn.cl.bid, knn.cl.uni)
  
  ###==== create line segments
  
  line.segments <- knn.cl.lines %>% select(cl.from, cl.to)
  nodes$cl <- as.numeric((as.character(nodes$cl)))
  line.segments <- left_join(line.segments, select(nodes, x, y, cl), by = c("cl.from" = "cl"))
  line.segments <- left_join(line.segments, select(nodes, x, y, cl), by = c("cl.to" = "cl"))
  
  colnames(line.segments) <- c("cl.from", "cl.to", "x.from", "y.from", "x.to", "y.to")
  
  line.segments <- data.frame(line.segments,
                              freq.from = knn.cl.lines$Freq,
                              freq.to = knn.cl.lines$Freq.to,
                              frac.from = knn.cl.lines$frac,
                              frac.to =  knn.cl.lines$frac.to,
                              node.pt.from =  knn.cl.lines$node.pt.from,
                              node.pt.to = knn.cl.lines$node.pt.to)
  
  
  ##from points to native coords
  line.segments$node.size.from <- line.segments$node.pt.from/10
  line.segments$node.size.to <- line.segments$node.pt.to/10
  
  
  #change to max_size*line.segments$frac.from to make consistent for all nodes
  
  #line.segments$line.width.from <- line.segments$node.size.from*line.segments$frac.from
  #line.segments$line.width.to <- line.segments$node.size.to*line.segments$frac.to
  
  line.segments$line.width.from <- line.segments$frac.from
  line.segments$line.width.to <- line.segments$frac.to
  
  line.segments$line.width.from <- 0.1
  line.segments$line.width.to <- 0.1
  
  
  ##max fraction to max point size
  #line.segments$line.width.from<- (line.segments$frac.from/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.from
  
  #line.segments$line.width.to<- (line.segments$frac.to/max(line.segments$frac.from, line.segments$frac.to))*line.segments$node.size.to
  
  
  ###=== create edges, exaggerated width
  
  line.segments$ex.line.from <-line.segments$line.width.from #true to frac
  line.segments$ex.line.to <-line.segments$line.width.to #true to frac
  
  line.segments$ex.line.from <- pmin((line.segments$line.width.from*exxageration),line.segments$node.size.from) #exxagerated width
  line.segments$ex.line.to <- pmin((line.segments$line.width.to*exxageration),line.segments$node.size.to) #exxagerated width
  
  
  line.segments <- na.omit(line.segments)
  line.segments = line.segments[(line.segments$frac.from>cutoff & line.segments$frac.to>cutoff),]
  print("calculating edges")
  
  allEdges <- lapply(1:nrow(line.segments), edgeMaker, len = 500, curved = curved, line.segments=line.segments)
  allEdges <- do.call(rbind, allEdges)  # a fine-grained path with bend
  
  
  groups <- unique(allEdges$Group)
  
  poly.Edges <- data.frame(x=numeric(), y=numeric(), Group=character(),w=numeric(),stringsAsFactors=FALSE)
  imax <- as.numeric(length(groups))
  
  interpolate_points_for_each_row <- function(row, n_points) {
    # Extract x and y coordinates for the current row
    x1 <- row["x1"]
    y1 <- row["y1"]
    x2 <- row["x2"]
    y2 <- row["y2"]
    w = row['w1']
    frac = row['frac1']
    # Linearly interpolate between the points
    x_interp <- seq(x1, x2, length.out = n_points + 1)
    y_interp <- seq(y1, y2, length.out = n_points + 1)
    
    # Combine x and y coordinates into a dataframe
    interpolated_points <- data.frame(
      x = x_interp[-length(x_interp)],  # Remove the last point to avoid duplication
      y = y_interp[-length(y_interp)],
      w=w,
      frac=frac
    )
    
    return(interpolated_points)
  }
  
  
  for(i in 1:imax) {
    #svMisc::progress(i)
    #svMisc::progress(i, progress.bar=TRUE)
    select.group <- groups[i]
    #print(select.group)
    select.edge <- allEdges[allEdges$Group %in% select.group,]
    
    x <- select.edge$x
    y <- select.edge$y
    w <- select.edge$fraction
    #frac = c(rep(line.segments[i,'frac.from'],250),rep(line.segments[i,'frac.to'],250)) 
    frac = seq(line.segments[i,'frac.from'], line.segments[i,'frac.to'], length.out = 500)
    
    N <- length(x)
    leftx <- numeric(N)
    lefty <- numeric(N)
    rightx <- numeric(N)
    righty <- numeric(N)
    w_all = numeric(N)
    frac_all = numeric(N)
    
    ## Start point
    perps <- perpStart(x[1:2], y[1:2], w[1]/2)
    leftx[1] <- perps[1, 1]
    lefty[1] <- perps[1, 2]
    rightx[1] <- perps[2, 1]
    righty[1] <- perps[2, 2]
    w_all[1] = w[1]
    frac_all[1] = frac[1]
    
    ### mid points
    for (ii in 2:(N - 1)) {
      seq <- (ii - 1):(ii + 1)
      perps <- perpMid(as.numeric(x[seq]), as.numeric(y[seq]), w[ii]/2)
      leftx[ii] <- perps[1, 1]
      lefty[ii] <- perps[1, 2]
      rightx[ii] <- perps[2, 1]
      righty[ii] <- perps[2, 2]
      w_all[ii] = w[ii]
      frac_all[ii] = frac[ii]
    }
    ## Last control point
    perps <- perpEnd(x[(N-1):N], y[(N-1):N], w[N]/2)
    leftx[N] <- perps[1, 1]
    lefty[N] <- perps[1, 2]
    rightx[N] <- perps[2, 1]
    righty[N] <- perps[2, 2]
    w_all[N] = w[N]
    frac_all[N] = frac[N]
    
    
    #connect (leftx[i], lefty[i]),(rightx[i], righty[i]) 
    lineleft <- data.frame(x=leftx, y=lefty,w=w_all,frac=frac_all)
    lineright <- data.frame(x=rightx, y=righty,w=w_all,frac=frac_all)
    lines_all = cbind(lineleft, lineright)
    lines_all <- data.frame(lapply(lineleft, function(x) replicate(nrow(lineleft), NA)))
    lines_all$x = (lineleft$x + lineright$x)/2
    lines_all$y = (lineleft$y + lineright$y)/2
    lines_all$w = lineleft$w
    lines_all$frac = lineleft$frac
    #lines_all = lineleft
    #colnames(lines_all) = c('x1', 'y1', 'w1', 'frac1', 'x2', 'y2', 'w2', 'frac2')
    #interpolated_points <- apply(lines_all, 1, interpolate_points_for_each_row, n_points = 1)
    #interpolated_points_df <- do.call(rbind, interpolated_points)
    #lineright <- lineright[nrow(lineright):1, ]
    #lines.lr <- rbind(lineleft, lineright, interpolated_points_df)
    lines.lr = lines_all
    lines.lr$Group <- select.group
    poly.Edges <- rbind(poly.Edges,lines.lr)
    
    Sys.sleep(0.01)
    cat("\r", i, "of", imax)
    
  }
  
  
  
  
  labels <- nodes$cluster_label
  
  poly.Edges <- poly.Edges[order(poly.Edges$frac), ]
  rownames(line.segments) = paste(line.segments$cl.from, line.segments$cl.to, sep = '_')
  rownames(knn.cl.lines) = paste(knn.cl.lines$cl.from, knn.cl.lines$cl.to, sep = '_')
  line.segments$cl.from.label = knn.cl.lines[rownames(line.segments), 'cl.from.label']
  line.segments$cl.to.label = knn.cl.lines[rownames(line.segments), 'cl.to.label']
  line.segments1 = line.segments[,c('cl.from.label', 'cl.to.label', 'frac.from', 'frac.to', 'node.size.from', 'node.size.to')]
  write.csv(line.segments1, paste(dir,'color_interp1/th',cutoff,'_frac.csv', sep = ''))
  ####plot edges
  p.edges <- ggplot(poly.Edges, aes(x=x, y=y, group=Group, color = w)) + scale_color_gradientn(colours = rainbow(5)) + geom_point() + theme_void()
  
  #p.edges <- ggplot(poly.Edges, aes(x=x, y=y, group=Group, color = w)) + scale_color_gradientn(colours = rainbow(5))
  #p.edges <- p.edges +geom_polygon(aes(x=x, y=y, group=Group)) + theme_void()
  #p.edges
  
  #geom_polygon(data=poly.Edges,
  #             alpha=0.2,
  #             aes(x=x, y=y, group=Group))+
  
  
  pdf(paste(dir,'color_interp1/th',cutoff,'_one_line.pdf',sep=''), 5.76, 6.86)
  print(ggplot() +     geom_point(data=poly.Edges, aes(x=x, y=y, color = frac)) +
          scale_color_gradient(low = 'white', high = 'black') + 
          #scale_color_gradientn(colours = rainbow(5)) + 
          new_scale_colour() + 
          geom_point(data=nodes,alpha=0.7,shape=19, aes(x=x, y=y,size=cluster_size), color = 'grey') + 
          #scale_color_discrete(nodes$cluster_color) + 
          scale_size_area(trans="sqrt", max_size=max_size, breaks = c(5000,10000,15000,20000)) +
          geom_text_repel(data=nodes,aes(x=x,y=y,label=labels),size = label.size) +
          theme_void())
  dev.off()
}

