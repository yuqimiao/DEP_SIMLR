library(plot.matrix)
library(patchwork)
library(reshape2)
mds_list_calc = function(sim1, data_num = 3, dim = 2){
  mds_list = lapply(sim1[1:data_num], function(x){
    dist = dist2(x,x)
    mds = cmdscale(dist,dim)
    return(mds)
  })
  return(mds_list)
}

mds_visual = function(sim = NA,data = NA, dist = NA, true_label = NA,
                      name=NA, data_num = 3, lim = 5000,
                      text_border = F, text_x = NA, text_y = NA,
                      dim = 2){
  if(!is.na(sim)){
    if(is.na(name)){
      name = deparse(substitute(sim))
    }
    true_label = sim$truelabel
    mds = as_tibble(list.rbind(mds_list_calc(sim, data_num,dim = dim))) %>%
      mutate(id = seq_along(V1))%>%
      mutate(data_type = rep(1:data_num, each = nrow(sim$data1)),
             truth = factor(rep(true_label,data_num)))

  }else{
    if(is.na(name)){
      name = deparse(substitute(data))
    }
    if(is.na(dist)){
      dist = dist2(data,data)
    }
    # data should have feature as columns
    mds = as_tibble(cmdscale(dist,dim)) %>%
      mutate(id = seq_along(V1))%>%
      mutate(truth = factor(true_label)) %>%
      mutate(data_type = "data")

  }


  g1 = mds %>%
    ggplot(aes(x = V1, y = V2, color = truth))+
    geom_point()+
    facet_grid(data_type~.)+
    xlim(-lim,lim)+
    ylim(-lim,lim)+
    theme(legend.position = "right")+
    ggtitle(paste("MDS plot for",name, sep = " "))
  g_all = g1
  if(dim == 4){
    g2 = mds %>%
      ggplot(aes(x = V3, y = V4, color = truth))+
      geom_point()+
      facet_grid(data_type~.)+
      xlim(-lim,lim)+
      ylim(-lim,lim)+
      theme(legend.position = "right")+
      ggtitle(paste("MDS plot for",name, sep = " "))
    g_all = g_all+g2
  }

  if(text_border){
    g1 = g1+
      ggrepel::geom_text_repel(data = subset(mds,
                                             (V1> text_x[1] & V1< text_x[2])&(V2> text_y[1] & V2< text_y[2])),
                                aes(x = V1, y = V2,label = id))
    border_id = subset(mds,(V1> text_x[1] & V1< text_x[2])&(V2> text_y[1] & V2< text_y[2])) %>% pull(id)
    print(paste(border_id, collapse = ","))
  }
  return(g_all)
}


heat_nohier = function(s, name = NA, diag = NULL){
  tmp = as.matrix(s)
  if(!is.null(diag)){
    diag(tmp) = diag
  }

  if(is.na(name)){
    name = deparse(substitute(s))
  }
  heatmap(tmp, Colv = NA, Rowv = NA, scale = "column", main = name)
}

heatmap_gg = function(mat, name, discrete = F){
  if(is.null(colnames(mat))){
    rownames(mat) = 1:dim(mat)[1]
    colnames(mat) = 1:dim(mat)[2]
  }
  dat_melt = as_tibble(melt(mat))
  if(!discrete){
    g = ggplot(data=dat_melt,
               aes(x = Var1,y=Var2, fill = value)) +
      geom_tile() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "bottom")+
      ggtitle(name) +
      scale_fill_gradient(high="red", low = "white")
  }else{
    g = ggplot(data=dat_melt,
               aes(x = Var1,y=Var2, fill = factor(value))) +
      geom_tile() +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            legend.position = "bottom")+
      ggtitle(name) +
      scale_fill_discrete(high="red", low = "white")
  }
  g
}


# input: list of clusters, list(method=cluster_label)
# output: colorbar plot
color_bar_cl = function(cluster_ls){
  tib = as_tibble(cluster_ls) %>%
    mutate(id = names(cluster_ls[[1]])) %>%
    pivot_longer(cols = c(part_cimlr, cimlr, SNF),
                 names_to = "method",
                 values_to = "cluster") %>%
    mutate(cluster = as.factor(cluster),
           method = factor(method, levels = c("part_cimlr","cimlr","SNF"))) %>%
    arrange(method,cluster)
  tib = tib %>%
    mutate(id = factor(id, levels = tib$id[1:length(cluster_ls[[1]])]))


  # Stacked
  g = ggplot(tib, aes(fill=cluster, x=id, y=method)) +
    geom_tile()+
    theme(axis.text.x = element_text(size = 3, angle = 45, hjust = 0.8))
  return(g)
}

feature_heatmap = function(mat = matrix(rnorm(50000), nrow = 500),
                           column_id_order = NA,
                           row_id_order = NA,
                           discrete = F,
                           xlab = "sample",
                           ylab = "feature",
                           name = "? feature heatmap"){
  if(is.null(colnames(mat))){
    rownames(mat) = 1:dim(mat)[1]
    colnames(mat) = 1:dim(mat)[2]
  }
  # dat_melt = as_tibble(melt(mat))
  dat_melt = as_tibble(mat) %>%
    mutate(rows = rownames(mat)) %>%
    pivot_longer(cols = colnames(mat),
                 names_to = "columns",
                 values_to = "value")

  if(!is.na(column_id_order)){
    dat_melt = dat_melt %>%
      mutate(columns = factor(columns, levels = column_id_order))
  }

  if(!is.na(row_id_order)){
    dat_melt = dat_melt %>%
      mutate(rows = factor(rows, levels = row_id_order))
  }

  g = ggplot(data=dat_melt,
             aes(x = rows,y=columns, fill = value)) +
    geom_tile() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "right")+
    labs(title = name,
         x = xlab, y = ylab)

  if(!discrete){
    g = g +
      scale_fill_gradient(high="red", low = "white")
  }else{
    g = g +
      scale_fill_discrete(high="red", low = "white")
  }
  return(g)
}



