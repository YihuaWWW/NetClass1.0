#' Microbiome Subnetwork Classification and Central Bacterial Identification
#'
#' @param community A microbial relative abundance matrix. Each row represents a microbe(species or genus) and each column represents a sample.
#' @param phylum Corresponding relationship of different taxa.
#' @param algorithm Community partitioning algorithm, default is "walktrap". Other optional parameters are "clusters", "edge.betweenness", "label.propagation", "leading.eigenvector", "fastgreedy" or "multilevel"
#' @param nsub The number of candidate subnetworks. Default is 9.
#' @param SubM Key subnetwork ID selected according to the evaluation results of "netcore.R". Default is NA.
#' @param organism The name of the microbe selected according to the evaluation results of "nodescore.R". Default is NA.
#'
#' @return This function returns a list with the following components:
#'   \item{corr_net.pdf}{Correlation network visualization results.}
#'   \item{corrnet_character.csv}{Correlation network topological properties table.}
#'   \item{class_subnetwork.csv}{The result of community partition.}
#'   \item{size_tax_subnetwork.csv}{The size of the subnetworks.}
#'   \item{netclass_allsubnet.pdf}{Community partition visualization results.}
#'   \item{subnet_character.csv}{Candidate subnetworks topological properties table}
#'   \item{selected_subnet.pdf}{Selected subnetwork visualization results}
#'   \item{nodes_character_selectedsubnet.csv}{Topological properties of the nodes in the selected subnetwork.}
#'   \item{centralnode_net.pdf}{Central node network visualization results.}
#'
#' @references Yihua Wang (2023): NetClass: A noncontrol-dependent approach for microbiome subnetwork classification and central bacterial identification.
#' @export
comm <- function(community,
                     phylum,
                     algorithm ="walktrap",
                     nsub=9,
                     SubM =NA,
                     organism = NA) {
  ###########################################检查安装包###########################################
  # 需要屏蔽消息的代码
  suppressMessages({
    # 要执行的R表达式
    # 检查并加载Hmisc包
    if (!require("Hmisc")) {
      install.packages("Hmisc")
      library("Hmisc")
    }

    # 检查并加载lattice包
    if (!require("lattice")) {
      install.packages("lattice")
      library("lattice")
    }

    # 检查并加载survival包
    if (!require("survival")) {
      install.packages("survival")
      library("survival")
    }

    # 检查并加载Formula包
    if (!require("Formula")) {
      install.packages("Formula")
      library("Formula")
    }

    # 检查并加载ggplot2包
    if (!require("ggplot2")) {
      install.packages("ggplot2")
      library("ggplot2")
    }

    # 检查并加载igraph包
    if (!require("igraph")) {
      install.packages("igraph")
      library("igraph")
    }

    # 检查并加载dplyr包
    if (!require("dplyr")) {
      install.packages("dplyr")
      library("dplyr")
    }
  })
  suppressWarnings({

    ####################构建相关性网络#####
    #读取物种与环境因子数据
    rownames(community)=community[,1]
    community=community[,-1]
    com=t(community)
    data=com
    #data=com[c(169:189),]#T
    #data=com[c(148:168),]#S
    #data=com[c(1:21),]#B
    #data=com[c(22:42),]#GCE
    #data=com[c(85:105),]#PCE
    data<-data[,colMeans(data == 0)<1]#delete the column with all 0
    #计算相关性矩阵并筛选p<0.05,rcorr>=0.7的数据，r可以自己设置
    corr=rcorr(data, type="spearman")
    rcorr=corr$r #提取相关系数
    pcorr=corr$P #提取检验结果p值
    #注意要去掉自相关（对角线上的数据）#
    for (i in 1:nrow(rcorr)) {
      for (j in 1:nrow(rcorr)) {
        if (i!=j) {
          if (pcorr[i, j]>0.05) {
            rcorr[i, j]=0
          }
        }
        if (rcorr[i, j]>-0.8 & rcorr[i, j]<0.8) {
          rcorr[i, j]=0
        }
      }
    }
    #构建相关网络
    g1=graph.adjacency(rcorr, mode="undirected", weighted=T, diag=F)
    g2=graph.adjacency(abs(rcorr), mode="undirected", weighted=T, diag=F)
    #删除孤立点
    g1 <- delete.vertices(g1, names(degree(g1)[degree(g1) == 0]))
    g2 <- delete.vertices(g2, names(degree(g2)[degree(g2) == 0]))
    ######plot network######
    g<-g1
    vg<-data.frame("species"=names(V(g)))
    data2_sub<-data[,vg$species]
    m=length(colnames(data2_sub))

    size1=numeric(m)
    for (i in 1:m) {
      size1[i]=sum(data[,i])
    }
    size1=(10000*size1)^0.25
    size=size1
    myshape=rep("circle", m)
    weight=E(g)$weight

    vg<-data.frame("species"=names(V(g)))
    library(dplyr)
    vg_phylum<-inner_join(vg,phylum,by="species")
    unique_phylum<-data.frame("phylum"=unique(vg_phylum$phylum))
    color_num = length(unique_phylum$phylum):1
    #create a color palette of the same size as the number of vertices.
    jet.colors <- colorRampPalette( rainbow( length( unique(color_num) ) ) )
    color_spectrum <- jet.colors( length( unique(color_num ) ) )
    #and over here we map the pallete to the order of values on vertices

    phylum_col<-data.frame("phylum"=unique(vg_phylum$phylum),"color"=color_spectrum)
    vg_color<-full_join(vg_phylum,phylum_col,by="phylum")
    V(g)$color = vg_color$color
    # 检查文件夹是否已经存在
    if (!file.exists("./data_test/result")) {
      # 如果不存在，则递归创建文件夹
      dir.create("./data_test/result", recursive = TRUE)
    } else {
      # 如果已经存在，则打印一条消息
      # cat("文件夹'./data_test/result'已经存在\n")
    }

    # 在文件夹中创建一个新文件
    file_path <- "./data_test/result/new_file.txt"
    if (!file.exists(file_path)) {
      file.create(file_path)
      cat("在文件夹'./data_test/result'中创建了新文件 'new_file.txt'\n")
    } else {
      # cat("文件 'new_file.txt' 已经存在\n")
    }

    par(mar=c(0,0,0,0))
    plot(g, vertex.size=size*2,vertex.label=NA,vertex.color= V(g)$color,
         vertex.shape=myshape, edge.width=3*weight^2,margin= rep(0, 6))
    pdf("./data_test/result/corr_net.pdf")
    par(mar=c(0,0,0,0))
    plot(g, vertex.size=size*2,vertex.label=NA,vertex.color= V(g)$color,
         vertex.shape=myshape, edge.width=3*weight^2,margin= rep(0, 6))
    dev.off()

    #####

    ######网络拓扑性质汇总#####
    t_data<-t(data)
    t_data2<-data.frame("species"=rownames(t_data),t_data)
    rcorr_name<-data.frame("species"=rownames(rcorr))
    t_data3<-inner_join(rcorr_name,t_data2,by="species")
    rownames(t_data3)<-t_data3[,1]
    t_data3<-t_data3[,-1]
    data_rcorr<-t(t_data3)

    #节点数量（number of nodes）和边数量（number of edges）
    nodes_num <- length(V(g2))
    edges_num <- length(E(g2))

    #平均度（average degree）
    #average_degree <- mean(degree(g2))
    #或者，2x边数量/节点数量
    average_degree <- 2*edges_num/nodes_num

    #节点和边连通度（connectivity）
    nodes_connectivity <- vertex.connectivity(g2)
    edges_connectivity <- edge.connectivity(g2)

    #平均路径长度（average path length）
    average_path_length <- average.path.length(g2, directed = FALSE)

    #网络直径（diameter）
    graph_diameter <- diameter(g2, directed = FALSE)

    #图密度（density）
    graph_density <- graph.density(g2)

    #聚类系数（clustering coefficient）
    clustering_coefficient <- transitivity(g2)

    #介数中心性（betweenness centralization)
    betweenness_centralization <- centralization.betweenness(g2)$centralization

    #度中心性（degree centralization）
    degree_centralization <- centralization.degree(g2)$centralization

    #模块性（modularity
    fc <- cluster_fast_greedy(g2)
    modularity <- modularity(g2, membership(fc))

    #选择部分做个汇总输出
    igraph_character <- data.frame(
      nodes_num,    #节点数量（number of nodes）
      edges_num,    #边数量（number of edges）
      graph_diameter,    #网络直径（diameter）
      graph_density,    #图密度（density）
      average_degree,    #平均度（average degree)
      nodes_connectivity,    #节点连通度（vertex connectivity）
      edges_connectivity,    #边连通度（edges connectivity）
      clustering_coefficient,    #聚类系数（clustering coefficient）
      modularity,    #模块性（modularity）
      average_path_length    #平均路径长度（average path length）
    )

    #igraph_character_T<-igraph_character
    #igraph_character_S<-igraph_character
    #igraph_character_B<-igraph_character
    #igraph_character_GCE<-igraph_character
    #igraph_character_PCE<-igraph_character
    write.csv(igraph_character,file="./data_test/result/corrnet_character.csv")#保存构建的相关性网路的拓扑性质



    ######接下来进行子群划分，不同算法具体如下：
    #随机游走
    #其中weights为边的权重，这里即为相关性大小，由于要计算加权概率，负的连接是会有歧义的，
    #所以这里使用g2；step为随机游走的步长，步长越长聚类越粗糙
    #sub2=walktrap.community(g2, weights=E(g2)$weight, step=4)
    #其它备选算法
    #sub2=clusters(g1)#基于点连接的社群发现
    #sub2=edge.betweenness.community(g2, weight=E(g2)$weight, directed=F)#基于度中心性的子群分割
    #sub2=label.propagation.community(g2, weights=V(g2)$weight)#标签传播
    #sub2=leading.eigenvector.community(g2, weights=V(g2)$weight)#特征值
    #sub2=fastgreedy.community(g2, weights=V(g2)$weight)#贪心算法
    #sub2=multilevel.community(g2, weights=V(g2)$weight)#多层次聚类
    # 根据算法类型应用社群发现算法
    sub2 <- switch(algorithm,
                   "walktrap" = walktrap.community(g2, weights=E(g2)$weight, step=4),
                   "clusters" = clusters(g1),
                   "edge.betweenness" = edge.betweenness.community(g2, weight=E(g2)$weight, directed=F),
                   "label.propagation" = label.propagation.community(g2, weights=V(g2)$weight),
                   "leading.eigenvector" = leading.eigenvector.community(g2, weights=V(g2)$weight),
                   "fastgreedy" = fastgreedy.community(g2, weights=V(g2)$weight),
                   "multilevel" = multilevel.community(g2, weights=V(g2)$weight)
    )

    sg2<-data.frame("species"=sub2$names,"sg"=sub2$membership)
    subgroup=vector("list",length(unique(sub2$membership)))
    for(i in 1:length(unique(sub2$membership))){subgroup[[i]] = as.character(sg2[sg2[,2]== i, 1])}

    size_subgroup<-vector(mode="numeric",length=0)
    for(i in 1:length(unique(sub2$membership)))
    {
      size_subgroup[i]<-length(subgroup[[i]])
    }
    size_subgroup2<-data.frame("subgroup"=c(1:length(unique(sub2$membership))),"size"=size_subgroup)
    size_subgroup3<-size_subgroup2[order(-size_subgroup2[,2]),]#选择节点数排名前面的子网络

    write.csv(sg2,file = "./data_test/result/class_subnetwork.csv")##保存划分结果
    write.csv(size_subgroup3,file = "./data_test/result/size_tax_subnetwork.csv")#保存子群大小

    ############子群可视化：############
    data2<-data[,sg2[,1]]#根据列名选择指定列
    rcorr2<-rcorr[sg2[,1],sg2[,1]]

    m=length(colnames(data2))
    t=length(colnames(rcorr2))
    size1=numeric(m)
    for (i in 1:m) {
      size1[i]=sum(data2[,i])
    }
    size1=(10000*size1)^0.25
    size=size1
    myshape=rep("circle", m)
    weight=E(g1)$weight


    vg2<-data.frame("species"=names(V(g2)))
    vg2_phylum<-inner_join(vg2,phylum,by="species")
    unique_phylum<-data.frame("phylum"=unique(vg2_phylum$phylum))
    color_num = length(unique_phylum$phylum):1
    #create a color palette of the same size as the number of vertices.
    jet.colors <- colorRampPalette( rainbow( length( unique(color_num) ) ) )
    color_spectrum <- jet.colors( length( unique(color_num ) ) )
    #and over here we map the pallete to the order of values on vertices

    phylum_col<-data.frame("phylum"=unique(vg2_phylum$phylum),"color"=color_spectrum)
    vg2_color<-full_join(vg2_phylum,phylum_col,by="phylum")
    vg2_color<-inner_join(sg2,vg2_color,by="species")
    V(g2)$color = vg2_color$color

    par(mar=c(0,0,0,0))
    plot(sub2,g2,layout=layout.fruchterman.reingold, vertex.size=1.5*size, vertex.color= vg2_color$color,
         edge.width=4*weight^2,vertex.frame.color=NA,margin= rep(0,6),vertex.label=NA)

    #Display the graph and the legend.
    pdf("./data_test/result/netclass_allsubnet.pdf")
    par(mar=c(0,0,0,0))
    plot(sub2,g2,layout=layout.fruchterman.reingold, vertex.size=1.5*size, vertex.color= vg2_color$color,
         edge.width=4*weight^2,vertex.frame.color=NA,margin= rep(0,6),vertex.label=NA)
    dev.off()


    ######子网络的性质##########
    igraph_character_all<-data.frame(matrix(,nrow=nsub,ncol=10))##nrow=选取子网络的个数；ncol=选取拓扑性质的个数
    for(i in 1:nsub)#选取子网络的个数
    {
      nodes<-V(g2)[sub2$membership==size_subgroup3$subgroup[i]]#T:2,3,22,1,31,20,19,21
      g<-induced_subgraph(g2,nodes)
      nodes_num <- length(V(g))
      edges_num <- length(E(g))
      average_degree <- mean(degree(g))
      nodes_connectivity <- vertex.connectivity(g)
      edges_connectivity <- edge.connectivity(g)
      average_path_length <- average.path.length(g, directed = FALSE)
      graph_diameter <- diameter(g, directed = FALSE)
      graph_density <- graph.density(g)
      clustering_coefficient <- transitivity(g)
      fc <- cluster_fast_greedy(g)
      modularity <- modularity(g, membership(fc))
      igraph_character <- data.frame(
        nodes_num,    #节点数量（number of nodes）
        edges_num,    #边数量（number of edges）
        graph_diameter,    #网络直径（diameter）
        graph_density,    #图密度（density）
        average_degree,    #平均度（average degree)
        edges_connectivity,    #边连通度（edges connectivity）
        nodes_connectivity,    #节点连通度（vertex connectivity）
        clustering_coefficient,    #聚类系数（clustering coefficient）
        modularity,    #模块性（modularity）
        average_path_length    #平均路径长度（average path length）
      )
      igraph_character_all[i,]<-igraph_character[1,]
      colnames(igraph_character_all)<-colnames(igraph_character)
    }
    rownames(igraph_character_all)<-size_subgroup3$subgroup[1:nsub]#选取子网络的个数

    write.csv(igraph_character_all,file = "./data_test/result/subnet_character.csv")
    ############对每个community可以得到子图#######
    if (is.na(SubM)) {

    } else {
      # 执行当zWL不是NA时的操作
      nodes<-V(g2)[sub2$membership==SubM]#T:2,3,22,1,31,20,19,21#参数：这里需要自己选
      g<-induced_subgraph(g2,nodes)
      vg<-data.frame("species"=names(V(g)))

      data2_sub<-data[,vg$species]
      m=length(colnames(data2_sub))

      size1=numeric(m)
      for (i in 1:m) {
        size1[i]=sum(data2[,i])
      }
      size1=(10000*size1)^0.25
      size=size1
      myshape=rep("circle", m)
      weight=E(g)$weight

      vg<-data.frame("species"=names(V(g)))
      vg_phylum<-inner_join(vg,phylum,by="species")
      unique_phylum<-data.frame("phylum"=unique(vg_phylum$phylum))
      color_num = length(unique_phylum$phylum):1
      #create a color palette of the same size as the number of vertices.
      jet.colors <- colorRampPalette( rainbow( length( unique(color_num) ) ) )
      color_spectrum <- jet.colors( length( unique(color_num ) ) )
      #and over here we map the pallete to the order of values on vertices

      phylum_col<-data.frame("phylum"=unique(vg_phylum$phylum),"color"=color_spectrum)
      vg_color<-full_join(vg_phylum,phylum_col,by="phylum")
      V(g)$color = vg_color$color

      par(mar=c(0,0,0,0))
      plot(g, vertex.size=size*10,vertex.label=V(g)$label,vertex.color= V(g2)$color,
           vertex.shape=myshape, edge.width=3*weight^2,margin= rep(0, 6))

      pdf("./data_test/result/selected_subnet.pdf")
      par(mar=c(0,0,0,0))
      plot(g, vertex.size=size*10,vertex.label=V(g)$label,vertex.color= V(g2)$color,
           vertex.shape=myshape, edge.width=3*weight^2,margin= rep(0, 6))
      dev.off()

      ###########节点性质
      if (is.na(SubM)) {

      } else {

        nodes<-V(g2)[sub2$membership==SubM]#T可以自己选哪个子网络
        g<-induced_subgraph(g2,nodes)

        E(g)$corr <- E(g)$weight
        E(g)$weight <- abs(E(g)$weight)

        ##节点特征
        #节点度（Degree）
        V(g)$degree <- degree(g)
        #加权度（Weighted degree）
        V(g)$weight_degree <- strength(g)
        #接近中心性（Closeness centrality）
        V(g)$closeness_centrality <- closeness(g)
        #介数中心性（Betweenness centrality）
        V(g)$betweenness_centrality <- betweenness(g)
        #特征向量中心性（Eigenvector centrality）
        V(g)$eigenvector_centrality <- evcent(g)$vector

        #输出列表
        node_character <- data.frame(
          node_id = V(g)$name,
          degree = V(g)$degree,
          weight_degree = V(g)$weight_degree,
          betweenness_centrality = V(g)$betweenness_centrality,
          closeness_centrality = V(g)$closeness_centrality,
          eigenvector_centrality = V(g)$eigenvector_centrality)

        write.csv(node_character,file="./data_test/result/nodes_character_selectedsubnet.csv")##保存所选子网络节点拓扑性质

        if (is.na(organism)) {

        } else {

          e<-E(g2)[[.inc(organism )]]#T，参数：选一个菌，需要自己选
          hub<-subgraph.edges(g2, e, delete.vertices = TRUE)

          t_data<-t(data)
          t_data2<-data.frame("species"=rownames(t_data),t_data)
          vg<-data.frame("species"=names(V(hub)))
          t_data3<-inner_join(t_data2,vg,by="species")
          rownames(t_data3)<-t_data3[,1]
          t_data3<-t_data3[,-1]
          data2<-t(t_data3)
          m=length(colnames(data2))

          size1=numeric(m)
          for (i in 1:m) {
            size1[i]=sum(data2[,i])
          }
          size1=(10000*size1)^0.25
          size=size1
          myshape=rep("circle", m)
          weight=E(hub)$weight
          vg<-data.frame("species"=names(V(hub)))
          vg_phylum<-inner_join(vg,phylum,by="species")
          unique_phylum<-data.frame("phylum"=unique(vg_phylum$phylum))
          color_num = length(unique_phylum$phylum):1
          #create a color palette of the same size as the number of vertices.
          jet.colors <- colorRampPalette( rainbow( length( unique(color_num) ) ) )
          color_spectrum <- jet.colors( length( unique(color_num ) ) )
          #and over here we map the pallete to the order of values on vertices

          phylum_col<-data.frame("phylum"=unique(vg_phylum$phylum),"color"=color_spectrum)
          vg_color<-full_join(vg_phylum,phylum_col,by="phylum")
          V(hub)$color = vg_color$color

          par(mar=c(0,0,0,0))
          plot(hub, vertex.size=size*6,vertex.color= V(hub)$color,
               vertex.shape=myshape, edge.width=3*weight^2,margin= rep(0, 6),layout=layout_in_circle)

          pdf("./data_test/result/centralnode_net.pdf")
          par(mar=c(0,0,0,0))
          plot(hub, vertex.size=size*6,vertex.color= V(hub)$color,
               vertex.shape=myshape, edge.width=3*weight^2,margin= rep(0, 6),layout=layout_in_circle)
          dev.off()

        }

      }

    }



  })






}
