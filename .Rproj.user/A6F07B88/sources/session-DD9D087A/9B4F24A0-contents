
#' Evaluating and scoring central node microbes
#'
#' @param df1 the topological attributes of the chosen subnetwork nodes (microbes) using five network centrality measures
#'
#' @return This function returns a list with the following components:
#'   \item{nodes_score.csv}{the topological attributes scores and rank of the chosen subnetwork nodes (microbes) using five network centrality measures.}
#' @export
#'
NodeScore<-function(df1){
  ###########################################检查安装包###########################################
  # 需要屏蔽消息的代码
  suppressMessages({
    # 要执行的R表达式
    # 检查并加载readxl包
    if (!require("readxl")) {
      install.packages("readxl")
      library("readxl")
    }

    # 检查并加载ggplot2包
    if (!require("ggplot2")) {
      install.packages("ggplot2")
      library("ggplot2")
    }

    # 检查并加载dplyr包
    if (!require("dplyr")) {
      install.packages("dplyr")
      library("dplyr")
    }
  })
  data1 <- df1[, 2:ncol(df1)]
  # 找到每列的最小值和最大值
  min_vals <- apply(data1, 2, min)
  max_vals <- apply(data1, 2, max)

  # 正向标准化
  data1_std <- t(apply(data1, 1, function(x) {
    (x - min_vals) / (max_vals - min_vals)
  }))


  m <- nrow(data1_std)
  n <- ncol(data1_std)
  data1_value <- as.matrix(data1_std)
  k <- 1 / log(m)

  pij <- apply(data1_value, 2, function(x) x / sum(x))

  test <- pij * log(pij)
  test[is.na(test)] <- 0
  ej <- -k * colSums(test)
  wi <- (1 - ej) / sum(1 - ej)

  R_result <- data.frame(matrix(ncol = ncol(df1)-1, nrow = nrow(df1)))

  for (i in 1:ncol(data1_std)) {
    X <- colnames(data1_std)[i]
    R_result[paste("X", i, ":", X, sep = "")] <- data1_std[, i]
    R_result[paste("R", i, ":", X, sep = "")] <- dense_rank(data1_std[, i])
  }

  R_result<-R_result[,-c(1:ncol(df1)-1)]

  # 选取奇数列乘以wi变量
  odd<-R_result[, seq(2, ncol(R_result), 2)]
  selected_columns <- sapply(1:ncol(odd), function(i) odd[, i] * wi[i])
  selected_columns<-as.data.frame(selected_columns)
  # 按行求和
  sum_row <- rowSums(selected_columns)
  # 计算RSR列的值
  R_result$RSR <- sum_row / n
  R_result$RSR_Rank <- rank(-R_result$RSR)

  #绘制RSR分布表
  RSR <- R_result$RSR
  RSR_RANK_DICT <- data.frame(rank = rank(RSR), index = RSR)
  Distribution <- data.frame(index = unique(RSR))
  Distribution <- arrange(Distribution, index)
  Distribution$f <- table(RSR)
  Distribution$cumsum_f <- cumsum(Distribution$f)
  possible_ranks <- data.frame(index = unique(Distribution$index))
  Distribution$average_rank <- unique(inner_join(possible_ranks, RSR_RANK_DICT,by="index")$rank)
  Distribution$average_rank_percent <- Distribution$average_rank / m

  #最后一项用 修正: 将最后一行最后一列的元素赋值为1 - 1 / (4 * n)
  Distribution[nrow(Distribution), ncol(Distribution)] <- 1 - 1 / (4 * n)

  Distribution$Probit <- 5 + qnorm(Distribution[, ncol(Distribution)])
  rownames(Distribution)<-Distribution$index
  Distribution<-Distribution[,-1]

  r0 <- coef(lm(rownames(Distribution) ~ Probit, data = Distribution))
  summary(lm(rownames(Distribution) ~ Probit, data = Distribution))

  R_result$Probit <- sapply(R_result$RSR, function(item) Distribution[as.character(item), "Probit"])
  R_result$RSR_Regression <- predict(lm(rownames(Distribution) ~ Probit, data = Distribution), newdata = R_result)
  threshold <- predict(lm(rownames(Distribution) ~ Probit, data = Distribution), newdata = data.frame(Probit = c(2, 4, 6, 8)))
  R_result$Level <- cut(R_result$RSR_Regression, threshold, labels = seq(length(threshold) - 1, 1))#较高的数值映射到较小的标签

  node_id <- df1["node_id"]
  R_result2 <- cbind(R_result, node_id)
  R_result3 <- merge(Distribution, R_result2, by = "Probit", all = TRUE)
  R_result4 <- merge(df1, R_result3, by = "node_id", all = TRUE)

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
  write.csv(R_result4, "./data_test/result/nodes_score.csv", row.names = FALSE)

}


