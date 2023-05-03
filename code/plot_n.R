library(purrr)
library(ggplot2)
load("C:/Users/Yanhang/Desktop/code/result/res_n_normal.rda")
t <- 8
res2 <- list()
for (i in 1:t) {
  temp <- res[[i]]
  temp[1, ] <- temp[1, ] - 50
  temp[5, ] <- temp[5, ] - 5
  temp <- temp[-c(2, 3, 7), ]
  res2[[i]] <- temp[c(1, 3, 2, 4), ]
}

data <- data.frame(matrix(0, nrow = 5*t*4, ncol = 4))
data[, 2] <- rep(seq(300, 300+(t-1)*100, 100), times = 5*4)
data[, 3] <- rep(c("SGLasso", "GBridge", "GEL", "cMCP", "ADSIHT"), each = 4*t)
data[, 4] <- rep(rep(c("SE", "GSE", "MCC", "EE"), each = t), times = 5)
for (i in 1:5) {
  for (j in 1:4) {
    data[((i-1)*t*4+(j-1)*t+1):((i-1)*t*4+j*t), 1] <- unlist(map(res2, function(x) x[j, i]))
  }
}
colnames(data) <- c("value", "Sample", "Methods", "Metric")
data$Metric <- factor(data$Metric, c("SE", "GSE", "MCC", "EE"))
data$Methods <- factor(data$Methods, c("SGLasso", "GBridge", "GEL", "cMCP", "ADSIHT"))
data <- data[-c(1:t), ]
p1 <- ggplot(data, aes(x = Sample, y = value, group = Methods, color = Methods, linetype = Methods, shape = Methods))+geom_line()+geom_point(size = 2)+scale_shape_manual(values=c(15:18, 4))+guides(fill = guide_legend(nrow = 1))+facet_wrap(~Metric, scale = "free", ncol = 1)+
  scale_linetype_manual(values=c("solid", "twodash", "dashed", "twodash", "longdash"))+
  theme_bw()+theme(legend.position = "bottom", legend.key.size = unit(20, "pt"), legend.box.spacing = unit(0, 'pt'))+ylab("")+xlab("Sample Size")+scale_x_continuous(breaks = data$Sample)

load("C:/Users/Yanhang/Desktop/code/result/res_n.rda")
t <- 8
res2 <- list()
for (i in 1:t) {
  temp <- res[[i]]
  temp[1, ] <- temp[1, ] - 50
  temp[5, ] <- temp[5, ] - 5
  temp <- temp[-c(2, 3, 7), ]
  res2[[i]] <- temp[c(1, 3, 2, 4), ]
}

data <- data.frame(matrix(0, nrow = 5*t*4, ncol = 4))
data[, 2] <- rep(seq(300, 300+(t-1)*100, 100), times = 5*4)
data[, 3] <- rep(c("SGLasso", "GBridge", "GEL", "cMCP", "ADSIHT"), each = 4*t)
data[, 4] <- rep(rep(c("SE", "GSE", "MCC", "EE"), each = t), times = 5)
for (i in 1:5) {
  for (j in 1:4) {
    data[((i-1)*t*4+(j-1)*t+1):((i-1)*t*4+j*t), 1] <- unlist(map(res2, function(x) x[j, i]))
  }
}
colnames(data) <- c("value", "Sample", "Methods", "Metric")
data$Metric <- factor(data$Metric, c("SE", "GSE", "MCC", "EE"))
data$Methods <- factor(data$Methods, c("SGLasso", "GBridge", "GEL", "cMCP", "ADSIHT"))
data <- data[-c(1:t), ]
p2 <- ggplot(data, aes(x = Sample, y = value, group = Methods, color = Methods, linetype = Methods, shape = Methods))+geom_line()+geom_point(size = 2)+scale_shape_manual(values=c(15:18, 4))+guides(fill = guide_legend(nrow = 1))+facet_wrap(~Metric, scale = "free", ncol = 1)+
  scale_linetype_manual(values=c("solid", "twodash", "dashed", "twodash", "longdash"))+
  theme_bw()+theme(legend.position = "bottom", legend.key.size = unit(20, "pt"), legend.box.spacing = unit(0, 'pt'))+ylab("")+xlab("Sample Size")+scale_x_continuous(breaks = data$Sample)


library(ggpubr)
p3 <- ggpubr::ggarrange(p2, p1,  labels = c("A", "B"), ncol = 2, nrow = 1, common.legend = T, legend = "bottom", font.label = list(size = 15), legend.grob = get_legend(p1))
p3
ggsave(file = "C:\\Users\\Yanhang\\Desktop\\n.pdf", dpi = 600, p3, height = 10, width = 8)

