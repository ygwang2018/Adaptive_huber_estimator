library(ggplot2)
library(MASS)
install.packages("gridExtra")   # 如果你没安装过
library(gridExtra)

# Plasma dataframe

dataset1<-read.csv("C:\\Users\\PC\\Desktop\\Polonomial_block_descen_algorithm\\application\\plasma_data.csv",header=F)
y1<-dataset1$V11
df_plasma <- y1
df_plasma <- data.frame(plasma = y1)


# Boston dataframe
df_boston <- data.frame(medv = Boston$medv)

library(ggplot2)
library(MASS)
library(gridExtra)

# Plasma data
df_plasma <- data.frame(plasma = y1)

# Boston data
df_boston <- data.frame(medv = Boston$medv)

p1 <- ggplot(df_plasma, aes(x = plasma)) +
  geom_histogram(bins = 20, fill = "#008080", color = "white", alpha = 0.85) +
  geom_density(color = "black", size = 0.9) +
  labs(title = "(a) Histogram of Plasma Beta-Carotene level",
       x = "Plasma Beta-Carotene Level (ng/ml)",
       y = "Frequency") +
  theme_minimal(base_size = 14)

p2 <- ggplot(df_boston, aes(x = medv)) +
  geom_histogram(bins = 20, fill = "#800080", color = "white", alpha = 0.85) +
  geom_density(color = "black", size = 0.9) +
  labs(title = "(b) Histogram of House Price (medv)",
       x = "Median Value of Owner-Occupied Homes (medv)",
       y = "Frequency") +
  theme_minimal(base_size = 14)

grid.arrange(p1, p2, ncol = 2)
library(ggplot2)
library(gridExtra)

# 打开 EPS 设备
postscript("plasma_boston.eps", 
           width = 10, height = 5, 
           horizontal = FALSE, 
           onefile = FALSE, 
           paper = "special")

# 绘图
grid.arrange(p1, p2, ncol = 2)

# 关闭设备
dev.off()
