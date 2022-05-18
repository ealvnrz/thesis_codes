library(reshape2)
library(ggplot2)
library(latex2exp)
library(tidyverse)


setwd("~")

Beta2_Clayton_shape_05_nu_1<-read.table("./sim_data/nu1_res_s05.txt", quote="\"", comment.char="")
Beta2_Clayton_shape_05_nu_2<-read.table("./sim_data/nu2_res_s05.txt", quote="\"", comment.char="")
Beta2_Clayton_shape_05_nu_6<-read.table("./sim_data/nu6_res_s05.txt", quote="\"", comment.char="")

Beta2_Clayton_shape_15_nu_1<-read.table("./sim_data/nu1_res_s15.txt", quote="\"", comment.char="")
Beta2_Clayton_shape_15_nu_2<-read.table("./sim_data/nu2_res_s15.txt", quote="\"", comment.char="")
Beta2_Clayton_shape_15_nu_6<-read.table("./sim_data/nu6_res_s15.txt", quote="\"", comment.char="")

Beta2_Clayton_shape_25_nu_1<-read.table("./sim_data/nu1_res_s25.txt", quote="\"", comment.char="")
Beta2_Clayton_shape_25_nu_2<-read.table("./sim_data/nu2_res_s25.txt", quote="\"", comment.char="")
Beta2_Clayton_shape_25_nu_6<-read.table("./sim_data/nu6_res_s25.txt", quote="\"", comment.char="")

mean=0.2
mean1=-0.2
scale=0.2
shape=1.5
true_values=c(mean,mean1,scale,shape)

# Data
data1<-Beta2_Clayton_shape_15_nu_1
data2<-Beta2_Clayton_shape_15_nu_2
data3<-Beta2_Clayton_shape_15_nu_6

# Bias
dat=data3
bias=sapply(dat, mean)-true_values
round(bias,4)
mse=sapply(dat,var)+bias^2
round(mse,4)

##################################################################
##                           Boxplots                           ##
##################################################################

### mean
d1<-data1[,1]-mean(data1[,1])
d2<-data2[,1]-mean(data2[,1])
d3<-data3[,1]-mean(data3[,1])
p<-data.frame(d1,d2,d3)
names(p)<-c("1","2","6")
df<-melt(p)
q1 <- ggplot(df, aes(variable, value)) + geom_boxplot() + theme_bw() + labs(x = TeX(r'($\nu$)'), y = "") +
  theme(plot.title = element_text(hjust = 0.5))

### mean1
d1<-data1[,2]-mean(data1[,2])
d2<-data2[,2]-mean(data2[,2])
d3<-data3[,2]-mean(data3[,2])
p2<-data.frame(d1,d2,d3)
names(p2)<-c("1","2","6")
df2<-melt(p2)
q2<- ggplot(df2, aes(variable, value)) + geom_boxplot() + theme_bw() + labs(x = TeX(r'($\nu$)'), y = "") +
  theme(plot.title = element_text(hjust = 0.5))

### scale
d1<-data1[,3]-mean(data1[,3])
d2<-data2[,3]-mean(data2[,3])
d3<-data3[,3]-mean(data3[,3])
p3<-data.frame(d1,d2,d3)
names(p3)<-c("1","2","6")
df3<-melt(p3)
q3<- ggplot(df3, aes(variable, value)) + geom_boxplot() + theme_bw() + labs(x = TeX(r'($\nu$)'), y = "") +
  theme(plot.title = element_text(hjust = 0.5))

# shape
d1<-data1[,4]-mean(data1[,4])
d2<-data2[,4]-mean(data2[,4])
d3<-data3[,4]-mean(data3[,4])
p4<-data.frame(d1,d2,d3)
names(p4)<-c("1","2","3")
df4<-melt(p4)
q4<- ggplot(df4, aes(variable, value)) + geom_boxplot() + theme_bw() + labs(x = TeX(r'($\nu$)'), y = "") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot 
library(gridExtra)
grid.arrange(q1,q2,q3,q4,ncol=2,nrow=2)

