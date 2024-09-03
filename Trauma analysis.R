# 数据准备及预处理--------------------------------
data_value=read.csv("C:/Users/Yuki/OneDrive/代码程序/R/immue analysis/数据/数据2/data_train.csv",
                    header = TRUE)
data_value$ISS.levelcha=0
data_value$ISS.levelcha[which(data_value$ISS.level==2)] <- "Ⅰ"
data_value$ISS.levelcha[which(data_value$ISS.level==3)] <- "Ⅱ"
data_value$ISS.levelcha[which(data_value$ISS.level==4)] <- "Ⅲ"
data_value$ISS.levelcha[which(data_value$ISS.level==1)] <- "Ctrl"
data_value$ISS.levelcha=as.factor(data_value$ISS.levelcha)

data_test=read.csv("C:/Users/Yuki/OneDrive/代码程序/R/immue analysis/数据/data_test.csv",
                    header = TRUE)
data_test$ISS.levelcha=0
data_test$ISS.levelcha[which(data_test$ISS.level==2)] <- "Ⅰ"
data_test$ISS.levelcha[which(data_test$ISS.level==3)] <- "Ⅱ"
data_test$ISS.levelcha[which(data_test$ISS.level==4)] <- "Ⅲ"
data_test$ISS.levelcha[which(data_test$ISS.level==1)] <- "Ctrl"
data_test$ISS.levelcha=as.factor(data_test$ISS.levelcha)

#数据基线统计--------------------------

data1=subset(data_value, ISS.levelcha=='Ⅰ', select=Gender:ISS.levelcha)
data2=subset(data_value, ISS.levelcha=='Ⅱ', select=Gender:ISS.levelcha)
data3=subset(data_value, ISS.levelcha=='Ⅲ', select=Gender:ISS.levelcha)
data4=subset(data_value, ISS.levelcha=='Ctrl', select=Gender:ISS.levelcha)

table(data1$Gender)
mean(data1$Age)
max(data1$Age)
min(data1$Age)

mystatistic <- function(x, na.omit=FALSE){
  if(na.omit)
    x <- x[!is.na(x)]  #如果x是非空值，则赋值给x
  m <- mean(x)
  n <- length(x)
  s <- sd(x)
  skew <- sum((x-m)^3/s^3)/n
  kurt <- sum((x-m)^4/s^4)/n-3
  return(c(n=n, mean=m, sd=s, skew=skew, kurt=kurt))
}

vars <- c("IL.6","PCT","CRP","LMR","PLR","NLR","Alb","WBC","Platelet","Neutrophil","Lymphocyte","Monocyte",
          "RDW.CV","RDW.SD","pH","PCO2","PO2","Lactic.acid","SaO2" )
vars <- c("IL.6","PCT","CRP","LMR","PLR","NLR","Alb","WBC","Platelet","Neutrophil","Lymphocyte","Monocyte",
          "RDW.CV","RDW.SD" )

datastati=round(sapply(data4[vars], mystatistic),2)
datastati=t(datastati)
datastati=data.frame(datastati)
result=as.data.frame(paste(datastati$mean,datastati$sd,sep="±"))

write.csv(result,"result.csv",fileEncoding = "GB18030", row.names = T)

#数据相关性
###################
data=data_value[,c(-1,-19)]        

install.packages('Hmisc')
library(Hmisc)
df_rcorr<-rcorr(as.matrix(data),type = "spearman")

install.packages('corrplot')
library(corrplot)
corrplot(df_rcorr$r,#数据
         type="lower",#可选择展示方式，"full", "lower", "upper"
         tl.col ="black",#文本颜色
         #tl.srt = 90#标签旋转
         tl.pos="lt"
)
corrplot(df_rcorr$r,add=TRUE, type="upper", method="number",order="original",diag=F,tl.pos="n", cl.pos="n")


install.packages('PerformanceAnalytics')
library(PerformanceAnalytics)
chart.Correlation(data_value[,1:15],histogram = T,pch="+")


library(PerformanceAnalytics)

mychart.Correlation <- function (R, histogram = TRUE, method = c("pearson",
                                                                 "kendall","spearman"), ...)
{
  x = checkData(R, method = "matrix")
  if (missing(method))
    method = method[1]
  cormeth <- method
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs",
                        method = cormeth, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor))
      cex <- 3.5   #字体大小
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***",
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * 0.8)
    text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "light grey", probability = TRUE, axes = FALSE,
         main = "", breaks = "FD")   #绘制对角线柱状图底部x分布图
    lines(density(x, na.rm = TRUE), col = 2, lwd = 2)
    rug(x)  #rug()是在坐标轴上标出元素出现的频数
  }
  
  if (histogram)
    pairs(x, gap = 0, lower.panel = mypanel.smooth,  upper.panel = panel.cor,
          diag.panel = hist.panel)
  # pairs(x, gap = 0, lower.panel = mypanel.smooth,upper.panel = panel.cor)
  else pairs(x, gap = 0, lower.panel = mypanel.smooth, upper.panel = panel.cor)
}
windows()
mychart.Correlation(data,histogram = T,pch="+")

#相关性和弦图------------------------------------------------------------------
install.packages("circlize")
install.packages("corrplot")
install.packages("psych")
install.packages("polycor")
install.packages("vcd")
install.packages("circlize")
install.packages("ltm")
library(polycor)
library(vcd)
library(circlize)
library(corrplot)
library(psych)
library(circlize)
library(ltm)
#数据准备
data=rbind.data.frame(data1,data2,data3,data4)
data=data[,c(-17,-18)]
data$ISS.level=as.factor(data$ISS.level)
data=as.data.frame(data)
# 计算创伤类别与每个指标的相关性

cor_trauma <- sapply(levels(data$ISS.level), function(t) {
  sapply(1:16, function(i) {
    biserial.cor(data[[i]], as.numeric(data$ISS.level == t), use = "complete.obs")
  })
})

colnames(cor_trauma) <- c("Ctrl","Ⅰ","Ⅱ","Ⅲ" )
rownames(cor_trauma) <- c("Gender","Age","IL.6","PCT","CRP","LMR","PLR","NLR","Alb","WBC","Platelet","Neutrophil","Lymphocyte","Monocyte",
                          "RDW.CV","RDW.SD" )

# 将相关性矩阵转换为长格式
cor_long <- as.data.frame(as.table(cor_trauma))

# 过滤掉弱相关（可以设置一个阈值，比如0.3）
cor_long <- cor_long[abs(cor_long$Freq) > 0.3,]

# 检查并清理 factors 向量
factors <- unique(c(cor_long$Var1, cor_long$Var2))
factors <- factors[!is.na(factors) & factors != ""]

# 确保 factors 没有空值或缺失值
print(factors)

# 生成颜色向量
colors <- rainbow(length(factors))

# 绘制和弦图
circos.clear()
circos.par(gap.degree = 10)
circos.initialize(factors = factors, xlim = c(0, 1))


# 添加和弦图轨道
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

# 添加和弦
for (i in 1:nrow(cor_long)) {
  color <- ifelse(cor_long$Freq[i] < 0, "#D6604DFF", "#4393C3FF")
  width <- abs(cor_long$Freq[i]) * 8  # 根据相关性调整线宽
  circos.link(cor_long$Var1[i], 0.5, cor_long$Var2[i], 0.5, col = color, lwd = width)
}


#箱型图---------------------------------------------------------------------------
install.packages('tidyverse')
install.packages('cowplot')
install.packages('ggthemes')
install.packages('scales')
install.packages('ggsci')
install.packages('ggpubr')
install.packages('hrbrthemes')


library(tidyverse)
library(hrbrthemes)
library(viridis)
library(cowplot)
library(ggthemes)
library(ggplot2)

require(cowplot)
require(tidyverse)
require(ggplot2)
require(ggsci)
require(ggpubr)


data=rbind.data.frame(data4,data1,data2,data3)
data=data[,c(-17)]
data$ISS.levelcha=as.factor(data$ISS.levelcha)
data=rename(data,"ISS.level"="ISS.levelcha")
#箱型图1-4
my_comparisons = list( c("Ctrl", "Ⅰ"), c("Ctrl", "Ⅱ"), c("Ctrl", "Ⅲ") )
box1=ggboxplot(data, x = "ISS.level", y = "IL.6",
          color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                ref.group = "Ctrl", size=7,label.y =c(300,400,500) )+theme(legend.position = "none")+ylab("IL-6")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,700)

box2=ggboxplot(data, x = "ISS.level", y = "PCT",
               color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(8,10,13) )+theme(legend.position = "none")+ylab("PCT")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,16)

box3=ggboxplot(data, x = "ISS.level", y = "CRP",
               color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(270,300,330) )+theme(legend.position = "none")+ylab("CRP")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,360)


box4=ggboxplot(data, x = "ISS.level", y = "Alb",
               color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(50,57,63) )+theme(legend.position = "none")+ylab("Alb")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(25,70)


f=plot_grid(box1, box2, box3, box4, rel_heights = c(1, 1), align = 'h' )
windows()
plot (f)




#箱型图5-8
box5=ggboxplot(data, x = "ISS.level", y = "WBC",
               color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(36,41,46) )+theme(legend.position = "none")+ylab("WBC")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,50)

box6=ggboxplot(data, x = "ISS.level", y = "PLR",
               color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(480,580,670) )+theme(legend.position = "none")+ylab("PLR")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,800)

box7=ggboxplot(data, x = "ISS.level", y = "NLR",
               color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(38,45,52) )+theme(legend.position = "none")+ylab("NLR")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,60)


box8=ggboxplot(data, x = "ISS.level", y = "LMR",
               color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(11.5,14,17) )+theme(legend.position = "none")+ylab("LMR")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,20)

f=plot_grid(box5, box6, box7, box8,align = 'h' )
windows()
plot (f)

#箱型图9-12
box9=ggboxplot(data, x = "ISS.level", y = "Platelet",
               color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(490,550,630) )+theme(legend.position = "none")+ylab("Platelet")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(50,700)

box10=ggboxplot(data, x = "ISS.level", y = "Neutrophil",
                color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(32,39,44) )+theme(legend.position = "none")+ylab("Neutrophil")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,50)

box11=ggboxplot(data, x = "ISS.level", y = "Lymphocyte",
                color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(3.8,4.8,5.5) )+theme(legend.position = "none")+ylab("Lymphocyte")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,6.5)


box12=ggboxplot(data, x = "ISS.level", y = "Monocyte",
                color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(2.4,2.9,3.2) )+theme(legend.position = "none")+ylab("Monocyte")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(0,3.5)


f=plot_grid(box9, box10, box11, box12, align = 'h' )
windows()
plot (f)

#箱型图13-14
box13=ggboxplot(data, x = "ISS.level", y = "RDW.CV",
                color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(20,21,22) )+theme(legend.position = "none")+ylab("RDW-CV")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(10,24)

box14=ggboxplot(data, x = "ISS.level", y = "RDW.SD",
                color = "ISS.level", palette = "jama",width = 0.4,lwd= 2,bxp.errorbar = TRUE,bxp.errorbar.width=0.4)+
  stat_summary(fun.y = mean, geom = "point", shape = 18, size=3, color = "red")+
  theme_bw(base_family = "Times_New_Roman")+stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test",
                                                               ref.group = "Ctrl", size=7,label.y =c(55,58,62))+theme(legend.position = "none")+ylab("RDW-SD")+xlab("ISS level")+
  theme(axis.text = element_text(size = 30),axis.title = element_text(size = 35),panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"))+ylim(36,68)


f=plot_grid(box13, box14)
windows()
plot (f)


# Lasso回归--------------------------------------------------------------------------------------

install.packages('glmnet')
install.packages('lars')
library(lars) 
library(glmnet)
# 示例数据准备

dat=read.csv("C:/Users/hyt95/OneDrive/代码程序/R/immue analysis/数据/data_train.csv",header = TRUE)
dat=dat[,2:18]

y<-as.numeric(dat[,17])
x<-as.matrix(dat[,1:16])


lasso_model = glmnet(x, y, family="multinomial", alpha=1, intercept = T) #这里alpha=1为LASSO回归，如果等于0就是岭回归
print(lasso_model)

plot(lasso_model, xvar="lambda", label=TRUE)


#交叉验证
windows()
set.seed(1234)
cv_model=cv.glmnet(x,y,alpha=1, nfolds = 10)
plot(cv_model)


#选择最佳模型
lambda_min=cv_model$lambda.min
print(lambda_min) 
lambda_1se=cv_model$lambda.1se
print(lambda_1se)

best_model <- glmnet(x, y, alpha = 1, lambda = cv_model$lambda.min)
coef(best_model)


#可视化图形输出
library(tidyverse)
library(broom)

tidy_df <- broom::tidy(lasso_model)
tidy_cvdf <- broom::tidy(cv_model)
tidy_df=tidy_df[2887:3630,]
library(ggplot2)
library(RColorBrewer)


mypalette <- c(brewer.pal(10,"Paired"),brewer.pal(11,"BrBG"),brewer.pal(11,"Spectral"))
p1=ggplot(tidy_df, aes(lambda, estimate, group = term, color = term)) +
  geom_line(size=1.2)+
  geom_hline(yintercept = 0)+
  scale_x_log10(name = "Log(λ)")+
  ylab("Coefficients")+
  scale_color_manual(name="Variable",values = mypalette)+ theme_bw(base_family = "Times_New_Roman")+
   theme(axis.text = element_text(size = 20),axis.title = element_text(size = 25),
         panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
         legend.text=element_text(size=14),legend.title=element_text(size=15))
 
 

  p2=ggplot(tidy_cvdf)+ geom_line(aes(lambda,estimate),colour="black")+geom_point(aes(lambda,estimate),colour="red")+
    geom_errorbar(data = tidy_cvdf, aes(x=lambda,ymin=conf.low,ymax=conf.high))+
    scale_x_log10(name = "Log(λ)")+
  ylab("Mean-Squared Error")+geom_vline(xintercept=c(lambda_1se,lambda_min),linetype="dotted",size=2)+theme_bw(base_family = "Times_New_Roman")+
    theme(axis.text = element_text(size = 20),axis.title = element_text(size = 25),panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
  
library(cowplot)
  f1=plot_grid(p1, p2, rel_heights = c(2, 1), labels = c("a","b","c","d"),label_size = 30, label_fontfamily = "Times_New_Roman",align = 'h' )
  windows()
  plot (f1)




#根据lambda值，确定那些变量应该被保留
coef_cv=coef(cv_model,         #拟合的Lasso模型
            s=lambda_min) #s为lambda大小，其越大表示模型正
coef_cv

#根据既定lambda下有意义的变量提取出来
coef_cv=as.matrix(coef_cv)#转为矩阵
coef_cv=data.frame(coef_cv)#转为数据框
coef_cv

coef_cv$OR=exp(coef_cv$s1)#计算每一个变量的OR

nonzero_vars=rownames(coef_cv[coef_cv$OR!=1,])#OR不为1就是筛选的变量

nonzero_vars=nonzero_vars[2:13]#去除截距项
nonzero_vars<- na.omit(nonzero_vars)

#构建筛选出的lasso矩阵
lasso_data= data.frame(y,data[,c(nonzero_vars)])
write.csv(lasso_data,"lasso_data.csv")                          

  
#==========数据划分================================= 
 train=data_value[,c(-1,-2,-11,-12,-15,-18)]
 test=data_test 
 xtest=test[,1:12]
 ytest=test[,13]  
  
#数据平衡————————————————————————————————————————————————————  
  library(UBL)
  
  table(train$ISS.levelcha)
  #set.seed(1234)
  balanced_data <- SmoteClassif(ISS.levelcha ~ ., train, C.perc = list(ctrl=1,Ⅰ=0.6,Ⅱ=1.2,Ⅲ=2.1))
  
  table(balanced_data$ISS.levelcha)
  
#==============Workflow==============    
  install.packages('tidymodels') 
  install.packages('bonsai') 
  install.packages('parsnip') 
  install.packages('discrim') 
  install.packages('baguette') 
  library(tidymodels)
  library(bonsai)
  library(parsnip)
  library(discrim)
  library(baguette) 
 
  install.packages('xgboost') 
  install.packages('naivebayes') 
  install.packages('kknn') 
  install.packages('ranger') 
  install.packages('kernlab') 
  install.packages('lightgbm') 
  library(xgboost)  
  library(naivebayes) 
  library(kknn) 
  library(ranger) 
  library(kernlab) 
  library(lightgbm) 
 
  install.packages('stacks') 
  install.packages('finetune') 
  install.packages('rules') 
  library(stacks)
  library(finetune)
  library(rules)
  library(RColorBrewer)


  #### 
  bagtree_spec <-
    bag_tree() %>%
    set_engine('rpart') %>%
    set_mode('classification')
  
  xgboost_spec <-
    boost_tree() %>%
    set_engine('xgboost') %>%
    set_mode('classification')
  
    lightgbm_spec <- boost_tree() %>%
    set_engine('lightgbm') %>%
    set_mode('classification')
  
 
  cart_spec <-
    decision_tree() %>%
    set_engine('rpart') %>%
    set_mode('classification')
  
  
  naivebayes_spec <-
    naive_Bayes() %>%
    set_engine('naivebayes')
  
  kknn_spec <-
    nearest_neighbor() %>%
    set_engine('kknn') %>%
    set_mode('classification')
  
  randomforest_spec <-
    rand_forest() %>%
    set_engine('ranger') %>%
    set_mode('classification')
  
  svm_spec <-
    svm_linear() %>%
    set_engine('kernlab') %>%
    set_mode('classification')
  

  # 10 折交叉验证
  
 # folds <- createFolds(factor(train$Injury.degree), k = 10, list = FALSE)
  
  folds <- vfold_cv(balanced_data,v=10,strata = ISS.levelcha)
 
  # 设置recipes
  rec <- recipe(ISS.levelcha~.,data = balanced_data)
  
  wf <- workflow_set(preproc = list(basic=rec),
                     models = list(bagtree=bagtree_spec,
                                   xgboost=xgboost_spec,
                                   dtree=cart_spec,
                                   naivebayes=naivebayes_spec,
                                   knn=kknn_spec,
                                   rf=randomforest_spec,
                                   svm=svm_spec,
                                   lightgbm=lightgbm_spec
                     ))
  wf
 
   wf <- wf %>% 
    workflow_map("fit_resamples",
                 resamples=folds) 
   
 wf$wflow_id=c("Bagged tree" ,"XGBoost","Decision tree","Native Bayes","KNN","RF", "SVM",  "LightGBM" )
 
 Mset=rank_results(wf,rank_metric = "roc_auc") %>% 
    filter(.metric=="roc_auc") 

  mypalette <- c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF") 
  Mset$wflow_id <- factor(Mset$wflow_id, levels = c("LightGBM", "RF", "XGBoost", "KNN", "Bagged tree", "SVM", "Native Bayes", "Decision tree"))
  ggplot(Mset,aes(rank,mean,color = wflow_id)) +
    scale_color_manual(values = mypalette)+ geom_boxplot(size=1.7,width=0.7)+
    geom_errorbar(aes(x=rank, ymin=mean-std_err, ymax=mean+std_err),size=0.8,width=0.4)+xlab("Model rank")+
    ylab("Accuracy rate")+geom_text(aes(x=rank,y = mean - 1/90,label = wflow_id),size=5.5,angle = 0, hjust = 0.5)+
    theme_bw(base_family = "Times_New_Roman")+ 
    theme(axis.text = element_text(size = 20),axis.title = element_text(size = 25),
          panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
          legend.text=element_text(size=14),legend.title=element_text(size=15))+
    theme(legend.position = "none")+scale_y_continuous(breaks=seq(0.55,1.05,0.02),labels = scales::percent)+
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6,7, 8))
 
 
  rf_wf<- 
    wf %>% 
    extract_workflow("LightGBM")
  rf_wf
  
  rf_wf_fit <- 
    rf_wf %>% 
    finalize_workflow(tibble(prod_degree = 1)) %>% 
    fit(data = balanced_data)
  
  rf_wf_fit
  
#混淆矩阵----------------------------------------------------------------- 
  pred_wd2 <- rf_wf_fit %>% predict(xtest,type="class")
  combin_data <-  test %>% 
    select(ISS.levelcha) %>% 
    bind_cols(pred_wd2)
  cm= combin_data %>%
    conf_mat(ISS.levelcha,.pred_class) 
  cm
  
  install.packages('cvms')
  install.packages('tidyverse')
  library(tidyverse)
  library(cvms)
  
  # 假设 'cm' 是一个混淆矩阵对象，需要转换为数据框以便使用geom_text
  cm_data <- as.data.frame(cm$table)
  
  # 生成热图并添加数字
  library(ggplot2)
  install.packages('ggfortify')
  library(ggfortify)
  
  ggplot(cm_data, aes(Truth,Prediction, fill= Freq)) +
    geom_tile() + geom_text(aes(label=Freq), size = 8) +
    scale_fill_gradient(low="white", high="#009194") +
    labs(x = "Reference",y = "Prediction") +
   # theme_minimal(base_family = "Times New Roman") + 
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 25),
     # legend.position = "none"
    )
    
#AUC曲线-----------------------------------------------------------------
   
   pred_wd <- rf_wf_fit %>% predict(xtest,type="prob")
  combin_data1 <-  test %>% 
    select(ISS.levelcha) %>% 
    bind_cols(pred_wd)
  combin_data1=cbind(combin_data, pred_wd)
    
  #计算auc
  # 定义计算每个分类AUC的函数
  library(pROC)
  calculate_auc <- function(data, label) {
    roc_obj <- roc(data$ISS.levelcha == label, data[[paste0(".pred_", label)]])
    auc(roc_obj)
  }
  labels <- c("Ⅰ", "Ⅱ", "Ⅲ", "Ctrl")
  aucs <- sapply(labels, function(label) calculate_auc(combin_data1, label))
  names(aucs) <- labels 
  
  sensitivt=sensitivity(combin_data1,ISS.levelcha,.pred_class)
  sepecify=specificity(combin_data1,ISS.levelcha,.pred_class)
  acu=accuracy(combin_data1,ISS.levelcha,.pred_class)
  auc=roc_auc(combin_data1,ISS.levelcha,.pred_Ⅰ:.pred_Ctrl)
  
 
  roc_curve(combin_data1,ISS.levelcha,.pred_Ⅰ:.pred_Ctrl) %>% 
    autoplot()
  roc= roc_curve(combin_data1,ISS.levelcha,.pred_Ⅰ:.pred_Ctrl)
  
 install.packages("ggalt")
 install.packages("colorspace")
  library(ggalt)
  library(colorspace)
  library(ggplot2)
 mypalette <- c( "#019875FF", "#99B898FF", "#C0392BFF", "#96281BFF") 
 ggplot(roc,aes(x = 1-specificity, y = sensitivity, color=.level)) +
    stat_xspline(spline_shape = 1,size=2)+geom_point(size=1)+
      geom_abline(slope = 1, intercept = 0, color = "grey10",linetype = 2) +
   scale_color_manual(values = mypalette)+# theme_bw()+theme(
      #strip.background = element_rect(
         #fill = "white"))+
    scale_x_continuous(name  = "1-Specificity",limits = c(0,1))+
    scale_y_continuous(name = "Sensitivity",limits = c(0,1))+ 
    theme(plot.title = element_text(hjust = 0.5,size = 24),
          axis.text=element_text(size=17),
          axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24)) 
  
  
  #模型解释可视化----------------------------------------------------------------------------
  install.packages('DALEXtra')  
  library(DALEXtra)
  wf_explain <- explain_tidymodels(rf_wf_fit,data=train[,1:12],y=train[,13],label = "LightGBM") 
  
  #局部模型解释提供有关单个观测值的预测的信息
  wf_bd <- predict_parts(wf_explain,
                          new_observation = train[18,])
   windows()
   plot(wf_bd,cex.axis = 3, cex.lab = 1.2, cex.main = 1.5, cex.sub = 1.2)
  
  
  #shap值查看各变量的对模型预测的贡献
   set.seed(1801)
   wf_shap <- 
     predict_parts(
       explainer = wf_explain, 
       new_observation = train[18,], 
       type = "shap",
       B = 20
     )
   windows()
   plot( wf_shap)
   
   
   #全局模型解释  
   set.seed(1804)
   vip_rf <- model_parts(wf_explain,type = "raw")
   windows()
   plot(vip_rf) 
   
   #局部模型解释来构建全局模型解释，比如部分依赖图（PDP）。
   
   library(DALEX)
   windows()
   pdp_rf <- model_profile(wf_explain,variables = "Neutrophil")
   plot(pdp_rf,geom = "profiles")
   
   
  