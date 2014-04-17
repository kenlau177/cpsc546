setwd("C://Users//Ken//Desktop//School//MSC1_Winter//cpsc546//Project//Visualizations")

rm(list=ls())

library(mlbench)
library(ggplot2)
library(gridExtra)
data(Ozone)

process_Ozone = function(x){
  out = x[,-9]
  names(out) = c("Month", "Month_Day", "Week_Day", "Max_O3",
                 "Pressure_Height", "Wind_Speed", "Humidity", 
                 "Temperature", "Inversion_Base_Height", 
                 "Pressure_Gradient", "Inversion_Base_Temperature", 
                 "Visibility")
  temp = out[,-match("Max_O3",colnames(out))]
  out = data.frame("Max_O3"=out$Max_O3, temp)  
  out = out[complete.cases(out),]
  return(out)
}
sqep_loss = function(x, ep){
  out = pmax(0, x^2-ep)
  return(out)
}

dat = process_Ozone(x=Ozone)
x1 = seq(-4, 4, by=.01)
x1df = data.frame(x1, "sqep"=sqep_loss(x1,ep=1))
ep = 30

#plot of mlr fit to all data obs vs pred plot
fit_lm = lm(Max_O3~., data=dat)
gg1.df = data.frame("obs"=dat$Max_O3, "pred"=fit_lm$fitted.values)
gg1 = ggplot(gg1.df, aes(x=obs, y=pred)) + 
  geom_point(size=2) + 
  geom_abline(intercept=0, slope=1, linetype="dashed") + 
  xlab("Observed") + ylab("Predicted") + xlim(c(0,45)) + 
  ylim(c(-15,45)) + 
  ggtitle("Obs vs Predicted for MLR Fit to All Data") + 
  theme(axis.title=element_text(size=16), 
    axis.text=element_text(size=16), title=element_text(size=16))
#ggsave("mlr_show_outlier_plot.pdf", gg, width=6, height=6)

#density plot observation data
gg2 = ggplot(dat, aes(x=Max_O3)) + geom_density(size=1.25) + 
      xlab("Maximum O3") + ylab("Density") + 
      ggtitle("Density Distribution of Maximum O3") + 
      theme(axis.title=element_text(size=16), 
       axis.text=element_text(size=16), title=element_text(size=16))
#ggsave("density_plot_Max_O3.pdf", gg, width=6, height=6)
pdf("outlier_examples_plots.pdf", width=12, height=6, compress=F)
grid.arrange(gg1, gg2, ncol=2)
dev.off()


#squared epsilon insensitive loss
gg = ggplot(x1df, aes(x=x1, y=sqep)) + geom_line(size=1.2) + 
      xlab("Residual") + ylab("Error") + 
      ggtitle("Squared Error Epsilon Insensitive") + 
      theme(axis.title=element_text(size=25), 
      axis.text=element_text(size=18), title=element_text(size=25))
ggsave("squared_error_epsilon_loss.pdf", gg, width=6, height=6)

#RSS vs RSS with epsilon insensitive
rho = seq(5, 16, by=.1)
rss1 = sapply(rho, function(x,dat){
        return(sum((x-dat$Max_O3)^2))}, dat=dat)
rss2 = sapply(rho, function(x,dat,ep){
        return(sum(sqep_loss(x-dat$Max_O3,ep)))}, dat=dat, ep=ep)
rsa = sapply(rho, function(x,dat){
        return(sum(abs(x-dat$Max_O3)))}, dat=dat)

rss1df = data.frame(rho, rss1)
rss2df = data.frame(rho, rss2)
rsadf = data.frame(rho, rsa)

gg1 = ggplot(rss1df, aes(x=rho, y=rss1)) + geom_line(size=1.2) +
        xlab("Constant Estimate") + ylab("Loss") + 
        ggtitle("Sum of Squared Errors") + 
        theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=16), title=element_text(size=16))
gg2 = ggplot(rss2df, aes(x=rho, y=rss2)) + geom_line(size=1.2) + 
        xlab("Constant Estimate") + ylab("Loss") + 
        ggtitle("SSE with Eps Insensitive") + 
        theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=16), title=element_text(size=16))
gg3 = ggplot(rsadf, aes(x=rho, y=rsa)) + geom_line(size=1.2) + 
        xlab("Constant Estimate") + ylab("Loss") + 
        ggtitle("Sum of Absolute Deviations") + 
        theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=16), title=element_text(size=16))
pdf("smoothness.pdf", width=12, height=4, compress=F)
grid.arrange(gg1, gg2, gg3, ncol=3)
dev.off()



