#load libraries
source("libraries.R")

#reading csv files

#aph5
{
aph_5<-read_excel("ESTARTAPH5.xlsx")
aph_data<-aph_5

LISTOFALL_ <- list()

for (k in 1:length(aph_data)) {
  if (is.null(aph_data[[k]])) {
    aph_data[[k]]<-NA
  }
  a <- as.data.frame(matrix(aph_5[[k]],nrow=1))
  LISTOFALL_<-append(LISTOFALL_,list(a))
}

ESTART <- rbind.fill(LISTOFALL_)
ESTART[ESTART == 0] <- NA

data_v1 <- data.frame(ESTART)
data_v1[is.na(data_v1)] <-1000 #substitute NA with 1000 so data manipulation is working

data_v2 <- data_v1-120
data_v2[,1]

data_v2_t<-t(data_v2)
data_v3<-(apply(data_v2_t,1,sort))

data_v4_<-data.frame(data_v3)
data_v4<-data.frame(t(data_v4_))
no_signal_aph5<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA

Detection <- numeric(nrow(data_v4))
for (k in 1:nrow(data_v4)){
  da <- data_v4[k,]
  Pos <- min(da[da > 0])
  Detection[k] <- min(Pos)
}
Detection_aph5<-Detection 

Negative <- numeric(nrow(data_v4))
for (k in 1:nrow(data_v4)){
  na <- data_v4[k,]
  Neg <- na[na <= 0]
  Negative[k] <- length(Neg)
}
Negative_aph5<-Negative
}

#total breaks
data_v5<-data_v4
data_v5[data_v5 == 880] <- NA
data_v5[data_v5 == 880] <- NA

#aph20
{
  aph_20<-read_excel("ESTARTAPH20.xlsx")
  aph_data<-aph_20
  
  LISTOFALL_ <- list()
  
  for (k in 1:length(aph_data)) {
    if (is.null(aph_data[[k]])) {
      aph_data[[k]]<-NA
    }
    a <- as.data.frame(matrix(aph_20[[k]],nrow=1))
    LISTOFALL_<-append(LISTOFALL_,list(a))
  }
  
  library(plyr)
  ESTART <- rbind.fill(LISTOFALL_)
  ESTART[ESTART == 0] <- NA
  
  data_v1 <- data.frame(ESTART)
  data_v1[is.na(data_v1)] <-1000 #substitute NA with 1000 so data manipulation is working
  
  data_v2 <- data_v1-120
  data_v2[,1]
  
  data_v2_t<-t(data_v2)
  data_v3<-(apply(data_v2_t,1,sort))
  
  data_v4_<-data.frame(data_v3)
  data_v4<-data.frame(t(data_v4_))
  no_signal_aph20<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA
  
  Detection <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    da <- data_v4[k,]
    Pos <- min(da[da > 0])
    Detection[k] <- min(Pos)
  }
  Detection_aph20<-Detection 
  
  Negative <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    na <- data_v4[k,]
    Neg <- na[na <= 0]
    Negative[k] <- length(Neg)
  }
  Negative_aph20<-Negative
}

#aph50
{
  aph_50<-read_excel("ESTARTAPH50.xlsx")
  aph_data<-aph_50
  
  LISTOFALL_ <- list()
  
  for (k in 1:length(aph_data)) {
    if (is.null(aph_data[[k]])) {
      aph_data[[k]]<-NA
    }
    a <- as.data.frame(matrix(aph_50[[k]],nrow=1))
    LISTOFALL_<-append(LISTOFALL_,list(a))
  }
  
  ESTART <- rbind.fill(LISTOFALL_)
  ESTART[ESTART == 0] <- NA
  
  data_v1 <- data.frame(ESTART)
  data_v1[is.na(data_v1)] <-1000 #substitute NA with 1000 so data manipulation is working
  
  data_v2 <- data_v1-120
  data_v2[,1]
  
  data_v2_t<-t(data_v2)
  data_v3<-(apply(data_v2_t,1,sort))
  
  data_v4_<-data.frame(data_v3)
  data_v4<-data.frame(t(data_v4_))
  no_signal_aph50<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA
  
  Detection <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    da <- data_v4[k,]
    Pos <- min(da[da > 0])
    Detection[k] <- min(Pos)
  }
  Detection_aph50<-Detection 
  
  Negative <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    na <- data_v4[k,]
    Neg <- na[na <= 0]
    Negative[k] <- length(Neg)
  }
  Negative_aph50<-Negative
}

#aph100
{
  aph_100<-read_excel("ESTARTAPH100.xlsx")
  aph_data<-aph_100
  
  LISTOFALL_ <- list()
  
  for (k in 1:length(aph_data)) {
    if (is.null(aph_data[[k]])) {
      aph_data[[k]]<-NA
    }
    a <- as.data.frame(matrix(aph_100[[k]],nrow=1))
    LISTOFALL_<-append(LISTOFALL_,list(a))
  }
  
  ESTART <- rbind.fill(LISTOFALL_)
  ESTART[ESTART == 0] <- NA
  
  data_v1 <- data.frame(ESTART)
  data_v1[is.na(data_v1)] <-1000 #substitute NA with 1000 so data manipulation is working
  
  data_v2 <- data_v1-120
  data_v2[,1]
  
  data_v2_t<-t(data_v2)
  data_v3<-(apply(data_v2_t,1,sort))
  
  data_v4_<-data.frame(data_v3)
  data_v4<-data.frame(t(data_v4_))
  no_signal_aph100<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA
  
  Detection <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    da <- data_v4[k,]
    Pos <- min(da[da > 0])
    Detection[k] <- min(Pos)
  }
  Detection_aph100<-Detection 
  
  Negative <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    na <- data_v4[k,]
    Neg <- na[na <= 0]
    Negative[k] <- length(Neg)
  }
  Negative_aph100<-Negative
}

Negative_gent20<-read_excel("E_facium_feno_df_negative.xlsx")
Detection_gent20<-read_excel("E_facium_feno_df_detection.xlsx")

{
#aph5
time_to_detect<-Detection_aph5
scen<-rep(c("aph5"), times=1001)
aph5<-data.frame(scen, time_to_detect)
#aph20
time_to_detect<-Detection_aph20
scen<-rep(c("aph20"), times=1001)
aph20<-data.frame(scen, time_to_detect)
#aph50
time_to_detect<-Detection_aph50
scen<-rep(c("aph50"), times=1001)
aph50<-data.frame(scen, time_to_detect)
#aph100
time_to_detect<-Detection_aph100
scen<-rep(c("aph100"), times=1001)
aph100<-data.frame(scen, time_to_detect)
#gent20
time_to_detect<-Detection_gent20$Detection
scen<-rep(c("gent20"), times=1001)
gent20<-data.frame(scen, time_to_detect)
head(gent20)

}
aph_time_to_detect<-rbind(aph5, aph20, aph50, aph100, gent20)
#aph_time_to_detect<-rbind(aph5, aph20, aph50, aph100)
head(aph_time_to_detect)
{
  #aph5
  Negative<-Negative_aph5
  scen<-rep(c("aph5"), times=1001)
  aph5<-data.frame(scen, Negative)
  aph5_<-aph5[order(aph5$Negative,decreasing=FALSE),]
  tail(aph5_)
  aph5_a<-aph5_[-(992:1001),]
  summary(aph5_a)
  #aph20
  Negative<-Negative_aph20
  scen<-rep(c("aph20"), times=1001)
  aph20<-data.frame(scen, Negative)
  aph20_<-aph20[order(aph20$Negative,decreasing=FALSE),]
  tail(aph20_)
  aph20_a<-aph20_[-(992:1001),]
  summary(aph20_a)

  #aph50
  Negative<-Negative_aph50
  scen<-rep(c("aph50"), times=1001)
  aph50<-data.frame(scen, Negative)
  aph50_<-aph50[order(aph50$Negative,decreasing=FALSE),]
  tail(aph50_)
  aph50_a<-aph50_[-(992:1001),]
  summary(aph50_a)
  
  #aph100
  Negative<-Negative_aph100
  scen<-rep(c("aph100"), times=1001)
  aph100<-data.frame(scen, Negative)
  aph100_<-aph100[order(aph100$Negative,decreasing=FALSE),]
  tail(aph100_)
  aph100_a<-aph100_[-(992:1001),]
  summary(aph100_a)
  
  #gent20  
  Negative<-Negative_gent20
  scen<-rep(c("gent"), times=1000)
  gent20<-data.frame(scen, Negative)
  gent20_<-gent20[order(gent20$Negative,decreasing=FALSE),]
  tail(gent20_)
  gent20_a<-gent20_[-(992:1001),]
  summary(gent20_a)
  

}
aph_negative<-rbind(aph5, aph20, aph50, aph100, gent20)
aph_negative_a<-rbind(aph5_a, aph20_a, aph50_a, aph100_a, gent20_a)


mynamestheme <- theme(
  panel.background = element_rect(fill = 'white', color = 'grey'),
  panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
  axis.title = element_text(size = (15), colour = "steelblue4"),
  axis.text = element_text(colour = "cornflowerblue", size = (12))
)
#false negative
library(ggplot2)
b<-ggplot(data=aph_negative_a, aes(x=factor(scen, levels = c("aph5", "aph20", "aph50", "aph100", "gent")),y=Negative))+
  geom_violin()+
  geom_jitter(width = 0.1, height = 0.1, color="brown", size=0.4, alpha=0.9)

c<-b + ylim(0, 50)+labs(x="", y="# false positive signals")+
  ggtitle("")    

d<-c+ scale_x_discrete(labels=c("aph5"="mg_5", "aph20"="mg_20", "aph50"="mg_50", "aph100"="mg_100", "gent"="feno_20"))

d+ theme_bw()

plot_false_positive_aph<-d+mynamestheme
plot_false_positive_aph

#time to detect
aph_time_to_detect["time_to_detect"][aph_time_to_detect["time_to_detect"] == 880] <- 150

b<-ggplot(data=aph_time_to_detect, aes(x=factor(scen, levels = c("aph5", "aph20", "aph50", "aph100", "gent20")),y=time_to_detect))+
  #geom_violin(adjust=1, bw=0.05)+
  geom_violin()+
  geom_jitter(width = 0.1, height = 0.1, color="brown", size=0.4, alpha=0.9)

c<-b + ylim(0, 180)+labs(x="", y="# time to detect")+
  ggtitle("")    

d <-c + annotate("text", x = 1.5, y = 145, label = "no detection",color="black", size=5)

d+ theme_bw()

e <- d + scale_y_continuous(breaks = c(0,12,24,36,48,60,72,84,96,108,120))

e
plot_time_detect_aph<-e+mynamestheme
plot_time_detect_aph


