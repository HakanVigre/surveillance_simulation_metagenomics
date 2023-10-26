#load libraries
source("libraries.R")

#reading xls files generated from dbest analysis

#tet5
{
tet_5<-read_excel("EStart5_Tetra.xlsx")
tet_data<-tet_5

LISTOFALL_ <- list()

for (k in 1:length(tet_data)) {
  if (is.null(tet_data[[k]])) {
    tet_data[[k]]<-NA
  }
  a <- as.data.frame(matrix(tet_5[[k]],nrow=1))
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
no_signal_tet5<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA

Detection <- numeric(nrow(data_v4))
for (k in 1:nrow(data_v4)){
  da <- data_v4[k,]
  Pos <- min(da[da > 0])
  Detection[k] <- min(Pos)
}
Detection_tet5<-Detection 

Negative <- numeric(nrow(data_v4))
for (k in 1:nrow(data_v4)){
  na <- data_v4[k,]
  Neg <- na[na <= 0]
  Negative[k] <- length(Neg)
}
Negative_tet5<-Negative
}

#tet20
{
  tet_20<-read_excel("EStart20_Tetra.xlsx")
  tet_data<-tet_20
  
  LISTOFALL_ <- list()
  
  for (k in 1:length(tet_data)) {
    if (is.null(tet_data[[k]])) {
      tet_data[[k]]<-NA
    }
    a <- as.data.frame(matrix(tet_20[[k]],nrow=1))
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
  no_signal_tet20<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA
  
  Detection <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    da <- data_v4[k,]
    Pos <- min(da[da > 0])
    Detection[k] <- min(Pos)
  }
  Detection_tet20<-Detection 
  
  Negative <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    na <- data_v4[k,]
    Neg <- na[na <= 0]
    Negative[k] <- length(Neg)
  }
  Negative_tet20<-Negative
}

#tet50
{
  tet_50<-read_excel("EStart50_Tetra.xlsx")
  tet_data<-tet_50
  
  LISTOFALL_ <- list()
  
  for (k in 1:length(tet_data)) {
    if (is.null(tet_data[[k]])) {
      tet_data[[k]]<-NA
    }
    a <- as.data.frame(matrix(tet_50[[k]],nrow=1))
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
  no_signal_tet50<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA
  
  Detection <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    da <- data_v4[k,]
    Pos <- min(da[da > 0])
    Detection[k] <- min(Pos)
  }
  Detection_tet50<-Detection 
  
  Negative <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    na <- data_v4[k,]
    Neg <- na[na <= 0]
    Negative[k] <- length(Neg)
  }
  Negative_tet50<-Negative
}

#tet100
{

  tet_100<-read_excel("EStart100_Tetra.xlsx")
  tet_data<-tet_100
  
  LISTOFALL_ <- list()
  
  for (k in 1:length(tet_data)) {
    if (is.null(tet_data[[k]])) {
      tet_data[[k]]<-NA
    }
    a <- as.data.frame(matrix(tet_100[[k]],nrow=1))
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
  no_signal_tet100<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA
  
  Detection <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    da <- data_v4[k,]
    Pos <- min(da[da > 0])
    Detection[k] <- min(Pos)
  }
  Detection_tet100<-Detection 
  
  Negative <- numeric(nrow(data_v4))
  for (k in 1:nrow(data_v4)){
    na <- data_v4[k,]
    Neg <- na[na <= 0]
    Negative[k] <- length(Neg)
  }
  Negative_tet100<-Negative
}

Negative_tet20_feno<-read_excel("E_coli_feno_df_negative.xlsx")
Detection_tet20_feno<-read_excel("E_coli_feno_df_detection.xlsx")

{
#tet5
time_to_detect<-Detection_tet5
scen<-rep(c("tet5"), times=1001)
tet5<-data.frame(scen, time_to_detect)
tet5_<-tet5[order(tet5$time_to_detect,decreasing=FALSE),]
tet5_a<-tet5_[-(1:10),]
summary(tet5_a)

#tet20
time_to_detect<-Detection_tet20
scen<-rep(c("tet20"), times=1001)
tet20<-data.frame(scen, time_to_detect)
tet20_<-tet20[order(tet20$time_to_detect,decreasing=FALSE),]
tet20_a<-tet20_[-(1:10),]
summary(tet20_a)
#tet50
time_to_detect<-Detection_tet50
scen<-rep(c("tet50"), times=1001)
tet50<-data.frame(scen, time_to_detect)
tet50_<-tet50[order(tet50$time_to_detect,decreasing=FALSE),]
tet50_a<-tet50_[-(1:10),]
summary(tet50_a)
#tet100
time_to_detect<-Detection_tet100
scen<-rep(c("tet100"), times=1000)
tet100<-data.frame(scen, time_to_detect)
tet100_<-tet100[order(tet100$time_to_detect,decreasing=FALSE),]
tet100_a<-tet100_[-(1:10),]
summary(tet100_a)
#tetracykline phenotype
time_to_detect<-Detection_tet20_feno$Detection
scen<-rep(c("tet_feno"), times=1000)
tet20_f<-data.frame(scen, time_to_detect)
tet20_f_<-tet20_f[order(tet20_f$time_to_detect,decreasing=FALSE),]
tet20_f_a<-tet20_f_[-(1:10),]
summary(tet20_f_a)

}
tet_time_to_detect<-rbind(tet5, tet20, tet50, tet100, tet20_f)
tet_time_to_detect_a<-rbind(tet5_a, tet20_a, tet50_a, tet100_a, tet20_f_a)

{
  #tet5
  Negative<-Negative_tet5
  scen<-rep(c("tet5"), times=1001)
  tet5<-data.frame(scen, Negative)
  tet5_<-tet5[order(tet5$Negative,decreasing=FALSE),]
  tail(tet5_)
  tet5_a<-tet5_[-(992:1001),]
  summary(tet5_a)
  #tet20
  Negative<-Negative_tet20
  scen<-rep(c("tet20"), times=1001)
  tet20<-data.frame(scen, Negative)
  tet20_<-tet20[order(tet20$Negative,decreasing=FALSE),]
  tail(tet20_)
  tet20_a<-tet20_[-(992:1001),]
  summary(tet20_a)
  #tet50
  Negative<-Negative_tet50
  scen<-rep(c("tet50"), times=1001)
  tet50<-data.frame(scen, Negative)
  tet50_<-tet50[order(tet50$Negative,decreasing=FALSE),]
  tail(tet50_)
  tet50_a<-tet50_[-(992:1001),]
  summary(tet50_a)
  #tet100
  Negative<-Negative_tet100
  scen<-rep(c("tet100"), times=1000)
  tet100<-data.frame(scen, Negative)
  tet100_<-tet100[order(tet100$Negative,decreasing=FALSE),]
  tail(tet100_)
  tet100_a<-tet100_[-(992:1001),]
  summary(tet100_a)
  #tet20_pheno  
  Negative<-Negative_tet20_feno
  scen<-rep(c("tet_feno"), times=1000)
  tet20_f<-data.frame(scen, Negative)
  tet20_f_<-tet20_f[order(tet20_f$Negative,decreasing=FALSE),]
  tail(tet20_f_)
  tet20_f_a<-tet20_f_[-(992:1001),]
  summary(tet20_f_a)
  
}
tet_negative<-rbind(tet5, tet20, tet50, tet100, tet20_f)
tet_negative_a<-rbind(tet5_a, tet20_a, tet50_a, tet100_a, tet20_f_a)


#false negative
b<-ggplot(data=tet_negative_a, aes(x=factor(scen, levels = c("tet5", "tet20", "tet50", "tet100", "tet_feno")),y=Negative))+
  geom_violin(draw_quantiles = c(0.5))+
  geom_jitter(width = 0.1, height = 0.1, color="brown", size=0.4, alpha=0.9)

c<-b + ylim(0, 40)+labs(x="method and number of samples", y="# false positive signals")+
  ggtitle("")    

d<-c+ scale_x_discrete(labels=c("tet5"="mg_5", "tet20"="mg_20", "tet50"="mg_50", "tet100"="mg_100", "tet_feno"="feno_20"))
  
f<-d+ scale_x_discrete(labels=c("tet5"="mg_5", "tet20"="mg_20", "tet50"="mg_50", "tet100"="mg_100", "tet_feno"="pheno_20"))

mynamestheme <- theme(
  panel.background = element_rect(fill = 'white', color = 'grey'),
  panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
  axis.title = element_text(size = (15), colour = "steelblue4"),
  axis.text = element_text(colour = "cornflowerblue", size = (10))
)

plot_no_detect<-f+mynamestheme
plot_no_detect

#time to detect
tet_time_to_detect["time_to_detect"][tet_time_to_detect["time_to_detect"] == 880] <- 150
tet_time_to_detect_a["time_to_detect"][tet_time_to_detect_a["time_to_detect"] == 880] <- 150

b<-ggplot(data=tet_time_to_detect_a, aes(x=factor(scen, levels = c("tet5", "tet20", "tet50", "tet100", "tet_feno")),y=time_to_detect))+
  geom_violin()+
  geom_jitter(width = 0.1, height = 0.1, color="brown", size=0.4, alpha=0.9)


c<-b + ylim(0, 180)+labs(x="", y="# months until detection")+
  ggtitle("")    

d <-c + annotate("text", x = 1.5, y = 145, label = "no detection",color="black", size=5)

e <- d + scale_y_continuous(breaks = c(0,12,24,36,48,60,72,84,96,108,120))
f<-e+ scale_x_discrete(labels=c("tet5"="mg_5", "tet20"="mg_20", "tet50"="mg_50", "tet100"="mg_100", "tet_feno"="pheno_20"))

plot_time_detect<-f+mynamestheme
plot_time_detect




