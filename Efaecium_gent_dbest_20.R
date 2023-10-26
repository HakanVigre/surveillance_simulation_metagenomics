#gentamycin resistance i Enterococcus faecalis/faecium

#load libraries
source("libraries.R")

#generating data
par(mar = c(1, 1, 1, 1))

{
  prop_res<-0.06#prop gentamycin Entercoccus faecalis in pigs 2021 from DANMAP
  n_samples_month<-20 #samples of isolates
  months<-120
  prop_res_1_120<-rep(c(prop_res), times=months)
  
  #5% true increase per year
  rate_months<-0.05/12
  month_serie<-seq(1:120)
  rate_sum<-(1+rate_months)^month_serie
  prop_res_121_240<-prop_res*rate_sum
  #generating data for month 1-240
  prop_240<-c(prop_res_1_120,prop_res_121_240)
  
  #series
  #series_n <- 1001#running 1001 series and analyse them in DBEST takes TIME!
  series_n<-2
  set.seed(123)
  t<-rbinom(size=n_samples_month,n=240,prob=prop_240)
  observed <- (matrix(nrow=240,ncol=series_n))
  for (i in 1:series_n){
    A<-rbinom(size=n_samples_month,n=240,prob=prop_240)
    B<-(A/n_samples_month)*100 #formatting into %
    observed[,i]<-B
  }
  observed[,1]
}

###break point analysis####
#Paramters detecting changes >=5%.
#same settings as for metagenom
set.seed(123)
duration = 3
first = prop_res*100#make extreme large so no abrupt changes?
second = prop_res*100#make extreme large so no abrupt changes?
change = prop_res*100*0.05 #5% change of initial proportion

library(DBEST)
feno20 <- list()
for (k in 1:series_n) {
  feno20[[k]] <-  observed[,k]
}

EStart20 <- list()
Esig20 <- list()
z <- length(feno20) # To chanage
for (g in 1:series_n) {
  xdata <- feno20[[g]] # To  chanage
  trytetst <- try(DEX <- DBEST(data= xdata, data.type="non-cyclical",
                               seasonality=12, algorithm="change detection",
                               first.level.shift=first,
                               second.level.shift=second, duration=duration, 
                               distance.threshold="default", alpha=0.05, change.magnitude=change))
  EStart20[[g]] <- DEX$Start
  if (class(trytetst)=="try-error") {
    EStart20[[g]]<-NULL
  }  
  Esig20[[g]] <- DEX$Significance
  if (class(trytetst)=="try-error") {
    Esig20[[g]]<-NULL 
  }
  print(g)
}

feno20<-Map('*',EStart20,Esig20)

LISTOFALL_ <- list()

for (k in 1:length(feno20)) {
  if (is.null(feno20[[k]])) {
    feno20[[k]]<-NA
  }
  a <- as.data.frame(matrix(feno20[[k]],nrow=1))
  LISTOFALL_<-append(LISTOFALL_,list(a))
}

ESTARTFENO20 <- rbind.fill(LISTOFALL_)

ESTARTFENO20[ESTARTFENO20 == 0] <- NA

data_v1 <- data.frame(ESTARTFENO20)
data_v1[is.na(data_v1)] <-1000 #substitute NA with 1000 so data manipulation is working

data_v2 <- data_v1-120

data_v3<-t(apply(data_v2,1,sort))

data_v4<-data.frame(data_v3)
#data_v4 is a fully sorted data frame.

no_signal<-data_v4$X1==880 #X1 = 880 as lowest value is 1000-120 equivalent to NA
summary(no_signal) # RESULT 1 number of time series that have no signal at ALL

#Positives ## RESULT 2 time to first detection AFTER start to increased
Detection <- numeric(nrow(data_v4))
for (k in 1:nrow(data_v4)){
  da <- data_v4[k,]
  Pos <- min(da[da > 0])
  Detection[k] <- min(Pos)
}
Detection

E_faecium_feno_df_detection<-data.frame(Detection)
#if the # is removed the excel sheet in the repository will be overwritten
#write.xlsx(E_faecium_feno_df_detection, "E_faecium_feno_df_detection.xlsx", rowNames=FALSE)


##########################################
#Negatives ## the number of false positives before start to increased
Negative <- numeric(nrow(data_v4))
for (k in 1:nrow(data_v4)){
  na <- data_v4[k,]
  Neg <- na[na <= 0]
  Negative[k] <- length(Neg)
}
Negative

E_facium_feno_df_negative<-data.frame(Negative)
#if the # is removed the excel sheet in the repository will be overwritten
#write.xlsx(E_facium_feno_df_negative, "E_facium_feno_df_negative.xlsx", rowNames=FALSE)

