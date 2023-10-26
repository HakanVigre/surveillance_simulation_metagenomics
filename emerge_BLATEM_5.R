#Code used for fitting the concentration, simulating data and plotting 
remove (list = objects() )

#load libraries
source("libraries.R")

reads_seqdepth_vetII_samples <- read_excel("reads_seqdepth_vetII_samples.xlsx")
reads_seqdepth_vetII_samples$blaTEM_CPM

Data_blaTEM <- reads_seqdepth_vetII_samples$blaTEM_CPM
Data_blaTEM <- Data_blaTEM[!Data_blaTEM %in% 0] #To remove the zero's
fln_blaTEM <- fitdist(Data_blaTEM, "lnorm")
fit_blaTEM <- summary(fln_blaTEM)
seqdepth <- reads_seqdepth_vetII_samples$TotalFragmentsN




set.seed(123)

n_Farm <- c(1,1,1,2,2,2,2,3,3,3,4,5,5,6,7,8,9,10,12,14,16,18,21,24,27,32,36,
            42,48,55,63,72,83,95,109,125,144,165,190,218,250,287,329,378,434,
            498,572,657,754,865,993,1140,1309,1503,1725,
            1981,2274,2610,2997,3440)

prevFarm<-n_Farm/4000
prev_pig<-0.1#asuming 10% of pigs in infected farms are carrying the gene
prev_month<-prev_pig*prevFarm

genePrev_CPM <-rlnorm(n=1, meanlog = fit_blaTEM$estimate[1], sdlog = fit_blaTEM$estimate[2])#From the Data



n_samples_month<-5 #samples of isolates
months<-60

#series
#generating the CPM in the technical sample of pooled samples 
set.seed(1234)
series_n <- 1001
meanConc_blaTEM <- data.frame(matrix(nrow=60,ncol=series_n)) 
for (k in 1:series_n) {#
  C <-  prev_month 
  conc_gene<-vector("numeric")
  for (i in 1:60) {
    pigs<-vector("numeric")
    for (j in 1:5){ # Change to number of sample
      inf_pig <- (rbinom(n = 1, size = 1, prob=C[i]))
      konc_pig <- rlnorm(n = 1,meanlog = fln_blaTEM$estimate[1], sdlog = fln_blaTEM$estimate[2])
      frag_pig <- inf_pig*konc_pig
      pigs[j] <- frag_pig
      conc_gene_ <- mean(pigs)
    }
    conc_gene[i] <- conc_gene_
  }
  meanConc_blaTEM[,k] <- conc_gene
}
meanConc_blaTEM #<- meanConc/1e6

#sequencing
library(truncnorm)
set.seed(1234)
blaTEM_count_5 <- matrix(0,ncol=series_n,nrow=60)
for (k in 1:series_n){
  vec=seq(from=1,by=1,length.out = 60)
  #Change to depth based on the money you are using
  depth_ <-rtruncnorm(vec,a=167e6*0.8,b=167e6*1.2,mean=167e6,sd=(167e6-167e6*0.8)/3)
  #20depth_ <-rtruncnorm(vec,a=149e6*0.8,b=149e6*1.2,mean=149e6,sd=(149e6-149e6*0.8)/3)
  #50depth_ <-rtruncnorm(vec,a=114e6*0.8,b=114e6*1.2,mean=114e6,sd=(114e6-114e6*0.8)/3)
  #100depth_ <-rtruncnorm(vec,a=54e6*0.8,b=54e6*1.2,mean=54e6,sd=(54e6-54e6*0.8)/3)
  seqdepth <- round(as.numeric(depth_), 0)
  A <- meanConc_blaTEM[,k]/1e6
  B <- seqdepth
  CPM_Beta <- numeric(60)
  for (m in 1:60) {
    X_Beta <- rbinom(n=1,size=B[m],prob=A[m]) #stochasticity
    #CPM_Beta[m] <- (X_Beta)/B[m]*1e6#CPM
    CPM_Beta[m] <- (X_Beta)#counts
  }
  blaTEM_count_5[,k] <- CPM_Beta
  if(k-1 %% 1){
    print(k)
  }
}

blaTEM_count_5
t<-t(blaTEM_count_5)
d<-which(t > 0, arr.ind=TRUE)  

d2<-data.frame(d)
head(d2)
d3 <- d2[!duplicated(d2$row), ]

d4<-(d3$col)

hist(d4)
dim(d4)

#adding values to vector 
a<-1001-969

a1<-rep(c(150), times=a)
time_to_detect<-c(d4,a1)
hist(time_to_detect)
time_to_detect<-data.frame(time_to_detect)
library("openxlsx")
#write.xlsx(time_to_detect, file = "blatem_time_to_detect_5s.xlsx",rowNames = FALSE)


