#Code used for fitting the concentration, simulating data and plotting 
remove (list = objects() )

#load libraries
source("libraries.R")

#A: estimating presens of aph genes"

#reading true data
reads_seqdepth_vetII_samples <- read_excel("reads_seqdepth_vetII_samples.xlsx")

aph_2_ia2<-reads_seqdepth_vetII_samples$`aph(2__)-Ia_2_AP009486`
aph_2_ia3<-reads_seqdepth_vetII_samples$`aph(2__)-Ia_3_AJ536195`
aac_6_aph2<-reads_seqdepth_vetII_samples$`aac(6_)-aph(2__)_1_M13771`

aph_2<-aph_2_ia2+aph_2_ia3+aac_6_aph2
as.numeric(aph_2)
aph_2_pos <- rep(0,length(aph_2)) 
aph_2_pos[aph_2 > 0] <- 1
prev_aph_2<-mean(aph_2_pos)


#Calcualting CPM - counts per million
seq_depth<-reads_seqdepth_vetII_samples$TotalFragmentsN
CPM_aph_2<-aph_2/seq_depth*1000000
CPM_aph_2_<-CPM_aph_2[CPM_aph_2 > 0]

#fitting distributions to positive
plotdist(CPM_aph_2_,histo=TRUE,demp=TRUE)

fw <- fitdist(CPM_aph_2_, "weibull") 
fg <- fitdist(CPM_aph_2_, "gamma") 
fln <- fitdist(CPM_aph_2_, "lnorm") 
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fw, fln, fg), legendtext = plot.legend) 
qqcomp(list(fw, fln, fg), legendtext = plot.legend) 
cdfcomp(list(fw, fln, fg), legendtext = plot.legend) 
ppcomp(list(fw, fln, fg), legendtext = plot.legend)

gofstat(list(fln, fg, fw))

fln <- fitdist(CPM_aph_2_, "lnorm") 
summary(fln)


#B simulating true occurence of CPM aph genes
#first 120 months no change, after that 5% increase per year
y <- seq(1,120,1)
d <- (1)^y
s <- seq(121,240,1)
c <- (1+0.05/12)^(s-120)
increase_prop <- c(d,c)

#Interested the prevalence ## We increased the proportion of pigs carrying the gene
#whereas the conc-dist is the same in pigs carrying the gene
prev_month <- prev_aph_2 * increase_prop # probabilities 
true_data_<-mean(rlnorm(n = 1000000,meanlog = fln$estimate[1], sdlog = fln$estimate[2]))
true_data<-prev_month*true_data_

#C - sampling and sequencing results
#results based on X1 n_sample per pool and X2 seq_depth#

vec=seq(from=1,by=1,length.out = 240)
series_n<-11 #number of series. skal være 1001

  #n_sample<-5
  #n_sample<-20
  #n_sample<-50
  n_sample<-100
  

#simulating CPM in sample given 5 samples
set.seed(1234)
meanConc_aph1 <- data.frame(matrix(nrow=240,ncol=series_n)) 
for (k in 1:series_n) {#
  C <-  prev_month 
  conc_gene<-vector("numeric")
  for (i in 1:240) {
    pigs<-vector("numeric")
    for (j in 1:n_sample){ # Change to number of sample
      inf_pig <- (rbinom(n = 1, size = 1, prob=C[i]))
      konc_pig <- rlnorm(n = 1,meanlog = fln$estimate[1], sdlog = fln$estimate[2])
      frag_pig <- inf_pig*konc_pig
      pigs[j] <- frag_pig
      conc_gene_ <- mean(pigs)
    }
    conc_gene[i] <- conc_gene_
  }
  meanConc_aph1[,k] <- conc_gene
}
meanConc_aph1 #<- meanConc/1e6

#generate time series for sample concentration
CPM_sample_aph <- ts(meanConc_aph1[,1:series_n], start = 2000, frequency = 12)

#simulating the CPM measured in metagenomics
set.seed(123)
SS_series_n_aph <- matrix(0,ncol=series_n,nrow=240)
for (k in 1:series_n){
  seq_depth<-rtruncnorm(vec,a=54e6*0.8,b=54e6*1.2,mean=54e6,sd=(54e6-54e6*0.8)/3)
  #seq_depth_5<-rtruncnorm(vec,a=167e6*0.8,b=167e6*1.2,mean=167e6,sd=(167e6-167e6*0.8)/3)
  #seq_depth_20<-rtruncnorm(vec,a=149e6*0.8,b=149e6*1.2,mean=149e6,sd=(149e6-149e6*0.8)/3)
  #seq_depth_50<-rtruncnorm(vec,a=114e6*0.8,b=114e6*1.2,mean=114e6,sd=(114e6-114e6*0.8)/3)
  #seq_depth_100<-rtruncnorm(vec,a=54e6*0.8,b=54e6*1.2,mean=54e6,sd=(54e6-54e6*0.8)/3)
  
  seqdepth <- round(as.numeric(seq_depth), 0)
  A <- meanConc_aph1[,k]/1e6
  B <- seqdepth
  CPM_aph <- numeric(240)
  for (m in 1:240) {
    X_aph <- rbinom(n=1,size=B[m],prob=A[m])
    CPM_aph[m] <- (X_aph)/B[m]*1e6
  }
  SS_series_n_aph[,k] <- CPM_aph
  if(k-1 %% 1){
    print(k)
  }
}
SS_series_n_aph

#Converting to ts

CPMaph_n <- ts(SS_series_n_aph[,10], start = 2000, frequency = 12)



#-----------------------------------------------------------------
# setting parameters for detecting changes greater than a threshold of for example change.magnitude=0.020


#Paramters
percent = 0.1
True_Mean = mean(CPM_aph_2)
duration = 6
first = True_Mean/2
second = True_Mean


plot_aph<-DBEST(data= CPMaph_n, data.type="non-cyclical",
                seasonality=12, algorithm="change detection",
                first.level.shift=first,
                second.level.shift=second, duration=duration, 
                distance.threshold="default", alpha=0.05, change.magnitude=percent*True_Mean,plot="fig1")

plot(plot_aph,figure=1, title="n")
