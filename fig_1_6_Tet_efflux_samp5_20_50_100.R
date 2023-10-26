#Code used for fitting the concentration, simulating data and plotting 
remove (list = objects() )

#load libraries
source("libraries.R")
#A: estimating presens of efflux genes"

#reading true data
reads_seqdepth_vetII_samples <- read_excel("reads_seqdepth_vetII_samples.xlsx")

#selecting the genes for tet efflux
tet_efflux<-reads_seqdepth_vetII_samples[c(359,362:363,366:383,435:438)]
seq_depth<-reads_seqdepth_vetII_samples[c(2)]
tet_efflux_sum<-c(rowSums(tet_efflux, na.rm = TRUE))

#Calcualting CPM - counts per million
CPM_tet_efflux<-((tet_efflux_sum/seq_depth)*1000000)
#fitting log normal to positive values of CPM of tet efflux genes
CPM_tet_efflux_<-as.numeric(unlist(CPM_tet_efflux))
fln <- fitdist(CPM_tet_efflux_, "lnorm") 
summary(fln)

#B simulating true occurence of CPM tet efflux genes
#first 120 months no change, after that 5% increase per year
y <- seq(1,120,1)
d <- (1)^y
s <- seq(121,240,1)
c <- (1+0.05/12)^(s-120)
increase_conc <- c(d,c)
#The prevalence is always 100%
prev_month <- 1 # probabilities 

#C - sampling and sequencing results
#results based on X1 n_sample per pool and X2 seq_depth#

vec=seq(from=1,by=1,length.out = 240)
series_n<-11 #number of series. skal være 1001

{#5 samples
  n_sample<-5
  #seq_depth<-rtruncnorm(vec,a=167e6*0.8,b=167e6*1.2,mean=167e6,sd=(167e6-167e6*0.8)/3)
}
# {#20 samples
#   n_sample<-20
#   seq_depth<-rtruncnorm(vec,a=149e6*0.8,b=149e6*1.2,mean=149e6,sd=(149e6-149e6*0.8)/3)
# }
# {#50 samples
#   n_sample<-50
#   seq_depth<-rtruncnorm(vec,a=114e6*0.8,b=114e6*1.2,mean=114e6,sd=(114e6-114e6*0.8)/3)
# }
# {#100 samples
#   n_sample<-100
#   seq_depth<-rtruncnorm(vec,a=54e6*0.8,b=54e6*1.2,mean=54e6,sd=(54e6-54e6*0.8)/3)
# }


#simulating CPM in sample given n samples
set.seed(1234)
meanConc_Tetra1 <- data.frame(matrix(nrow=240,ncol=series_n)) 
for (k in 1:series_n) {#
  conc_gene<-vector("numeric")
  for (i in 1:240) {
    pigs<-vector("numeric")
    for (j in 1:n_sample){
      inf_pig <- 1 #pigs are always caring a gene 
      conc_pig <- rlnorm(n = 1,meanlog = fln$estimate[1], sdlog = fln$estimate[2])
      frag_pig <- inf_pig*conc_pig
      pigs[j] <- frag_pig
      conc_gene_ <- mean(pigs)
    }
    conc_gene[i] <- conc_gene_
  }
  meanConc_Tetra1[,k] <- conc_gene
}

meanConc_Tetra1<-meanConc_Tetra1*increase_conc

#generate time series for sample concentration
CPM_sample_Tetra <- ts(meanConc_Tetra1[,1:series_n], start = 2000, frequency = 12)

#simulating the CPM measured in metagenomics
set.seed(1234)
SS_series_n_Tetra <- matrix(0,ncol=series_n,nrow=240)
for (k in 1:series_n){
    seq_depth<-rtruncnorm(vec,a=167e6*0.8,b=167e6*1.2,mean=167e6,sd=(167e6-167e6*0.8)/3)
    #seq_depth_5<-rtruncnorm(vec,a=167e6*0.8,b=167e6*1.2,mean=167e6,sd=(167e6-167e6*0.8)/3)
    #seq_depth_20<-rtruncnorm(vec,a=149e6*0.8,b=149e6*1.2,mean=149e6,sd=(149e6-149e6*0.8)/3)
    #seq_depth_50<-rtruncnorm(vec,a=114e6*0.8,b=114e6*1.2,mean=114e6,sd=(114e6-114e6*0.8)/3)
    #seq_depth_100<-rtruncnorm(vec,a=54e6*0.8,b=54e6*1.2,mean=54e6,sd=(54e6-54e6*0.8)/3)
  seqdepth <- round(as.numeric(seq_depth), 0)
  A <- meanConc_Tetra1[,k]/1e6
  B <- seqdepth
  CPM_Tetra <- numeric(240)
  for (m in 1:240) {
    X_Tetra <- rbinom(n=1,size=B[m],prob=A[m])
    CPM_Tetra[m] <- (X_Tetra)/B[m]*1e6
  }
  SS_series_n_Tetra[,k] <- CPM_Tetra
  if(k-1 %% 1){
    print(k)
  }
}
SS_series_n_Tetra

#Converting to ts
CPM_Tetra <- list(SS_series_n_Tetra[,1:series_n])
for (e in 1:length(CPM_Tetra[[1]])) {
  CPMTetra <- ts(SS_series_n_Tetra[,1:series_n], start = 2000, frequency = 12)
}

#data management before plotting
a<-CPMTetra[,11]#selecting time series 11 to plot
b<-CPM_sample_Tetra[,11]#selecting time series 11 to plot
data_plot<-data.frame(a, b)
data_ts<-ts(data=data_plot, start=c(2000, 1), end=c(2019, 12), frequency=12)

CPM_obs<-as.vector(data_plot$a)
CPM_sample<-as.vector(data_plot$b)
dates_1 <- as.yearmon(time(data_ts))

#generating true CPM over time 
true_data<-mean(CPM_tet_efflux_)
true_data_<-true_data*increase_conc

#creating a new data object with dates
data_dates<-data.frame(dates_1, CPM_obs, CPM_sample, true_data_)

######plotting fig 1 in paper - Tetracykline ##############################
a_plot<-ggplot(data_dates)+
  geom_line(aes(x=dates_1, y=CPM_obs),col="green",linetype="solid", linewidth=1)+
  geom_line(aes(x=dates_1+0.1, y=CPM_sample),col="brown",linetype="solid", linewidth=1)+
  geom_line(aes(x=dates_1, y = true_data_), linewidth=1, color = "blue")

mynamestheme <- theme(
  panel.background = element_rect(fill = 'white', color = 'grey'),
  panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
  axis.title = element_text(size = (15), colour = "steelblue4"),
  axis.text = element_text(colour = "cornflowerblue", size = (10))
)

b_plot<-a_plot+mynamestheme

b_plot + ylim(15, 65)+labs(x="date", y="CPM")

####plotting forecast - fig 6#####################

CPMTetra[,11] %>%
  ets(model="AAN", alpha = NULL, phi = NULL) %>% 
  #("A" additive error and trend type N" no season)
  forecast(h=60, level=c(90)) %>%
  ggplot2::autoplot(ts.colour='green',ts.size=0.5)+ xlab("date") + ylab("CPM") +
  theme_bw()+ ylim(20,75)
 

