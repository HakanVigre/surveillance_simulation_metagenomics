#Code for simulating 1001 time series and analysing using DBEST package 

#fitting CPM distribution####
remove (list = objects() )
#load libraries
source("libraries.R")

reads_seqdepth_vetII_samples <- read_excel("reads_seqdepth_vetII_samples.xlsx")

tet_efflux<-reads_seqdepth_vetII_samples[c(359,362:363,366:383,435:438)]
seq_depth<-reads_seqdepth_vetII_samples[c(2)]

tet_efflux_sum<-c(rowSums(tet_efflux, na.rm = TRUE))
dim(tet_efflux_sum)
dim(seq_depth)

a<-data.frame(tet_efflux,seq_depth)
CPM_tet_efflux<-((tet_efflux_sum/seq_depth)*1000000)
dim(CPM_tet_efflux)

CPM_tet_efflux_<-as.numeric(unlist(CPM_tet_efflux))

#fitting log normal to positive
fln <- fitdist(CPM_tet_efflux_, "lnorm") 
summary(fln)

y <- seq(1,120,1)
d <- (1)^y
s <- seq(121,240,1)
c <- (1+0.05/12)^(s-120)
increase_conc <- c(d,c)
#The prevalence is always 100%=1
prev_month <- 1 # probabilities 

set.seed(1234)
#series_n <- 1001#running 1001 series and analyse them in DBEST takes TIME!
series_n <- 2
n_sample<-100
meanConc_Tetra1 <- data.frame(matrix(nrow=240,ncol=series_n)) 
for (k in 1:series_n) {#
  conc_gene<-vector("numeric")
  for (i in 1:240) {
    pigs<-vector("numeric")
    for (j in 1:n_sample){
      inf_pig <- 1 #pigs are always caring agene
      konc_pig <- rlnorm(n = 1,meanlog = fln$estimate[1], sdlog = fln$estimate[2])
      frag_pig <- inf_pig*konc_pig
      pigs[j] <- frag_pig
      conc_gene_ <- mean(pigs)
    }
    conc_gene[i] <- conc_gene_
  }
  meanConc_Tetra1[,k] <- conc_gene
}

meanConc_Tetra1<-meanConc_Tetra1*increase_conc

#simulating the CPM measured in metagenomics
set.seed(1234)
SS_series_n_Tetra <- matrix(0,ncol=series_n,nrow=240)
for (k in 1:series_n){
  vec=seq(from=1,by=1,length.out = 240)
  depth_ <-rtruncnorm(vec,a=54e6*0.8,b=54e6*1.2,mean=54e6,sd=(54e6-54e6*0.8)/3)
  seqdepth <- round(as.numeric(depth_), 0)
  A <- meanConc_Tetra1[,k]/1e6
  B <- seqdepth
  CPM_Tetra <- numeric(240)
  for (m in 1:240) {
    X_Tetra <- rbinom(n=1,size=B[m],prob=A[m]) #stochasticity
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
  CPMTetra_n <- ts(SS_series_n_Tetra[,1:series_n], start = 2000, frequency = 12)
}
#CPMTetra_n

#-----------------------------------------------------------------
# setting parameters for detecting changes greater than a threshold of for example change.magnitude=0.020

Tetra_n <- list()
for (k in 1:series_n) {
  Tetra_n[[k]] <-  CPMTetra_n[,k]
}


#Paramters
percent = 0.1
True_Mean = mean(CPM_tet_efflux_)
duration = 6
first = True_Mean/2
second = True_Mean

EStart_n <- list() # month  of the break point
Esig_n <- list() # significant of the break point
#EChange <- list()
z <- length(Tetra_n) # To  chanage
for (g in 1:series_n) {
  xdata <- Tetra_n[[g]] # To  chanage
  trytetst <- try(DEX <- DBEST(data= xdata, data.type="non-cyclical",
                               seasonality=12, algorithm="change detection",
                               first.level.shift=first,
                               second.level.shift=second, duration=duration, 
                               distance.threshold="default", alpha=0.05, change.magnitude=percent*True_Mean))#, plot="fig1"
  EStart_n[[g]] <- DEX$Start
  if (class(trytetst)=="try-error") {
    EStart_n[[g]]<-NULL
  }
  Esig_n[[g]] <- DEX$Significance
  if (class(trytetst)=="try-error") {
    Esig_n[[g]]<-NULL
  }
  print(g)
}

f_n<-Map('*',EStart_n,Esig_n) # change name

# For Break points
LISTOFALL <- list()

for (k in 1:length(f_n)) {
  if (is.null(f_n[[k]])) {
    f_n[[k]]<-NA
  }
  a <- as.data.frame(matrix(f_n[[k]],nrow=1))
  LISTOFALL<-append(LISTOFALL,list(a))
}

ESTARTTetra_n <- rbind.fill(LISTOFALL)
write.xlsx(ESTARTTetra_n, file = "EStart100_Tetra.xlsx",rowNames = FALSE)

