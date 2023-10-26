#Code used for fitting the concentration, simulating data and plotting 
remove (list = objects() )

#load libraries
source("libraries.R")

reads_seqdepth_vetII_samples <- read_excel("reads_seqdepth_vetII_samples.xlsx")

aph_2_ia2<-reads_seqdepth_vetII_samples$`aph(2__)-Ia_2_AP009486`
aph_2_ia3<-reads_seqdepth_vetII_samples$`aph(2__)-Ia_3_AJ536195`
aac_6_aph2<-reads_seqdepth_vetII_samples$`aac(6_)-aph(2__)_1_M13771`

aph_2<-aph_2_ia2+aph_2_ia3+aac_6_aph2
as.numeric(aph_2)
aph_2_pos <- rep(0,length(aph_2)) 
aph_2_pos[aph_2 > 0] <- 1
prev_aph_2<-mean(aph_2_pos)

#fitting log normal to positive
seq_depth<-reads_seqdepth_vetII_samples$TotalFragmentsN
CPM_aph_2<-aph_2/seq_depth*10^6
CPM_aph_2_<-CPM_aph_2[CPM_aph_2 > 0]
fln <- fitdist(CPM_aph_2_, "lnorm") 
summary(fln)

#increase_prop # will be the same vector
y <- seq(1,120,1)
d <- (1)^y
s <- seq(121,240,1)
c <- (1+0.05/12)^(s-120)
increase_prop <- c(d,c)

#Interested the prevalence ## We increased the proportion of prevalence
#whereas the conc is the same
prev_month <- prev_aph_2 * increase_prop # probabilities 

set.seed(1234)
#series_n <- 1001 #running 1001 series and analyse them in DBEST takes TIME!
series_n <- 2
n_sample<-50

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
meanConc_aph1 

set.seed(1234)
SS_series_n_aph <- matrix(0,ncol=series_n,nrow=240)
for (k in 1:series_n){
  vec=seq(from=1,by=1,length.out = 240)
  #Change to depth based on the money you are using
  depth_ <-rtruncnorm(vec,a=114e6*0.8,b=114e6*1.2,mean=114e6,sd=(114e6-114e6*0.8)/3)
  seqdepth <- round(as.numeric(depth_), 0)
  A <- meanConc_aph1[,k]/1e6
  B <- seqdepth
  CPM_aph <- numeric(240)
  for (m in 1:240) {
    X_aph <- rbinom(n=1,size=B[m],prob=A[m]) #stochasticity
    CPM_aph[m] <- (X_aph)/B[m]*1e6
  }
  SS_series_n_aph[,k] <- CPM_aph
  if(k-1 %% 1){
    print(k)
  }
}
SS_series_n_aph

#Converting to ts

CPM_aph <- list(SS_series_n_aph[,1:series_n])
for (e in 1:length(CPM_aph[[1]])) {
  CPMaph_n <- ts(SS_series_n_aph[,1:series_n], start = 2000, frequency = 12)
}
CPMaph_n

#-----------------------------------------------------------------
# setting parameters for detecting changes greater than a threshold of for example change.magnitude=0.020

aph_n <- list()
for (k in 1:series_n) {
  aph_n[[k]] <-  CPMaph_n[,k]
}

#Paramters
percent = 0.1
True_Mean = mean(CPM_aph_2)
duration = 6
first = True_Mean/2
second = True_Mean

EStart_n <- list() # month  of the break point
Esig_n <- list() # significant of the break point
#EChange <- list()
z <- length(aph_n) # To  chanage
for (g in 1:series_n) {
  xdata <- aph_n[[g]] # To  chanage
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


ESTARTaph_n <- rbind.fill(LISTOFALL)
#removing # will overwrite data in depository
#write.xlsx(ESTARTaph_n, file = "EStart50APH.xlsx",rowNames = FALSE)
