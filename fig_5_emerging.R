#load libraries
source("libraries.R")

#tem5
{blatem5<-read_excel("blatem_time_to_detect_5s.xlsx")
  tem_data<-blatem5
time_to_detect<-blatem5
scen<-rep(c("blatem5"), times=1001)
tem5<-data.frame(scen, time_to_detect)
}

#tem20
{blatem20<-read_excel("blatem_time_to_detect_20s.xlsx")
  tem_data<-blatem20
  time_to_detect<-blatem20
  scen<-rep(c("blatem20"), times=1001)
  tem20<-data.frame(scen, time_to_detect)
}

#tem50
{blatem50<-read_excel("blatem_time_to_detect_50s.xlsx")
  tem_data<-blatem50
  time_to_detect<-blatem50
  scen<-rep(c("blatem50"), times=1001)
  tem50<-data.frame(scen, time_to_detect)
}

#tem100
{blatem100<-read_excel("blatem_time_to_detect_100s.xlsx")
  tem_data<-blatem100
  time_to_detect<-blatem100
  scen<-rep(c("blatem100"), times=1001)
  tem100<-data.frame(scen, time_to_detect)
}

#esbl20
{esbl20<-read_excel("blatem_time_to_detect_ESBL_20s.xlsx")
  tem_data<-esbl20
  time_to_detect<-esbl20
  scen<-rep(c("ESBL20"), times=1001)
  ESBL20<-data.frame(scen, time_to_detect)
}


blatem_time_to_detect<-rbind(tem5, tem20, tem50, tem100, ESBL20)


mynamestheme <- theme(
  panel.background = element_rect(fill = 'white', color = 'grey'),
  panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
  axis.title = element_text(size = (15), colour = "steelblue4"),
  axis.text = element_text(colour = "cornflowerblue", size = (10))
)

b<-ggplot(data=blatem_time_to_detect, aes(x=factor(scen, levels = c("blatem5", "blatem20", "blatem50", "blatem100", "ESBL20")),y=time_to_detect))+
  geom_violin(draw_quantiles = c(0.5))+
  geom_jitter(width = 0.1, height = 0.1, color="brown", size=0.4, alpha=0.9)

c<-b + ylim(0, 60)+labs(x="", y="# time to detect")+
  ggtitle("")    

d <-c + annotate("text", x = 1.5, y = 145, label = "no detection",color="black", size=5)

e <- d + scale_y_continuous(breaks = c(0,12,24,36,48,60))

f<-e+ scale_x_discrete(labels=c("blatem5"="mg_5", "blatem20"="mg_20", "blatem50"="mg_50", "blatem100"="mg_100", "ESBL20"="pheno_20"))

plot_time_to_detect<-f+mynamestheme
plot_time_to_detect


#plot true prevalence farms

n_Farm <- c(1,1,1,2,2,2,2,3,3,3,4,5,5,6,7,8,9,10,12,14,16,18,21,24,27,32,36,
            42,48,55,63,72,83,95,109,125,144,165,190,218,250,287,329,378,434,
            498,572,657,754,865,993,1140,1309,1503,1725,
            1981,2274,2610,2997,3440)

prevFarm<-(n_Farm/4000)*100
prevFarm
month<-1:60

dat_prev<-data.frame(prevFarm, month)

b<-ggplot(data=dat_prev, aes(x=prevFarm,y=month))+
  geom_line()

c<-b + ylim(0, 60)+labs(x="farm prevalence %", y="# time after introduction")+
  ggtitle("")    

e <- c + scale_y_continuous(breaks = c(0,12,24,36,48,60))

plot_farm_prevalence_time_to_detect<-e+mynamestheme
plot_farm_prevalence_time_to_detect





