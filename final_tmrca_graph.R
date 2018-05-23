markerfile <- read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/ooa5k.marker.txt")
mfile <- cbind(markerfile$MarkerID, (markerfile$AlleleCount1/(markerfile$AlleleCount0+markerfile$AlleleCount1)))   #mfile store the frequency of all markers from the simulation, e.g. N3
#mfrare1 <- mfile[mfile[,2] <= 0.05, ]             #all markers with derived allele <= 5% (used when we want to consider only rare variants under a certain threshold)
#mfrare <- mfrare1[,1]              #all markers with derived allele <= 5%
mfrare <- mfile[,1] #used when all variants are considered (have to neutralize it when graphs/stats for rare variants)

a1<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_1.txt")      #true vs inferred TMRCA calculated thanks to patrick TMRCA algorithm (cf link in report)
a2<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_2.txt")
a3<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_3.txt")
a4<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_4.txt")
a5<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_5.txt")
a6<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_6.txt")
a7<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_7.txt")
a8<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_8.txt")
a9<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_9.txt")
a10<-read.table(header=TRUE, "/media/francoismallord/USB DISK/chr1_for_real_age/n3/n3_true_ccf_10.txt")
asum=rbind(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)    #all true vs inferred TMRCA from N3 gathered together
rm(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)     #to have more space to do the graph, better to suppress data when they become useless

e1a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)        #we only consider the variants previously selected in mfrare
e1a<-subset(e1a, as.character(e1a$Clock)=='M')
e1<-e1a$TMRCA         #real age
e1<-as.data.frame(e1)
f1<-(e1a$Shape/e1a$Rate)*2*ne #estimated age
f1<-as.data.frame(f1)

e2a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)
e2a<-subset(e2a, as.character(e2a$Clock)=='R')
e2<-e2a$TMRCA         #real age
f2<-(e2a$Shape/e2a$Rate)*2*ne #estimated age
e2<-as.data.frame(e2)
f2<-as.data.frame(f2)

e3a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)
e3a<-subset(e3a, as.character(e3a$Clock)=='C')
e3<-e3a$TMRCA         #real age
f3<-(e3a$Shape/e3a$Rate)*2*ne #estimated age
e3<-as.data.frame(e3)
f3<-as.data.frame(f3)

e4a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)
e4a<-subset(e4a, as.character(e4a$Clock)=='M')
e4a<-subset(e4a, as.numeric(e4a$Shared)==1) #conc pairs
e4<-e4a$TMRCA         #real age
f4<-(e4a$Shape/e4a$Rate)*2*ne #estimated age
e4<-as.data.frame(e4)
f4<-as.data.frame(f4)

e5a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)
e5a<-subset(e5a, as.character(e5a$Clock)=='R')
e5a<-subset(e5a, as.numeric(e5a$Shared)==1) #conc pairs
e5<-e5a$TMRCA         #real age
f5<-(e5a$Shape/e5a$Rate)*2*ne #estimated age
e5<-as.data.frame(e5)
f5<-as.data.frame(f5)

e6a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)
e6a<-subset(e6a, as.character(e6a$Clock)=='C')
e6a<-subset(e6a, as.numeric(e6a$Shared)==1) #conc pairs
e6<-e6a$TMRCA         #real age
f6<-(e6a$Shape/e6a$Rate)*2*ne #estimated age
e6<-as.data.frame(e6)
f6<-as.data.frame(f6)

e7a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)
e7a<-subset(e7a, as.character(e7a$Clock)=='M')
e7a<-subset(e7a, as.numeric(e7a$Shared)==0) #disc pairs
e7<-e7a$TMRCA         #real age
f7<-(e7a$Shape/e7a$Rate)*2*ne #estimated age
e7<-as.data.frame(e7)
f7<-as.data.frame(f7)

e8a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)
e8a<-subset(e8a, as.character(e8a$Clock)=='R')
e8a<-subset(e8a, as.numeric(e8a$Shared)==0) #disc pairs
e8<-e8a$TMRCA         #real age
f8<-(e8a$Shape/e8a$Rate)*2*ne #estimated age
e8<-as.data.frame(e8)
f8<-as.data.frame(f8)

e9a<-subset(asum, as.numeric(asum$MarkerID) %in% mfrare)
e9a<-subset(e9a, as.character(e9a$Clock)=='C')
e9a<-subset(e9a, as.numeric(e9a$Shared)==0) #disc pairs
e9<-e9a$TMRCA         #real age
f9<-(e9a$Shape/e9a$Rate)*2*ne #estimated age
e9<-as.data.frame(e9)
f9<-as.data.frame(f9)

names(e1)<-"V1"
names(e2)<-"V1"
names(e3)<-"V1"
names(e4)<-"V1"
names(e5)<-"V1"
names(e6)<-"V1"
names(e7)<-"V1"
names(e8)<-"V1"
names(e9)<-"V1"
names(f1)<-"V1"
names(f2)<-"V1"
names(f3)<-"V1"
names(f4)<-"V1"
names(f5)<-"V1"
names(f6)<-"V1"
names(f7)<-"V1"
names(f8)<-"V1"
names(f9)<-"V1"

e=rbind(e1, e2, e3, e4, e5, e6, e7, e8, e9)
f=rbind(f1, f2, f3, f4, f5, f6, f7, f8, f9)
Xe=c(nrow(e1), nrow(e2), nrow(e3), nrow(e4), nrow(e5), nrow(e6), nrow(e7), nrow(e8), nrow(e9))        #to find the total number of pairs by type and clock model
Xf=c(nrow(f1), nrow(f2), nrow(f3), nrow(f4), nrow(f5), nrow(f6), nrow(f7), nrow(f8), nrow(f9))        #to find the total number of pairs by type and clock model


RMSL10E=function(a,b) {
  n=length(a)
  a=log10(a)
  b=log10(b)
  sqrt(sum((a-b)**2/n))
}

#number=vector(mode = "numeric", length=length(e[[1]]))
exp=vector(mode = "character", length=length(e[[1]]))
clock=vector(mode = "character", length=length(e[[1]]))
exp1=vector(mode = "character", length=9)
clock1=vector(mode = "character", length=9)
r2=vector(mode = "numeric", length=9)
rh=vector(mode = "numeric", length=9)
rmsle=vector(mode = "numeric", length=9)

i=1
k=1
frontiers=c(0, length(e1[[1]]), length(e2[[1]]), length(e3[[1]]), length(e4[[1]]), length(e5[[1]]), length(e6[[1]]), length(e7[[1]]), length(e8[[1]]), length(e9[[1]]))
ranges=vector(mode = "numeric", length=10)
for (m in 1:10) {
  ranges[m]=sum(frontiers[1:m])
}
type_pairs=c("All pairs", "All pairs", "All pairs", "Concordant pairs", "Concordant pairs", "Concordant pairs", "Discordant pairs", "Discordant pairs", "Discordant pairs")
type_clock=c("Mutation clock", "Recombination clock", "Combined clock", "Mutation clock", "Recombination clock", "Combined clock", "Mutation clock", "Recombination clock", "Combined clock")

rm(e1, e2, e3, e4, e5, e6, e7, e8, e9)
rm(f1, f2, f3, f4, f5, f6, f7, f8, f9)
rm(e1a, e2a, e3a, e4a, e5a, e6a, e7a, e8a, e9a)  #suppress data that have become useless

while (k<=length(e[[1]])) {
  if (ranges[i]<k & k<=ranges[i+1]) {
    exp[k]=type_pairs[i]
    clock[k]=type_clock[i]
    if (k==ranges[i]+1) {
      exp1[i]=type_pairs[i]
      clock1[i]=type_clock[i]
      b=ranges[i]+1
      c=ranges[i+1]
      b0=cor(log(e$V1[b:c]), log(f$V1[b:c]))**2
      r2[i]=as.numeric(format(b0, digits=3))
      b1=cor.test(~e$V1[b:c]+f$V1[b:c], cbind(e,f), method="spearman", continuity=FALSE, conf.level=0.95)
      b1=b1$estimate     
      rh[i]=as.numeric(format(b1, digits=3))
      b2=RMSL10E(e$V1[b:c], f$V1[b:c])
      rmsle[i]=as.numeric(format(b2, digits=3))
    } else if (k==ranges[i+1]) {
      i=i+1
    }
    
  }
  k=k+1
}

az=cbind(e, f, exp, clock)   # the dataset used for plotting TMRCA graph


#install.packages("ggplot2")
#install.packages("cowplot")
#install.packages("Metrics")
#install.packages("ggthemes")
library(ggthemes)
library(Metrics)
library(ggplot2)
library(cowplot)
library(broom)
library(dplyr)

exp=exp1
clock=clock1
dataz=cbind(exp, clock, r2, rh, rmsle)
datazz=as.data.frame(dataz)           #the dataset of the metrics of the TMRCA panels

rm(markerfile, mfile, b, b0, b1, b2, c, clock, clock1, exp, exp1, frontiers, i, k, m, ne, r2, ranges, rh, rmsle, type_clock, type_pairs, Xe, Xf)
rm(asum, dataz, e, f, mfrare)     #we suppress all that isn't usefull to plot the graph

g=ggplot(az, aes(x=az[[1]], y=az[[2]])) +
  facet_grid(exp~clock) +
  geom_raster(aes(az[[1]], az[[2]], fill = (..count..)/tapply(..count..,..PANEL..,max)[..PANEL..]), position = "identity", stat = "bin2d", binwidth = c(0.1, 0.1)) +
  geom_tile(aes(az[[1]], az[[2]], fill = (..count..)/tapply(..count..,..PANEL..,max)[..PANEL..]), position = "identity", stat = "bin2d", binwidth =   c(0.1, 0.1)) +
  geom_abline(intercept = c(0,0), slope = 1, alpha = 1/2) +
  geom_text(data=datazz, inherit.aes=FALSE, size=3, aes(x=3000, y=3, label=paste("\n rho =", rh, "\n RMSLE =", rmsle)))+
  coord_cartesian(xlim = c(0.5, 99000), ylim = c(0.5, 99000), expand = F) +
  scale_x_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
  scale_y_log10(breaks = 10^(0:6), labels = trimws(format(10^(0:6), big.mark = ',', scientific = F))) +
  scale_fill_gradientn(colours = c("ivory","sandybrown", "sienna2", "sienna3", "sienna4", "chocolate4"), na.value = "white", limits = c(-0.05,1.05), breaks = c(0,0.25,0.5,0.75,1)) +
  theme_few() +
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 9),
        panel.border = element_rect(fill = NA, colour = "grey20", size = 2/3),
        panel.background = element_rect(fill = NA, colour = "grey20", size = 2/3),
        strip.text.x = element_text(face = "bold", size = 11),
        strip.text.y = element_text(size = 13),
        legend.justification=c(0,1), legend.position=c(1 - (3*0.985)/3, (1*0.98)/1),
        legend.background = element_rect(fill = "grey90", colour = "grey70", size = 1/4),
        axis.ticks = element_line(colour = "grey20"),
        axis.text.x = element_text(hjust = 0),
        axis.title = element_text(size = 11),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.25, "cm"),
        legend.margin = margin(-1,1,0.5,0.5, "mm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank()) +
  ylab("Estimated allele age (generations)") + xlab("True allele age (generations)")

