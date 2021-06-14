
library(ggplot2)
library(plyr)
library(dplyr)
library(nlme)
library(faraway)

###Crozet

load("CrozetFin.Rdata")

# proba starting depred 
g1<-g1low<-g1up<-NA
g1<-resage[41:56,1]
g1low<-resage[41:56,3]
g1up<-resage[41:56,4]
# gamma removal entry probability 
gamma<-resage[58,1]
gammalow<-resage[58,3]
gammaup<-resage[58,4]

# survival 1 juvenile 11 non depred 12 depred
S11<-resage[34,1]
S11low<-resage[34,3]
S11up<-resage[34,4]
S12<-resage[35,1]
S12low<-resage[35,3]
S12up<-resage[35,4]
S11
S12
# survival 2 adults 21 non depred 22 depred
S21<-resage[36,1]
S21low<-resage[36,3]
S21up<-resage[36,4]
S22<-resage[37,1]
S22low<-resage[37,3]
S22up<-resage[37,4]
S21
S22
# q observed from the coast
q<-resage[58,1]
qlow<-resage[58,3]
qup<-resage[58,4]
q
# p observed from boat
p<-resage[59,1]
plow<-resage[59,3]
pup<-resage[59,4]
p


wC1<-1/resage[1:16,2]
wD1<-1/resage[17:32,2]
low_ND1=resage[17:32,3]
up_ND1=resage[17:32,4]
low_NC1=resage[1:16,3]
up_NC1=resage[1:16,4]
NC1=resage[1:16,1]
ND1=resage[17:32,1]
Year=c(1:16)
Year<-seq(2003,2018)
Year2<-seq(2003,2018)
ad=acf(ND1)
ac=acf(NC1)


D1=data.frame(ND1,Year,wD1,low_ND1,up_ND1)
D1=data.frame(ND1,Year2,wD1,low_ND1,up_ND1)
C1=data.frame(NC1,Year,wC1,low_NC1,up_NC1)
C1=data.frame(NC1,Year2,wC1,low_NC1,up_NC1)
model1_D<-gnls(ND1~SSlogis(Year,Asym,xmid,scal), data=D1,
               weights =~wD1, correlation=corAR1(0.75))
model2_D<-gls(ND1~Year, data=D1,weights =~wD1, correlation=corAR1(0.75))
summary(model2_D)
# Generalized least squares fit by REML
# Model: ND ~ Year 
# Data: D1 
# AIC      BIC    logLik
# 119.6246 122.1808 -55.81228
# 
# Correlation Structure: AR(1)
# Formula: ~1 
# Parameter estimate(s):
#   Phi 
# 0.7921113 
# Variance function:
#   Structure: fixed weights
# Formula: ~wD 
# 
# Coefficients:
#   Value Std.Error  t-value p-value
# (Intercept) 39.77671 11.918365 3.337430  0.0049
# Year         3.90215  1.276349 3.057274  0.0085
# 
# Correlation: 
#   (Intr)
# Year -0.791
# 
# Standardized residuals:
#   Min         Q1        Med         Q3        Max 
# -1.3486443 -0.1493441  0.3522264  0.8619780  1.2595474 
# 
# Residual standard error: 15.50982 
# Degrees of freedom: 16 total; 14 residual
model3_D<-gls(ND1~1, data=D1,weights =~wD1, correlation=corAR1(0.75))
AIC(model1_D,model2_D,model3_D)
# df      AIC
# model1_D  5 106.2884
# model2_D  4 119.6246
# model3_D  3 126.6177

model1_C<-gnls(NC1~SSlogis(Year,Asym,xmid,scal), data=C1,
               weights =~wC1, correlation=corAR1(0.75))
model2_C<-gls(NC1~Year, data=C1,weights =~wC1, correlation=corAR1(0.75))
summary(model2_C)
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept) 12.729022 1.1892132 10.703734       0
# Year        -0.855108 0.1441423 -5.932384       
model3_C<-gls(NC1~1, data=C1,weights =~wC1, correlation=corAR1(0.75))
AIC(model1_C,model2_C,model3_C)
# df      AIC
# model1_C  5 75.21294
# model2_C  4 80.16151
# model3_C  3 84.32412
plot(D1$Year,C1$NC1, type='p',col='darkblue')
curve(predict(model1_C,newdata=data.frame(Year=x)),add=TRUE, col='orange')
plot(D1$Year,D1$ND1, type='p',col='darkblue')
curve(predict(model1_D,newdata=data.frame(Year=x)),add=TRUE, col='orange')



ggplot(C1, aes(x=Year, y=NC1)) +
  geom_errorbar(aes( ymin = low_NC1, ymax = up_NC1), width=0.4)+
 geom_point(aes(x=Year, y=NC1))+
  scale_x_continuous(breaks = round(seq(2003, 2018, by = 3),1)) +
  ggtitle("(A)")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(D1, aes(x=Year, y=ND1)) +
  geom_errorbar(aes( ymin = low_ND1, ymax = up_ND1), width=0.4)+
 geom_point(aes(x=Year, y=ND1))+
  scale_x_continuous(breaks = round(seq(2003, 2018, by = 3),1)) +
  ggtitle("(B)")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

## Type C
tt<-predict(model1_D)
dd<-data.frame(Year2=Year2, pred=tt)
D1<-merge(D1,dd)

CC1<-merge(C1,D1)
ggplot(CC1, aes(x=Year, y=NC1)) +
  geom_errorbar(aes( ymin = low_ND1, ymax = up_ND1), width=0.4)+
  geom_point(aes(x=Year, y=ND1), shape=21, fill = "lightgrey",
             color = "black", size = 2)+
  geom_errorbar(aes( ymin = low_NC1, ymax = up_NC1), width=0.4)+
  geom_point(aes(x=Year, y=NC1), shape=21, fill = "black",
             color = "black", size = 2)+
  geom_line(aes(x=Year2, y=pred),color="lightgrey")+
  scale_x_continuous(breaks = round(seq(2003, 2018, by = 3),1)) +
  ylab("Number of individuals")+
  scale_y_continuous(breaks = round(seq(0, 100, by = 10),1), limits=c(0,100)) +
  ggtitle("(A)")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

 

####Type D


load("DFin.Rdata")

# proba starting depred 
g2<-g2low<-g2up<-NA
g2<-resageD[22:37,1]
g2low<-resageD[22:37,3]
g2up<-resageD[22:37,4]
# survival 1 juveniles
S1<-resageD[18,1]
S1low<-resageD[18,3]
S1up<-resageD[18,4]
# survival 2 adults
S2<-resageD[19,1]
S2low<-resageD[19,3]
S2up<-resageD[19,4]
# p observed from boat
p<-resageD[39,1]
plow<-resageD[39,3]
pup<-resageD[39,4]



wD2<-1/resageD[1:16,2]

low_ND2=resageD[1:16,3]
up_ND2=resageD[1:16,4]
ND2=resageD[1:16,1]
Year=c(1:16)
Year<-seq(2003,2018)
Year2<-seq(2009,2018)
ad=acf(ND2)
ac=acf(NC2)

DD1=data.frame(ND2,Year,wD2)
DD1=data.frame(ND2=ND2[7:16],Year2=Year2,wd2=wD2[7:16], low_ND2=low_ND2[7:16],up_ND2=up_ND2[7:16])
model1_D<-gnls(ND2~SSlogis(Year2,Asym,xmid,scal), data=DD1,
               weights =~wd2, correlation=corAR1(0.75))
model2_D<-gls(ND2~Year2, data=DD1,weights =~wd2, correlation=corAR1(0.75))
summary(model2_D)
model3_D<-gls(ND2~1, data=DD1,weights =~wd2, correlation=corAR1(0.75))
AIC(model1_D,model2_D,model3_D)
# df      AIC
# model1_D  5 59.56215
# model2_D  4 50.81241
# model3_D  3 58.99627

tt<-predict(model2_D)
dd<-data.frame(Year2=Year2, pred=tt)
DD1<-merge(DD1,dd)

ggplot(DD1, aes(x=Year2, y=ND2)) +
  geom_errorbar(aes( ymin = low_ND2, ymax = up_ND2), width=0.4)+
  geom_point(aes(x=Year2, y=ND2), color="lightgrey")+
  geom_line(aes(x=Year2, y=pred),color="lightgrey")+
  ylab("Number of individuals")+
  xlab("Year")+
  scale_x_continuous(breaks = round(seq(2003, 2018, by = 3),1),limits=c(2003,2018)) +
  scale_y_continuous(breaks = round(seq(0, 100, by = 10),1), limits=c(0,100)) +
  ggtitle("(B)")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


ggplot(CC1, aes(x=Year, y=NC2)) +
  geom_errorbar(aes( ymin = ND2-sd_ND2, ymax = ND2+sd_ND2), width=0.4)+
  geom_point(aes(x=Year, y=ND2), shape=21, fill = "lightgrey",
             color = "black", size = 2)+
  geom_errorbar(aes( ymin = NC2-sd_NC2, ymax = NC2+sd_NC2), width=0.4)+
  geom_point(aes(x=Year, y=NC2), shape=21, fill = "black",
             color = "black", size = 2)+
  scale_x_continuous(breaks = round(seq(2003, 2018, by = 3),1)) +
  scale_y_continuous(breaks = round(seq(0, 100, by = 10),1), limits=c(0,100)) +
  ggtitle("(B)")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


CD=data.frame(NC2,Year,wC2, ND2, wD2)
CD2<-CD
CD2$NC2[CD2$Year<2009]<-NA
CD2$ND2[CD2$Year<2009]<-NA
CD2$wC2[CD2$Year<2009]<-NA
CD2$wD2[CD2$Year<2009]<-NA


ggplot(CD2, aes(x=Year, y=NC2)) +
  geom_errorbar(aes( ymin = ND2-wD2, ymax = ND2+wD2), width=0.4)+
  geom_point(aes(x=Year, y=ND2), shape=21, fill = "lightgrey",
             color = "black", size = 2)+
  geom_errorbar(aes( ymin = NC2-wC2, ymax = NC2+wC2), width=0.4)+
  geom_point(aes(x=Year, y=NC2), shape=21, fill = "black",
             color = "black", size = 2)+
  scale_x_continuous(breaks = round(seq(2003, 2018, by = 3),1)) +
  scale_y_continuous(breaks = round(seq(0, 100, by = 10),1), limits=c(0,100)) +
  ggtitle("(B)")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


###All together pop###

CCD<-data.frame(Year=CD$Year,NC1=C1$NC1, ND1=D1$ND1, wD1=D1$wD1, wC1=C1$wC1, 
                NC2=CD2$NC2, ND2=CD2$ND2, wD2=CD2$wD2, wC2=CD2$wC2)

ggplot(CCD, aes(x=Year, y=NC1)) +
  geom_errorbar(aes( ymin = ND1-wD1, ymax = ND1+wD1), 
                width=0.4, color="black")+
  geom_point(aes(x=Year, y=ND1), shape=25, fill = "black", color="black",
              size = 2)+
  geom_errorbar(aes( ymin = NC1-wC1, ymax = NC1+wC1),
                width=0.4, color="#999999", linetype=1)+
  geom_point(aes(x=Year, y=NC1), shape=21, fill = "#999999", color="#999999",
              size = 2)+
  geom_errorbar(aes( ymin = ND2-wD2, ymax = ND2+wD2),width=0.4, linetype=3)+
  geom_point(aes(x=Year, y=ND2), shape=25, fill = "black",
             color = "black", size = 2)+
  geom_errorbar(aes( ymin = NC2-wC2, ymax = NC2+wC2),
                width=0.4,color="#999999", linetype=3)+
  geom_point(aes(x=Year, y=NC2), shape=21, fill = "#999999",
             color = "#999999", size = 2)+
  scale_x_continuous(breaks = round(seq(2003, 2018, by = 3),1)) +
  scale_y_continuous(breaks = round(seq(-10, 100, by = 10),1), limits=c(-5,100)) +
  ggtitle("B")+
  ylab("")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

###All together proba depred###

gg<-data.frame(Year=seq(2003,2018),g1=g1, g1low=g1low,g1up=g1up, 
                g2=g2, g2low=g2low, g2up=g2up)

gg$g2[gg$Year<2009]<-NA
gg$g2low[gg$Year<2009]<-NA
gg$g2up[gg$Year<2009]<-NA

ggplot(gg, aes(x=Year, y=g1)) +
  geom_errorbar(aes( ymin = g1low, ymax = g1up), 
                width=0.4, color="black")+
  geom_point(aes(x=Year, y=g1), shape=21, fill = "black", color="black",
             size = 2)+
  geom_errorbar(aes( ymin = g2low, ymax = g2up),
                width=0.4, linetype=1, color="#999999")+
  geom_point(aes(x=Year, y=g2), shape=21, fill = "#999999",
             color = "#999999", size = 2)+
  scale_x_continuous(breaks = round(seq(2003, 2018, by = 3),1)) +
  scale_y_continuous(breaks = round(seq(0, 1.1, by = 0.2),1), limits=c(0,1.1)) +
  ggtitle("")+
  ylab("Probability of starting depredation")+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(gg, aes(x=Year, y=g1)) +
  geom_errorbar(aes( ymin = g2low, ymax = g2up),
                width=0.4, linetype=1, color="black")+
  geom_point(aes(x=Year, y=g2), shape=21, fill = "black",
             color = "#999999", size = 2)+