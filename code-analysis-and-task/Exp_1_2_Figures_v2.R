# K. Garner - 11/09/17 Code plots data points for each cue (validity x spatial-certainty x spatial-value ), 
# model predictions = dashed lines, and BFs for winning model relative to next 5.
###### EXPERIMENT 1 PLOT
#------------------------------------------------------------------------------------------------------------------
rm(list = ls())
library(wesanderson)
library(cowplot)
library(ggplot2)
library(plyr)
library(dplyr)


# manual instructions
# 1. SETWD TO SOURCE FILE LOCATION


#setwd("~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/EXP_3_CUE_PROB_OUT/exp3")
load("EXP1/ANALYSIS/EXP1_ANALYSIS_MIXDMDLS_BIAS.R")
# PLOT

eb_length = .1
c0 = -(.5*log2(.5))-(.5*log2(.5)) # defined according to Prinzmetal et al, 2015, JEP:HPP
xs = c( - .98*log2(.98)  - .02*log2(.02) ,
        - .92*log2(.92)  - .08*log2(.08) ,
        -  .9*log2(.9)   -  .1*log2(.1) ,
        -  .8*log2(.8)   -  .2*log2(.2) ,
        -  .6*log2(.6)   -  .4*log2(.4)  )
xs = c0 - xs

####### PLOT USING GGPLOT2
figxs = rep(xs, each = 4)
sum.dat.all.plot$xs = figxs
# get 1 std error of the mean for plotting the error bars
# decided against the 95% within subject confidence intervals as would have to do one
# for each comparison
sum.dat.win.se = ddply( sum.dat.err, .(value, fixprob, valid), summarise,
                        se = sd(mu)/sqrt(length(mu)))
sum.dat.all.plot = inner_join(sum.dat.all.plot, sum.dat.win.se, by=c("value", "fixprob", "valid"))

l_wd = 1
f_size=8
alph = 0.5
this.legend.title = "Value"
rt.plot <- ggplot(sum.dat.all.plot, aes(x=xs, y=mean, color=as.factor(r1))) +
  ylim(0.5,0.9) +
  geom_point(size=1.5) + geom_errorbar(aes(x=xs, ymin=mean-se, ymax=mean+se)) + 
  geom_line(aes(x=xs, y=predict, color=as.factor(r1)), linetype=1, lwd=l_wd, alpha=alph) +
  ylab("RT") + xlab("Spatial Certainty (bits)") +
  facet_wrap(~valid) + 
  scale_color_manual(this.legend.title, values=c(wes_palette("Royal1")[1], wes_palette("Royal1")[2]), 
                     guide=guide_legend(direction = "horizontal")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=f_size),
        axis.text.x = element_text(size=f_size),
        axis.text.y = element_text(size=f_size),
        legend.position=c(0.05, 0.8)) 
ggsave("e1_rt_data.bmp", plot=rt.plot, width=10, height=5, units="cm") 
###### NOW ACCURACY
acc.dat.win.err = ddply( acc.dat.err, .(value, fixprob, valid), summarise,
                         se = sd(acc)/sqrt(length(acc)))
acc.dat.plot = inner_join(acc.dat.plot, acc.dat.win.err, by=c("value", "fixprob", "valid"))
acc.dat.plot$xs = figxs

acc.plot <- ggplot(acc.dat.plot, aes(x=xs, y=acc, color=as.factor(r1))) +
  ylim(0.5,1) +
  geom_point(size=1.5) + geom_errorbar(aes(x=xs, ymin=acc-se, ymax=acc+se)) + 
  geom_line(aes(x=xs, y=predict, color=as.factor(r1)), linetype=1, lwd=l_wd, alpha=alph) +
  ylab("Acc") + xlab("Spatial Certainty (bits)") +
  facet_wrap(~valid) + 
  scale_color_manual(this.legend.title, values=c(wes_palette("Royal1")[1], wes_palette("Royal1")[2])) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=f_size),
        axis.text.x = element_text(size=f_size),
        axis.text.y = element_text(size=f_size),
        legend.position="none")
ggsave("e1_acc_data.bmp", plot=acc.plot, width=10, height=5, units="cm") 

# BAYES RT BARPLOT
# ORDER FOR BARPLOT
dot_size = 2
this.legend.title = "Model Group"
BF = c(rt.bf.for.plot$BF) # adding a zero to place between the two groups of models
names = c(as.character(rt.bf.for.plot$names))
pnames=c("iii", "i", "iv", "ii", "v")
lower = c(rt.bf.for.plot$lower) # adding a zero to place between the two groups of models
upper = c(rt.bf.for.plot$upper) # adding a zero to place between the two groups of models
bf.rt.for.plot = data.frame(BF, names, lower, upper, pnames)
bar.order = c(2, 4, 1, 3, 5) #c(2, 4, 6, 1, 3, 5)
bf.rt.for.plot = bf.rt.for.plot[bar.order, ]
bf.rt.for.plot$grp = c("add", "add", "exp", "exp", "exp")
bf.rt.for.plot$x = seq(1, 5, 1)
bf.rt.for.plot$x = as.factor(bf.rt.for.plot$x)
cols = c(wes_palette("Royal1")[4], wes_palette("Royal1")[3])

bf.mars = c(0.2, 0.2, 0.2, 0.2)
bf.rt <- ggplot(bf.rt.for.plot, aes(x, BF-1, fill=grp)) +
          geom_bar(stat="identity") +
          geom_errorbar(aes(ymin=lower-1, ymax=upper-1), width=.2) +
          coord_flip() +
          scale_fill_manual(this.legend.title, values=cols) +
          scale_y_continuous(breaks=seq(0, 19, by=2), labels=seq(0, 19, by=2)+1) +
          scale_x_discrete(breaks=as.factor(seq(1,5,1)), labels=pnames[bar.order]) +
          geom_hline(yintercept = 0, linetype="dotted", color="black") +
          geom_hline(yintercept = 2, linetype="dashed", color="black") +
          xlab("Alt Models") + ylab("") +
          theme(legend.position="none", 
                text=element_text(size=f_size),
                plot.title = element_text(size=f_size),
                axis.text.x = element_text(size=f_size),
                axis.text.y = element_text(size=f_size)) 
ggsave("e1_rt_bf.bmp", plot=bf.rt, width=5, height=5, units="cm")
#ggtitle("Winning Model: v*sc+iv") +         
# BAYES ACC BARPLOT         
BF = c(acc.bf.for.plot$BF)
names = c(as.character(acc.bf.for.plot$names))
pnames=c("i", "ii", "iii", "iv", "v")
lower = c( acc.bf.for.plot$lower) 
upper = c(acc.bf.for.plot$upper) 
bf.acc.for.plot = data.frame(BF, names, lower, upper)
bar.order = c(1, 2, 3, 4, 5)
bf.acc.for.plot = bf.acc.for.plot[bar.order, ]
bf.acc.for.plot$grp = c("add", "add", "exp", "exp", "exp")
bf.acc.for.plot$x = seq(1, 5, 1)
bf.acc.for.plot$x = as.factor(bf.acc.for.plot$x)
cols = c(wes_palette("Royal1")[4], wes_palette("Royal1")[3])

bf.acc <- ggplot(bf.acc.for.plot, aes(x, BF-1, fill=grp)) +
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=lower-1, ymax=upper-1), width=0.2) +
  coord_flip() +
  scale_fill_manual(this.legend.title, values=cols) +
  scale_y_continuous(breaks=seq(0, 19, by=2), labels=seq(0, 19, by=2)+1) +
  scale_x_discrete(breaks=as.factor(seq(1,5,1)), labels=pnames[bar.order]) +
  geom_hline(yintercept = 0, linetype="dotted", color="black") +
  geom_hline(yintercept = 2, linetype="dashed", color="black") +
  xlab("") + ylab("BF (P[Win]/P[Alt])") +
  theme(legend.position="none", 
        text=element_text(size=f_size),
        axis.text.x = element_text(size=f_size),
        axis.text.y = element_text(size=f_size),
        plot.title = element_text(size=f_size)) 
ggsave("e1_acc_bf.bmp", plot=bf.acc, width=5, height=5, units="cm")





###### EXPERIMENT 2 PLOT
#------------------------------------------------------------------------------------------------------------------
rm(list = ls())
# manual instructions
# 1. SETWD TO SOURCE FILE LOCATION
load("EXP2/ANALYSIS/EXP2_ANALYSIS_MIXDMDLS_BIAS.R")

eb_length = .1
c0 = -(.5*log2(.5))-(.5*log2(.5))
xs = c(-  .6*log2(.6)  -  .4*log2(.4),
       -  .8*log2(.8)  -  .2*log2(.2) )
xs = c0 - xs
####### PLOT USING GGPLOT2
figxs = rep(xs, each = 8)
sum.dat.all.plot$xs = figxs
# get 1 std error of the mean for plotting the error bars
# decided against the 95% within subject confidence intervals as would have to do one
# for each comparison
sum.dat.win.se = ddply( sum.dat.err, .(value, fixprob, valid, rew_cond), summarise,
                        se = sd(mu)/sqrt(length(mu)))
sum.dat.all.plot = inner_join(sum.dat.all.plot, sum.dat.win.se, by=c("value", "fixprob", "valid", "rew_cond"))

l_wd = 1
f_size=8
alph = 0.5
this.legend.title = "Value"
rt.plot <- ggplot(sum.dat.all.plot, aes(x=xs, y=mean, color=as.factor(r1*100))) +
  ylim(0.5,0.8) +
  geom_point(size=1.5) + geom_errorbar(aes(x=xs, ymin=mean-se, ymax=mean+se), width = 0.1, size=0.5) + 
  geom_line(aes(x=xs, y=predict, color=as.factor(r1*100)), linetype=1, lwd=l_wd, alpha=alph) +
  ylab("RT") + xlab("Spatial Certainty (bits)") +
  facet_wrap(~rew_cond*valid, nrow=1) + 
  scale_color_manual(this.legend.title, values=c(wes_palette("Royal1")[1], wes_palette("Royal1")[2]), 
                     guide=guide_legend(direction = "horizontal")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=f_size),
        axis.text.x = element_text(size=f_size),
        axis.text.y = element_text(size=f_size),
        legend.position=c(0.01, 0.1)) 
ggsave("e2_rt_data.bmp", plot=rt.plot, width=12, height=6, units="cm") 


###### NOW ACCURACY
acc.dat.win.err = ddply( acc.dat.err, .(value, fixprob, valid, rew_cond), summarise,
                         se = sd(acc)/sqrt(length(acc)))
acc.dat.plot = inner_join(acc.dat.plot, acc.dat.win.err, by=c("value", "fixprob", "valid", "rew_cond"))
acc.dat.plot$xs = figxs

acc.plot <- ggplot(acc.dat.plot, aes(x=xs, y=acc, color=as.factor(r1*100))) +
  ylim(0.65,1) +
  geom_point(size=1.5) + geom_errorbar(aes(x=xs, ymin=acc-se, ymax=acc+se), width = 0.1, size=0.5) + 
  geom_line(aes(x=xs, y=predict, color=as.factor(r1*100)), linetype=1, lwd=l_wd, alpha=alph) +
  ylab("Acc") + xlab("Spatial Certainty (bits)") +
  facet_wrap(~rew_cond*valid, nrow=1) + 
  scale_color_manual(this.legend.title, values=c(wes_palette("Royal1")[1], wes_palette("Royal1")[2]), 
                     guide=guide_legend(direction = "horizontal")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(size=f_size),
        axis.text.x = element_text(size=f_size),
        axis.text.y = element_text(size=f_size),
        legend.position=c(0.01, 0.1)) 
ggsave("e2_acc_data.bmp", plot=acc.plot, width=12, height=6, units="cm") 


#### BAYES FACTOR
# BAYES RT BARPLOT
# ORDER FOR BARPLOT
dot_size = 2
this.legend.title = "Model Group"
BF = c(rt.bf.for.plot$BF) # adding a zero to place between the two groups of models
names = c(as.character(rt.bf.for.plot$names))
pnames=c("v", "i", "ii", "iii", "iv")
lower = c(rt.bf.for.plot$lower) # adding a zero to place between the two groups of models
upper = c(rt.bf.for.plot$upper) # adding a zero to place between the two groups of models
bf.rt.for.plot = data.frame(BF, names, lower, upper, pnames)
bar.order = c(2, 3, 4, 5, 1)
bf.rt.for.plot = bf.rt.for.plot[bar.order, ]
bf.rt.for.plot$grp = c("add", "add", "add", "add", "exp")
bf.rt.for.plot$x = seq(1, 5, 1)
bf.rt.for.plot$x = as.factor(bf.rt.for.plot$x)
cols = c(wes_palette("Royal1")[4], wes_palette("Royal1")[3])

bf.rt <- ggplot(bf.rt.for.plot, aes(x, BF-1, fill=grp)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=lower-1, ymax=upper-1), width=.2) +
  coord_flip() +
  scale_fill_manual(this.legend.title, values=cols) +
  scale_y_continuous(breaks=seq(0, 19, by=2), labels=seq(0, 19, by=2)+1) +
  scale_x_discrete(breaks=as.factor(seq(1,5,1)), labels=pnames[bar.order]) +
  geom_hline(yintercept = 0, linetype="dotted", color="black") +
  geom_hline(yintercept = 2, linetype="dashed", color="black") +
  xlab("Alt Models") + ylab("") +
  theme(legend.position="none", 
        text=element_text(size=f_size),
        plot.title = element_text(size=f_size),
        axis.text.x = element_text(size=f_size),
        axis.text.y = element_text(size=f_size)) 
ggsave("e2_rt_bf.bmp", plot=bf.rt, width=5, height=5, units="cm")

# ACCURACY
BF = c(acc.bf.for.plot$BF)
names = c(as.character(acc.bf.for.plot$names))
pnames=c("i", "ii", "iii", "iv", "v", "vi", "vii")
lower = c( acc.bf.for.plot$lower) 
upper = c(acc.bf.for.plot$upper) 
bf.acc.for.plot = data.frame(BF, names, lower, upper)
bar.order = c(1, 2, 3, 4, 5, 6, 7)
bf.acc.for.plot = bf.acc.for.plot[bar.order, ]
bf.acc.for.plot$grp = c("add", "add", "add", "add", "add", "exp", "exp")
bf.acc.for.plot$x = seq(1, 7, 1)
bf.acc.for.plot$x = as.factor(bf.acc.for.plot$x)
cols = c(wes_palette("Royal1")[4], wes_palette("Royal1")[3])

bf.acc <- ggplot(bf.acc.for.plot, aes(x, BF-1, fill=grp)) +
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=lower-1, ymax=upper-1), width=0.2) +
  coord_flip() +
  scale_fill_manual(this.legend.title, values=cols) +
  scale_y_continuous(breaks=seq(0, 5, by=2), labels=seq(0, 5, by=2)+1) +
  scale_x_discrete(breaks=as.factor(seq(1,7,1)), labels=pnames[bar.order]) +
  geom_hline(yintercept = 0, linetype="dotted", color="black") +
  geom_hline(yintercept = 2, linetype="dashed", color="black") +
  xlab("") + ylab("BF (P[Win]/P[Alt])") +
  theme(legend.position="none", 
        text=element_text(size=f_size),
        axis.text.x = element_text(size=f_size),
        axis.text.y = element_text(size=f_size),
        plot.title = element_text(size=f_size)) 
ggsave("e2_acc_bf.bmp", plot=bf.acc, width=5, height=5, units="cm")








