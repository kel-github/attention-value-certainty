# Analysis code for Garner & Raymond - Distinct selection mechanisms for when
# predictions and rewards guide visual selective attention
# https://docs.google.com/document/d/1j_9irYF9LbegaxuxnWxSFdmkLe8nSrm6AsWl-ShVCYI/edit?usp=sharing

# (c) K. Garner, University of Birmingham, Sept 2017
# if you have any questions please email getkellygarner@gmail.com

# for each experiment I
# load the pre-cleaned data
# summarise data for plotting
# conduct for RT and accuracy separately
# - an NHST mixed models analysis (did this out of curiosity to see if it
# corroborates the Bayesian approach - it does)
# - a Bayesian analysis to find evidence for the winning model
# use the winning model to predict data for plots
# save to the workspace

# EXP 1 -
#####################################################################################

rm(list=ls())
setwd("~/ADDBIAS_REPOS") # set this to own directory

### tried to use packrat but error with library creation. will come back to this and update
# in the future - sorry world!
# packrat::init(getwd()) # preserve package versions (start R from the project directory and
# get directed to the package library for this project)
library(wesanderson)
#library(lm.beta)
library(lme4)
library(BayesFactor)
library(plyr)
#________________________________________________________________________________________________
load("EXP1/ANALYSIS/EXP1_ANALYSIS_MIXDMDLS_BIAS.R")
# NOTE: this loads the workspace that has the outputs that this code was written to 
# produce using the dataset commented out below.
# load("EXP1/ANALYSIS/exp3_clean_BS_fv_18_03_17.R")

# DATA PREP
# ____________________________________________________________________________________________
#### first, trim the data to remove correct RTs < or > 2.5 std devs from the mu, then calculate regressors
data = dat
data = data[data$cor_resp == data$resp, ]
data = na.omit(data)
tmp <- by(data.frame(rt=data$RT, sub=data$sub, fixprob=data$fixprob, 
                     value=data$value, valid=data$valid, cresp=data$cor_resp), 
                     list(data$loc_prob, data$value, data$session), trim) # remove RTs > or < 2.5 stdevs from the mu
data = as.data.frame(do.call(rbind, tmp))
rm(tmp)
data$prob[data$fixprob == "98:02" & data$valid == "valid"] = .98
data$prob[data$fixprob == "98:02" & data$valid == "invalid"] = .02
data$prob[data$fixprob == "94:06" & data$valid == "valid"] = .94
data$prob[data$fixprob == "94:06" & data$valid == "invalid"] = .06
data$prob[data$fixprob == "90:10" & data$valid == "valid"] = .9
data$prob[data$fixprob == "90:10" & data$valid == "invalid"] = .1
data$prob[data$fixprob == "80:20" & data$valid == "valid"] = .8
data$prob[data$fixprob == "80:20" & data$valid == "invalid"] = .2
data$prob[data$fixprob == "60:40" & data$valid == "valid"] = .6
data$prob[data$fixprob == "60:40" & data$valid == "invalid"] = .4

data$val_num[data$value == "high"] = 50
data$val_num[data$value == "low"] = 1
data$rel_num_comp[data$value == "high"] = 1
data$rel_num_comp[data$value == "low"] = 50

data$cue[data$valid == "valid" & data$value == "high"] = "h2l"
data$cue[data$valid == "invalid" & data$value == "low"] = "h2l"
data$cue[data$valid == "valid" & data$value == "low"] = "l2h"
data$cue[data$valid == "invalid" & data$value == "high"] = "l2h"

# SUMMARISE
# __________________________________________________________________________________________________

###### RTs
sum.dat.all = ddply(data, .(sub, fixprob, value, valid), summarise,
                mu = mean(rt),
                r1 = val_num[1],
                r2 = rel_num_comp[1])

# PLOT AVERAGE DATA
sum.dat.all.plot = ddply(sum.dat.all, .(fixprob, value, valid), summarise,
                         mean = mean(mu),
                         N = length(mu) )
# calc confidence intervals
sum.dat.sub = ddply(sum.dat.all, .(sub), summarise,
                         mean = mean(mu))
sum.dat.all.gmu = mean(sum.dat.sub$mean)
sum.dat.err = sum.dat.all
for (i in levels(sum.dat.err$sub)) sum.dat.err$mu[sum.dat.err$sub == i] = sum.dat.err$mu[sum.dat.err$sub == i] - sum.dat.sub$mean[sum.dat.sub$sub == i] + sum.dat.all.gmu
# crit.t = 2.086
sum.dat.cis = ddply( sum.dat.err, .(value, fixprob, valid), summarise,
                     ci = (sd(mu)/sqrt(length(mu))) * ( 9/8 ) * 2.086
                     )
###### ACC
acc.dat.all = ddply(dat, .( sub, fixprob, value, valid ), summarise, 
                N=length(resp), 
                acc = sum(cor_resp == resp, na.rm=T)/N)
acc.dat.sub = ddply(dat, .(sub), summarise,
                N = length(resp),
                acc = sum(cor_resp == resp, na.rm=T)/N)
acc.dat.gmu = mean(acc.dat.sub$acc)
acc.dat.plot = ddply(dat, .( fixprob, value, valid ), summarise, 
                     N=length(resp), 
                     acc = sum(cor_resp == resp, na.rm=T)/N)
acc.dat.err = acc.dat.all
for (i in levels(acc.dat.err)) acc.dat.err$acc[acc.dat.err$sub == i] = acc.dat.err$acc[acc.dat.err$sub == i] - acc.dat.sub$acc[acc.dat.sub$sub == i] + acc.dat.all.gmu
acc.dat.cis = ddply( acc.dat.err, .(value, fixprob, valid), summarise,
                     ci = (sd(acc)/sqrt(length(acc))) * ( 9/8 ) * 2.086
                    )

# ANALYSE RT AND ACCURACY DATA
# ______________________________________________________________________________________________________________________________________
##### NOTE - I PLAYED AROUD WITH NHST w LINEAR MIXED MODELS. NOT REPORTED IN PAPER. JUST WANTED TO 
# SEE IF IT CORROBORATED THE BAYES PICTURE (IT DOES)
sum.dat.all$info_gain = rep(xs, each = 4, times = 21)
### 1 = check each variable contributes
rt.full = lmer( mu ~ valid + info_gain + r1 + info_gain:valid + info_gain:r1 + (1|sub), REML = FALSE, data = sum.dat.all)
# knock out valid
rt.valid = lmer( mu ~ info_gain + r1 + info_gain:valid + info_gain:r1 + (1|sub), REML = FALSE, data = sum.dat.all)
# knock out info
rt.valinf = lmer( mu ~ valid + info_gain + r1 + info_gain:valid  +  (1|sub), REML = FALSE, data = sum.dat.all)
rt.infgain = lmer( mu ~ valid + info_gain + r1 + (1|sub), REML = FALSE, data = sum.dat.all)
rt.meinfgain = lmer( mu ~ valid + r1 + (1|sub), REML = FALSE, data = sum.dat.all)
# knock out reward
rt.rew.int = lmer( mu ~ valid + info_gain + info_gain:valid + (1|sub), REML = FALSE, data = sum.dat.all)

e1.rt.a1 = anova(rt.full, rt.valid) # **** me valid
e1.rt.a2 = anova(rt.full, rt.valinf )
e1.rt.a3 = anova( rt.valinf, rt.infgain) ### info * valid
e1.rt.a4 = anova(rt.infgain, rt.meinfgain)
e1.rt.a5 = anova(rt.valinf, rt.rew.int) # **** me value
e1.rt.ps = c(2.2e-16, 0.3409, 0.002196, 0.1387, 2.977e-12  )
e1.rt.p.adj = p.adjust(e1.rt.ps, method = "fdr")
e1.rt.p.win = which(e1.rt.p.adj < .05)
summary(rt.full)

# BAYESIAN ANALYSIS (reported in paper)
rt.all.mods = generalTestBF( mu ~ valid * info_gain * r1 + sub, data = sum.dat.all,
                             whichRandom = "sub", neverExclude="^sub$")
rt.all.mods = recompute(rt.all.mods, iterations = 500000) # careful running this - takes a long time, reduce iterations if you 
# want to run it yourself in a speedier fashion
rt.top = head(rt.all.mods)
# comparing winning model to all models at top that included an interaction between value and either validity or info gain
rt.bf2 = rt.top[1]/rt.top[2] # valid + info_gain + valid:info_gain + r1 + info_gain:r1 + sub
rt.bf3 = rt.top[1]/rt.top[3] # valid + r1 + sub 
rt.bf4 = rt.top[1]/rt.top[4] # valid + info_gain + valid:info_gain + r1 + valid:r1 + sub 
rt.bf5 = rt.top[1]/rt.top[5] # valid + info_gain + r1 + sub   
rt.bf6 = rt.top[1]/rt.top[6] # valid + info_gain + valid:info_gain + r1 + valid:r1 + info_gain:r1 + sub

rt.bf.for.plot = data.frame(  BF = c(3.814597, 4.84, 5.189638, 13.53112, 17.23548),
                              upper = c(3.814597 + (3.814597*.0015), 4.83534 + ( 4.83534*.004 ), 5.189638 + (5.189638*.0053), 13.53112 + (13.53112 * .0057), 17.23548 + (17.23548*.0072)),
                              lower = c(3.814597 - (3.814597*.0015),4.83534 - ( 4.83534*.004 ), 5.189638 - (5.189638*.0053), 13.53112 - (13.53112 * .0057), 17.23548 - (17.23548*.0072)),
                              names = c( "v*c + c*va", "v + va", "v*c + v*va", "v + c + va", "v*c + v*va + c*va" ))

# NOW DO ACCURACY DATA 
acc.dat.all$info_gain = rep(xs, each = 4, times = 21)
acc.dat.all$r1[acc.dat.all$value == "high"] = 50
acc.dat.all$r1[acc.dat.all$value == "low"] = 1

# AGAIN, NHST APPROACH FOR FUNSIES
acc.full = lmer( acc ~ valid + info_gain + r1 + info_gain:valid + info_gain:r1 + (1|sub), REML = FALSE, data = acc.dat.all)
# knock out valid
acc.valid = lmer( acc ~ info_gain + r1 + info_gain:valid + info_gain:r1 + (1|sub), REML = FALSE, data = acc.dat.all)
# knock out info
acc.valinf = lmer( acc ~ valid + info_gain + r1 + info_gain:valid  +  (1|sub), REML = FALSE, data = acc.dat.all)
acc.infgain = lmer( acc ~ valid + info_gain + r1 + (1|sub), REML = FALSE, data = acc.dat.all)
acc.meinfgain = lmer( acc ~ valid + r1 + (1|sub), REML = FALSE, data = acc.dat.all)
# knock out reward
acc.rew.int = lmer( acc ~ valid + info_gain + info_gain:valid + (1|sub), REML = FALSE, data = acc.dat.all)

e1.acc.a1 = anova( acc.full, acc.valid) # ********
e1.acc.a2 = anova(acc.full, acc.valinf )
e1.acc.a3 = anova( acc.valinf, acc.infgain)
e1.acc.a4 = anova( acc.infgain, acc.meinfgain)
e1.acc.a5 = anova( acc.valinf, acc.rew.int) # ******
e1.acc.ps = c( 8.29e-14, 0.3209, 0.07941, 0.1635, 0.0004885)
e1.acc.p.adj = p.adjust(e1.acc.ps, method = "fdr")
e1.acc.p.win = which(e1.acc.p.adj < .05)
summary( acc.full )

# BAYESIAN ANALYSIS - in paper
acc.all.mods = generalTestBF( acc ~ valid * info_gain * r1 + sub, data = acc.dat.all,
                             whichRandom = "sub", neverExclude="^sub$")
acc.all.mods = recompute(acc.all.mods, iterations = 500000)
acc.top = head(acc.all.mods) # win = valid + r1 + sub  
acc.bf2 = acc.top[1]/acc.top[2] # valid + info_gain + r1 + sub    
acc.bf3 = acc.top[1]/acc.top[3] # valid + info_gain + valid:info_gain + r1 + sub  
acc.bf4 = acc.top[1]/acc.top[4] # valid + r1 + valid:r1 + sub 
acc.bf5 = acc.top[1]/acc.top[5] # valid + info_gain + r1 + info_gain:r1 + sub
acc.bf6 = acc.top[1]/acc.top[6] # valid + info_gain + valid:info_gain + r1 + info_gain:r1 + sub

acc.bf.for.plot = data.frame(  BF = c(2.57, 3.187716, 5.694131, 8.452145 , 9.173843),
                               upper = c(2.57 + (2.57*.0068), 3.187716 + ( 3.187716*.0078 ), 5.694131 + (5.694131*.0066), 8.452145 + (8.452145 * .0086), 9.173843 + (9.173843*.0095)),
                               lower = c(2.57 - (2.57*.0068), 3.187716 - ( 3.187716*.0078 ), 5.694131 - (5.694131*.0066), 8.452145 - (8.452145 * .0086), 9.173843 - (9.173843*.0095)),
                               names = c( "v + c + va", "v*c + va", "v*va", "v + c*va", "v*c + c*va"))


################################################################################################
########### PLOT DATA PAPER ####################################################################
################################################################################################
# 1 - get data together
# 2 - predict based on winning model
# JUST PREDICTING FIXED EFFECTS FOR PLOT
win.mod.rt = lm( mu ~ valid + info_gain + valid:info_gain + r1, data = sum.dat.all) # top model
sum.dat.all.plot$info_gain = rep(xs, each = 4)
sum.dat.all.plot$r1 = rep(c(50,1), each = 2, times = 5)
sum.dat.all.plot$predict = predict( win.mod.rt, sum.dat.all.plot )

win.mod.acc = lm( acc ~ valid + r1, data = acc.dat.all) # top model
acc.dat.plot$info_gain = rep(xs, each = 4)
acc.dat.plot$r1 = rep(c(50,1), each = 2, times = 5)
acc.dat.plot$predict = predict( win.mod.acc, acc.dat.plot)

# save.image("EXP1/ANALYSIS/EXP1_ANALYSIS_MIXDMDLS_BIAS.R") # run this line if you want to save new 
# stuff

################# EXPERIMENT ONE ANALYSIS COMPLETE - HUZZAH!

# EXP 2
#________________________________________________________________________________________________
rm(list=ls())
library(wesanderson)
library(lm.beta)
library(lme4)
library(BayesFactor)
library(plyr)

# EXP 2 
# IS AN ADDITIVE OR AN INTERACTIVE MODEL BETTER TO ACCOUNT FOR DATA?
setwd("~/ADDBIAS_REPOS") # set this to own directory
##### trim functions - remove RTs 2.5 standard deviations above and below the mean
source("EXP2/ANALYSIS/trim_functions.R")
# NOTE: this loads the workspace that has the outputs that this code was written to 
# produce using the dataset commented out below.
# load("EXP2/ANALYSIS/exp2_clean_BS_v1_28_02_17")
load("EXP2/ANALYSIS/EXP2_ANALYSIS_MIXDMDLS_BIAS.R")

data = dat
### trim rts
tmp <- by(data.frame(rt=data$RT, sub=data$sub, loc_prob=data$loc_prob, 
                     value=data$value, rew_cond=data$rew_cond, valid = data$valid,
                     fixprob = data$fixprob), 
          list(data$sub, data$loc_prob, data$value, data$rew_cond), trim) # remove RTs > or < 2.5 stdevs from the mu
data = as.data.frame(do.call(rbind, tmp))
data = data[!is.na(data$rt),] # remove na's
data$loc_prob_fact = data$loc_prob
data$loc_prob = varhandle::unfactor(data$loc_prob)
data$val_num[data$value == "high"] = 50
data$val_num[data$value == "low"] = 1
data$rel_num_comp[data$value == "high"] = 1
data$rel_num_comp[data$value == "low"] = 50

data$prob[data$fixprob == "80:20" & data$valid == "valid"] = .8
data$prob[data$fixprob == "80:20" & data$valid == "invalid"] = .2
data$prob[data$fixprob == "60:40" & data$valid == "valid"] = .6
data$prob[data$fixprob == "60:40" & data$valid == "invalid"] = .4

data$cue[data$valid == "valid" & data$value == "high"] = "h2l"
data$cue[data$valid == "invalid" & data$value == "low"] = "h2l"
data$cue[data$valid == "valid" & data$value == "low"] = "l2h"
data$cue[data$valid == "invalid" & data$value == "high"] = "l2h"

# SUMMARISE AND GET DATA TO PLOT
# __________________________________________________________________________________________________
###### RTs
sum.dat.all = ddply(data, .(sub, fixprob, value, valid, rew_cond), summarise,
                    mu = mean(rt),
                    r1 = val_num[1],
                    r2 = rel_num_comp[1])

# DATA FOR PLOTTING
sum.dat.all.plot = ddply(sum.dat.all, .(fixprob, value, valid, rew_cond), summarise,
                         mean = mean(mu),
                         N = length(mu) )
sum.dat.sub = ddply(sum.dat.all, .(sub), summarise,
                    mean = mean(mu))
sum.dat.all.gmu = mean(sum.dat.sub$mean)
sum.dat.err = sum.dat.all
for (i in levels(sum.dat.err$sub)) sum.dat.err$mu[sum.dat.err$sub == i] = sum.dat.err$mu[sum.dat.err$sub == i] - sum.dat.sub$mean[sum.dat.sub$sub == i] + sum.dat.all.gmu
# crit.t = 2.086
sum.dat.cis = ddply( sum.dat.err, .(value, fixprob, valid, rew_cond), summarise,
                     ci = (sd(mu)/sqrt(length(mu))) * ( 3/2 ) * 2.086
)
###### ACC
acc.dat.all = ddply(dat, .( sub, fixprob, value, valid, rew_cond ), summarise, 
                    N=length(resp), 
                    acc = sum(cor_resp == resp, na.rm=T)/N)
acc.dat.sub = ddply(dat, .(sub), summarise,
                    N = length(resp),
                    acc = sum(cor_resp == resp, na.rm=T)/N)
acc.dat.gmu = mean(acc.dat.sub$acc)
acc.dat.plot = ddply(dat, .( fixprob, value, valid, rew_cond ), summarise, 
                     N=length(resp), 
                     acc = sum(cor_resp == resp, na.rm=T)/N)
acc.dat.err = acc.dat.all
for (i in levels(acc.dat.err)) acc.dat.err$acc[acc.dat.err$sub == i] = acc.dat.err$acc[acc.dat.err$sub == i] - acc.dat.sub$acc[acc.dat.sub$sub == i] + acc.dat.all.gmu
acc.dat.cis = ddply( acc.dat.err, .(value, fixprob, valid, rew_cond), summarise,
                     ci = (sd(acc)/sqrt(length(acc))) * ( 3/2 ) * 2.086
)

# ANALYSE RT AND ACCURACY DATA;
# ______________________________________________________________________________________________________________________________________
sum.dat.all$info_gain = rep(xs, each = 8, times = 26)
### 1 = check each variable contributes
rt.full = lmer( mu ~ valid + r1 + info_gain + rew_cond + valid:info_gain + info_gain:r1 + info_gain:rew_cond + r1:rew_cond + info_gain:r1:rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)
### drop interactions using info gain
rt.d3way = lmer( mu ~ valid + r1 + info_gain + rew_cond + valid:info_gain + info_gain:r1 + info_gain:rew_cond + r1:rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)
rt.dinf.gain = lmer( mu ~ valid + r1 + info_gain + rew_cond + valid:info_gain + info_gain:r1 + r1:rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)
rt.dinf.r1 = lmer( mu ~ valid + r1 + info_gain + rew_cond + valid:info_gain + r1:rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)
rt.dinf.val = lmer( mu ~ valid + r1 + info_gain + rew_cond +  r1:rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)
rt.dme.inf = lmer( mu ~ valid + r1 + rew_cond +  r1:rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)

#### now drop validity
rt.dval = lmer( mu ~ r1 + rew_cond +  r1:rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)

#### drop rew by rew_cond int
rt.dr1.rewcond = lmer( mu ~ valid + r1 + rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)

#### drop main effect of rew cond
rt.drewcond = lmer( mu ~ valid + r1  + (1|sub), REML = FALSE, data = sum.dat.all)

#### drop main effect of value
rt.drewr1 = lmer( mu ~ valid +  rew_cond + (1|sub), REML = FALSE, data = sum.dat.all)

e2.rt.a1 = anova(rt.full, rt.d3way) # info gain does not interact with rew_cond and reward
e2.rt.a2 = anova(rt.d3way, rt.dinf.gain) # info gain does not interact with rew_cond - decay v fixed
e2.rt.a3 = anova(rt.dinf.gain, rt.dinf.r1) # info gain does not interact with value
e2.rt.a4 = anova(rt.dinf.r1, rt.dinf.val ) # info gain does not interact with validity
e2.rt.a5 = anova(rt.dinf.val, rt.dme.inf ) # no main effect of info gain
e2.rt.a6 = anova(rt.dme.inf, rt.dval) # main effect validity ********
e2.rt.a7 = anova(rt.dme.inf, rt.dr1.rewcond) # no interaction with reward cond
e2.rt.a8 = anova(rt.dr1.rewcond, rt.drewcond) # main effect of reward condition ******
e2.rt.a9 = anova(rt.dr1.rewcond, rt.drewr1) # main effect of loc-value***
e2.rt.p = c( 0.9221, 0.09444, 0.9073, 0.09309, 0.9344, 3.635e-10, 0.4477, 2.406e-09, 0.0004157)
e2.rt.p = p.adjust(e2.rt.p, method = "fdr")
e2.rt.win = which(e2.rt.p < .05)

##### BAYES ANALYSIS
rt.all.mods = generalTestBF( mu ~ valid * info_gain * r1 * rew_cond + sub, data = sum.dat.all,
                             whichRandom = "sub", neverExclude="^sub$" )
rt.all.mods = recompute(rt.all.mods, iterations = 500000)
rt.top = head(rt.all.mods) # valid + r1 + rew_cond + sub 
rt.bf2 = rt.top[1]/rt.top[2] # valid + r1 + valid:r1 + rew_cond + sub
rt.bf3 = rt.top[1]/rt.top[3] # valid + r1 + rew_cond + r1:rew_cond + sub
rt.bf4 = rt.top[1]/rt.top[4] # valid + info_gain + r1 + rew_cond + sub 
rt.bf5 = rt.top[1]/rt.top[5] # valid + r1 + rew_cond + valid:rew_cond + sub 
rt.bf6 = rt.top[1]/rt.top[6] # valid + info_gain + valid:info_gain + r1 + rew_cond + sub 

rt.bf.for.plot = data.frame(  BF = c(4.645585, 5.031692, 6.63436, 6.642919, 8.969172),
                              lower = c(4.645585 - (4.645585*.024), 5.031692 - (5.031692*.02), 6.63436 - (6.63436*.018), 6.642919 - (6.642919*.019), 8.969172 - (8.969172*.018)),
                              upper = c(4.645585 + (4.645585*.024), 5.031692 + (5.031692*.02), 6.63436 + (6.63436*.018), 6.642919 + (6.642919*.019), 8.969172 + (8.969172*.018)),
                              names = c("v*va + rc", "v + va*rc", "v + c + va + rc", "v*rc + val", "v*c + va + rc"))

# accuracy model to report values in paper
acc.report = lmer( acc ~ valid*info_gain + rew_cond + r1 + (1|sub), REML = FALSE, data = acc.dat.all)
# now accuracy
acc.dat.all$info_gain = rep(xs, each = 8, times = 26)
acc.dat.all$r1[acc.dat.all$value == "high"] = 50
acc.dat.all$r1[acc.dat.all$value == "low"] = 1
acc.full = lmer( acc ~ valid + r1 + info_gain + rew_cond + valid:info_gain + info_gain:r1 + info_gain:rew_cond + r1:rew_cond + info_gain:r1:rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)
### drop interactions using info gain
acc.d3way = lmer( acc ~ valid + r1 + info_gain + rew_cond + valid:info_gain + info_gain:r1 + info_gain:rew_cond + r1:rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)
acc.dinf.gain = lmer( acc ~ valid + r1 + info_gain + rew_cond + valid:info_gain + info_gain:r1 + r1:rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)
acc.dinf.r1 = lmer( acc ~ valid + r1 + info_gain + rew_cond + valid:info_gain + r1:rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)
acc.dinf.val = lmer( acc ~ valid + r1 + info_gain + rew_cond +  r1:rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)
acc.dme.inf = lmer( acc ~ valid + r1 + rew_cond +  r1:rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)

#### now drop validity
acc.dval = lmer( acc ~ r1 + rew_cond +  r1:rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)

#### drop rew by rew_cond int
acc.dr1.rewcond = lmer( acc ~ valid + r1 + rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)

#### drop main effect of rew cond
acc.drewcond = lmer( acc ~ valid + r1  + (1|sub), REML = FALSE, data = acc.dat.all)

#### drop main effect of value
acc.drewr1 = lmer( acc ~ valid +  rew_cond + (1|sub), REML = FALSE, data = acc.dat.all)

e2.acc.a1 = anova(acc.full, acc.d3way) # info gain does not interact with rew_cond and reward
e2.acc.a2 = anova(acc.d3way, acc.dinf.gain) # info gain does not interact with rew_cond - decay v fixed
e2.acc.a3 = anova(acc.dinf.gain, acc.dinf.r1) # info gain does not interact with value
e2.acc.a4 = anova(acc.dinf.r1, acc.dinf.val ) # info gain DOES with validity
e2.acc.a5 = anova(acc.dinf.val, acc.dme.inf ) # MAIN EFFECT OF INFO GAIN****
e2.acc.a6 = anova(acc.dme.inf, acc.dval) # main effect validity****
e2.acc.a7 = anova(acc.dme.inf, acc.dr1.rewcond) # no interaction with reward cond
e2.acc.a8 = anova(acc.dr1.rewcond, acc.drewcond) # main effect of reward condition****
e2.acc.a9 = anova(acc.dr1.rewcond, acc.drewr1) # main effect of loc-value
e2.acc.p = c(  0.9911, 0.4061, 0.999, 0.03348, 0.007604, 8.222e-09, 0.9816, 1.031e-08, 0.03825 )

e2.acc.p.adj = p.adjust(e2.acc.p, method = "fdr")
e2.acc.p.win = which( e2.acc.p.adj < .05 )

acc.all.mods = generalTestBF( acc ~ valid * info_gain * r1 * rew_cond + sub, data = acc.dat.all,
                             whichRandom = "sub", neverExclude="^sub$" )
acc.all.mods = recompute(acc.all.mods, iterations = 500000)
acc.top = head(acc.all.mods) # valid + info_gain + valid:info_gain + r1 + rew_cond + sub
acc.bf2 = acc.top[1]/acc.top[2] # valid + info_gain + valid:info_gain + rew_cond + sub
acc.bf3 = acc.top[1]/acc.top[3] # valid + info_gain + r1 + rew_cond + sub 
acc.bf4 = acc.top[1]/acc.top[4] # valid + info_gain + rew_cond + sub  
acc.bf5 = acc.top[1]/acc.top[5] # valid + info_gain + valid:info_gain + r1 + valid:r1 + rew_cond + sub
acc.bf6 = acc.top[1]/acc.top[6] # valid + info_gain + valid:info_gain + r1 + rew_cond + info_gain:rew_cond + sub

### doing specific models because evidence not strong
bf_int1 = lmBF( acc ~ valid + info_gain + r1 + rew_cond + valid:info_gain + info_gain:r1:rew_cond + sub,
                data = acc.dat.all, whichRandom = "sub" )
bf_int1 = recompute(bf_int1, iterations = 500000)
acc.bf7 = acc.top[1]/bf_int1

bf_int2 = lmBF( acc ~ valid + info_gain + r1 + rew_cond + valid:info_gain + info_gain:r1:rew_cond:valid + sub,
                data = acc.dat.all, whichRandom = "sub" )
bf_int2 = recompute(bf_int2, iterations = 500000)
acc.bf8 = acc.top[1]/bf_int2

acc.bf.for.plot = data.frame(  BF = c(1.561972, 1.624413, 2.070495, 2.646015, 3.292799, 4.395263, 3.935153),
                              lower = c(1.561972 - (1.561972*.0075), 1.624413 - (1.624413*.0075), 2.070495 - (2.070495*.016), 2.646015 - (2.646015*.092), 3.292799 - (3.292799*.092), 4.395263 - (4.395263*.053), 3.935153 - (3.935153*.05)),
                              upper = c(1.561972 + (1.561972*.0075), 1.624413 + (1.624413*.0075), 2.070495 + (2.070495*.016), 2.646015 + (2.646015*.092), 3.292799 + (3.292799*.092), 4.395263 + (4.395263*.053), 3.935153 + (3.935153*.05)),
                              names = c("v*c + rc", "v + c + va + rc", "v + c + rc", "v*c + v*va + rc", "v*c + c*rc + va", "v*c + c*va*rc", "v*c + vc*va*rc*v" ))


################################################################################################
########### PREDICTED DATA FOR PLOTS ####################################################################
################################################################################################
# 2 PLOTS - 1 FOR RT AND 1 FOR ACCURACY
# ROW 1 = RT: HIGH VALID & LOW INVALID, LOW VALID & HIGH INVALID, BAYES FACTORS - STATIC CONDITION
# ROW 2 = SAME BUT DECAY CONDITION
# 1 - get data together
# 2 - predict based on winning model
# JUST PREDICTING FIXED EFFECTS
win.mod.rt = lm( mu ~ valid + r1 + rew_cond, data = sum.dat.all)
sum.dat.all.plot$info_gain = rep(xs, each = 8)
sum.dat.all.plot$r1 = rep(c(50,1), each = 4, times = 2)
sum.dat.all.plot$predict = predict( win.mod.rt, sum.dat.all.plot )

win.mod.acc = lm( acc ~ rew_cond + valid + info_gain + valid:info_gain + r1, data = acc.dat.all)
acc.dat.plot$info_gain = rep(xs, each = 8)
acc.dat.plot$r1 = rep(c(50,1), each = 4, times = 2)
acc.dat.plot$predict = predict( win.mod.acc, acc.dat.plot)

# save.image("EXP2_ANALYSIS_MIXDMDLS_BIAS.R") uncomment and use to save new stuff

