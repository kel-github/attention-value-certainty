# load data and working environment for exp 2
#rm(list=ls())

#########################################################################################################
############################################# pckges and wd #############################################
#########################################################################################################
require(wesanderson)
require(plyr)
#require(ez)
require(ggplot2)
#require(pwr)
require(reshape2)
###### set working directory

##### trim functions - remove RTs 2.5 standard deviations above and below the mean
#source("~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/functions/trim_functions.R")

#########################################################################################################
###################################### read in dat          #############################################
#########################################################################################################
dat <- read.csv('KG_exp2_1_CueProb_subs501to526.csv',header = T, na.strings = c("NaN"))
# 20800
dat$RT[which(is.na(dat$resp))] = NA
# 608
dat$sub <- as.factor(dat$sub)
dat$block <- as.factor(dat$block)

dat$rew_cond[ dat$fixprob == 1 | dat$fixprob == 3 ] = 1
dat$rew_cond[ dat$fixprob == 2 | dat$fixprob == 4 ] = 2
dat$rew_cond <- as.factor(dat$rew_cond)
levels(dat$rew_cond) <- c("static", "decay")

dat$fixprob[ dat$fixprob == 2 ] = 1
dat$fixprob[ dat$fixprob > 2 ] = 2
dat$fixprob <- as.factor(dat$fixprob)
levels(dat$fixprob) <- c( ".4/.6",".2/.8" )

dat$leftval <- as.factor(dat$leftval)
dat$rightval <- as.factor(dat$rightval)
levels(dat$leftval) <- c("high","low")
levels(dat$rightval) <- c("high","low")
dat$cuedir <- as.factor(dat$cuedir)
levels(dat$cuedir) <- c("left","right")
dat$tgtloc <- as.factor(dat$tgtloc)
levels(dat$tgtloc) <- c("left","right")
dat$cor_resp <- as.factor(dat$cor_resp)
levels(dat$cor_resp) <- c("h","n")
dat$resp <- as.factor(dat$resp)
levels(dat$resp) <- c("h","n")

########### create h v l factor
dat$rel[dat$leftval == "high" & dat$cuedir == "left"] <- 1
dat$rel[dat$rightval == "high" & dat$cuedir == "right"] <- 1
dat$rel[dat$leftval == "low" & dat$cuedir == "left"] <- 2
dat$rel[dat$rightval == "low" & dat$cuedir == "right"] <- 2
dat$rel <- as.factor(dat$rel)
levels(dat$rel) <- c("cue_high","cue_low")
dat$valid[dat$cuedir ==  dat$tgtloc] <- 1
dat$valid[dat$cuedir !=  dat$tgtloc] <- 2
dat$valid <- as.factor(dat$valid)
levels(dat$valid) <- c("valid","invalid")

dat$loc_prob[dat$fixprob == ".4/.6" & dat$valid == "invalid"] = .4
dat$loc_prob[dat$fixprob == ".2/.8" & dat$valid == "invalid"] = .2
dat$loc_prob[dat$fixprob == ".4/.6" & dat$valid == "valid"] = .6
dat$loc_prob[dat$fixprob == ".2/.8" & dat$valid == "valid"] = .8

dat$loc_num = dat$loc_prob
dat$loc_prob <- as.factor(dat$loc_prob)
levels(dat$loc_prob) <- c( ".2",".4",".6",".8"  )

dat$value[dat$loc_prob == ".2" & dat$rel == "cue_low"] = 1
dat$value[dat$loc_prob == ".4" & dat$rel == "cue_low"] = 1
dat$value[dat$loc_prob == ".6" & dat$rel == "cue_high"] = 1
dat$value[dat$loc_prob == ".8" & dat$rel == "cue_high"] = 1
dat$value[dat$loc_prob == ".2" & dat$rel == "cue_high"] = 2
dat$value[dat$loc_prob == ".4" & dat$rel == "cue_high"] = 2
dat$value[dat$loc_prob == ".6" & dat$rel == "cue_low"] = 2
dat$value[dat$loc_prob == ".8" & dat$rel == "cue_low"] = 2

dat$value <- as.factor(dat$value)
levels(dat$value) = c( "high", "low" )



