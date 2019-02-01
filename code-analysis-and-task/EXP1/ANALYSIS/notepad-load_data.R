# load data and working environment for exp 3
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
dat <- read.csv('KG_exp3_1_CueProb_subs301to323.csv',header = T, na.strings = c("NaN"))
dat$sub <- as.factor(dat$sub)
dat$block <- as.factor(dat$block)
dat$fixprob <- as.factor(dat$fixprob)
levels(dat$fixprob) <- c("98:02","94:06","90:10","80:20", "60:40")
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

# and session factor
dat$session = rep(c(1,2,3,4), each = 500, times = length(levels(dat$sub)))
dat$session = as.factor(dat$session)

dat$loc_prob[dat$fixprob == "98:02" & dat$valid == "invalid"] = .1
dat$loc_prob[dat$fixprob == "94:06" & dat$valid == "invalid"] = .1
dat$loc_prob[dat$fixprob == "90:10" & dat$valid == "invalid"] = .1
dat$loc_prob[dat$fixprob == "80:20" & dat$valid == "invalid"] = .2
dat$loc_prob[dat$fixprob == "60:40" & dat$valid == "invalid"] = .4
dat$loc_prob[dat$fixprob == "60:40" & dat$valid == "valid"] = .6
dat$loc_prob[dat$fixprob == "80:20" & dat$valid == "valid"] = .8
dat$loc_prob[dat$fixprob == "90:10" & dat$valid == "valid"] = .9
dat$loc_prob[dat$fixprob == "94:06" & dat$valid == "valid"] = .92
dat$loc_prob[dat$fixprob == "98:02" & dat$valid == "valid"] = .96
dat$loc_num = dat$loc_prob
dat$loc_prob <- as.factor(dat$loc_prob)
levels(dat$loc_prob) <- c(".1",".2",".4",".6",".8",".9",".92",".96")

dat$value[dat$loc_prob == ".1" & dat$rel == "cue_low"] = 1
dat$value[dat$loc_prob == ".1" & dat$rel == "cue_low"] = 1
dat$value[dat$loc_prob == ".1" & dat$rel == "cue_low"] = 1
dat$value[dat$loc_prob == ".2" & dat$rel == "cue_low"] = 1
dat$value[dat$loc_prob == ".4" & dat$rel == "cue_low"] = 1
dat$value[dat$loc_prob == ".6" & dat$rel == "cue_high"] = 1
dat$value[dat$loc_prob == ".8" & dat$rel == "cue_high"] = 1
dat$value[dat$loc_prob == ".9" & dat$rel == "cue_high"] = 1
dat$value[dat$loc_prob == ".92" & dat$rel == "cue_high"] = 1
dat$value[dat$loc_prob == ".96" & dat$rel == "cue_high"] = 1


dat$value[dat$loc_prob == ".1" & dat$rel == "cue_high"] = 2
dat$value[dat$loc_prob == ".1" & dat$rel == "cue_high"] = 2
dat$value[dat$loc_prob == ".1" & dat$rel == "cue_high"] = 2
dat$value[dat$loc_prob == ".2" & dat$rel == "cue_high"] = 2
dat$value[dat$loc_prob == ".4" & dat$rel == "cue_high"] = 2
dat$value[dat$loc_prob == ".6" & dat$rel == "cue_low"] = 2
dat$value[dat$loc_prob == ".8" & dat$rel == "cue_low"] = 2
dat$value[dat$loc_prob == ".9" & dat$rel == "cue_low"] = 2
dat$value[dat$loc_prob == ".92" & dat$rel == "cue_low"] = 2
dat$value[dat$loc_prob == ".96" & dat$rel == "cue_low"] = 2

dat$value <- as.factor(dat$value)
levels(dat$value) <- c("high", "low")

dat$trial_num = rep(c(1:2000), times = length(unique(dat$sub)))

#### add previous trial's validity

get.prev.val <- function(data){
  data$prev_val <- NA
    for (i in 2:max(data$trial_num)){
      data$prev_val[i] = data$valid[i-1] 
    }
  return(data)
}

tmp = by(dat, dat$sub, get.prev.val)
dat = as.data.frame(do.call(rbind, tmp))
rm(tmp)
dat$prev_val <- as.factor(dat$prev_val)
levels(dat$prev_val) = c("valid", "invalid")


