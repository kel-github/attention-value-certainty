# K. Garner - 11/09/17 Code plots data points for each cue (validity x spatial-certainty x spatial-value ), 
# model predictions = dashed lines, and BFs for winning model relative to next 5.
###### EXPERIMENT 1 PLOT
#------------------------------------------------------------------------------------------------------------------
rm(list = ls())
setwd("~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/EXP_3_CUE_PROB_OUT/exp3")
load("EXP1_ANALYSIS_MIXDMDLS.R")
# PLOT
svg(filename="~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/EXP_3_CUE_PROB_OUT/Figures/EXP1_RT_ACCs.svg",
    width = 12, height = 6) 
par(mfrow = c(1, 2))
eb_length = .1
c0 = -(.5*log2(.5))-(.5*log2(.5))
xs = c( - .98*log2(.98)  - .02*log2(.02) ,
        - .92*log2(.92)  - .08*log2(.08) ,
        -  .9*log2(.9)  - .1*log2(.1) ,
        -  .8*log2(.8)  -  .2*log2(.2) ,
        -  .6*log2(.6)  -  .4*log2(.4)  )
xs = c0 - xs

svg(filename="~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/EXP_3_CUE_PROB_OUT/Figures/EXP1_RT_ACC_MDL_FITS.svg",
    width = 8, height = 8) 
par(mfrow = c(2,3))
#### RT ####
ylims = c(min(sum.dat.all.plot$mean)-.15, max(sum.dat.all.plot$mean)+.15)
plot(xs, sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[2], xlab = "information gain in bits", ylab = "RT", ylim = ylims, 
     main = "valid cue", axes = FALSE)
points(xs, c(sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[1] )
axis(side = 1)
axis(side = 2)
legend(0.1, 0.55, c("50", "1"), pch = 19, lwd = c(3,3), col = c(wes_palette("Royal1")[2], wes_palette("Royal1")[1]), bty = "n")
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid"] - (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "valid"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid"] + (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "valid"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid"] - (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "valid"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid"] + (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "valid"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])

plot(xs, sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[1], xlab = "information gain in bits", ylab = "RT", ylim = ylims, 
     main = "invalid cue", axes = FALSE)
points(xs, c(sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[2] )
axis(side = 1)
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid"] - (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "invalid"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid"] + (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "invalid"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid"] - (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "invalid"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid"] + (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "invalid"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])

barcenters = barplot(rt.bf.for.plot$BF, col = wes_palette("Royal1")[4], main = "~ v*c + va", 
                     names.arg = rt.bf.for.plot$names, ylab = "BF", las = 2, ylim = c(0, 20))
arrows(barcenters, rt.bf.for.plot$lower, barcenters, rt.bf.for.plot$upper, lwd = 1.5, angle = 90, code = 3, length = eb_length)
abline(h=3, lty = 4, lwd = 3)

#### ACCURACY ####
ylims = c(min(acc.dat.plot$acc)-.1, max(acc.dat.plot$acc)+.1)
plot(xs, c(acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid"]), 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[2], xlab = "information gain in bits", ylab = "Accuracy", ylim = ylims, 
     axes = FALSE)
points(xs, c(acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[1] )
axis(side = 1)
axis(side = 2)
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid"] - (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "valid"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid"] + (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "valid"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid"] - (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "valid"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid"] + (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "valid"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
points(xs, acc.dat.plot$predict[ acc.dat.plot$value == "high" &  acc.dat.plot$valid == "valid"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])
points(xs,  acc.dat.plot$predict[ acc.dat.plot$value == "low" &  acc.dat.plot$valid == "valid"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])

plot(xs, c(acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid"]), 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[1], xlab = "information gain in bits", ylab = "Accuracy", ylim = ylims, 
     axes = FALSE)
points(xs, c(acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[2] )
axis(side = 1)
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid"] - (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "invalid"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid"] + (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "invalid"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid"] - (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "invalid"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid"] + (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "invalid"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
points(xs, acc.dat.plot$predict[ acc.dat.plot$value == "low" &  acc.dat.plot$valid == "invalid"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])
points(xs,  acc.dat.plot$predict[ acc.dat.plot$value == "high" &  acc.dat.plot$valid == "invalid"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])

barcenters = barplot(acc.bf.for.plot$BF, col = wes_palette("Royal1")[4], main = "~ v + va", 
                     names.arg = acc.bf.for.plot$names, ylab = "BF", las = 2, ylim = c(0,10))
arrows(barcenters, acc.bf.for.plot$lower, barcenters, acc.bf.for.plot$upper, lwd = 1.5, angle = 90, code = 3, length = eb_length)
abline(h=3, lty = 4, lwd = 3)

dev.off()

###### EXPERIMENT 2 PLOT
#------------------------------------------------------------------------------------------------------------------
rm(list = ls())
setwd("~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/EXP_3_CUE_PROB_OUT/exp5")
load("EXP2_ANALYSIS_MIXDMDLS_BIAS.R")

# 3 - plot data as dots, winning model as lines
svg(filename="~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/EXP_3_CUE_PROB_OUT/Figures/EXP2_RT_ACC_MDL_FITS.svg",
    width = 12, height = 8)
par(mfrow = c(2,5))
ylims = c(min(sum.dat.all.plot$mean)-.15, max(sum.dat.all.plot$mean)+.15)
plot(xs, sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "static"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[2], xlab = "information gain in bits", ylab = "RT", ylim = ylims, xlim = c(-.029, .35),
     main = "valid", axes = FALSE)
points(xs, c(sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "static"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[1] )
axis(side = 1)
axis(side = 2)
legend(0.1, 0.55, c("50", "1"), pch = 19, lwd = c(3,3), col = c(wes_palette("Royal1")[2], wes_palette("Royal1")[1]), bty = "n")
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "static"] - (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "valid" & sum.dat.all.plot$rew_cond == "static"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "static"] + (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "valid" & sum.dat.all.plot$rew_cond == "static"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "static"] - (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "valid" & sum.dat.all.plot$rew_cond == "static"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "static"] + (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "valid" & sum.dat.all.plot$rew_cond == "static"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "static"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "static"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])

plot(xs, sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[1], xlab = "information gain in bits", ylab = "RT", ylim = ylims, xlim = c(-.029, .35),
     main = "invalid", axes = FALSE)
points(xs, c(sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[2] )
axis(side = 1)

arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"] - (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"] + (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"] - (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"] + (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "static"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])


plot(xs, sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[2], xlab = "information gain in bits", ylab = "RT", ylim = ylims, xlim = c(-.029, .35),
     axes = FALSE, main = "valid")
points(xs, c(sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[1] )
axis(side = 1)
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"] - (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"] + (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"] - (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"] + (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "valid" & sum.dat.all.plot$rew_cond == "decay"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])

plot(xs, sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[1], xlab = "information gain in bits", ylab = "RT", ylim = ylims, xlim = c(-.029, .35),
     axes = FALSE, main = "invalid")
points(xs, c(sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[2] )
axis(side = 1)

arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"] - (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"] + (sum.dat.cis$ci[ sum.dat.cis$value == "low" & sum.dat.cis$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
arrows( x0= xs,
        y0= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"] - (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"]/2),
        x1= xs,
        y1= sum.dat.all.plot$mean[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"] + (sum.dat.cis$ci[ sum.dat.cis$value == "high" & sum.dat.cis$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "low" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])
points(xs, sum.dat.all.plot$predict[sum.dat.all.plot$value == "high" & sum.dat.all.plot$valid == "invalid" & sum.dat.all.plot$rew_cond == "decay"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])

barcenters = barplot(rt.bf.for.plot$BF, col = wes_palette("Royal1")[4], main = "~ v + va + rc", 
                     names.arg = rt.bf.for.plot$names, ylab = "BF", las = 2, ylim = c(0, 10))
arrows(barcenters, rt.bf.for.plot$lower, barcenters, rt.bf.for.plot$upper, lwd = 1.5, angle = 90, code = 3, length = eb_length)
abline(h=3, lty = 4, lwd = 3)
### ACCURACY
ylims = c(min(acc.dat.plot$acc)-.15, max(acc.dat.plot$acc)+.15)
plot(xs, acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "static"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[2], xlab = "information gain in bits", ylab = "Accuracy", ylim = ylims, xlim = c(-.029, .35),
     axes = FALSE)
points(xs, c(acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "static"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[1] )
axis(side = 1)
axis(side = 2)

arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "static"] - (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "valid" & acc.dat.plot$rew_cond == "static"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "static"] + (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "valid" & acc.dat.plot$rew_cond == "static"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "static"] - (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "valid" & acc.dat.plot$rew_cond == "static"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "static"] + (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "valid" & acc.dat.plot$rew_cond == "static"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
points(xs, acc.dat.plot$predict[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "static"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])
points(xs, acc.dat.plot$predict[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "static"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])

plot(xs, acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "static"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[1], xlab = "information gain in bits", ylab = "Accuracy", ylim = ylims, xlim = c(-.029, .35),
     axes = FALSE)
points(xs, c(acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "static"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[2] )
axis(side = 1)

arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "static"] - (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "invalid" & acc.dat.plot$rew_cond == "static"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "static"] + (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "invalid" & acc.dat.plot$rew_cond == "static"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "static"] - (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "invalid" & acc.dat.plot$rew_cond == "static"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "static"] + (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "invalid" & acc.dat.plot$rew_cond == "static"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
points(xs, acc.dat.plot$predict[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "static"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])
points(xs, acc.dat.plot$predict[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "static"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])


plot(xs, acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "decay"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[2], xlab = "information gain in bits", ylab = "Accuracy", ylim = ylims, xlim = c(-.029, .35),
     axes = FALSE)
points(xs, c(acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "decay"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[1] )
axis(side = 1)
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "decay"] - (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "valid" & acc.dat.plot$rew_cond == "decay"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "decay"] + (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "valid" & acc.dat.plot$rew_cond == "decay"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "decay"] - (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "valid" & acc.dat.plot$rew_cond == "decay"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "decay"] + (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "valid" & acc.dat.plot$rew_cond == "decay"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
points(xs, acc.dat.plot$predict[acc.dat.plot$value == "high" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "decay"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])
points(xs, acc.dat.plot$predict[acc.dat.plot$value == "low" & acc.dat.plot$valid == "valid" & acc.dat.plot$rew_cond == "decay"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])

plot(xs, acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "decay"], 
     pch = 19, cex = 2.5, col = wes_palette("Royal1")[1], xlab = "information gain in bits", ylab = "Accuracy", ylim = ylims, xlim = c(-.029, .35),
     axes = FALSE)
points(xs, c(acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "decay"]), 
       pch = 19, cex = 2.5, col = wes_palette("Royal1")[2] )
axis(side = 1)

arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "decay"] - (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "invalid" & acc.dat.plot$rew_cond == "decay"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "decay"] + (acc.dat.cis$ci[ acc.dat.cis$value == "low" & acc.dat.cis$valid == "invalid" & acc.dat.plot$rew_cond == "decay"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[1])
arrows( x0= xs,
        y0= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "decay"] - (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "invalid" & acc.dat.plot$rew_cond == "decay"]/2),
        x1= xs,
        y1= acc.dat.plot$acc[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "decay"] + (acc.dat.cis$ci[ acc.dat.cis$value == "high" & acc.dat.cis$valid == "invalid" & acc.dat.plot$rew_cond == "decay"]/2),
        angle=90, code=3, length = eb_length, col = wes_palette("Royal1")[2])
points(xs, acc.dat.plot$predict[acc.dat.plot$value == "low" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "decay"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[1])
points(xs, acc.dat.plot$predict[acc.dat.plot$value == "high" & acc.dat.plot$valid == "invalid" & acc.dat.plot$rew_cond == "decay"], lwd = 2.5,
       type = "l", lty = 4, col = wes_palette("Royal1")[2])

barcenters = barplot(acc.bf.for.plot$BF, col = wes_palette("Royal1")[4], main = "~ v*c + va + rc", 
                     names.arg = acc.bf.for.plot$names, ylab = "BF", las = 2, ylim = c(0, 5))
arrows(barcenters, acc.bf.for.plot$lower, barcenters, acc.bf.for.plot$upper, lwd = 1.5, angle = 90, code = 3, length = .01)
abline(h=3, lty = 4, lwd = 3)
dev.off()





