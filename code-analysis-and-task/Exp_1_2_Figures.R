# K. Garner - 11/09/17 Code plots data points for each cue (validity x spatial-certainty x spatial-value ), 
# model predictions = dashed lines, and BFs for winning model relative to next 5.
###### EXPERIMENT 1 PLOT
#------------------------------------------------------------------------------------------------------------------
library(wesanderson)
rm(list = ls())
# manual instructions
# 1. SETWD TO SOURCE FILE LOCATION
# 2. run code @ bottom that defines plot.dvs function

#setwd("~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/EXP_3_CUE_PROB_OUT/exp3")
load("EXP1/ANALYSIS/EXP1_ANALYSIS_MIXDMDLS_BIAS.R")
# PLOT
pdf(file=paste(getwd(), '/Fig2_Exp1_Results.pdf', sep=""),
    width = 5, height = 5) 

eb_length = .1
c0 = -(.5*log2(.5))-(.5*log2(.5)) # defined according to Prinzmetal et al, 2015, JEP:HPP
xs = c( - .98*log2(.98)  - .02*log2(.02) ,
        - .92*log2(.92)  - .08*log2(.08) ,
        -  .9*log2(.9)   -  .1*log2(.1) ,
        -  .8*log2(.8)   -  .2*log2(.2) ,
        -  .6*log2(.6)   -  .4*log2(.4)  )
xs = c0 - xs
rows = c(2,3)
omas = c(2, 4, 2, 2)  # bottom, left, top, right
par(mfrow = rows,    
    oma = omas)

#### PLOT THE BEHAVIOURAL DATA (MU RTs -/+ 95% CIs, MODEL FIT = DOTTED LINES)
# note: this calls a function defined at the bottom of the script
plot.dvs(xs, sum.dat.all.plot, sum.dat.cis, xlab = "SC", ylab = "RT", ylims=c(.5, 1), xlims = c(0, .9), mars = c(2, 2, 3, 1) + 0.1,
         dv = "mean", dat_a_value = "high", dat_b_value = "low", cue_validity = "valid", legend_on = 0, leg.x = 0.5, leg.y = 0.55, 
         dot_size = 2, fit_line = 3, title = "valid cue", x_axis = 1, y_axis = 1)
plot.dvs(xs, sum.dat.all.plot, sum.dat.cis, xlab = " ", ylab = " ", ylims=c(.5, 1), xlims = c(0, .9), mars = c(2, 2, 2, 1) + 0.1,
         dv = "mean", dat_a_value = "high", dat_b_value = "low", cue_validity = "invalid", legend_on = 1, leg.x = 0.2, leg.y = 0.65, 
         legend_low = 1, legend_high = 50, dot_size = 2, fit_line = 3, title = "invalid cue", x_axis = 1, y_axis = 0)

# BAYES BARPLOT
# ORDER FOR BARPLOT
dot_size = 2
par(mar = c(2, 5, 5, 1) + 0.1)
tmp.dat = c(rt.bf.for.plot$BF, 0) # adding a zero to place between the two groups of models
tmp.names = c(as.character(rt.bf.for.plot$names), "")
tmp.lower = c( rt.bf.for.plot$lower, 0) # adding a zero to place between the two groups of models
tmp.upper = c(rt.bf.for.plot$upper, 0) # adding a zero to place between the two groups of models
bar.order = c(2, 4, 6, 1, 3, 5)
cols = c(wes_palette("Royal1")[4], wes_palette("Royal1")[4], wes_palette("Royal1")[3], wes_palette("Royal1")[3], wes_palette("Royal1")[3], wes_palette("Royal1")[3])
barcenters = barplot(tmp.dat[bar.order], col = cols, 
                     ylab = "", las = 2, ylim = c(0, 20), 
                     axes = FALSE)
#names.arg = tmp.names[bar.order], 
arrows(barcenters, tmp.lower[bar.order], barcenters, tmp.upper[bar.order], lwd = 1.5, angle = 90, code = 3, length = eb_length)

abline(h=3, lty = 1, lwd = 2)
abline(h=1, lty = 1, lwd = 2)
axis(side = 1, tck = 0, labels = FALSE)
axis(side = 2, tck = 0, las = 2, cex.axis = dot_size)
mtext(side = 2, "BF (IV + SC)", line = 3, cex = dot_size)
#mtext("~ v*c + va (winning model)", side=1, line = 2, cex.main = dot_size)

###### NOW ACCURACY
plot.dvs(xs, acc.dat.plot, acc.dat.cis, xlab = " ", ylab = "Acc", ylims = c(0.6, 1), xlims = c(0, .9), mars = c(2, 2, 3, 1) + 0.1,
         dv = "acc", dat_a_value = "high", dat_b_value = "low", cue_validity = "valid", legend_on = 0, leg.x = 0.5, leg.y = 0.55, 
         dot_size = 2, fit_line = 3, title = " ", x_axis = 1, y_axis = 1)
plot.dvs(xs, acc.dat.plot, acc.dat.cis, xlab = " ", ylab = " ", ylims = c(0.6, 1), xlims = c(0, .9), mars = c(2, 2, 2, 1) + 0.1,
         dv = "acc", dat_a_value = "high", dat_b_value = "low", cue_validity = "invalid", legend_on = 0, 
         dot_size = 2, fit_line = 3, title = " ", x_axis = 1, y_axis =0 )

par(mar = c(2, 5, 5, 1) + 0.1)
tmp.dat = c(acc.bf.for.plot$BF, 0)
tmp.names = c(as.character(acc.bf.for.plot$names), "")
tmp.lower = c( acc.bf.for.plot$lower, 0) 
tmp.upper = c(acc.bf.for.plot$upper, 0) 
bar.order = c(1, 2, 6, 3, 4, 5)
cols = c(wes_palette("Royal1")[4], wes_palette("Royal1")[4], wes_palette("Royal1")[3], wes_palette("Royal1")[3], wes_palette("Royal1")[3], wes_palette("Royal1")[3])
barcenters = barplot(tmp.dat[ bar.order ], col = cols, 
                     names.arg = tmp.names[ bar.order ], ylab = "", las = 2, ylim = c(0,10),
                     axes = FALSE)
arrows(barcenters, tmp.lower[ bar.order ], barcenters, tmp.upper[ bar.order ], lwd = 1.5, angle = 90, code = 3, length = eb_length)
abline(h=3, lty = 1, lwd = 2)
abline(h=1, lty = 1, lwd = 2)
#axis(side = 1, tck = 0, labels = FALSE)
axis(side = 2, tck = 0, las = 2, cex.axis = dot_size)
mtext(side = 2, "BF (valid + IV)", line = 3, cex = dot_size)
#title("~ v + va (winning model)", line = -2, cex.main = dot_size)

dev.off()


###### EXPERIMENT 2 PLOT
#------------------------------------------------------------------------------------------------------------------
rm(list = ls())
# manual instructions
# 1. SETWD TO SOURCE FILE LOCATION
# 2. run code @ bottom that defines plot.dvs function
load("EXP2/ANALYSIS/EXP2_ANALYSIS_MIXDMDLS_BIAS.R")

svg(filename="~/Dropbox/BHAMPROJECTS/RelValue_StudyProgramme/Analysis/EXP_3_CUE_PROB_OUT/Figures/EXP2_RT_ACCs.svg",
    width = 9, height = 9) 

eb_length = .1
c0 = -(.5*log2(.5))-(.5*log2(.5))
xs = c( -  .8*log2(.8)  -  .2*log2(.2) ,
        -  .6*log2(.6)  -  .4*log2(.4)  )
xs = c0 - xs
rows = c(2,5)
omas = c(3, 3, 3, 3) + 0.1 # bottom, left, top, right
par(mfrow = rows,    
    oma = omas)
rt.ylims =c(min(sum.dat.all.plot$mean)-.15, max(sum.dat.all.plot$mean)+.15)
acc.ylims = c(min(acc.dat.plot$acc)-.15, max(acc.dat.plot$acc)+.1)
xs = c(xs[2], xs[1])
plot.dvs(xs, sum.dat.all.plot[sum.dat.all.plot$rew_cond == "static", ], sum.dat.cis[sum.dat.all.plot$rew_cond == "static",], 
         xlab = "SC", ylab = "RT", xlims = c(0, .3), rt.ylims, mars = c(2, 3, 2, 1) + 0.1,
         dv = "mean", dat_a_value = "high", dat_b_value = "low", cue_validity = "valid", legend_on = 0, 
         dot_size = 2, fit_line = 3, title = "valid cue", x_axis = 1, y_axis = 1)
plot.dvs(xs, sum.dat.all.plot[sum.dat.all.plot$rew_cond == "static", ], sum.dat.cis[sum.dat.all.plot$rew_cond == "static",], 
         xlab = " ", ylab = " ", xlims = c(0, .3), rt.ylims, mars = c(2, 2, 2, 1) + 0.1, 
         dv = "mean", dat_a_value = "high", dat_b_value = "low", cue_validity = "invalid", legend_on = 0, 
         dot_size = 2, fit_line = 3, title = "invalid cue", x_axis = 1, y_axis = 0)
plot.dvs(xs, sum.dat.all.plot[sum.dat.all.plot$rew_cond == "decay", ], sum.dat.cis[sum.dat.all.plot$rew_cond == "static",], 
         xlab = " ", ylab = " ", xlims = c(0, .3), rt.ylims, mars = c(2, 2, 2, 1) + 0.1, 
         dv = "mean", dat_a_value = "high", dat_b_value = "low", cue_validity = "valid", legend_on = 0, 
         dot_size = 2, fit_line = 3, title = "valid cue", x_axis = 1, y_axis = 0)
plot.dvs(xs, sum.dat.all.plot[sum.dat.all.plot$rew_cond == "decay", ], sum.dat.cis[sum.dat.all.plot$rew_cond == "static",], 
         xlab = " ", ylab = " ", xlims = c(0, .3), rt.ylims, mars = c(2, 2, 2, 1) + 0.1, 
         dv = "mean", dat_a_value = "high", dat_b_value = "low", cue_validity = "invalid", legend_on = 1, 
         leg.x = -0.01, leg.y = 0.55, legend_low = 1, legend_high = 50, dot_size = 2, fit_line = 3, 
         title = "invalid cue", x_axis = 1, y_axis = 0)

#### BAYES FACTOR
# BAYES BARPLOT
# ORDER FOR BARPLOT
dot_size = 2
par(mar = c(2, 3, 2, 1) + 0.1)
tmp.dat = c(rt.bf.for.plot$BF, 0)
tmp.names = c(as.character(rt.bf.for.plot$names), "")
tmp.lower = c( rt.bf.for.plot$lower, 0) 
tmp.upper = c( rt.bf.for.plot$upper, 0) 
bar.order = c(2, 3, 4, 5, 6, 1)
cols = c(wes_palette("Royal1")[4], wes_palette("Royal1")[4], wes_palette("Royal1")[4], wes_palette("Royal1")[4], wes_palette("Royal1")[3], wes_palette("Royal1")[3])
barcenters = barplot(tmp.dat[bar.order], col = cols, 
                     names.arg = tmp.names[bar.order], ylab = "", las = 2, ylim = c(0, 10), 
                     axes = FALSE)
arrows(barcenters, tmp.lower[bar.order], barcenters, tmp.upper[bar.order], lwd = 1.5, angle = 90, code = 3, length = eb_length*0.5)
abline(h=3, lty = 1, lwd = 2)
abline(h=1, lty = 1, lwd = 2)
axis(side = 1, tck = 0, labels = FALSE)
axis(side = 2, tck = 0, las = 2, cex.axis = dot_size)
mtext(side = 2, "BF (valid + IV)", line = 2, cex = 1.5)
#title("~ v + va + rc", line = -2, cex.main = dot_size)

plot.dvs(xs, acc.dat.plot[acc.dat.plot$rew_cond == "static", ], acc.dat.cis[acc.dat.cis$rew_cond == "static",], 
         xlab = " ", ylab = "Acc", xlims = c(0, .3), acc.ylims, mars = c(2, 3, 2, 1) + 0.1,
         dv = "acc", dat_a_value = "high", dat_b_value = "low", cue_validity = "valid", legend_on = 0, 
         dot_size = 2, fit_line = 3, title = "", x_axis = 1, y_axis = 1)
plot.dvs(xs, acc.dat.plot[acc.dat.plot$rew_cond == "static", ], acc.dat.cis[acc.dat.cis$rew_cond == "static",], 
         xlab = " ", ylab = " ", xlims = c(0, .3), acc.ylims, mars = c(2, 2, 2, 1) + 0.1,
         dv = "acc", dat_a_value = "high", dat_b_value = "low", cue_validity = "invalid", legend_on = 0, 
         dot_size = 2, fit_line = 3, title = "", x_axis = 1, y_axis = 0)
plot.dvs(xs, acc.dat.plot[acc.dat.plot$rew_cond == "decay", ], acc.dat.cis[acc.dat.cis$rew_cond == "decay",], 
         xlab = " ", ylab = " ", xlims = c(0, .3), acc.ylims, mars = c(2, 2, 2, 1) + 0.1,
         dv = "acc", dat_a_value = "high", dat_b_value = "low", cue_validity = "valid", legend_on = 0, 
         dot_size = 2, fit_line = 3, title = "", x_axis = 1, y_axis = 0)
plot.dvs(xs, acc.dat.plot[acc.dat.plot$rew_cond == "decay", ], acc.dat.cis[acc.dat.cis$rew_cond == "decay",], 
         xlab = " ", ylab = " ", xlims = c(0, .3), acc.ylims, mars = c(2, 2, 2, 1) + 0.1,
         dv = "acc", dat_a_value = "high", dat_b_value = "low", cue_validity = "invalid", legend_on = 0, 
         dot_size = 2, fit_line = 3, title = "", x_axis = 1, y_axis = 0)

par(mar = c(2, 3, 2, 1) + 0.1)
tmp.dat = c(acc.bf.for.plot$BF, 0)
tmp.names = c(as.character(acc.bf.for.plot$names), "")
tmp.lower = c( acc.bf.for.plot$lower, 0) 
tmp.upper = c( acc.bf.for.plot$upper, 0) 
bar.order = c(1, 2, 3, 4, 5, 8, 6, 7)
cols = c(wes_palette("Royal1")[4], wes_palette("Royal1")[4], wes_palette("Royal1")[4], wes_palette("Royal1")[4], 
         wes_palette("Royal1")[4], wes_palette("Royal1")[3], wes_palette("Royal1")[3], wes_palette("Royal1")[3])
barcenters = barplot(tmp.dat[bar.order], col = cols, 
                     names.arg = tmp.names[bar.order], ylab = "", las = 2, ylim = c(0, 5), 
                     axes = FALSE)
arrows(barcenters, tmp.lower[bar.order], barcenters, tmp.upper[bar.order], lwd = 1.5, angle = 90, code = 3, length = eb_length*0.5)
abline(h=3, lty = 1, lwd = 2)
abline(h=1, lty = 1, lwd = 2)
axis(side = 1, tck = 0, labels = FALSE)
axis(side = 2, tck = 0, las = 2, cex.axis = dot_size)
mtext(side = 2, "BF (IV + SC)", line = 2, cex =1.5)
# title("~ v + va + rc", line = -2, cex.main = dot_size)
dev.off()



######### CALL THIS FUNCTION TO PLOT DVS
plot.dvs <- function(xs, dat.df.plot, dat.cis, xlab, ylab, xlims, ylims, mars, dv, dat_a_value, dat_b_value, cue_validity, legend_on, leg.x, leg.y, legend_low, legend_high, dot_size, fit_line, title, x_axis, y_axis){
            # set margins for fig
            par(mar = mars) 
            # get dvs
            if (dv == "mean") {
                ylims = ylims
                dat_a = dat.df.plot$mean[dat.df.plot$value == dat_a_value & dat.df.plot$valid == cue_validity]
                dat_b = dat.df.plot$mean[dat.df.plot$value == dat_b_value & dat.df.plot$valid == cue_validity]
                ci_a = dat.cis$ci[ dat.cis$value == dat_a_value & dat.cis$valid == cue_validity]/2
                ci_b = dat.cis$ci[ dat.cis$value == dat_b_value & dat.cis$valid == cue_validity]/2
                fit_a = dat.df.plot$predict[dat.df.plot$value == dat_a_value & dat.df.plot$valid == cue_validity]
                fit_b = dat.df.plot$predict[dat.df.plot$value == dat_b_value & dat.df.plot$valid == cue_validity]
            } else if (dv == "acc"){
                ylims = ylims
                dat_a = dat.df.plot$acc[dat.df.plot$value == dat_a_value & dat.df.plot$valid == cue_validity]
                dat_b = dat.df.plot$acc[dat.df.plot$value == dat_b_value & dat.df.plot$valid == cue_validity]
                ci_a = dat.cis$ci[ dat.cis$value == dat_a_value & dat.cis$valid == cue_validity]/2
                ci_b = dat.cis$ci[ dat.cis$value == dat_b_value & dat.cis$valid == cue_validity]/2
                fit_a = dat.df.plot$predict[ dat.df.plot$value == dat_a_value &  dat.df.plot$valid == cue_validity ]
                fit_b = dat.df.plot$predict[ dat.df.plot$value == dat_b_value &  dat.df.plot$valid == cue_validity ]
            }
                
            if (dat_a_value == "high"){
              palette = c(wes_palette("Royal1")[2], wes_palette("Royal1")[1])
            } else {
              palette = c(wes_palette("Royal1")[1], wes_palette("Royal1")[2])
            }
                
            plot(xs, dat_a, pch = 19, cex = dot_size, col = palette[1], 
                 xlab = " ", ylab = " ", ylim = ylims, xlim = xlims, axes = FALSE,
                 cex.axis = dot_size)
            points( xs, dat_b, pch = 19, cex = dot_size, col = palette[2] )
            if (x_axis == 1) {
              axis(side = 1, tck = 0, cex.axis = dot_size)
              mtext(side = 1, xlab, line = 4, cex = dot_size)
            }
            if (y_axis == 1) {
              axis(side = 2, tck = 0, las = 2, cex.axis = dot_size)
              mtext(side = 2, ylab, line = 4, cex = dot_size)
            }
            arrows( x0= xs, y0= dat_a + ci_a, x1= xs, y1= dat_a - ci_a,
                    angle=90, code=3, length = eb_length, col = palette[1], lwd = dot_size)
            arrows( x0= xs, y0= dat_b + ci_b, x1= xs, y1 = dat_b - ci_b,
                    angle=90, code=3, length = eb_length, col = palette[2], lwd = dot_size)
            points(xs, fit_a, lwd = dot_size, type = "l", lty = fit_line, col = palette[1])
            points(xs, fit_b, lwd = dot_size, type = "l", lty = fit_line, col = palette[2])
            title(title, line = -2, cex.main = dot_size)
            if (legend_on == 1) {
              legend(leg.x, leg.y, c(legend_high, legend_low), pch = 19, lwd = c(3,3), 
                     col = c(wes_palette("Royal1")[2], wes_palette("Royal1")[1]), 
                     bty = "n", title = "IV", cex = dot_size)
            }
}

            












