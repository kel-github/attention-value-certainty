rm(list=ls())
library(wesanderson)
time = seq(0,1, .01)
decay_rate = -4
high_decay = 500*exp(decay_rate*time)
low_decay = 100*exp(decay_rate*time)

bmp(filename = "decay_manip.bmp",
    width = 5, height = 5, units = "cm", res=256)
par(ps = 8, cex = 1, cex.axis = 1, cex.lab = 1, mar=c(3, 3, 1, 1))
plot(time, high_decay, ylim=c(0, 500), xlab = "Time", ylab = "Points", type="l", col=wes_palette("Royal1")[2], lwd = 1.5, bty="n",
     yaxt="n", xaxt="n")
points(time, low_decay, type="l", col=wes_palette("Royal1")[1], lwd = 1.5, bty="n")
legend(0.2, 500, c("5000", "100"), bty="n", lty=1, col=c(wes_palette("Royal1")[2], wes_palette("Royal1")[1]), cex=1)
title(ylab="Points", line=2)
title(xlab="Time (s)", line=2)
axis(1, col = 1, col.ticks = NA)
axis(2, col = 1, col.ticks= NA)
dev.off()

