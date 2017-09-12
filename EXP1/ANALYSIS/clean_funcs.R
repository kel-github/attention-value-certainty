# block switch cleaning__________________________________________________________________
clean.out.block.switches <- function(data,ntrials){
  
  # clean data
  data$rnum = seq(1, length(data$sub), 1)
  cdata = data[data$cor_resp == data$resp, ]
  tmp <- by(data.frame(rt=cdata$RT, sub=cdata$sub, fixprob=cdata$fixprob,
                       loc_prob=cdata$loc_prob, value=cdata$value,
                       rnum=cdata$rnum),
            list(cdata$loc_prob, cdata$value), trim)
  cdata = as.data.frame(do.call(rbind, tmp))
  cdata <- cdata[order(cdata$rnum), ]
  rm(tmp)
  
  
  ntrials = ntrials
  # step 1 - mark first 30 trials after switch
  trials = length(cdata$rt)
  mark.trials = array(data=0, dim=trials)
  for (x in 2:trials){
    if (cdata$fixprob[x-1] != cdata$fixprob[x]) mark.trials[c(x:(x+ntrials-1))] = 1
  }
  cdata$mark = mark.trials
  cdata$mark <- as.factor(cdata$mark)
  
  rt = cdata$rt[cdata$mark == "1"]
  n_bs = length(rt)/ntrials
  block = rep(c(1:ntrials), times = n_bs)
  block = as.factor(block)
  tmp = as.data.frame(cbind(rt, block))
  tmp2 = ddply(tmp, .(block), summarise,
               mean = mean(rt, na.rm = TRUE))
  tmp = tmp2; rm(tmp2)
  # now have data for regression, run piecewise models along x, 
  # moving break point along x
  x = as.numeric(tmp$block)
  breaks = x[which(x >= 1 & x <= 25)] #
  mse <- numeric(length(unique(breaks))) # to hold MSE
  for(i in 1:length(breaks)){
    tmp.mod <- lm(tmp$mean ~ x*(x < breaks[i]) + x*(x>=breaks[i]))
    mse[i] <- summary(tmp.mod)[6]
  }
  mse <- as.numeric(mse)  
  # get (1st) lowest breakpoint error
  brk_pnt = breaks[which(mse == min(mse))][1]
  # run final model
  final.mod <- lm(tmp$mean ~ x*(x<brk_pnt) + x*(x>brk_pnt))
  summary(final.mod)
  
  # plot participant data and save breakpoint
  plot(x, tmp$mean, pch = 19, col = wesanderson::wes_palette("FantasticFox")[5],
       ylab="mean RT", xlab="bin")
  curve((final.mod$coefficients[1] + final.mod$coefficients[3]) + (final.mod$coefficients[2] + final.mod$coefficients[5])*x, 
        add=T, from=1, to=brk_pnt)
  curve((final.mod$coefficients[1] + final.mod$coefficients[4]) + final.mod$coefficients[2]*x, 
        add = T, from=brk_pnt, to=max(x))
  abline(v = brk_pnt, lty = 3)
  
  # cut out data from pre-breakpoint
  row = brk_pnt
  return(row)
}
