############################################## FUNCTIONS ##############################################
############ creating a function to trim rts - which will use tapply
trim <- function(data){
  mu <- mean(data$rt[data$cor_resp == data$resp],na.rm=TRUE)
  sigma <- sd(data$rt[data$cor_resp == data$resp],na.rm=TRUE)
  upper <- mu + (2.5*sigma)
  lower <- mu - (2.5*sigma)
  data$rt[(data$rt[data$cor_resp == data$resp] < lower)] <- NA
  data$rt[(data$rt[data$cor_resp == data$resp] > upper)] <- NA
  
  return(data)
}

########## and scale
scale.rts <- function(data){
  min <- min(data$rt,na.rm=TRUE)
  max <- max(data$rt,na.rm=TRUE)
  data$rt <- (data$rt - min)/(max - min)
  data
}

