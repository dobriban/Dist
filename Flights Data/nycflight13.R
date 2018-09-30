

writePath <-
  "C:/Dropbox/Projects/Distributed/Experiments/Flights Data"
setwd(writePath)

#load packages
if (!require("pacman"))
  install.packages("pacman")

pacman::p_load(dplyr,
               nycflights13,
               ggplot2,
               pracma,
               caret,
               corpcor,
               install = TRUE,
               update = FALSE)

#toBibtex(citation("nycflights13"))

#create data frame
dat = flights # omit the missing values
# join the separate datasets: weather, planes, airlines
dat = left_join(dat, weather)
dat = left_join(dat, planes, by = "tailnum")
dat = left_join(dat, airlines, by = "carrier")
# the following are the variables with real values
variate = c(
  "arr_delay",
  "month",
  "day",
  "dep_time",
  "sched_dep_time",
  "dep_delay",
  "arr_time",
  "sched_arr_time",
  "air_time",
  "distance",
  "hour",
  "temp",
  "dewp",
  "humid",
  "wind_dir",
  "wind_speed",
  "wind_gust",
  "precip",
  "pressure",
  "visib",
  "year.y",
  "seats"
)
dat = dat[, variate]
dat = na.omit(dat)
p = dim(dat)[2]

#remove highly correlated variables
sc = cor(dat)
thr = 0.8
largec = abs(sc) > 0.8
for (i in 1:(p - 1)) {
  for (j in (i + 1):p) {
    if (largec[i, j] == 1)
    {
      print(c(i, j))
    }
  }
}
dat.full = dat
vr = c(4, 5, 6, 9, 12, 16) #vars to remove
dat = dat[, -vr]


m = dim(dat)[1]  # total number of rows
p = dim(dat)[2] - 1
#int = ones(m,1)
#dat = cbind(dat, int)
set.seed(4280)

n.arr = c(100,500,1000,3000,10000,20000,30000)
for (n in n.arr) {
  #n = 1000
  #p = 21
  idx.all = sample(m, 2 * n)
  idx = idx.all[1:n]
  idx.test = idx.all[(n + 1):(2 * n)]
  Y = as.matrix(dat[idx, 1]) # the first column is the response
  X = as.matrix(dat[idx, 2:(p + 1)]) # the 2-22 columns are the features
  Yt = as.matrix(dat[idx.test, 1]) # the first column is the response
  Xt = as.matrix(dat[idx.test, 2:(p + 1)]) # the 2-22 columns are the features
  
  
  #Compute OLS
  P  = pseudoinverse(X)
  #beta.full = solve(t(X) %*% X) %*% t(X) %*% Y
  beta.full = P %*% Y
  err = Yt - Xt %*% beta.full
  Of = norm(err) ^ 2 / n
  
  #Compute distributed
  k.max = floor(n / (p))-1
  #k.max = 10
  Od = array(0, dim = c(k.max, 1))
  k.arr = (1:k.max)
  for (k in k.arr) {
    RSS.d.track = rep(0, k)
    MSE.d = rep(0, k)
    beta.arr = array(0, dim = c(p, k))
    beta.d = array(0, dim = c(p, 1))
    kf = createFolds(Y, k)
    i = 1
    for (index in kf) {
      X.d = X[index,]
      Y.d = Y[index]
      P.d  = pseudoinverse(X.d)
      beta.sub = P.d %*% Y.d
      #beta.sub = solve(t(X.d) %*% X.d) %*% t(X.d) %*% Y.d
      beta.arr[, i] = beta.sub
      #MSE.d[i] = sum(diag(solve(t(X.d) %*% X.d)))
      MSE.d[i] = sum(diag(P.d %*% t(P.d)))
      beta.d = beta.d + beta.sub / MSE.d[i]
      i = i + 1
    }
    
    IV = sum(1 / MSE.d)
    beta.d = beta.d / IV
    Od[k] = norm(Yt - Xt %*% beta.d) ^ 2 / n
    
  }
  gamma = p / n
  OE = Of / Od
  
  #OE.theo = 1/(1+(k.arr-1)*gamma^2/(1-k.arr*gamma))
  OE.theo = 1 / (1 - gamma) * 1 / (1 + 1 / (1 / gamma - k.arr))
  
  
  #plot
  fin = paste("oe_flights_k=", max(k.arr), "_n=", n, ".png", sep = "")
  png(file = fin,
      width = 800,
      height = 600)
  
  par(cex = "2")
  plot(OE, xlab = "k", ylab = "OE", 
       ylim = c(0,max(OE, OE.theo) + 0.1),
       cex.axis = 1.5, cex.lab = 1.5, pch = 15,
       main=paste("n = ",  n)
       )
  points(OE.theo, col = 'red', pch = 16)
  legend(
    0.5,
    legend = c("Empirical", "Theory"),
    col = c("black", "red"),
    lty = 1:1,
    lwd = 4,
    cex = 1.8,
    bty = "n"
  )
  
  #dev.copy(png,fin)
  dev.off()
}
