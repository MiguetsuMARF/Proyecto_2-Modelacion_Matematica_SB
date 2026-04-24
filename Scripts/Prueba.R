source("Scripts/Grind.R")

library(deSolve)
par(mfrow = c(1,1))
model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
    dx = alfa+TGFb*(x^n/(BAMBI^n + x^n))- x*(beta + (mu1*y))
    dy = gamma*x - y*((mu2*SnoN) + delta)
    return(list(c(dx, dy)))
  })
}
p <-  c(alfa = 3, TGFb = 10, BAMBI = 0.2, 
        beta = 1, mu1 = 0.5, n = 6, gamma = 1, mu2 = 0.4, SnoN = 1, delta = 1) 

s <- c(x=0,y=0)

run(tmax = 10, tstep = 0.01, state = s, parms = p, odes = model)


run(tmax = 200, tstep = 0.1, state = s, parms = p, odes = model)

plane(xmax=10, ymax = 10)
mid <- newton(s,plot=T)
low <- newton(c(x=10,y=0),plot=T)
hig <- newton(c(x=0,y=20),plot=T)

# parametro n
continue(state=low, parms = p, odes=model, x="n", step=0.001, 
         xmin=0,xmax=10,y="x", ymin=0, ymax=25)
continue(state=mid,parms = p,odes=model, x="n",  step=0.001, 
         xmin=0, xmax=50,y="x", ymin=0, ymax=25, add = TRUE)
continue(state=hig, parms = p, odes=model, x="n", step=0.001, 
         xmin=0, xmax=50,y="x", ymin=0, ymax=25, add = TRUE)

continue(state=low, parms = p, odes=model, x="n", step=0.001, 
         xmin=0,xmax=10,y="y", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="n",  step=0.001, 
         xmin=0, xmax=10,y="y", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="n", step=0.001, 
         xmin=0, xmax=10,y="y", ymin=0, ymax=50, add = TRUE)

# parametro alfa
continue(state=low, parms = p, odes=model, x="alfa", step=0.001, 
         xmin=0,xmax=100,y="y", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="alfa",  step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="alfa", step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)

continue(state=low, parms = p, odes=model, x="alfa", step=0.001, 
         xmin=0,xmax=100,y="x", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="alfa",  step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="alfa", step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)

# parametro TGFb
continue(state=low, parms = p, odes=model, x="TGFb", step=0.001, 
         xmin=0,xmax=100,y="y", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="TGFb",  step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="TGFb", step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)

continue(state=low, parms = p, odes=model, x="TGFb", step=0.001, 
         xmin=0,xmax=100,y="x", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="TGFb",  step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="TGFb", step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)

# parametro BAMBI
continue(state=low, parms = p, odes=model, x="BAMBI", step=0.001, 
         xmin=0,xmax=10,y="y", ymin=0, ymax=6)
continue(state=mid,parms = p,odes=model, x="BAMBI",  step=0.001, 
         xmin=0, xmax=10,y="y", ymin=0, ymax=6, add = TRUE)
continue(state=hig, parms = p, odes=model, x="BAMBI", step=0.001, 
         xmin=0, xmax=10,y="y", ymin=0, ymax=6, add = TRUE)

continue(state=low, parms = p, odes=model, x="BAMBI", step=0.001, 
         xmin=0,xmax=10,y="x", ymin=0, ymax=6)
continue(state=mid,parms = p,odes=model, x="BAMBI",  step=0.001, 
         xmin=0, xmax=10,y="x", ymin=0, ymax=6, add = TRUE)
continue(state=hig, parms = p, odes=model, x="BAMBI", step=0.001, 
         xmin=0, xmax=10,y="x", ymin=0, ymax=6, add = TRUE)

p[3] <- 3.3
par(mfrow = c(1,1))
plane(xmax=10, ymax = 4)
mid <- newton(s,plot=T)
low <- newton(c(x=2,y=1),plot=T)
hig <- newton(c(x=5,y=4),plot=T)

ini_1 <- c(x=2,y=1)
ini_2 <- c(x=5,y=4)

par(mfrow = c(3,2))
for (BAMBI2 in c(0.1, 3.3, 10)) { 
  p[3] <-  BAMBI2
  out1 <- ode(y = ini_1, times = seq(from = 0, to = 10, by = 0.01), func = model, parms = p)
  
  plot(out1[,1], out1[,2],type = "l", ylim=c(0,5),
       col="red", xlab = "Time", ylab = "X(t)", 
       main = paste("BAMBI=", toString(BAMBI2), sep=" "))
  
  out2 <- ode(y = ini_2, times = seq(from = 0, to = 10, by = 0.01), func = model, parms = p)
  lines(out2[,1], out2[,2],type = "l", col="blue")
  legend=c(paste("I.C.=", toString(ini_1), sep=" "), paste("I.C.=", toString(ini_2), sep=" "))
  
  poli.flowField <- flowField(
    model, xlim = c(0, 10), ylim = c(0, 10), 
    parameters = p, points = 20, add = FALSE)
  poli.trajectory <- trajectory(
    model, y0 = ini_1, tlim = c(0,20), 
    parameters = p, col = "blue")
  poli.trajectory <- trajectory(
    model, y0 = ini_2, tlim = c(0,20), 
    parameters = p, col = "red")
}

par(mfrow = c(1,3))
for (BAMBI2 in c(0.1, 3.3, 10)){
  p[3] <-  BAMBI2
  flowField(model, xlim = c(0, 10), ylim = c(0, 10), parameters = p, points = 20, add = FALSE)
  s <- c(x=0,y=0)
  mid <- newton(s,plot=T)
  low <- newton(c(x=2,y=1),plot=T)
  hig <- newton(c(x=5,y=4),plot=T)
  
  for (ii in seq(1,20,1) ){
    r1=runif(1); r2=runif(1); r3=runif(1);
    
    if (r1<0.5){
      ini = c(as.numeric(r2<0.5) + rnorm(1,sample(1:7)), r3 + rnorm(1,sample(1:7)))
    } else {
      ini = c(r3, as.numeric(r2<0.5) + 0.1^r1)
    }
    
    trajectory(model, y0 = ini, tlim = c(0,20), parameters = p, col = "blue")
  }
}

par(mfrow = c(1,1))
p[3] <-  0.2

# Parametro beta
continue(state=low, parms = p, odes=model, x="beta", step=0.001, 
         xmin=0,xmax=50,y="y", ymin=0, ymax=4)
continue(state=mid,parms = p,odes=model, x="beta",  step=0.001, 
         xmin=0, xmax=50,y="y", ymin=0, ymax=4, add = TRUE)
continue(state=hig, parms = p, odes=model, x="beta", step=0.001, 
         xmin=0, xmax=50,y="y", ymin=0, ymax=4, add = TRUE)

continue(state=low, parms = p, odes=model, x="beta", step=0.001, 
         xmin=0,xmax=50,y="x", ymin=0, ymax=5)
continue(state=mid,parms = p,odes=model, x="beta",  step=0.001, 
         xmin=0, xmax=50,y="x", ymin=0, ymax=5, add = TRUE)
continue(state=hig, parms = p, odes=model, x="beta", step=0.001, 
         xmin=0, xmax=50,y="x", ymin=0, ymax=5, add = TRUE)

p[4] <-  35
par(mfrow = c(1,1))
plane(xmax=0.5, ymax = 20)
mid <- newton(s,plot=T)
low <- newton(c(x=0.2,y=1),plot=T)
hig <- newton(c(x=0.5,y=4),plot=T)



ini_1 <- c(x=0,y=0)
ini_2 <- c(x=0.5,y=4)

par(mfrow = c(3,2))
for (i in 1:3) {
  a <- c(5,0.5,0.5)
  b<- c(0.1, 35, 50)
  beta2 <- b[i]
  out1 <- ode(y = ini_1, times = seq(from = 0, to = 10, by = 0.01), func = model, parms = p)
  
  plot(out1[,1], out1[,2],type = "l", ylim=c(0,a[i]),
       col="red", xlab = "Time", ylab = "X(t)", 
       main = paste("beta=", toString(beta2), sep=" "))
  
  out2 <- ode(y = ini_2, times = seq(from = 0, to = 10, by = 0.01), func = model, parms = p)
  lines(out2[,1], out2[,2],type = "l", col="blue")
  legend=c(paste("I.C.=", toString(ini_1), sep=" "), paste("I.C.=", toString(ini_2), sep=" "))
  
  poli.flowField <- flowField(
    model, xlim = c(0, a[i]), ylim = c(0, 10), 
    parameters = p, points = 20, add = FALSE)
  poli.trajectory <- trajectory(
    model, y0 = ini_1, tlim = c(0,20), 
    parameters = p, col = "blue")
  poli.trajectory <- trajectory(
    model, y0 = ini_2, tlim = c(0,20), 
    parameters = p, col = "red")
}

par(mfrow = c(1,3))
for (i in 1:3) {
  a <- c(10,0.5,0.5)
  b<- c(0.1, 35, 50)
  beta2 <- b[i]
  p[4] <-  beta2
  flowField(model, xlim = c(0, a[i]), ylim = c(0, 10), parameters = p, points = 20, add = FALSE)
  s <- c(x=0,y=0)
  mid <- newton(s,plot=T)
  low <- newton(c(x=0.2,y=1),plot=T)
  hig <- newton(c(x=0.5,y=4),plot=T)
  
  for (ii in seq(1,20,1) ){
    r1=runif(1); r2=runif(1); r3=runif(1);
    
    if (r1<0.5){
      ini = c(as.numeric(r2<0.5) + rnorm(1,sample(1:7)), r3 + rnorm(1,sample(1:7)))
    } else {
      ini = c(r3, as.numeric(r2<0.5) + 0.1^r1)
    }
    
    trajectory(model, y0 = ini, tlim = c(0,20), parameters = p, col = "blue")
  }
}

par(mfrow = c(1,1))
p[4] <-  1


continue(state=low, parms = p, odes=model, x="mu1", step=0.001, 
         xmin=0,xmax=100,y="y", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="mu1",  step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="mu1", step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)

continue(state=low, parms = p, odes=model, x="mu1", step=0.001, 
         xmin=0,xmax=100,y="x", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="mu1",  step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="mu1", step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)


continue(state=low, parms = p, odes=model, x="gamma", step=0.001, 
         xmin=0,xmax=100,y="y", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="gamma",  step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="gamma", step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)

continue(state=low, parms = p, odes=model, x="gamma", step=0.001, 
         xmin=0,xmax=100,y="x", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="gamma",  step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="gamma", step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)


continue(state=low, parms = p, odes=model, x="mu2", step=0.001, 
         xmin=0,xmax=100,y="y", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="mu2",  step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="mu2", step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)

continue(state=low, parms = p, odes=model, x="mu2", step=0.001, 
         xmin=0,xmax=100,y="x", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="mu2",  step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="mu2", step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)


continue(state=low, parms = p, odes=model, x="SnoN", step=0.001, 
         xmin=0,xmax=100,y="y", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="SnoN",  step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="SnoN", step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)

continue(state=low, parms = p, odes=model, x="SnoN", step=0.001, 
         xmin=0,xmax=100,y="x", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="SnoN",  step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="SnoN", step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)


continue(state=low, parms = p, odes=model, x="delta", step=0.001, 
         xmin=0,xmax=100,y="y", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="delta",  step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="delta", step=0.001, 
         xmin=0, xmax=100,y="y", ymin=0, ymax=50, add = TRUE)

continue(state=low, parms = p, odes=model, x="delta", step=0.001, 
         xmin=0,xmax=100,y="x", ymin=0, ymax=50)
continue(state=mid,parms = p,odes=model, x="delta",  step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)
continue(state=hig, parms = p, odes=model, x="delta", step=0.001, 
         xmin=0, xmax=100,y="x", ymin=0, ymax=50, add = TRUE)




