setwd("E:/Documents/POMP/Parus Major/")
#library(reshape2)
library(pomp)
library(plyr)
library(repr)
library(plotly)

loc <- url("http://kingaa.github.io/short-course/intro/parus.csv")
dat <- read.csv(loc)

parus<-pomp(dat,time='year',t0=1959)
p <- plot_ly(dat,x=~year,y=~pop, type = "scatter",mode='lines+markers',colors = "Dark2")
p

stochStep <- Csnippet("
  e = rnorm(0,sigma);
  N = r*N*exp(-N+e);
                    ")

pomp(parus,rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
     paramnames=c("r","sigma"),statenames=c("N","e")) -> parus