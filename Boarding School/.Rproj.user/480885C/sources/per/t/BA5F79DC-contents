#### SEIR Stochastic Simulations

#### SEIR MODEL Inference with Demography, Seasonality and Immigration #####


#### Clear the environment
rm(list=ls())
### Libraries
library(reshape2)
library(repr)
library(ggplot2)
library(pomp)
library(plyr)
library(dplyr)
library(xlsx)
set.seed(594709947L)


### read.table("https://kingaa.github.io/sbied/stochsim/bsflu_data.txt") -> bsflu


data_file_csv=file.path(getwd(), "Data","Boarding.csv")
reported=read.table(data_file_csv,sep=",",header=TRUE)

rproc <- Csnippet("
                  double N = 764;
                  double Lambda = (Beta*I)/N;
                  double p1 = rbinom(S,1-exp(-Lambda*dt));
                  double p2 = rbinom(I,1-exp(-Gamma*dt));
                  double p3 = rbinom(C,1-exp(-Delta*dt));
                  
                  S  -= p1;
                  I  += p1 - p2;
                  C  += p2 - p3;
                  R  += N - S - I - C;
                  ")

init <- Csnippet("
                 S = 763;
                 I = 1;
                 C = 0;
                 R = 0;
                 ")

dmeas <- Csnippet("
           ///       lik = dpois(flu,C,give_log);
                     double tau=0.1;
                     double nbr=1/tau/tau;
                     lik = dnbinom(flu,nbr,1/(1+C/nbr),give_log);
                  ")

rmeas <- Csnippet("
          ///        flu = rpois(C);
                     double tau=0.1;
                     double nbr=1/tau/tau;
                     flu = rnbinom(nbr,1/(1+C/nbr));

                  ")

reported %>%
  pomp(time="day",t0=0,
       rprocess=euler.sim(rproc,delta.t=1/5),
       initializer=init,
       rmeasure=rmeas,
       dmeasure=dmeas,
       paramnames=c("Beta","Gamma","Delta"),
       statenames=c("S","I","C","R")
  ) -> boarding



params <- c(Beta=3.2,Gamma=1.5,Delta=0.5)
pf <- pfilter(boarding,params=params,Np=10000,save.states=TRUE)
print(logLik(pf))



