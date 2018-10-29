#### SEIR MODEL with Demography, Seasonality and Immigration #####


#### Clear the environment
rm(list=ls())
### Libraries
library(reshape2)
library(repr)
library(ggplot2)
library(pomp)
library(plyr)
library(dplyr)


##### Fixed Parameters
Mu       <- 1/70     ## Birth rate (years)
Nu       <- 1/70     ## Death rate (years)
Beta0    <- 450      ## B0 for seasonality
Beta1    <- 0.05     ## B1 for seasonality
Sigma    <- 365/8    ## Average latent period
Gamma    <- 365/14   ## Average recovery period
immig    <- 5        ## Immigration rate (per year)
S0       <- 100000    ## inital susceptible population
E0       <- 40       ## inital exposed population
I0       <- 100       ## inital infectious population
H0       <- 0        ## Case coounts: integral(Gamma*I)
N0       <- 1e6      ## Inital population

dt       <- 7/365    ## Time step size (in years)
n_years  <- 60       ## analysis period (years) 
t_i       <- -7/365   ## Start time

R0 <- N0-S0-E0-I0
####### Deterministic Skeleton
SEIR.ode <- Csnippet("
  double Beta = beta0*(1+beta1*sin(2*M_PI*t));
  DS = mu*N + i - Beta*S*I/N - nu*S;
  DE = Beta*S*I/N - sigma*E - nu*E;
  DI = sigma*E - gamma*I - nu*I;
  DR = gamma*I-nu*R;
")


### Event increment S,E,I,R,N,C

###### Initial condition
init <- Csnippet("
  S = S_0;
  E = E_0;
  I = I_0;
  R = N - S_0 - E_0 - I_0;
")


###### Pomp object
pomp(data=data.frame(time=seq(0,n_years,by=dt),data=NA),
     times="time",t0=t_i,
     skeleton=vectorfield(SEIR.ode),
     initializer=init,
     statenames=c("S","E","I","R"),
     paramnames=c("beta0","beta1","gamma","mu","nu","sigma","i","N","S_0","E_0","I_0")
) -> SEIR.model


#### Assigning parameters
params1 <- c(beta0=Beta0,beta1=Beta1,gamma=Gamma,mu=Mu,nu=Nu,sigma=Sigma,i=immig,N=N0,S_0=S0,E_0=E0,I_0=I0)

###### Deterministic Trajectory
trac <- trajectory(SEIR.model,params=params1,as.data.frame=TRUE)


#### Obtaining Cae counts (C) from I
C <- data.frame(Gamma*trac['I']*(7/365)) ### Case counts
colnames(C) <- "C"
df <- data.frame(trac[c('time','I')],C)

##### Plotting Case counts (C) vs I
df <- melt(df ,  id.vars = 'time')
fig<-ggplot(df)+geom_point(aes(time, value, group=time, colour=variable))
print(fig)




df2 <- data.frame(trac[c('time','S','E','I')])
df2 <- melt(df2 ,  id.vars = 'time')
fig3 <-ggplot(df2)+geom_point(aes(time, value, group=time, colour=variable))
print(fig3)



init2 <- Csnippet("
  S = S_0;
  E = E_0;
  I = I_0;
  R = N_0 - S_0 - E_0 - I_0;
  N=  N_0;
  H = H_0;
")


dt=1/365
n_year=30
pomp(data=data.frame(time=seq(0,n_years,by=dt),data=NA),
     times="time",t0=t_i,
     #skeleton=vectorfield(SEIR.ode),
     rprocess=gillespie.sim(
       v=cbind(
         birth=c(1,0,0,0,1,0),
         immigration=c(1,0,0,0,1,0),
         sdeath=c(-1,0,0,0,-1,0),
         exposure=c(-1,1,0,0,0,0),
         infection=c(0,-1,1,0,0,0),
         edeath=c(0,-1,0,0,-1,0),
         ideath=c(0,0,-1,0,-1,0),
         recovery=c(0,0,-1,1,0,1),
         rdeath=c(0,0,0,-1,-1,0)
       ),
       rate.fun=function(j,x,t,params,...) {
         POP <- sum(x[c("S","E","I","R")])
         BETA <- params["beta0"]*(1+params["beta1"]*sin(2*pi*t))
         switch(j,
                params["mu"]*POP,         ## Birth
                params["i"],              ## Immigration 
                params["nu"]*x["S"],      ##Sdeath
                BETA*x["S"]*x["I"]/POP,   ## Exposure
                params["sigma"]*x["E"],   ## Infection
                params["nu"]*x["E"],      ##Edeath
                params["nu"]*x["I"],      ##Ideath
                params["gamma"]*x["I"],   ##Recovery
                params["nu"]*x["R"]       ##Rdeath
         )
       }
     ),
     #rmeasure=Csnippet("reports=rbinom(H,rho);"),
     paramnames=c("beta0","beta1","gamma","mu","nu","sigma","i","S_0","E_0","I_0","R_0","N_0","H_0"),
     statenames=c("S","E","I","R","N","H"),
     zeronames="H",
     initializer=init2
) -> SEIR.Stochstic


params2 <- c(beta0=Beta0,beta1=Beta1,gamma=Gamma,mu=Mu,nu=Nu,sigma=Sigma,i=immig,S_0=S0,E_0=E0,I_0=I0,R_0=R0,N_0=N0,H_0=H0)


#pomp(SEIR.Stochstic,initializer=init2) -> SEIR.Stochstic


simulate(SEIR.Stochstic,params=params2,nsim=5,states=TRUE,obs=FALSE,as.data.frame=TRUE) -> sims
fig2<-ggplot(sims)+geom_point(aes(time, I, group=time, colour=sim))
print(fig2)

### Figure out why simulated infections does not make a return, like the trajectory??

