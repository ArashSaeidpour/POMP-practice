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
Iota     <- 1.5      ## mean number of visiting infectives at anytime
rho      <- 0.7      ## reporting probability
tol      <- 1e-18    ## reporting variance


#### Endemic equilbrium  
Ro       <- (Beta0*Sigma)/((Mu+Gamma)*(Mu+Sigma))
S_eq     <- 1/Ro
E_eq     <- (Mu*(Mu+Gamma)*(Ro-1))/(Beta0*Sigma)
I_eq     <- (Mu/Beta0)*(Ro-1)

#### Total population
N0      <- 1e6      ## Total Population


#### Initial condition from endemic equil.
S0      <- as.integer(S_eq * N0)       ## inital susceptible population
E0      <- as.integer(E_eq * N0)       ## inital exposed population
I0      <- as.integer(I_eq * N0)       ## inital infectious population
R0      <- N0-(S0+E0+I0)

#### Time and setep size parameters
dt       <- 1/365.25    ## Time step size (in years)
n_years  <- 65        ## analysis period (years) 
t_i      <- -1/365.25   ## Start time




####### Deterministic Skeleton
SEIR.ode <- Csnippet("
                     double Beta = beta0*(1+beta1*sin(2*M_PI*t));
                     DS = nearbyint(mu*N  - Beta*S*(I+iota)/N - nu*S);
                     DE = nearbyint(Beta*S*(I+iota)/N - sigma*E - nu*E);
                     DI = nearbyint(sigma*E - gamma*I - nu*I);
                     DR = nearbyint(gamma*I-nu*R);
                     DC = nearbyint(DR);
                     
                     ")


### Event increment S,E,I,R,N,C

###### Initial condition
init <- Csnippet("
                 S = S_0;
                 E = E_0;
                 I = I_0;
                 R = N - S_0 - E_0 - I_0;
                 C=0;
                 ")


###### Pomp object
pomp(data=data.frame(time=seq(0,n_years,by=dt),data=NA),
     times="time",t0=t_i,
     skeleton=vectorfield(SEIR.ode),
     initializer=init,
     statenames=c("S","E","I","R","C"),
     zeronames="C",
     paramnames=c("beta0","beta1","gamma","mu","nu","sigma","iota","N","S_0","E_0","I_0")
) -> SEIR.model


#### Assigning parameters
params1 <- c(beta0=Beta0,beta1=Beta1,gamma=Gamma,mu=Mu,nu=Nu,sigma=Sigma,iota=Iota,N=N0,S_0=S0,E_0=E0,I_0=I0)

###### Deterministic Trajectory
trac <- trajectory(SEIR.model,params=params1,as.data.frame=TRUE)


df <- data.frame(trac[c('time','I','C')])

##### Plotting Case counts (C) vs I
df <- melt(df ,  id.vars = 'time')
fig<-ggplot(df)+geom_point(aes(time, value, group=time, colour=variable))
print(fig)




#df2 <- data.frame(trac[c('time','S','E','I')])
#df2 <- melt(df2 ,  id.vars = 'time')
#fig3 <-ggplot(df2)+geom_point(aes(time, value, group=time, colour=variable))
#print(fig3)