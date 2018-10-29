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



data_file_txt=file.path(getwd(), "Data","incidence.txt")
reported=read.table(data_file_txt)



##### Fixed Parameters
mu       <- 1/70     ## Birth rate (years)
nu       <- 1/70     ## Death rate (years)
beta0    <- 450      ## B0 for seasonality
beta1    <- 0.05     ## B1 for seasonality
sigma    <- 365/8    ## Average latent period
gamma    <- 365/14   ## Average recovery period
iota     <- 1.5      ## mean number of visiting infectives at anytime
rho      <- 0.15      ## reporting probability
#tau      <- 8.0      ## reporting overdispersion


#### Endemic equilbrium  
Ro       <- (beta0*sigma)/((mu+gamma)*(mu+sigma))
S_eq     <- 1/Ro
E_eq     <- (mu*(mu+gamma)*(Ro-1))/(beta0*sigma)
I_eq     <- (mu/beta0)*(Ro-1)

#### Total population
N_0      <- 1e6      ## Total Population


#### Initial condition from endemic equil.
S_0      <- as.integer(S_eq * N_0)       ## inital susceptible population
E_0      <- as.integer(E_eq * N_0)       ## inital exposed population
I_0      <- as.integer(I_eq * N_0)       ## inital infectious population
R_0      <- N_0-(S_0+E_0+I_0)

#### Time and setep size parameters
dt       <- 1/365.25    ## Time step size (in years)
n_years  <- 65        ## analysis period (years) 
t_i      <- -1/365.25   ## Start time



rproc <- Csnippet("
                  double Beta,POP,dPOP,r_birth;
                  double rate[6],trans[6];
                  const double pi=3.1415926535897932384626;
                  
                  POP  =  S+E+I+R;
                  Beta =  Beta0*(1+Beta1*sin(2*pi*t));
                  
                  
                  /// S rates
                  r_birth=Mu*POP;                   // Birth
                  rate[0]=(Beta*(I+Iota))/POP;      // Force of infection
                  rate[1]=Nu;                       // Death S
                  
                  /// E rates
                  rate[2]=Sigma;                    // Infection
                  rate[3]=Nu;                       // Death E
                  
                  /// I rates
                  rate[4]=Gamma;                    // Recovery
                  rate[5]=Nu;                       // Death I
                  
                  
                  // transitions between classes
                  dPOP = rpois(r_birth * dt);
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,E,&rate[2],dt,&trans[2]);
                  reulermultinom(2,I,&rate[4],dt,&trans[4]);
                  
                  S  += dPOP - trans[0] - trans[1];                       ///  DS/DT = mu*N  - Beta*S*(I+i)/N - nu*S
                  E  += trans[0] - trans[2] - trans[3];                   ///  DE/DT = Beta*S*(I+i)/N - sigma*E - nu*E
                  I  += trans[2] - trans[4] - trans[5];                   ///  DI/DT = sigma*E - gamma*I - nu*I 
                  R  = POP - S - E - I;
                  C  = trans[4];                                          ///  True Incidenc
")




########### Initial condition
initlz <- Csnippet("
                   S  = S0;
                   E  = E0;
                   I  = I0;
                   R  = R0;
                   C  = 0;
                   ")



## ----dmeasure------------------------------------------------------------
dmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
")




toEst <- Csnippet("
                  TBeta0  = log(Beta0);
                  TBeta1  = log(Beta1);
                  TGamma  = log(Gamma);
                  TSigma  = log(Sigma);
                  TIota   = log(Iota);
                  to_log_barycentric (&TS0, &S0, 4);
                  ")

fromEst <- Csnippet("
                    TBeta0  = exp(Beta0);
                    TBeta1  = exp(Beta1);
                    TGamma  = exp(Gamma);
                    TSigma  = exp(Sigma);
                    TIota   = exp(Iota);
                    from_log_barycentric (&TS0, &S0, 4);
                    ")


reported %>%
pomp(t0=with(reported,days[1]-days[2]),
     times="time",
     rprocess=euler.sim(rproc,delta.t=dt),
     rmeasure=rmeas,
     toEstimationScale=toEst,
     fromEstimationScale=fromEst,
     initializer=initlz,
     zeronames=c("C"),
     statenames=c("S","E","I","R","C"),
     paramnames=c("Mu","Nu","Sigma","Gamma","Beta0","Beta1","Iota","Rho",
                   "S0","E0","I0","R0")
) -> SEIR.StochasticEst


theta <- c(mu,nu)
names(theta) <- c("Mu","Nu")

estpars <- setdiff(names(theta),c("Mu","Nu"))






