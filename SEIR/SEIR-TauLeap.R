#### SEIR deterministic trajectory and stochastic simulations


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
library(xlsx)
set.seed(594709947L)


##### Fixed Parameters
mu       <- 1/70     ## Birth rate (years)
nu       <- 1/70     ## Death rate (years)
beta0    <- 450      ## B0 for seasonality
beta1    <- 0.05     ## B1 for seasonality
sigma    <- 365/8    ## Average latent period
gamma    <- 365/14   ## Average recovery period
iota     <- 1.5      ## mean number of visiting infectives at anytime
rho      <- 0.15      ## reporting probability
tau      <- 8.0      ## reporting overdispersion


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
  C  = trans[4];                                          ///  True Incidence
  CR = rbinom(C,Rho);                                    /// Reported Incidence

  

  if (CR > 0.0) {
                  CR = nearbyint(CR);
                  } else {
                  CR = 0.0;
                  }
                                           
")

########### Initial condition
initlz <- Csnippet("
  S  = S0;
  E  = E0;
  I  = I0;
  R  = R0;
  C  = 0;
  CR = 0;
")



pomp(
     data=data.frame(time=seq(0,n_years,by=dt),data=NA),
    
     times="time",t0=t_i,
     rprocess=euler.sim(rproc,delta.t=dt),
#     rmeasure=rmeas,
     initializer=initlz,
     zeronames=c("C","CR"),
     statenames=c("S","E","I","R","C","CR"),
     paramnames=c("Mu","Nu","Sigma","Gamma","Beta0","Beta1","Iota","Rho","Tau",
                  "S0","E0","I0","R0")
) -> SEIR.StochasticSim

params <- c(Beta0=beta0,Beta1=beta1,Gamma=gamma,Mu=mu,Nu=nu,Sigma=sigma,Iota=iota,S0=S_0,E0=E_0,I0=I_0,R0=R_0,Rho=rho,Tau=tau)
simulate(params=params,SEIR.StochasticSim,nsim=10,states=TRUE,obs=FALSE,as.data.frame=TRUE) -> sims


########### Deterministic trajectory
################
#########################
####################################
####### Deterministic Skeleton
SEIR.ode <- Csnippet("
                     double Beta = Beta0*(1+Beta1*sin(2*M_PI*t));
                     double POP  =  S+E+I+R;
                     DS = nearbyint(Mu*POP  - Beta*S*(I+Iota)/POP - Nu*S);
                     DE = nearbyint(Beta*S*(I+Iota)/POP - Sigma*E - Nu*E);
                     DI = nearbyint(Sigma*E - Gamma*I - Nu*I);
                     DR = nearbyint(Gamma*I-Nu*R);
                     DC = nearbyint(Gamma*I);
                     ")


### Event increment S,E,I,R,C

###### Initial condition
init <- Csnippet("
                 S = S0;
                 E = E0;
                 I = I0;
                 R = R0;
                 C=0;
                 ")


###### Pomp object
pomp(data=data.frame(time=seq(0,n_years,by=dt),data=NA),
     times="time",t0=t_i,
     skeleton=vectorfield(SEIR.ode),
     initializer=init,
     statenames=c("S","E","I","R","C"),
     zeronames="C",
     paramnames=c("Beta0","Beta1","Gamma","Mu","Nu","Sigma","Iota","S0","E0","I0","R0")
) -> SEIR.model


#### Assigning parameters
params1 <- c(Beta0=beta0,Beta1=beta1,Gamma=gamma,Mu=mu,Nu=nu,Sigma=sigma,Iota=iota,S0=S_0,E0=E_0,I0=I_0,R0=R_0)

###### Deterministic Trajectory
trac_float <- trajectory(SEIR.model,params=params1,as.data.frame=TRUE)
trac <- trac_float
trac[,c(1,2,3,4,5)] <-round(trac[,c(1,2,3,4,5)],0)

#################################
############################################
#######################################################

start_year=50

days   <- data.frame(365.25*sims['time']) ### days
colnames(days) <- "days"
sims <- data.frame(sims,days)
df <- data.frame(sims,trac['C'])
colnames(df)[10] <- "C_deter"

df3 <- filter(df,days>start_year*365)



sims_monthly <- filter(df,row_number() %% 30==1)

#df1 <- subset(sims_monthly, select=c(C,C_deter,days,sim))
df1 <- subset(df, select=c(C,days,sim))
df2 <- melt(df1, id=c("days","sim"))


fig <- ggplot(df2, aes(x=days, y=value,group=sim,color=variable))+
geom_line(size=0.4)
fig + facet_wrap(~sim,ncol = 2)


plot_dir="E:/Documents/POMP/SEIR/Plots/"
fname='Stochastic sims.pdf'
# ggsave(filename=fname,path=plot_dir,device="pdf",
#        scale = 1, width = 6, height = 9, units = "in",
#        dpi = 300, limitsize = TRUE)



sim1 <- filter(df3, sim==1)
fig2 <- ggplot(sim1,aes(x=days))+
        geom_line(aes(y=C,colour='Stochastic Simulation'),size = 0.5)+
        geom_line(aes(y=C_deter,color='Deterministic Trajectory'),size = 1)+
        labs(x= "Time (years)",y="True Incidence")+
        scale_x_continuous(breaks=c(50*365,55*365,60*365,65*365),labels=c(50,55,60,65))
fig2


plot_dir=file.path(getwd(), "Plots")
fname='Deterministic trac vs stochastic - last 15yrs.pdf'
# ggsave(filename=fname,path=plot_dir,device="pdf",
#        scale = 1, width = 9, height = 6, units = "in",
#        dpi = 300, limitsize = TRUE)

sim1 <- filter(df, sim==1)
fig3 <- ggplot(sim1,aes(x=days))+
  geom_line(aes(y=C,colour='Stochastic Simulation'),size = 0.5)+
  geom_line(aes(y=C_deter,color='Deterministic Trajectory'),size = 1)+
  labs(x= "Time (years)",y="True Incidence")+
  scale_x_continuous(breaks=c(10*365,20*365,30*365,40*365,50*365,60*365),labels=c(10,20,30,40,50,60))
fig3


plot_dir=file.path(getwd(), "Plots")
fname='Deterministic trac vs stochastic - 65yrs.pdf'
# ggsave(filename=fname,path=plot_dir,device="pdf",
#        scale = 1, width = 9, height = 6, units = "in",
#        dpi = 300, limitsize = TRUE)

########### Stochastic parameter estimation
#################################
############################################
#######################################################

reported <- filter(sims,sim==1 & days>=50*365 )
reported <- select(reported,C,time,days)
colnames(reported) <- c("cases", "time","days")
reported['days']  <- reported['days']-reported[1,'days']
reported['time']  <- reported['time']-reported[1,'time']

data_file_txt=file.path(getwd(), "Data","incidence.txt")
data_file_xlsx=file.path(getwd(), "Data","incidence.xlsx")

# write.table(reported, data_file_txt, sep="\t")
# write.xlsx(reported, data_file_xlsx, sheetName = "Sheet1", 
#            col.names = TRUE, row.names = TRUE, append = FALSE)






