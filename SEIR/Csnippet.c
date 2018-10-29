/* pomp model file: Csnippet */
/*  Time: 2018-07-31 12:50:45.652 -0400 Salt: C9B4771D22B256144EE9360D */

#include <E:/Documents/R/win-library/3.5/pomp/include/pomp.h>
#include <R_ext/Rdynload.h>

 
#define Mu	(__p[__parindex[0]])
#define Nu	(__p[__parindex[1]])
#define Sigma	(__p[__parindex[2]])
#define Gamma	(__p[__parindex[3]])
#define Beta0	(__p[__parindex[4]])
#define Beta1	(__p[__parindex[5]])
#define Iota	(__p[__parindex[6]])
#define Rho	(__p[__parindex[7]])
#define Tau	(__p[__parindex[8]])
#define S0	(__p[__parindex[9]])
#define E0	(__p[__parindex[10]])
#define I0	(__p[__parindex[11]])
#define R0	(__p[__parindex[12]])
#define S	(__x[__stateindex[0]])
#define E	(__x[__stateindex[1]])
#define I	(__x[__stateindex[2]])
#define R	(__x[__stateindex[3]])
#define C	(__x[__stateindex[4]])
#define CR	(__x[__stateindex[5]])
#define data	(__y[__obsindex[0]])
#define DS	(__f[__stateindex[0]])
#define DE	(__f[__stateindex[1]])
#define DI	(__f[__stateindex[2]])
#define DR	(__f[__stateindex[3]])
#define DC	(__f[__stateindex[4]])
#define DCR	(__f[__stateindex[5]])
#define TMu	(__pt[__parindex[0]])
#define TNu	(__pt[__parindex[1]])
#define TSigma	(__pt[__parindex[2]])
#define TGamma	(__pt[__parindex[3]])
#define TBeta0	(__pt[__parindex[4]])
#define TBeta1	(__pt[__parindex[5]])
#define TIota	(__pt[__parindex[6]])
#define TRho	(__pt[__parindex[7]])
#define TTau	(__pt[__parindex[8]])
#define TS0	(__pt[__parindex[9]])
#define TE0	(__pt[__parindex[10]])
#define TI0	(__pt[__parindex[11]])
#define TR0	(__pt[__parindex[12]])
#define lik	(__lik[0])

void __pomp_initializer (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{

                   S  = S0;
                   E  = E0;
                   I  = I0;
                   R  = R0;
                   C  = 0;
                   CR = 0;
                    
}


void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)
{

                  if (C > 0.0) {
                  lik = rbinom(C,Rho);
                  } else {
                  lik = rbinom(C+1,Rho);
                  }
                   
}


void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{

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
                  
                   
}

#undef Mu
#undef Nu
#undef Sigma
#undef Gamma
#undef Beta0
#undef Beta1
#undef Iota
#undef Rho
#undef Tau
#undef S0
#undef E0
#undef I0
#undef R0
#undef S
#undef E
#undef I
#undef R
#undef C
#undef CR
#undef data
#undef DS
#undef DE
#undef DI
#undef DR
#undef DC
#undef DCR
#undef TMu
#undef TNu
#undef TSigma
#undef TGamma
#undef TBeta0
#undef TBeta1
#undef TIota
#undef TRho
#undef TTau
#undef TS0
#undef TE0
#undef TI0
#undef TR0

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {
  ++__pomp_load_stack;
}

void __pomp_load_stack_decr (int *val) {
  *val = --__pomp_load_stack;
}

void R_init_Csnippet (DllInfo *info)
{
R_RegisterCCallable("Csnippet", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("Csnippet", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("Csnippet", "__pomp_initializer", (DL_FUNC) __pomp_initializer);
R_RegisterCCallable("Csnippet", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("Csnippet", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
}
