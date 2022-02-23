
#ifndef _utilities_h_included_
#define _utilities_h_included_
#include "utilities.h"
#endif

#ifndef _math_h_included_
#define _math_h_included_
#include "math.h"
#endif

#ifndef _stdio_h_included_
#define _stdio_h_included_
#include <stdio.h>
#endif

class Astro
{

public:

//-------------- Astrocyte Parameters -----------------------
// All values are from De Pitta et al (2009) except 'aNa.'
static constexpr double aNa=20;     // Number of IP3 channels in cluster;  ; Nadkarni & Jung (2007)

static constexpr double ac0=2e+3;   // Total cell free Ca2+ concentration; nM

static constexpr double ac1=0.185;  // Ratio of ER volume to cytosol volume;  

static constexpr double av1=6e-3;        // Maximal IP3R flux; per ms
static constexpr double av2=0.11e-3;     // Maximal rate of Ca2+ leak from ER; per ms
static constexpr double av3=0.9;         // Maximal rate of SERCA uptake; nM per ms
static constexpr double ak3=0.1e+3;      // SERCA Ca2+ affinity; nM

static constexpr double ad1=0.13e+3;     // IP3 dissociation constexprant; nM
static constexpr double ad2=1.049e+3;    // Ca2+ inactivation dissociation constexprant; nM
static constexpr double ad3=0.9434e+3;   // IP3 dissociation constexprant; nM
static constexpr double ad5=0.08234e+3;  // Ca2+ activation dissociation constexprant; nM
static constexpr double aa2=0.2e-6;      // IP3R binding rate for Ca2+ Inhibition; per nM ms

// IP3 regulation parameters
// Following Hofer (2002), two processes mediate IP3 production.
// 1) the plc Beta  IP3 production process is activated thru a GPCR mechanism,
//    which requires an extracellular agonist.
// 2) the plc Delta IP3 production process depends on [Ca2+]i

             double av_plcb;         // Maximal rate of IP3 production by PLC Beta; nM per ms
static constexpr double av_plcd=0.05;    // Maximal rate of IP3 production by PLC Delta; nM per ms 

static constexpr double ak_plcd=1.5e+3;  // Inhibition constexprant of PLC Delta activity; nM
static constexpr double aK_plcd=0.1e+3;  // Ca2+ affinity of PLC Delta; nM
static constexpr double ar5p=0.05e-3;    // Maximal rate of degradation by IP-5P; per ms

static constexpr double aK_R=1.3e-3;  // Glutamate affinity of the receptor; mM
static constexpr double aK_P=10e-3;   // Ca2+/PKC-dependent inhibition factor; mM
static constexpr double aK_pi=0.6e+3; // Ca2+ affinity of PKC; nM
static constexpr double av_3K=2;      // Maximal rate of degradation by IP3-3K; nM per ms
static constexpr double aK_D=0.7e+3;  // Ca2+ affinity of IP3-3K; nM
static constexpr double aK3=1e+3;     // IP3 affinity of IP3-3K; nM

// Gliotransmitter release parameters
static constexpr double nva=12;        // Number of fusogenic SLMV inside astrocyte; 
                                   //       ; Malarkey & Parpura (2011)
static constexpr double gva=20;        // Glutamate concentration inside one SLMV; mM;  Montana et al (2006)
static constexpr double ak1=3.75e-6;   // Ca2+ association  rate for S1; per nM per ms; Bertram et al (1996)
static constexpr double ak2=2.5e-6;    // Ca2+ association  rate for S2; per nM per ms; Bertram et al (1996)
static constexpr double akk3=1.25e-05; // Ca2+ association  rate for S3; per nM per ms; Computationally determined
static constexpr double ak_1=4e-4;     // Ca2+ dissociation rate for S1; per ms;        Bertram et al (1996)
static constexpr double ak_2=1e-3;     // Ca2+ dissociation rate for S2; per ms;        Bertram et al (1996)
static constexpr double ak_3=10e-3;    // Ca2+ dissociation rate for S3; per ms;        Computationally determined

// Glutamate regulatory parameters
static constexpr double adegG=10;  // Glutamate clearance rate from the extra-synaptic cleft; per ms; Destexhe et al (1998)
static constexpr double atau_rec=800;  // SLMV recovery time constexprant; ms; Tsodyks & Markram (1997)
static constexpr double atau_inact=3;  // SLMV inactivation time constexprant; ms; Tsodyks & Markram (1997)


// Astro Vm depends on, guessing, Vrest versus leak and calcium currents (channels and pumps)
static constexpr double gL=0.1; 
static constexpr double gCa=10;  
static constexpr double VCa=80;    
static constexpr double VL=-60;
static constexpr double Vrest=-85;

// ??
static constexpr double R_in=0.7985e+8;    // Input resistance of dendrite spine; k-ohm; Calculated
static constexpr double tau_mem=50;        // Post-synaptic membrane time constexprant; ms; Tsodyks & Markram (1997)

//------ Astrocyte Variables -------------
double * ca;      // Calcium concentration
double * ax;      // IP3R gating variable
double * a_ip3;   // IP3 concentration
double * aO1;     // S1 with calcium bound 
double * aO2;     // S2 with calcium bound
double * aO3;     // S3 with calcium bound
double * aE_syn;  // Fraction of effective SLMV in extra-synaptic cleft 
double * aI_syn;  // Fraction of inactivated SLMV
double * aR_syn;  // Fraction of releasable SLMV
double * aG_syn;  // Glutamate concentration in extra-synaptic cleft
 
 

Astro();  
Astro(int tn);
    
void set(int tn); 
void astro_model(int i, double t, double deltaT, double tdr, double G_syn); 
};


Astro::Astro() { ; }

Astro::Astro(int tn) 
{    
    ca=    init_double(tn+2);  
    ax=    init_double(tn+2);
    a_ip3= init_double(tn+2);
    aO1=   init_double(tn+2);
    aO2=   init_double(tn+2); 
    aO3=   init_double(tn+2);  
    aE_syn= init_double(tn+2);   
    aI_syn= init_double(tn+2);
    aR_syn= init_double(tn+2);  
    aG_syn= init_double(tn+2);  
}

void Astro::set(int tn) 
{ 
    for(int i=0; i < tn+2; ++i)
    {
        ca[i]=0;  
        ax[i]=0;
        a_ip3[i]=0;
        aO1[i]=0;
        aO2[i]=0; 
        aO3[i]=0;  
        aE_syn[i]=0;   
        aI_syn[i]=0;
        aR_syn[i]=0;  
        aG_syn[i]=0;  
    }
    
    av_plcb=0.5;    // max ip3 production rate
            
    // initial conditions   
    ca[1]=100;    // Resting calcium concentration; nM
    
    ax[1]=0.5;    // Gating variable for IP3R 
    aO1[1]=0.48;  // 1st calcium binding site of SLMV
    aO2[1]=0.2;   // 2nd calcium binding site of SLMV
    aO3[1]=5e-4;  // 3rd calcium binding site of SLMV
    aE_syn[1]=0;  // Fraction of effective SLMV in extra-synaptic cleft
    aI_syn[1]=0;  // Fraction of inactivated SLMV
    aR_syn[1]=1;  // Fraction of releasable SLMV
    a_ip3[1]=160;    // Resting [IP3]; nM
    aG_syn[1]=1e-3;  // Basal glutamate in the extra-synaptic cleft; mM
}

void Astro::astro_model(int i, double t, double deltaT, double tdr, double G_syn) 
{ 
//  // Astrocyte Processes
// Time-constexprants for three binding sites of SLMV
double atau1=(ak1*ca[i]+ak_1);   // closure of O1; per ms
double atau2=(ak2*ca[i]+ak_2);   // closure of O2; per ms
double atau3=(akk3*ca[i]+ak_3);  // closure of O3; per ms

//IP3 production & degradation terms.     All are from De Pitta et al (2009)
double aplcb_ca= 1+(aK_P/aK_R)*(ca[i]/(ca[i]+aK_pi));  // Calcium-dependent inhibition of 'ap_glu'
double ap_glu=  av_plcb* pow( G_syn, 0.7) / (pow(G_syn, 0.7) + pow((aK_R*aplcb_ca),0.7));      // Agonist-dep IP3 production
double ap_plcd= av_plcd*(1/(1+(a_ip3[i]/ak_plcd)))*pow(ca[i],2)/(pow(ca[i],2)+pow(aK_plcd,2)); // Agonist-indep IP3 production
double ap_mapk= av_3K* (pow(ca[i],4)/(pow(ca[i],4) + pow(aK_D,4))) *(a_ip3[i]/(a_ip3[i]+aK3)); // IP3 deg due to IP3-3K

double ap_deg=ar5p*(a_ip3[i]);  // IP3 degradation due to IP-5P

  // IP3R gating variables
double aaq=aa2*ad2*(a_ip3[i]+ad1)/(a_ip3[i]+ad3);

double abq=aa2*ca[i];            // gating variables for IP3R

double aminf=a_ip3[i]/(a_ip3[i]+ad1);
double aninf=ca[i]/(ca[i]+ad5);  // Steaty-state functions governing IP3 & Ca2+ based closing of IP3R  
  
double auer=(ac0-ca[i])/ac1;     // ER Calcium concentration

// Box-Muller Algorithm (Fox, 1997) for noise-term in equation (12) of Tewari & Majumdar 2012.
double au1=rnd();   
double au2=rnd();                // uniformly distributed random variables

double aa=((aaq*(1-ax[i])-abq*ax[i]))/aNa;  // co-variance 
double dW=sqrt(-(2*deltaT*(aa)*log(au1)))*cos(2* M_PI *au2);  // independent Gaussian random number
  
// ax[i] must be in range [0,1].
if (ax[i] + dW >=0   &&   ax[i] + dW <=1) {
    ax[i+1]=deltaT*( aaq*(1-ax[i])-abq*ax[i] )+ax[i]+dW;
}
else {
    ax[i+1]=deltaT*( aaq*(1-ax[i])-abq*ax[i] )+ax[i];
}
 
//printf("%f %f  %f \n", ax[i], aaq, abq);
//ax[i]=1.0;   // 0.5 an .75 drops down to 60,  1.0 goes up and oscillates

double ajchan = ac1*av1*pow(aminf,3) * pow(aninf,3)* pow(ax[i],3)*(ca[i]-auer);  // Calcium flux through IP3R

double ajpump = av3*pow(ca[i],2) / ( pow(ak3,2)+pow(ca[i],2) );  // SERCA pump flux

double ajleak = ac1*av2*(ca[i]-auer);  // Calcium leak from cytosol into ER
  
// Astrocyte calcium clamp.  
// clamp on from t=0 to t=30 seconds.   
//     if t[i]>=0 && t[i]<=30e+3
//         ca[i+1]=ca[1];                 // Calcium clamped at resting concentration
//     elseif t[i]>30e+3 && t[i]<=600e+3

// resting calcium should be 100 nM ???   should not decay to ~60 nM

ca[i+1]= ca[i]+ deltaT*(-ajchan-ajpump-ajleak);  // Dynamic calcium concentration
     
// printf("ca[i+1]=%f  ca[i]=%f  ax=%f,  ajchan=%f,  ajpump=%f,  ajleak=%f \n", \
// ca[i+1], ca[i], ax[i], ajchan, ajpump, ajleak);
   
if ( ca[i+1] < 0 ||  isnan(ca[i+1])  ) 
{
   printf("ca[i+1]=%f  ca[i]=%f  ax=%f,  ajchan=%f,  ajpump=%f,  ajleak=%f \n", \
   ca[i+1], ca[i], ax[i], ajchan, ajpump, ajleak);
     
   exit(0);
}

a_ip3[i+1]=a_ip3[i]+deltaT*(ap_glu+ap_plcd-ap_mapk-ap_deg);  // IP3 concentration

aO1[i+1]=aO1[i]+deltaT*(ak1*ca[i]-(aO1[i]*atau1));     // Site 1 of SLMV with calcium bound
aO2[i+1]=aO2[i]+deltaT*(ak2*ca[i]-(aO2[i]*atau2));     // Site 2 of SLMV with calcium bound
aO3[i+1]=aO3[i]+deltaT*(akk3*ca[i]-(aO3[i]*atau3));    // Site 3 of SLMV with calcium bound

double aRRP=aO1[i]*aO2[i]*aO3[i];  // Releasable SLMVs due to Ca increase

// Fraction of releasable SLMVs
aR_syn[i+1]=aR_syn[i]+\
   deltaT*(((aI_syn[i])/atau_rec)-((heaviside(ca[i]-196.69)*aRRP)*aR_syn[i]));  
   
// Fraction of effective SLMVs
aE_syn[i+1]=aE_syn[i]+\
   deltaT*(((heaviside(ca[i]-196.69)*aRRP)*aR_syn[i])-(aE_syn[i]/atau_inact));  
   
aI_syn[i+1]=1-aR_syn[i+1]-aE_syn[i+1];  // Frac of inactivated SLMVs

// Glutamate concentration in the extra-synaptic cleft
aG_syn[i+1]=aG_syn[i]+deltaT*(nva*gva*aE_syn[i]-adegG*(aG_syn[i]));  
}

