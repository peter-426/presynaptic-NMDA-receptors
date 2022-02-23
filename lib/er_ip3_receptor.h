//! IP3 Receptor class.

/*!  IPRs reside in the ER membrane.

When IP3 binds to IP3 receptors, calcium is released from the ER into the cytosol.

*/

#ifndef _utilities_h_included_
#define _utilities_h_included_
#include "utilities.h"
#endif

#ifndef _math_h_included_
#define _math_h_included_
#include "math.h"
#endif

class IP3_Receptor
{
 
public:
 
// ER constants for ip3R 
static constexpr double c1=0.185;  // Ratio of Vol ER to Vol Neuron cell;dimless;   
                               // Shuai & Jung (2002)
                               
static constexpr double v1=10e-3;  // 30e-3;   // Max calcium flux thru IP3R;per ms;  
                               // Determined by simulation.
                                
                                
static constexpr double v2=0.2374e-3; // Leak of calcium from ER to cytosol;per ms;   
                                  // Determined by simulation.
                                  
static constexpr double v3=90;        // Max SERCA pump rate;nM/ms;    
                                  // Determined by simulation.
                                
static constexpr double k3=0.1e+3;    // MM constant for SERCA;nM;  Erler (2004)
 
// these don't match Tewari's for ip3R "related" SERCA pump
static constexpr double v_serca = 2e-4;  // nM/(cm^2 ms)  Gabbiani p203 
static constexpr double K_serca = 2000;  // nM,   2 uM    Gabbiani p203

 
 
//IP3R-kinetics, Shuai,Jung(2002)
static constexpr double d1=0.13e+3;   // IP3 dissociation constexprant;nM;
static constexpr double d2=1.049e+3;  // Inhibitory Ca2+ dissociation constexprant;  nM;   
                                  
static constexpr double d3=0.9434e+3; // IP3 dissociation constexpr;nM;
static constexpr double d5=0.08234e+3;// Activation Ca2+ dissociation constexprant; nM;   
                                  
static constexpr double a2=0.2e-6;    // Inhibitory Ca2+ binding constexprant;  per nM ms;   
                              

//IP3 generation parameters
static constexpr double v_glu=0.062;    // Maximal IP3 production rate from mGluRs;  
                                    // nM per ms;   Nadkarni & Jung (2008)
static constexpr double k_glu=0.78e-3;  // Glutamate conc at which v_glu is halved;  
                                    // mM;   Nadkarni & Jung (2008)
static constexpr double tau_ip3=1400;   // 0.14e-3;// IP3 degradation constexprant;  per ms;   
                                    // Nadkarni & Jung (2008)
static constexpr double np=0.7;         // 0.3;
////////////////////////////
 
 
double * p;      // IP3 concentration
double * q;      // IP3 gating variable



IP3_Receptor(int tn) {
 
    p=init_double(tn+2);  // IP3 concentration
    q=init_double(tn+2);  // IP3 gating variable

    p[1]=160;   // Resting [IP3];  nM
    q[1]=0.22;  // Gating variable for IP3R
}

 
double ip3R(int i, double deltaT, double ca,  double cer, double aGsyn) {
 
   // IP3R Kinetics:  steady states, fluxes, and IP3 production
   // ip3R gating
   double minf=p[i]/(p[i]+d1);  // Steaty-state: [IP3]  based closing of IP3R
   double ninf=ca/(ca +d5);     // Steaty-state: [Ca2+] based closing of IP3R
    
   double aq=a2*d2*(p[i]+d1)/(p[i]+d3); 
   double bq=a2*ca;  
    
   // flux between ER and cytoplasm wrt ER:  
   //      negative flux adds Ca to cytoplasm [Ca]
   //      positive flux adds Ca to ER [Ca]
   // jchan: Calcium flux through IP3R
   // jpump: serca pump flux
   // jleak: Calcium leak through ER into cytosol
   // p_glu: IP3 production due to extra-synaptic glutamate
   double jchan= c1 * v1 * pow(minf,3) * pow(ninf,3) * pow(q[i],3) * (ca - cer);  
   double jpump= 0; // v3*pow(ca,2)/( pow(k3,2)+pow(ca,2)  );                    // no c1 ??
   double jleak= 0; // c1*v2*(ca-cer);  
   double p_glu= v_glu*pow(aGsyn,np)/( pow(k_glu,np) + pow(aGsyn,np) );  
   
   // deltaT is applied to ip3 and gate, but not flux
   p[i+1]= p[i]+deltaT*( p_glu -  ( (p[i]-p[1])/tau_ip3) );  // p=Intracellular [IP3]
   q[i+1] =q[i]+deltaT*(aq*(1-q[i])-bq*q[i]);                // q: IP3R gating variable: fraction of ip3Rs activated  
   
   return (-jchan - jpump - jleak);
}

double serca(int i, double c) 
{
   return  v_serca *pow(c,2)/( pow(K_serca,2)+pow(c,2) );
}

double leak(int i, double c, double cer)
{
    return c1 * v2* (c - cer);
}
};





