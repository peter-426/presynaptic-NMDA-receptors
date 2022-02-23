//! Endoplasmic reticulum (ER)

#ifndef _er_ryr_receptor_h_included_
#define _er_ryr_receptor_h_included_
#include "er_ryr_receptor.h"
#endif

class ER {
 
public: 

static constexpr double c_rest_ER = 5e+6; //1.3e+6; // 5,000,000 nM  was 1.30e+6 nm in Tewari. 
                                      // Resting ER [Ca2+]

                                  
                                   // Might not use leak and pump fluxes for removal of Ca2+ because  
                                   // a single calcium decay time constexprant usually works well enough.
static constexpr double v2=0.2374e-3;  // Leak of calcium from ER to cytosol in ms   
static constexpr double v_serca =2e-4; // nM/cm^2 ms,   p203 Gabbiani                                
static constexpr double k_serca=2000;  // nM,  2 uM on  p203 Gabbiani


double * cer;    // ER Calcium concentration

RyR ryr; 


ER() { ;  }


ER(int tn) 
{    
   cer= init_double(tn+2);   // ER Calcium concentration
   
   ryr = RyR(tn+2);
}

void set(int tn)
{
    for(int i = 1; i < tn+1; ++i)
    {
       cer[i]=0; // ER Calcium concentration
    }
    
    cer[1]=c_rest_ER;
}
};