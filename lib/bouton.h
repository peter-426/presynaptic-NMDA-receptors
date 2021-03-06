//!  Presynaptic Bouton class
/*!
   The Bouton class includes an implementation of the Hodgkin-Huxley model of
   action potentials.
   
   This bouton has voltage gated calcium channels and presynaptic NMDA receptors
   (preNMDARs) on its plasma membrane, a vesicle object that includes at least one calcium sensor, 
   and an endoplasmic reticulum (ER) with ryanodine receptors (RyR). 
*/   


#ifndef _bouton_receptors_h_included_
#define _bouton_receptors_h_included_
#include "bouton_receptors.h"
#endif

#ifndef _vesicles_hill_h_included_
#define _vesicles_hill_h_included_
#include "vesicles_hill.h"
#endif

#ifndef _er_h_included_
#define _er_h_included_
#include "er.h"
#endif

#ifndef _utilities_h_included_
#define _utilities_h_included_
#include "utilities.h"
#endif

class Bouton 
{
public:
 
//  Bouton parameters
static constexpr double bouton_rad=0.3141;  // um;  bouton radius;  e-4 if  cm; 
double bouton_vol;                   // bouton vol;  L;   
                                     // 0.13 um^3 in Koester & Sakmann (2000), citing Shepherd & Harris (1998).
                                     // Direct measurement gave an average varicosity diameter of 0.40 +/- 0.13 um. 
                                     // By comparison, the calculated diameter derived from length and volume was 
                                     // 0.33 +/- 0.10 um.
                                     
double domain_vol;                   // domain volume (vol of local region at channel)
double bouton_sa;                    // Surface area of bouton;  cm^2;   Koester & Sakmann (2000)
 
// CA3 presynaptic	CA1 post synaptic: not used here, using original HH paramters (for squid!).
// gNa	   100	     100        // Assumption: we just need to simulate an action potential,                 
// gK(DR)  	8	    12          // the exact type is not important for this simulation.
// gK(A)	  17	      17
// gK(AHP)	 7	      2.7
// gK(C )	 36.6	    33
// gKtotal	 68.6	    64.7
// gLeak	   0.33    	0.33
// VNa	    50      	50
// VK 	   -80 	    -80
// Vleak	 -65 	    -65 
 
static constexpr double vr= -70;  // Resting membrane potential of bouton;   mV
static constexpr double gna=120;  // Sodium conductance density;  mS/cm^2;   HH (1952)
static constexpr double gk=36;    // Potassium conductance density; mS/cm^2; HH (1952)

static constexpr double gl= 0.27; // Leak conductance density;  mS/cm^2;  HH (1952)
static constexpr double vl=-59.4; // Rev potential for leak;    mV;       HH (1952)
static constexpr double vna=45;   // Rev potential for Sodium ion;  mV;   HH (1952)
static constexpr double vk=-82;   // Rev potential for Potassium ion; mV; HH (1952) 
    
/////// ER
static constexpr double c1=0.185;     // Ratio of Vol ER to Vol Neuron cell; ;   
                                  // Shuai & Jung (2002)
                                   
////////////////  Plasma membrane (PM)  ////////////////////////////////////////
static constexpr double Ip=0.225;     // Max PM Ca current; was was 0.4 uA per cm2;
                                  // determined thru computer simulations
                                
static constexpr double k_pump=100;   // MM constexprant for PM Ca; nM; Erler et al (2004)
static constexpr double v_leak=0.001022664392140;// Maximum leak of Ca2+ thru PM; 
                                             // per ms; determined thru sim by Tewari
                                  
static constexpr double c_rest_bouton = 100; // Resting [Ca] for bouton average nM
static constexpr double R=8.3144598;      // Real Gas constant;  J /mole /K 
                                      // R is equivalent to the Boltzmann constant, but expressed in units of energy.
// R is the constant of proportionality that relates the energy scale in physics to the temperature scale, 
// when a mole of particles at the stated temperature is being considered.
//
// R occurs in the ideal gas law:  PV = nRT where P is the absolute pressure (SI unit pascals), V is the volume of gas (SI unit cubic metres), 
// n is the amount of gas (SI unit moles), m is the mass (SI unit kilograms) contained in V, and T is the thermodynamic temperature (SI unit kelvins).

static constexpr double T=273.15+21;    // Temperature in the experiments of Perea & Araque (2007); Kelvin.
static constexpr double F=96485.33289;  // Faraday's Constant;  Coulombs/mole; http://physics.nist.gov/cgi-bin/cuu/Value?F

// A joule and a watt second are the same thing, watt seconds and joules are units of energy
// An electrical engineer would say a watt second is the energy represented by one ampere of 
// current through a one volt potential for one second. 

// Electrical units
// Volts:   a measure of potential; 
// amperes: a measure of current (rate of flow of electrons);
// coulomb: a measure of charge.

// 6.241 x 10^18  electrons per coulomb.  
// One coulomb flowing for one second equals one ampere of current. 

// 1 J = 1 C x 1 V

// time constant for [Ca] decay in milliseconds
static constexpr  double tau_dec = 100; // bouton avg; 238ms 1/2 decay (Wu 1994),  27ms (p138 Sterrat)                        
    
// Pre-synaptic Bouton Variables
double G_syn;    // Synaptic glutamate concentration
double * ca_local;     // Calcium concentration near vesicles
double * ca_global;    // Calcium concentration for bouton (average)

double * v;      // Membrane potential
double * m;      // Sodium channel activation 
double * h;      // Sodium channel inactivation
double * n;      // Potassium channel activation


double *  Ivgcc;
double *  IPump;
double *  ICa_leak;

// Reversal potential for Calcium ion determined through Nernst equation. 
// Assumes extracellular [Ca2+]=2 mM as in Perea & Araque (2007)
double vca;

// Receptors
double * Inmda;
double * Inmda_Ca;

double * ca_VGCC;
double * ca_PreNMDAR;
double * ca_PreNMDAR_mean;

double * ca_RyR;
double * ca_VGCC_RyR;

//double * ca_IP3R;

PreNMDAR nmdaR;

VGCC_bouton vgcc;



#ifdef Hill
Vesicle_Hill ves;
#endif

#ifdef Markov
Vesicle_Markov ves;
#endif 

#ifdef Markov6
Vesicle_Markov_6 ves;
#endif 

#ifdef Allosteric
Vesicle_Allosteric ves;
#endif 


ER er;       

Bouton() { ; }

Bouton(int tn, double v_ca) 
{    
    

    vca = v_ca;  
    bouton_vol= 0.13;  // 1e-3 * (4.0/3.0) * M_PI *pow(bouton_rad,3);  // Volume of bouton; liter; Koester & Sakmann (2000)
    domain_vol=bouton_vol/10000;   // hypothetical volume of the channel domain, volume at source of Ca2+ influx 
    
    // sphere V = 4/3 pi r^3 
    // circle A =     pi r^2
    
    // if 1 AZ per bouton and vol_bouton is 0.13 um^3 and vol_local is V=pi*r^2*h
    //
    // vol_bouton/vol_local == 0.13/(pi * 0.03^2 * 0.050) = 919               
    //                                               
    bouton_sa= 1;  // 4* M_PI *pow(bouton_rad,2);  //bouton surface area; cm^2 Koester & Sakmann (2000)
    //
    // Pre-synaptic Bouton Variables
    ca_local  =  init_double(tn);   // Calcium concentration (local,  for vesicles)
    ca_global =  init_double(tn);   // Calcium concentration (global, for bouton)
    
    v=   init_double(tn);   // Membrane potential
    m=   init_double(tn);   // Sodium channel activation 
    h=   init_double(tn);   // Sodium channel inactivation
    n=   init_double(tn);   // Potassium channel activation

    Ivgcc   = init_double(tn);
    IPump   = init_double(tn);
    ICa_leak= init_double(tn);
       
    // Receptors on Bouton
    Inmda     = init_double(tn);     
    Inmda_Ca  = init_double(tn); 
    
    // printf("\n bouton.h line 158, tn = %d \n", tn);
    nmdaR = PreNMDAR(tn, vca);
    
    ca_VGCC          = init_double(tn); 
    ca_PreNMDAR      = init_double(tn); 
    ca_PreNMDAR_mean = init_double(tn); 
    
    ca_RyR      = init_double(tn);
    ca_VGCC_RyR = init_double(tn);
       
    // ca_IP3R = init_double(tn); 
    
    vgcc = VGCC_bouton(tn, vca); 
    
    #ifdef Markov
    ves = Vesicle_Markov(tn);
    #endif
    
    #ifdef Markov6
    ves = Vesicle_Markov_6(tn);
    #endif
    
    #ifdef Allosteric
    ves = Vesicle_Allosteric(tn);
    #endif 
    
    #ifdef Hill
    ves = Vesicle_Hill(tn);
    #endif

    er = ER(tn);
};
  
  
void set(int tn)
{
    for(int i = 1; i < tn+1; ++i)
    {
       ca_local[i]  =0;   // Calcium concentration
       ca_global[i] =0;
       
       v[i]=0;   // Membrane potential
       m[i]=0;   // Sodium channel activation 
       h[i]=0;   // Sodium channel inactivation
       n[i]=0;   // Pottassium channel activation
   
       Ivgcc[i]=0;
       IPump[i]=0;
       ICa_leak[i]=0;
       
       Inmda[i]=0;     
       Inmda_Ca[i]=0; 
       
       ca_VGCC[i]=0;
       ca_PreNMDAR[i]=0;
       ca_RyR[i]=0;
    }
    
    // Initial conditions
    
    ca_local[1] =c_rest_bouton;  
    ca_global[1]=c_rest_bouton;          
    
    v[1]=-70;  // Resting membrane potential of bouton;  mV
    m[1]=0.1;  // Gating variable for sodium channel (activation)
    h[1]=0.6;  // Gating variable for sodium channel (inactivation)
    n[1]=0.3;  // Gating variable for potassium channel (activation)
 
    Inmda[1]=0;
    Inmda_Ca[1]=0;
};  
  
void bouton_model(int i, EX ex, double aG_syn, int AP5, int RY) 
{
    // Gating Variables
    // an Opening: K channel activation 
    // bn Closing: K channel activation
    // am Opening: activation of Na channel
    // bm Closing: activation of Na channel
    // ah Opening: inactivation of Na channel
    // aq Closing: inactivation of Na channel

    double an=0.01*((-v[i]-60)/(exp((-v[i]-60)/10)-1)); 
    double bn=0.125*exp((-v[i]-70)/80);                
    double am=0.1*((-v[i]-45)/(exp((-v[i]-45)/10)-1)); 
    double bm=4*exp((-v[i]-70)/18);      
    double ah=0.07*exp((-v[i]-70)/20);   
    double bh=1/(exp((-v[i]-40)/10)+1);   

    // Ionic Currents
    double I_Na   = gna* (v[i]-vna);    // Sodium current;     uA per cm^2
    double I_K    = gk * (v[i]-vk);     // Potassium current;  uA per cm^2
    double I_Leak = gl * (v[i]-vl);     // Leak current;       uA per cm^2
    
    m[i+1] =m[i]+ex.deltaT*(am*(1-m[i])-bm*m[i]);   // m: Sodium channel activation
    h[i+1] =h[i]+ex.deltaT*(ah*(1-h[i])-bh*h[i]);   // h: Sodium channel inactivation
    n[i+1] =n[i]+ex.deltaT*(an*(1-n[i])-bn*n[i]);   // n: Potassium channel activation
  
    v[i+1] =v[i]+ex.deltaT*\
               ( ex.Iapp[i] - (pow(m[i],3)*h[i]* I_Na + pow(n[i],4) * I_K + I_Leak) );  //  (m^3 * h * I_Na) + (n^4 * I_K) 
    
    
    // Ca2+ plasma membrane (PM) flux, using tau_decay instead of explicit pump and leak fluxes
    //
    IPump[i] =0;    // Ip*pow(ca_local[i],2)/(pow(ca_local[i],2)+pow(k_pump,2));  // PMCA; uA per cm^2
    ICa_leak[i] =0; // v_leak*(v[i]-vca);                                         // Calcium leak; uA per cm^2  
    //
    //  
    //
    Ivgcc[i] = vgcc.I_Ca(i, ex.deltaT, v[i], ca_VGCC[i]);   // calcium current due to a number (1?) of  VGCCs
    
    //
    double fluxRyR=0, fluxVGCC=0, fluxPreNMDAR=0;  // change in concentration due to these channels
    //
    // convert Ca current to Ca flux in nM.  Note: local [Ca] at channel may be >> global [Ca]

 
    // fluxVGCC = ((-Ivgcc[i]  * bouton_sa)/(2*F*bouton_vol)); 
    //
    // F:  Faraday's Constant;  Coulombs/mole: 9.648 A 10^4 C/mol,  Fall, Computational Cell Biology, 2002, p107
    //
    // Ca2+ generally enters cells through ion channels. 
    // Ionic currents are measured in pA (10^12 C/s) or nA (10^9 C/s), 
    // so an additional factor is needed to convert to umoles:
    //                jflux = - alpha*ICa
    //
    // where alpha =  1/(2FVi).
    //
    // As an example, if ICa is measured in pA, and cytosolic volume is measured
    // in um^3, then jflux will have units of M/ms.
    //
    //
    double number_of_VGCCs     = 33;
    double number_of_preNMDARs = 33;
    
    double surface_area = 0.03;  // not used.  // surface area of membrane where channels are located; A=pi*r^2: active zone = pi * 0.5^2 
    double bouton_volume= 100;   // used to estimate global [Ca] from local [Ca] 
    
    double Fmicro = F/1e6;   // because F is in moles (M), and this model uses micro-moles (uM)
    
    //Note: ACh in NMJ synaptic cleft: 4 x e-6 cm^2/s == 4 x e2 um^2/s == 4 x e5 nm^2/s == 4 x e2 nm^2/ms
    //      Assuming a 50 nm cleft and that t ~ x^2/2D, t = 3.1 us
    //
    double vgcc_distance  = 0.090; // distance from VGCC to vesicle;  .10 um == 100 nm
    double nmdaR_distance = 0.030; // distance from nmdaR to vesicle; um
    // double DCa=0.220;           // diffusion coefficient: 0.220 um^2 /ms; 
    double DCa=0.050;              // diffusion coefficient: 0.050 um^2/ms   // Nadkarni et al. 2010
    //                                  
    // Nadkarni 2012:  50 um^2/s   == 0.05  um^2/ms  
    //
    // See Sterratt (2011) p 138:  Using F etc. we get the rate of change in [Ca], i.e. the flux.
    //
    // Jcc = - (area * ICa) / (2*F*volume):  Jcc is the change in calcium concentration in the compartment.
    // 
    // Change in [Ca2+] at distance vgcc_distance form the VGCCs due to the influx of Ca2+ ions from the cluster of VGCCs.
    // Key point:  Surface area is 1 (enough area for one cluster of VGCCs); and here we divide by distance, not volume.
    fluxVGCC = (-Ivgcc[i] * 1)/(2 * Fmicro * DCa * vgcc_distance); 
    
    ca_VGCC[i+1] = ca_VGCC[i] + ex.deltaT * (fluxVGCC - ((ca_VGCC[i] - 0)/tau_dec));
     
    // Calcium influx from RyR.   The release of Ca into a confined space between
    // the cell membrane and the SR can result in a much higher local [Ca] than is
    // predicted by a spatially homogeneous model.   Keener and Sneyd 1998, p181.
    //
    if (RY==0)  // if no RyR blocker
    {                                                          // ex.beg_pad (padding) is in milliseconds
       double delay_time_steps = (ex.beg_pad + 1.0)/ex.deltaT; // delay is very short, 1 ms??
       
       if (i > delay_time_steps) 
       {
         fluxRyR = 1 * er.ryr.Jcicr(i, ca_local[(int)(i - delay_time_steps)], \
                                     er.cer[(int)(i - delay_time_steps)]); 
         
         ca_RyR[i+1]  = ca_RyR[i] + ex.deltaT*(fluxRyR - ((ca_RyR[i] - 0)/tau_dec)); 
       }
       else
       {
         fluxRyR=0;
         ca_RyR[i+1] = 0;
       }
    }
    else
    {
       ca_RyR[i+1] = 0; 
    }
    
    // adjust for volume of ER versus volume bouton, assumes ER is c1=0.185 of bouton cytosolic volume;                               
    //
    er.cer[i+1] = er.cer[i] - ex.deltaT*( (fluxRyR -  (( ca_RyR[i] - 0)/tau_dec ))/c1 );
    
    if (AP5 == 0)  
    {    
      double delay_time_steps = (ex.beg_pad+1.0)/ex.deltaT;   // delay is probably ~1 ms or 20 dt steps
     
      double temp = 0; 
      if (i > delay_time_steps) {
        int gluTimePoint = (int) i - delay_time_steps;
        temp = nmdaR.I_Ca(i, ex.deltaT, ves.G_syn[gluTimePoint], v[i], ca_PreNMDAR[i]);
      }
      
      Inmda[i]    = temp;           // for plotting calcium current
      Inmda_Ca[i] = temp;

      // Convert Ca current to Ca flux in nM, but this flux is based on Glu evoked by previous spike.
      // Assumption: it is not possible for the preNMDARs to facilitate vesicle release evoked by 
      // the first spike.

      fluxPreNMDAR = (-Inmda_Ca[i] * 1)/(2*Fmicro * DCa * nmdaR_distance); 

      ca_PreNMDAR[i+1] = ca_PreNMDAR[i] + ex.deltaT*(fluxPreNMDAR - ((ca_PreNMDAR[i] - 0)/tau_dec)); 

      ca_PreNMDAR_mean[i+1] += ca_PreNMDAR[i+1];  // calculate mean at end of simulation
    }
   
    
    // ============ divide local Ca fluxes by bouton_volume to approximate global [Ca] =============
                                                                                                                  
    ca_local[i+1] = ca_local[i] + ex.deltaT * \
                    (fluxVGCC + fluxRyR + fluxPreNMDAR  - (ca_local[i] - c_rest_bouton)/tau_dec);   
                                                               
                                  
    ca_global[i+1]= ca_global[i]+ex.deltaT * ( (number_of_VGCCs * fluxVGCC/bouton_volume) + \
                                              fluxRyR/bouton_volume + \
                                              (number_of_preNMDARs * fluxPreNMDAR/bouton_volume)  \
                 - (ca_global[i]-c_rest_bouton)/tau_dec);
        
    ves.release(i, ex, v[i], vr, ca_local[i], AP5);
 }
 
};    

 