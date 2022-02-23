
#ifndef _math_h_included_
#define _math_h_included_
#include "math.h"
#endif

class AMPA_astro
{
    public:
     //   units are mM  and  ms        
         
    static constexpr double a_r=1.1;  
    static constexpr double a_d=0.19;    //   Large and fast EPSCs
    /*   
      In typical cortical cells, the rise time is 0.4 to 0.8 milliseconds.
      Ex.  Tmax/(1.1 + 0.19) ==> 0.8 rise time, 5 ms decay time.
       
      AMPA receptors onto inhibitory interneurons are about twice as fast 
      in rise and fall times as those onto excitatory neurons.
      Real AMPA synapses show quite strong depression. That is, the 
      peak amplitude of the AMPA current decreases with each subsequent spike. 
    */
    static constexpr double V_ampa=0;  
    static constexpr double gAMPA=0.35e-6;  // conductance mS;   Destexhe 1998:  0.35 .. 1.0 nS per syn
 
    double s;
 
    AMPA_astro() {
        s=0;
    }
 
 
    //Ermentrout and Terman, Mathematical Foundations of Neuroscience,2010,p161
    double syn(double deltaT, double glu, double Vm)
    {   
        double Iampa = gAMPA * s * (Vm - V_ampa);
        s = s + deltaT*(a_r * glu * (1 - s) - a_d * s);
        return Iampa;
    }
};


/*
    Astrocytes:
    1) NMDARs are not blocked by Mg++, [Mg++] has almost no impact on V_T,
    2) resting membrane potential of ~ -85 mV, if KO of Kir4.1 
       inward rectifying potassium channel, then depolarized to ~ -30 mV
*/
class NMDA_astro
{
    public:
    
    //  units are mM and ms, max [glu] in synaptic cleft ~ 1mM
    
    static constexpr double a_r=0.072;  
    static constexpr double a_d=0.0066;    // Destexhe 1998: max gNMDA = 0.01 .. 0.6 nS per syn
    static constexpr double Enmda=0;   
    static constexpr double gNMDA=0.7e-6;  // Traub 1994, p382, excitatory neuron
                                       // Gerstner, p53, max gNMDA = 1.2 nS (granule cell)
    static constexpr double Mg=0.2;        // astrocyte effective conc of extracellular Mg in mM

    double s;
    
    NMDA_astro( ) 
    {
       s=0;
    }

    double syn(double deltaT, double glu, double Vastro)
    {
        double Mg=0.02;   // conc of extracellular Mg in mM
        
        //Magnesium block B is an instantaneous function of voltage.
        double B =  1/( 1 + exp(-0.062 * Vastro) * Mg/3.57 );  
        
        double Inmda = gNMDA * s * B * (Vastro - Enmda);  // add Inmda to Vd_dt
        s = s + a_r * glu * (1 - s) - a_d * s;
        return Inmda;
    }
};




// from Garbo 2009.
//
class P2X_astro
{
     // Based on Angelo Di Garbo, Dynamics of ...,  2009
     
     public:
     static constexpr double ATP_ex  =3.0;   // extracellular conc of ATP
     static constexpr double A_ATP   =0.05;  // the corresponding amplitude: mu A/cm**2 / mu M
     static constexpr double HATP_P2X=0.9;   // half-saturation const for ATP-evoked ionotropic Ca2+, muM
     static constexpr double kATP_P2X=0.08;  // max rate of ATP-evoked ionotropic Ca influx, muM per sec
     static constexpr double ATP_life=1.0;   // ??
     
     int state;
     
     P2X_astro() {   
        state=0;
     }
     P2X_astro(int n) {   
        state=n;
     }
        
        
    
    // Depends on ATP release due to pre-synaptic APs.
    // t:   current time
    // tdr: time of evoked vesicle release
    double syn(double t, double tdr) 
    {
        double I_ATP=0;
        
        if  ( (t - tdr) < ATP_life ) {
         
            I_ATP = A_ATP * ATP_ex;   //   Garbo 2009, p 367, BUT was for interneuron NOT astrocyte.    FIX???
        }
        else {
            I_ATP=0;
        }
        return  I_ATP;  
    }

    
    // returns rate of influx of extracellular Ca in mu M/s
    //
    double nuATP_P2X(double t, double tdr ) 
    {    
    
        //   HATP_P2X mu M, Half saturation constant for ATP evoked ionotropic 
        //                  calcium influx amplitude.

        double nuATP=0, G_ATP_ex=0;
        
        if  ( (t - tdr) < ATP_life ) {

            G_ATP_ex = pow(ATP_ex,1.4) / ( HATP_P2X + pow(ATP_ex,1.4)  );
            nuATP=kATP_P2X * G_ATP_ex;
        }

        return nuATP;
    }
};
        
// from Garbo 2009.
//
class P2Y_astro
{

    double nuPLC_beta(double k_ATP_P2Y, double XF)
    {                                     
        double nuPLC_b = k_ATP_P2Y * XF;
        return nuPLC_b;  
    }
}; 