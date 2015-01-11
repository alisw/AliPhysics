//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtResonance.hh
//
// Description:resonance-defining class
//
// Modification history:
//
//    NK     September 4, 1997         Module created
//
//------------------------------------------------------------------------

#ifndef EVTRESONANCE_HH
#define EVTRESONANCE_HH

#include "EvtGenBase/EvtVector4R.hh"

class EvtComplex;



class EvtResonance {
public:

    EvtResonance& operator = (const EvtResonance &);
    
    //constructor with all information about the resonance
    EvtResonance(const EvtVector4R& p4_p, const EvtVector4R& p4_d1, 
		 const EvtVector4R& p4_d2, 
		 double ampl = 0.0, double theta = 0.0, double gamma = 0.0, 
		 double bwm = 0.0, int spin = 0);

    //destructor
    virtual ~EvtResonance();
    
    //accessors
    //return 4-momenta of the particles involved
    inline const EvtVector4R& p4_p() { return _p4_p; }
    inline const EvtVector4R& p4_d1() { return _p4_d1; }
    inline const EvtVector4R& p4_d2() { return _p4_d2; }  
    

     //return amplitude
    inline double amplitude() { return _ampl; }  

    //return theta
    inline double theta() { return _theta; } 

    //return gamma
    inline double gamma() { return _gamma; } 

    //return bwm
    inline double bwm() { return _bwm; } 

    //return spin
    inline int spin() { return _spin; } 

    //calculate amplitude for this resonance
    EvtComplex resAmpl();
    
    //calculate relativistic Breit-Wigner amplitude for P-decays of scalars
    EvtComplex relBrWig(int i);
   
private:

    EvtVector4R _p4_p, _p4_d1, _p4_d2;
    double _ampl, _theta, _gamma, _bwm;
    int _spin;

}; 

#endif

