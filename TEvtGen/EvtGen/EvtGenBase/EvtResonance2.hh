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
// Module: EvtGen/EvtResonance2.hh
//
// Description:resonance-defining class
//
// Modification history:
//
//    lange   Nov 21, 2000         Module created
//
//------------------------------------------------------------------------

#ifndef EVTRESONANCE2_HH
#define EVTRESONANCE2_HH

#include "EvtGenBase/EvtVector4R.hh"

class EvtComplex;



class EvtResonance2 {
public:

    //operator
    EvtResonance2& operator = (const EvtResonance2 &);
    
    //constructor with all information about the resonance
    // invmass_angdenom chooses whether to use the resonance mass (false)
    // or the daughter invariant mass (true) for the denominators in
    // angular distributions

    EvtResonance2(const EvtVector4R& p4_p, const EvtVector4R& p4_d1, 
		  const EvtVector4R& p4_d2, 
		  double ampl = 0.0, double theta = 0.0, double gamma = 0.0, 
		  double bwm = 0.0, int spin = 0, bool invmass_angdenom = false);

    //destructor
    virtual ~EvtResonance2();

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
   
  private:

    EvtVector4R _p4_p, _p4_d1, _p4_d2;
    double _ampl, _theta, _gamma, _bwm;
    int _spin;
    bool _invmass_angdenom;

}; 

#endif

