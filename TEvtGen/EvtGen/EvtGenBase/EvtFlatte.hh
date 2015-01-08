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
// Module: EvtGen/EvtFlatte.hh
//
// Description:resonance-defining class
//
// Modification history:
//
//    ponyisi  18 Feb 2008  created
//
//------------------------------------------------------------------------

#ifndef EVTFLATTE_HH
#define EVTFLATTE_HH

#include "EvtGenBase/EvtVector4R.hh"
#include <vector>

using std::vector;

class EvtComplex;

// Helper class

class EvtFlatteParam {
public:
  EvtFlatteParam(double m1, double m2, double g):
    _m1(m1), _m2(m2), _g(g) {}
      
  inline double m1() const { return _m1; }
  inline double m2() const { return _m2; }
  inline double g() const { return _g; }

private:
  double _m1, _m2, _g;
};

//class declaration

class EvtFlatte {
public:

  //operator
  EvtFlatte& operator = (const EvtFlatte &);

  //constructor with all information about the resonance
  EvtFlatte(const EvtVector4R& p4_p, const EvtVector4R& p4_d1, 
	    const EvtVector4R& p4_d2, 
	    double ampl, double theta,
	    double mass, 
             vector<EvtFlatteParam>& params
	    //           double m1a = 0.0, double m1b = 0.0, double g1 = 0.0,
	    //           double m2a = 0.0, double m2b = 0.0, double g2 = 0.0
	    );

  //destructor
  virtual ~EvtFlatte();

  //accessors
  //return 4-momenta of the particles involved
  inline const EvtVector4R& p4_p() { return _p4_p; }
  inline const EvtVector4R& p4_d1() { return _p4_d1; }
  inline const EvtVector4R& p4_d2() { return _p4_d2; }  
    

  //return amplitude
  inline double amplitude() { return _ampl; }  

  //return theta
  inline double theta() { return _theta; } 

  //return bwm
  inline double mass() { return _mass; } 

  //functions

  //calculate amplitude for this resonance
  EvtComplex resAmpl();
   
private:

  inline EvtComplex sqrtCplx(double in) { return (in > 0) ? EvtComplex(sqrt(in), 0) : EvtComplex
					    (0, sqrt(-in)); }

  EvtVector4R _p4_p, _p4_d1, _p4_d2;
  double _ampl, _theta, _mass;
  vector<EvtFlatteParam> _params;
  //      double _m1a, _m1b, _g1;
  //      double _m2a, _m2b, _g2;
}; 

#endif

