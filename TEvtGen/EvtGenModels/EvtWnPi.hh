#ifndef EvtWNPI_HH
#define EvtWNPI_HH

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"

class EvtWnPi {
public:
  EvtVector4C WCurrent(EvtVector4R q1);
  EvtVector4C WCurrent(EvtVector4R q1, EvtVector4R q2);
  EvtVector4C WCurrent(EvtVector4R q1, EvtVector4R q2, EvtVector4R q3);
  EvtVector4C WCurrent(EvtVector4R q1, EvtVector4R q2, EvtVector4R q3, EvtVector4R q4, EvtVector4R q5);
protected:
  EvtVector4C JB(EvtVector4R q1, EvtVector4R q2, EvtVector4R q3, EvtVector4R q4, EvtVector4R q5); 
  EvtComplex BWa( EvtVector4R q);
  EvtComplex BWf( EvtVector4R q);
  EvtComplex BWr( EvtVector4R q);
  double pi3G(double Q2);
};

#endif
