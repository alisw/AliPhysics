//--------------------------------------------------------------------------
//
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtItgSimpsonIntegrator.hh
//
// Description:
//      Abstraction of a generic function for use in integration methods elsewhere
//      in this package. (Stolen and modified from 
//      the BaBar IntegrationUtils package - author: Phil Strother).
//
// Modification history:
//
//    Jane Tinslay                March 21, 2001       Module adapted for use in 
//                                                     EvtGen
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenModels/EvtItgSimpsonIntegrator.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

#include <math.h>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "EvtGenModels/EvtItgAbsFunction.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;


EvtItgSimpsonIntegrator::EvtItgSimpsonIntegrator(const EvtItgAbsFunction &theFunction, double precision, int maxLoop):
  EvtItgAbsIntegrator(theFunction),
  _precision(precision),
  _maxLoop(maxLoop)
{}


//--------------
// Destructor --
//--------------

EvtItgSimpsonIntegrator::~EvtItgSimpsonIntegrator()
{}

double
EvtItgSimpsonIntegrator::evaluateIt(double lower, double higher) const{
  
  // report(Severity::Info,"EvtGen")<<"in evaluate"<<endl;
  int j;
  double result(0.0);
  double s, st, ost(0.0);
  for (j=1;j<4;j++) {
    st = trapezoid(lower, higher, j, result);
    s = (4.0 * st - ost)/3.0;
    ost=st;
  }

  double olds(s);
  st = trapezoid(lower, higher, j, result);
  s = (4.0 * st - ost)/3.0;

  if (fabs(s - olds) < _precision*fabs(olds) || (s==0.0 && olds==0.0))     return s;
  
  ost=st;

  for (j=5;j<_maxLoop;j++){

    st = trapezoid(lower, higher, j, result);
    s = (4.0 * st - ost)/3.0;
    
    if (fabs(s - olds) < _precision*fabs(olds) || (s==0.0 && olds==0.0))    return s;
    olds=s;
    ost=st;
  }
  
  report(Severity::Error,"EvtGen") << "Severe error in EvtItgSimpsonIntegrator.  Failed to converge after loop with 2**"
		 << _maxLoop << " calls to the integrand in." << endl;
  
  return 0.0;
    
}
