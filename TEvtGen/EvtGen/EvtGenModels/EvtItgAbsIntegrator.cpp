//--------------------------------------------------------------------------
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtItgIntegrator.cc
//
// Description:
//      Simpson integrator (Stolen and modified from 
//      the BaBar IntegrationUtils package - author: Phil Strother).
//
// Modification history:
//
//    Jane Tinslay                March 21, 2001       Module adapted for use in 
//                                                     EvtGen
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenModels/EvtItgAbsIntegrator.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

#include <math.h>
#include <iostream>

#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtItgAbsFunction.hh"
using std::endl;


EvtItgAbsIntegrator::EvtItgAbsIntegrator(const EvtItgAbsFunction &theFunction):
  _myFunction(theFunction)
{}

EvtItgAbsIntegrator::~EvtItgAbsIntegrator()
{}
  
double 
EvtItgAbsIntegrator::normalisation() const {
  return evaluateIt(_myFunction.lowerRange(), _myFunction.upperRange());
}

double
EvtItgAbsIntegrator::evaluate(double lower, double upper) const{

  double newLower(lower), newUpper(upper);

  boundsCheck(newLower, newUpper);

  return evaluateIt(newLower, newUpper);
}

double 
EvtItgAbsIntegrator::trapezoid(double lower, double higher, int n, double &result) const {

  if (n==1) return 0.5*(higher-lower)*(_myFunction(lower) + _myFunction(higher));
  
  int it, j;
  
  for (it=1, j=1;j<n-1;j++) it <<=1;
  
  double itDouble(it);
  
  double sum(0.0);

  double deltaX((higher - lower)/itDouble);
  
  double x(lower + 0.5* deltaX);
    
  for (j=1;j<=it;j++){
    sum+=_myFunction(x);
    x+=deltaX;
  }
  
  result = 0.5*(result+(higher - lower)*sum/itDouble);

  return result;
}

void 
EvtItgAbsIntegrator::boundsCheck(double &lower, double &upper) const{

  if (lower < _myFunction.lowerRange() ) {
    report(Severity::Warning,"EvtGen") << "Warning in EvtItgAbsIntegrator::evaluate.  Lower bound " << lower << " of integral " 
		    << " is less than lower bound " << _myFunction.lowerRange() 
		    << " of function.  No contribution from this range will be counted." << endl;
    lower = _myFunction.lowerRange();
  }

  if (upper > _myFunction.upperRange() ) {
    report(Severity::Warning,"EvtGen") << "Warning in EvtItgAbsIntegrator::evaluate.  Upper bound " << upper << " of integral "
		    << " is greater than upper bound " << _myFunction.upperRange() 
		    << " of function.  No contribution from this range will be counted." << endl;  
    upper = _myFunction.upperRange();
  }

}
