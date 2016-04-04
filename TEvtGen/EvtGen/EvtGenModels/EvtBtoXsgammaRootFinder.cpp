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
// Module: EvtBtoXsgammaRootFinder.cc
//
// Description:
//      Root finders for EvtBtoXsgammaKagan module.
//
// Modification history:
//
//      Jane Tinslay       March 21, 2001       Module created
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenModels/EvtBtoXsgammaRootFinder.hh"
#include "EvtGenModels/EvtItgTwoCoeffFcn.hh"
#include "EvtGenModels/EvtItgSimpsonIntegrator.hh"
#include "EvtGenBase/EvtReport.hh"
#include <math.h>
using std::endl;

//-------------
// C Headers --
//-------------
extern "C" {
}

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

#define EVTITGROOTFINDER_MAXIT 100
#define EVTITGROOTFINDER_RELATIVEPRECISION 1.0e-16


EvtBtoXsgammaRootFinder::EvtBtoXsgammaRootFinder() {}

EvtBtoXsgammaRootFinder::~EvtBtoXsgammaRootFinder( )
{}

double
EvtBtoXsgammaRootFinder::GetRootSingleFunc(const EvtItgAbsFunction* theFunc, double functionValue, double lowerValue, double upperValue, double precision) {
  
  // Use the bisection to find the root.
  // Iterates until find root to the accuracy of precision

  double xLower = 0.0, xUpper = 0.0;
  double root=0;

  double f1 = theFunc->value(lowerValue) - functionValue;
  double f2 = theFunc->value(upperValue) - functionValue;

  if ( f1*f2 > 0.0 ) {
    report(Severity::Warning,"EvtGen") << "EvtBtoXsgammaRootFinder: No root in specified range !"<<endl;  
    return 0;
  }

  // Already have root
  if (fabs(f1) < precision) {
    root = lowerValue;
    return root;
  }
  if (fabs(f2) < precision) {
    root = upperValue;
    return root;
  }
  
  // Orient search so that f(xLower) < 0
  if (f1 < 0.0) {
    xLower = lowerValue;
    xUpper = upperValue;
  } else {
    xLower = upperValue;
    xUpper = lowerValue;
  }
  
  double rootGuess = 0.5*(lowerValue + upperValue);
  double dxold = fabs(upperValue - lowerValue);
  double dx = dxold;
  
  double f = theFunc->value(rootGuess) - functionValue;
  
  for (int j = 0; j< EVTITGROOTFINDER_MAXIT; j++) {
    
      dxold = dx;
      dx = 0.5*(xUpper-xLower);
      rootGuess = xLower+dx;

      // If change in root is negligible, take it as solution.
      if (fabs(xLower - rootGuess) < precision) {
	root = rootGuess;
	return root;
      }
      
      f = theFunc->value(rootGuess) - functionValue;
 
      if (f < 0.0) {
	xLower = rootGuess;
      } else {
	xUpper = rootGuess;
      }
      
  }
  
  report(Severity::Warning,"EvtGen") << "EvtBtoXsgammaRootFinder: Maximum number of iterations "
			   <<"in EvtBtoXsgammaRootFinder::foundRoot exceeded!" 
			   <<" Returning false."<<endl;
  return 0;
  
}

double
EvtBtoXsgammaRootFinder::GetGaussIntegFcnRoot(EvtItgAbsFunction *theFunc1, EvtItgAbsFunction *theFunc2, double integ1Precision, double integ2Precision, int maxLoop1, int maxLoop2, double integLower, double integUpper, double lowerValue, double upperValue, double precision) {
 
  // Use the bisection to find the root.
  // Iterates until find root to the accuracy of precision
  
  //Need to work with integrators
  EvtItgAbsIntegrator *func1Integ = new EvtItgSimpsonIntegrator(*theFunc1, integ1Precision, maxLoop1);
  EvtItgAbsIntegrator *func2Integ = new EvtItgSimpsonIntegrator(*theFunc2, integ2Precision, maxLoop2);
  
  
  //coefficient 1 of the integrators is the root to be found
  //need to set this to lower value to start off with
  theFunc1->setCoeff(1,0,lowerValue);
  theFunc2->setCoeff(1,0,lowerValue);
  
  double f1 = func1Integ->evaluate(integLower,integUpper) - theFunc2->getCoeff(1,2)*func2Integ->evaluate(integLower,integUpper);
  theFunc1->setCoeff(1,0,upperValue);
  theFunc2->setCoeff(1,0,upperValue);
  double f2 = func1Integ->evaluate(integLower,integUpper) - theFunc2->getCoeff(1,2)*func2Integ->evaluate(integLower,integUpper);
  
  double xLower = 0.0, xUpper = 0.0;
  double root=0;

  if ( f1*f2 > 0.0 ) {
    report(Severity::Warning,"EvtGen") << "EvtBtoXsgammaRootFinder: No root in specified range !"<<endl;  
    return false;
  }

  // Already have root
  if (fabs(f1) < precision) {
    root = lowerValue;
    return root;
  }
  if (fabs(f2) < precision) {
    root = upperValue;
    return root;
  }
  
  // Orient search so that f(xLower) < 0
  if (f1 < 0.0) {
    xLower = lowerValue;
    xUpper = upperValue;
  } else {
    xLower = upperValue;
    xUpper = lowerValue;
  }
  
  double rootGuess = 0.5*(lowerValue + upperValue);
  double dxold = fabs(upperValue - lowerValue);
  double dx = dxold;
  
  theFunc1->setCoeff(1,0,rootGuess);
  theFunc2->setCoeff(1,0,rootGuess);
  double f = func1Integ->evaluate(integLower,integUpper) - theFunc2->getCoeff(1,2)*func2Integ->evaluate(integLower,integUpper);
  
  for (int j = 0; j< EVTITGROOTFINDER_MAXIT; j++) {
    
    dxold = dx;
    dx = 0.5*(xUpper-xLower);
    rootGuess = xLower+dx;
    
    // If change in root is negligible, take it as solution.
    if (fabs(xLower - rootGuess) < precision) {
      root = rootGuess;
      return root;
    }
    
    theFunc1->setCoeff(1,0,rootGuess);
    theFunc2->setCoeff(1,0,rootGuess);
    f = func1Integ->evaluate(integLower,integUpper) - theFunc2->getCoeff(1,2)*func2Integ->evaluate(integLower,integUpper);
    
    if (f < 0.0) {
      xLower = rootGuess;
    } else {
      xUpper = rootGuess;
    }
    
  }
  
  report(Severity::Warning,"EvtGen") << "EvtBtoXsgammaRootFinder: Maximum number of iterations "
			   <<"in EvtBtoXsgammaRootFinder::foundRoot exceeded!" 
			   <<" Returning false."<<endl;
  return 0;
  
}


