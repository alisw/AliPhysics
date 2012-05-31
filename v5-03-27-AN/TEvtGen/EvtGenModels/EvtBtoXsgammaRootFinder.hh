//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2001      Brunel University
//
// Module: EvtGen/EvtBtoXsgammaRootFinder.hh
//
// Description:
//       Root finding algorithms using the bilear method. Basic structure
//       lifted from the BaBar IntegrationUtils root finding algorithm 
//       (author John Back).
//
// Modification history:
//
//       Jane Tinslay     March 21, 2001         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOXSGAMMAROOTFINDER_HH
#define EVTBTOXSGAMMAROOTFINDER_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//#include "EvtGenBase/EvtItgAbsFunction.hh"
class EvtItgAbsFunction;

class EvtBtoXsgammaRootFinder{

public:

  // Constructors

  EvtBtoXsgammaRootFinder();
  
  // Destructor
  virtual ~EvtBtoXsgammaRootFinder( );
 
  double GetRootSingleFunc(const EvtItgAbsFunction* theFunc, double functionValue, 
			   double lowerValue, double upperValue, double precision);

  double GetGaussIntegFcnRoot(EvtItgAbsFunction *theFunc1, EvtItgAbsFunction *theFunc2, 
			      double integ1Precision, double integ2Precision, 
			      int maxLoop1, int maxLoop2, double integLower, 
			      double integUpper, double lowerValue, double upperValue, 
			      double precision);

private:
  
  EvtBtoXsgammaRootFinder( const EvtBtoXsgammaRootFinder& );                // Copy Constructor
  EvtBtoXsgammaRootFinder& operator= ( const EvtBtoXsgammaRootFinder& );    // Assignment op
  
};

#endif


