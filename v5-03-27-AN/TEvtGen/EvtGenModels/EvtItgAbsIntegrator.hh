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
// Module: EvtItgAbsIntegrator.hh
//
// Description:
//      Abstraction of a generic integrator (Stolen and modified from 
//      the BaBar IntegrationUtils package - author: Phil Strother).
//
// Modification history:
//
//    Jane Tinslay                March 21, 2001       Module adapted for use in 
//                                                     EvtGen
//
//------------------------------------------------------------------------

#ifndef EVTITGABSINTEGRATOR_HH
#define EVTITGABSINTEGRATOR_HH


#include "EvtGenModels/EvtItgAbsFunction.hh"

class EvtItgAbsIntegrator {

public:
  
  EvtItgAbsIntegrator(const EvtItgAbsFunction &);
  
  virtual ~EvtItgAbsIntegrator( );

  double evaluate(double lower, double upper) const;
 
  double normalisation() const;

protected:

   double trapezoid(double lower, double higher, int n, 
		   double &result) const;
  
  virtual double evaluateIt(double lower, double higher) const=0;
  
  double myFunction(double x) const {return _myFunction(x);}
 
private:
  
  const EvtItgAbsFunction &_myFunction;

  void boundsCheck(double &, double &) const;

  // Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  EvtItgAbsIntegrator();
  EvtItgAbsIntegrator( const EvtItgAbsIntegrator& );                // Copy Constructor
  EvtItgAbsIntegrator& operator= ( const EvtItgAbsIntegrator& );    // Assignment op
 
};

#endif // EVTITGABSINTEGRATOR_HH
