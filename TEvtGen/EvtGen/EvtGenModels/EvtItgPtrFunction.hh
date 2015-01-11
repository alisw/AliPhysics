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
// Module: EvtItgPtrFunction.hh
//
// Description:
//      Class describing a function with one vector of coefficients. (Stolen and 
//      modified from the BaBar IntegrationUtils package - author: Phil Strother).
//
// Modification history:
//
//    Jane Tinslay                March 21, 2001       Module adapted for use in 
//                                                     EvtGen
//
//------------------------------------------------------------------------

#ifndef EVTITGPTRFUNCTION_HH
#define EVTITGPTRFUNCTION_HH

#include <vector>
#include "EvtGenModels/EvtItgAbsFunction.hh"

class EvtItgPtrFunction: public EvtItgAbsFunction {

public:

  EvtItgPtrFunction( double (*theFunction)(double, const std::vector<double> &),
		     double lowerRange, double upperRange, const std::vector<double> &coeffs1);
 
  virtual ~EvtItgPtrFunction( );

  virtual void setCoeff(int, int, double);
  virtual double getCoeff(int, int);

protected:
  
  virtual double myFunction(double x) const;
 
private:
 
  // Data members
  double (*_myFunction)(double x, const std::vector<double> & coeffs1);

  // Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  EvtItgPtrFunction( const EvtItgPtrFunction& );                //// Copy Constructor
  EvtItgPtrFunction& operator= ( const EvtItgPtrFunction& );    // Assignment op
  std::vector<double> _coeffs1;

};

#endif // EVTITGPTRFUNCTION_HH
