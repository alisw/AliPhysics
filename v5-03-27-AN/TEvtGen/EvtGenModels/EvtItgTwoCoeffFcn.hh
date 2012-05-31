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
// Module: EvtItgTwoCoeffFcn.hh
//
// Description:
//      Class describing a function with two vectors of coefficients. 
//
// Modification history:
//
//    Jane Tinslay                March 21, 2001       Module created
//
//------------------------------------------------------------------------

#ifndef EVTITTWOCOEFFFCN_HH
#define EVTITTWOCOEFFFCN_HH

#include <vector>

//-------------
// C Headers --
//-------------
extern "C" {
}

#include "EvtGenModels/EvtItgAbsFunction.hh"

class EvtItgTwoCoeffFcn: public EvtItgAbsFunction {

public:

  EvtItgTwoCoeffFcn( double (*theFunction)(double, const std::vector<double> &, const std::vector<double> &),
		     double lowerRange, double upperRange, const std::vector<double> &coeffs1, const std::vector<double> &coeffs2);

  virtual ~EvtItgTwoCoeffFcn( );

  virtual void setCoeff(int, int, double);
  virtual double getCoeff(int, int);


protected:
  
  virtual double myFunction(double x) const;

private:
 
  // Data members
  double (*_myFunction)(double x, const std::vector<double> & coeffs1, const std::vector<double> & coeffs2);

  // Note: if your class needs a copy constructor or an assignment operator, 
  // make one of the following public and implement it.
   EvtItgTwoCoeffFcn( const EvtItgTwoCoeffFcn& );                //// Copy Constructor
  EvtItgTwoCoeffFcn& operator= ( const EvtItgTwoCoeffFcn& );    // Assignment op
 
  std::vector<double> _coeffs1;
  std::vector<double> _coeffs2;

};

#endif // EvtITGTWOCOEFFFUNCTION_HH
