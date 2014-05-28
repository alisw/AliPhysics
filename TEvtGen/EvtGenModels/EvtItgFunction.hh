//--------------------------------------------------------------------------
//
// Environment:
//      This software was developed for the BaBar collaboration.  If you
//      use all or part of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      LBNL
//
//------------------------------------------------------------------------

#ifndef EVTITGFUNCTION_HH
#define EVTITGFUNCTION_HH

#include "EvtGenModels/EvtItgAbsFunction.hh"

/**
 *  Copyright (C) 1998 LBNL
 *  
 *  Generic function where the pointer to the function is available.
 *
 *  The function is taken as type pointer to function returning double and 
 *  taking a double (the abscissa) and a const RWTValVector<double> reference
 *  (the parameter values of the function) as arguments.
 *
 *  @see EvtItgFunctionEvtItgFunction
 *
 *  @version $Id: EvtItgFunction.hh,v 1.2 2009-03-16 16:34:00 robbep Exp $ 
 *
 *  @author Phil Strother       Originator
 */

class EvtItgFunction: public EvtItgAbsFunction {

public:

  // Constructors
  EvtItgFunction( double (*theFunction)(double),
		     double lowerRange, double upperRange);
 
 
  // Destructor
  virtual ~EvtItgFunction( );

  virtual void setCoeff(int, int, double) {};
  virtual double getCoeff(int, int) {return 0.0;};
 
protected:
  
  // Helper functions

  virtual double myFunction(double x) const;

private:
 
  // Data members
  double (*_myFunction)(double x);

  // Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
   EvtItgFunction( const EvtItgFunction& );                // Copy Constructor
  EvtItgFunction& operator= ( const EvtItgFunction& );    // Assignment op
};

#endif // EvtITGFUNCTION_HH
