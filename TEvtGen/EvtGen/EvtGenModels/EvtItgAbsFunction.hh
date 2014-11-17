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
// Module: EvtItgAbsFunction.hh
//
// Description:
//      Abstraction of a generic function for use in integration methods elsewhere
//      in this package. (Stolen and modified from the BaBar IntegrationUtils package 
//      - author: Phil Strother).
//
// Modification history:
//
//    Jane Tinslay                March 21, 2001       Module adapted for use in 
//                                                     EvtGen
//
//------------------------------------------------------------------------

#ifndef EVTITGABSFUNCTION_HH
#define EVTITGABSFUNCTION_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

class EvtItgAbsFunction {

public:

  // Constructors

  EvtItgAbsFunction(double lowerRange, double upperRange);

  // Destructor
  virtual ~EvtItgAbsFunction( );

  virtual double value( double x) const;

  virtual double operator()(double x) const;
  
  // Selectors (const)
  
  inline double upperRange() const {return _upperRange;}
  inline double lowerRange() const {return _lowerRange;}
  inline void   getRange(double &lower,double &upper) const { lower = _lowerRange; upper = _upperRange; } 
  virtual void setCoeff(int, int, double)=0;
  virtual double getCoeff(int, int)=0;
  
protected:
  
  virtual double myFunction(double x) const=0;
  void    setRange(double x1,double x2) { _lowerRange=x1; _upperRange=x2; };

private:
  
  double _upperRange;
  double _lowerRange;
 
  EvtItgAbsFunction( const EvtItgAbsFunction& );                // Copy Constructor
  EvtItgAbsFunction& operator= ( const EvtItgAbsFunction& );    // Assignment op

};

#endif // EVTITGABSFUNCTION_HH
