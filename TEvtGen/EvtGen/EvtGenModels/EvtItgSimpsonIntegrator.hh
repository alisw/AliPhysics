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
//      Simpson integrator (Stolen and modified from 
//      the BaBar IntegrationUtils package - author: Phil Strother).
//
// Modification history:
//
//    Jane Tinslay                March 21, 2001       Module adapted for use in 
//                                                     EvtGen
//
//------------------------------------------------------------------------

#ifndef EVTITGSIMPSONINTEGRATOR_HH
#define EVTITGSIMPSONINTEGRATOR_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

#include "EvtGenModels/EvtItgAbsIntegrator.hh"

class EvtItgSimpsonIntegrator: public EvtItgAbsIntegrator {

public:
  
  EvtItgSimpsonIntegrator(const EvtItgAbsFunction &, double precision=1.0e-5, int maxLoop=20);

  virtual ~EvtItgSimpsonIntegrator( );
  
protected:
  
  virtual double evaluateIt(double , double) const;
  
private:
  
  double _precision;
  double _maxLoop;

  EvtItgSimpsonIntegrator();
  EvtItgSimpsonIntegrator( const EvtItgSimpsonIntegrator& );                //// Copy Constructor
  EvtItgSimpsonIntegrator& operator= ( const EvtItgSimpsonIntegrator& );    // Assignment op
  
};



#endif // ITGSIMPSONINTEGRATOR_HH
