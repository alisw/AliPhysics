//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtSemiLeptonicFF.hh
// This is the base class for semileptonic form factor calculations.
//
// Description:
//
// Modification history:
//
//    DJL/RYD     September 30, 1997         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSEMILEPTONICFF_HH
#define EVTSEMILEPTONICFF_HH

#include "EvtGenBase/EvtId.hh"

class EvtSemiLeptonicFF{

public:
  virtual ~EvtSemiLeptonicFF( ) { } ;

  virtual void getscalarff( EvtId parent, EvtId daught,
			    double t, double mass, double *fpf,
			    double *fmf ) = 0;
  virtual void getvectorff( EvtId parent, EvtId daught,
			    double t, double mass, double *a1f,
			    double *a2f, double *vf, double *a0f ) = 0;
  virtual void gettensorff( EvtId parent, EvtId daught,
			    double t, double mass, double *a1f,
			    double *a2f, double *vf, double *a0f ) = 0;

  virtual void getbaryonff(EvtId parent, EvtId daught,
			   double t,
                           double m_meson,
                           double *f1v, 
                           double *f1a, 
                           double *f2v, 
			   double *f2a) = 0;

  virtual void getdiracff( EvtId parent, EvtId daught,
			   double q2, double mass, 
			   double *f1, double *f2, double *f3,
			   double *g1, double *g2, double *g3 ) = 0;

  virtual void getraritaff( EvtId parent, EvtId daught,
			    double q2, double mass, 
			    double *f1, double *f2, 
			    double *f3, double *f4,
			    double *g1, double *g2, 
			    double *g3, double *g4 ) = 0;

};

#endif
