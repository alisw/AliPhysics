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
// Module: EvtGen/EvtHQETFF.hh
//
// Description:
//
// Modification history:
//
//    DJL     April 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTHQETFF_HH
#define EVTHQETFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

class EvtHQETFF : public EvtSemiLeptonicFF {

public:
  EvtHQETFF(double hqetrho2, double hqetr1, double hqetr2, double hqetc=0.);
  EvtHQETFF(double hqetrho2,  double hqetc=0.);
  void getvectorff(EvtId parent,EvtId daught,
                       double t, double mass, double *a1f,
                       double *a2f, double *vf, double *a0f );

  void getscalarff(EvtId parent,EvtId daught,
		   double t, double mass, double *f0p, double *f0m);

  void gettensorff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getbaryonff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getdiracff(EvtId, EvtId, double, double, double*, double*,
                  double*, double*, double*, double*);

  void getraritaff(EvtId, EvtId, double, double, double*, double*, 
		   double*, double*, double*, double*, double*, double*);

private:
  double r1;
  double rho2;
  double r2;
  double c;
};

#endif

