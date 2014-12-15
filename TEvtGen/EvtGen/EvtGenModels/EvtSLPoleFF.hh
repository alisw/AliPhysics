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
// Module: EvtGen/EvtSLPoleFF.hh
//
// Description:Form factor routines for EvtSLPole
//
// Modification history:
//
//    DJL     April 23, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSLPOLEFF_HH
#define EVTSLPOLEFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtId.hh"

class EvtSLPoleFF : public EvtSemiLeptonicFF {

public:
  EvtSLPoleFF(int numarg, double *arglist);
  void getscalarff(EvtId parent,EvtId daught,
                       double t, double mass, double *fpf,
                       double *f0f );
  void getvectorff(EvtId parent,EvtId daught,
                       double t, double mass, double *a1f,
                       double *a2f, double *vf, double *a0f );
  void gettensorff(EvtId parent,EvtId daught,
                       double t, double mass, double *hf,
                       double *kf, double *bp, double *bm );

  void getbaryonff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getdiracff(EvtId, EvtId, double, double, double*, double*,
                  double*, double*, double*, double*);

  void getraritaff(EvtId, EvtId, double, double, double*, double*, 
		   double*, double*, double*, double*, double*, double*);

private:
   int numSLPoleargs;
   double SLPoleargs[16];

};

#endif



