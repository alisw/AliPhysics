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
// Module: EvtGen/EvtSLBKPoleFF.hh
//
// Description:Form factor routines for EvtSLBKPole,
//             according to Becirevic and Kaidalov(BK)
//
// Modification history:
//
//    liheng     October 20, 2005         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSLBKPOLEFF_HH//modified
#define EVTSLBKPOLEFF_HH//modified

#include "EvtGenBase/EvtSemiLeptonicFF.hh"//modified
#include "EvtGenBase/EvtId.hh"

class EvtSLBKPoleFF : public EvtSemiLeptonicFF {//modified

public:
  EvtSLBKPoleFF(int numarg, double *arglist);//modified
  void getscalarff(EvtId parent,EvtId daught,
                       double t, double mass, double *fpf,
                       double *f0f );
  void getvectorff(EvtId parent,EvtId daught,
       		       double t, double mass, double *a1f,
                       double *a2f, double *vf, double *a0f );
  void gettensorff(EvtId parent,EvtId daught,//need to be modified, but not yet
                       double t, double mass, double *hf,
                       double *kf, double *bp, double *bm );

  void getbaryonff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getdiracff(EvtId, EvtId, double, double, double*, double*,
                  double*, double*, double*, double*);

  void getraritaff(EvtId, EvtId, double, double, double*, double*, 
		   double*, double*, double*, double*, double*, double*);

private:
   int numSLBKPoleargs;//modified
   double SLBKPoleargs[16];//modified

};

#endif



