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
// Module: EvtGen/EvtISGWFF.hh
//
// Description:Form factor routines specific to EvtISGW
//
// Modification history:
//
//    DJL/RYD     September 25, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTISGWFF_HH
#define EVTISGWFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;


class EvtISGWFF : public EvtSemiLeptonicFF {

  void getscalarff( EvtId parent, EvtId daught,
                       double t, double mass, double *fpf,
                       double *f0f );
  void getvectorff( EvtId parent, EvtId daught,
                       double t, double mass, double *a1f,
                       double *a2f, double *vf, double *a0f );
  void gettensorff( EvtId parent, EvtId daught,
                       double t, double mass, double *hf,
                       double *kf, double *bpf, double *bmf );

  void getbaryonff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getdiracff(EvtId, EvtId, double, double, double*, double*,
                  double*, double*, double*, double*);

  void getraritaff(EvtId, EvtId, double, double, double*, double*, 
		   double*, double*, double*, double*, double*, double*);

    // getscalarff, getvectorff, and gettensorff call the
    // correct isgw form factor routine which computes 
    // form factors according to the ISGW paper.


  void EvtISGW1FF3S1( EvtId parent, EvtId daught, 
                       double t, double mass, double *ff, double *gf,
                       double *apf, double *amf);
  void EvtISGW1FF23S1( EvtId parent, EvtId daught,
                        double t, double mass,double *fpf, double *gpf,
                        double *app, double *apm);
  void EvtISGW1FF3P1( EvtId parent, EvtId daught, 
                       double t, double mass,double *lf, double *qf,
                      double *cpf, double *cmf);
  void EvtISGW1FF3P0( EvtId parent, EvtId daught, 
                       double t, double mass, double *upf, double *umf);
  void EvtISGW1FF1S0( EvtId parent, EvtId daught, 
		       double t, double mass,double *fpf, double *fmf);
  void EvtISGW1FF21S0( EvtId parent, EvtId daught,
                        double t, double mass, double *fppf, double *fpmf);
  void EvtISGW1FF3P2( EvtId parent, EvtId daught, 
                       double t, double mass, double *h, double *k,
                       double *bp, double *bm);
  void EvtISGW1FF1P1( EvtId parent, EvtId daught, 
                       double t, double mass, double *vf, double *rf,
                       double *spf, double *smf);

};

#endif

