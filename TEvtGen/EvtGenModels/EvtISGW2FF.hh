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
// Module: EvtGen/EvtISGW2FF.hh
//
// Description:Form factor routines specific to EvtISGW2
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTISGW2FF_HH
#define EVTISGW2FF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

class EvtISGW2FF : public EvtSemiLeptonicFF {

public:

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

private:

    // getscalarff, getvectorff, and gettensorff call the
    // correct isgw2 form factor routine which computes 
    // form factors according to the ISGW2 paper.


  void EvtISGW2FF3S1( EvtId parent, EvtId daught, 
                       double t, double mass, double *ff, double *gf,
                       double *apf, double *amf);
  void EvtISGW2FF23S1( EvtId parent, EvtId daught,
                        double t, double mass,double *fpf, double *gpf,
                        double *app, double *apm);
  void EvtISGW2FF3P1( EvtId parent, EvtId daught, 
                       double t, double mass,double *lf, double *qf,
                       double *cpf, double *cmf);
  void EvtISGW2FF3P0( EvtId parent, EvtId daught, 
                       double t, double mass, double *upf, double *umf);
  void EvtISGW2FF1S0( EvtId parent, EvtId daught, 
		       double t, double mass,double *fpf, double *fmf);
  void EvtISGW2FF21S0( EvtId parent, EvtId daught,
                        double t, double mass, double *fppf, double *fpmf);
  void EvtISGW2FF3P2( EvtId parent, EvtId daught, 
                       double t, double mass, double *h, double *k,
                       double *bp, double *bm);
  void EvtISGW2FF1P1( EvtId parent, EvtId daught, 
                       double t, double mass, double *rf, double *vf,
                       double *spf, double *smf);

  double EvtGetas( double mass );
  double EvtGetas( double mass,double mass1  );
  double EvtGetGammaji( double z );

};

#endif


