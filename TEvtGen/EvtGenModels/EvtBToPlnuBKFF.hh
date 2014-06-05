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
// Module: EvtGenModels/EvtBToPlnuBK.hh
//
// Description:   B->Xu l nu with BK (Becirevic-Kaidalov) parametrization
//                Xu is a pseudoscalar (pi_plus,pi0,eta or eta_prime)
//
// Modification history:
//
//    Martin Simard, U. de Montreal, 08/01/2007    Module created
//
//------------------------------------------------------------------------
#ifndef EVTBTOPLNUBKFF_HH
#define EVTBTOPLNUBKFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

class EvtBToPlnuBKFF : public EvtSemiLeptonicFF {

public:
  EvtBToPlnuBKFF(double alpha, double beta);
  
  void getscalarff(EvtId parent,EvtId daught,
		   double t, double mass, double *fp, double *f0);

  void getvectorff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void gettensorff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getbaryonff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getdiracff(EvtId, EvtId, double, double, double*, double*,
                  double*, double*, double*, double*);

  void getraritaff(EvtId, EvtId, double, double, double*, double*, 
		   double*, double*, double*, double*, double*, double*);

private:
  double _alpha;
  double _beta;
};

#endif

