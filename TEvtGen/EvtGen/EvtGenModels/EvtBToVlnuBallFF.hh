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
// Description:   B->Xu l nu with the Ball/Zwicky decay model
//                Xu is a vector (rho, rho0, omega)
//
// Modification history:
//
//    Wells Wulsin      2008 Aug 14         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOVLNUBALLFF_HH
#define EVTBTOVLNUBALLFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"

class EvtId;

class EvtBToVlnuBallFF : public EvtSemiLeptonicFF {
  
public:
  EvtBToVlnuBallFF(double    r2_A1, 
		   double mfit2_A1,
		   double    r1_A2,
		   double    r2_A2,
		   double mfit2_A2,
		   double    r1_V,
		   double    r2_V,
		   double mfit2_V);
  
  void getvectorff(EvtId parent,EvtId daught,
		   double t, double mass, double *a1f,
		   double *a2f, double *vf, double *a0f );


  void getscalarff(EvtId, EvtId, double, double, double*, 
		   double*);

  void gettensorff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getbaryonff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getdiracff(EvtId, EvtId, double, double, double*, double*,
                  double*, double*, double*, double*);

  void getraritaff(EvtId, EvtId, double, double, double*, double*, 
		   double*, double*, double*, double*, double*, double*);

private:
  
  double    _r2_A1; 
  double _mfit2_A1;
  double    _r1_A2;
  double    _r2_A2;
  double _mfit2_A2;
  double    _r1_V;
  double    _r2_V;
  double _mfit2_V;

};

#endif

