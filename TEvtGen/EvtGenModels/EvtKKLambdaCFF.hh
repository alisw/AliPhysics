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

#ifndef EVTKKLAMBDACFF_HH
#define EVTKKLAMBDACFF_HH

#include "EvtGenBase/EvtSemiLeptonicFF.hh"
#include "EvtGenBase/EvtId.hh"

class EvtKKLambdaCFF : public EvtSemiLeptonicFF {

public:
  EvtKKLambdaCFF(int numarg, double *arglist);

  void getscalarff(EvtId parent,EvtId daught,
		   double t, double mass, double *f0p, double *f0m);
  
  void getvectorff(EvtId, EvtId,
		   double, double, double*,
		   double*, double*, double* );
  
  void gettensorff(EvtId, EvtId, double, double, double*, 
		   double*, double*, double*);

  void getbaryonff(EvtId parent,EvtId daught,
		   double t, double mass, double *f1v,
		   double *f1a, double *f2v, double *f2a );

  void getdiracff(EvtId, EvtId, double, double, double*, double*,
                  double*, double*, double*, double*);

  void getraritaff(EvtId, EvtId, double, double, double*, double*, 
		   double*, double*, double*, double*, double*, double*);

private:
   int _nargs;
   double _args[2];

};

#endif



