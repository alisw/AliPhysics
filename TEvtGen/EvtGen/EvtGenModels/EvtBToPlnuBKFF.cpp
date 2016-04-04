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
// Module: EvtGenModels/EvtBToPlnuBK.cc
//
// Description: B->Xu l nu with BK (Becirevic-Kaidalov) parametrization
//              Xu is a pseudoscalar (pi_plus,pi0,eta or eta_prime)
//
// Modification history:
//
//    Martin Simard, U. de Montreal, 08/01/2007    Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBToPlnuBKFF.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtPDL.hh"
#include <math.h>
#include <stdlib.h>

EvtBToPlnuBKFF::EvtBToPlnuBKFF(double alpha, double beta) {

  _alpha = alpha;
  _beta = beta;
  
  return;
}


void EvtBToPlnuBKFF::getscalarff(EvtId parent,EvtId /*daught*/,
			    double t, double /*mass*/, double *fp, double *f0) {

  //Define mBstar
  EvtId Bplus = EvtPDL::getId("B+");
  EvtId Bminus = EvtPDL::getId("B-");
  double mBstar = EvtPDL::getMeanMass(EvtPDL::getId("B*0"));
  if(parent==Bplus || parent==Bminus) mBstar = EvtPDL::getMeanMass(EvtPDL::getId("B*+")); 
  double mBstar2=mBstar*mBstar;
  
  //Compute BK parametrization (t==q2)
  double fplus=1.0/((1.0-t/mBstar2)*(1.0-_alpha*t/mBstar2));
  double fzero=1.0/(1.0-t/(mBstar2*_beta));

  *fp=fplus;
  *f0=fzero;
  
  return;
}


void EvtBToPlnuBKFF::getvectorff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getvectorff in EvtBToPlnuBKFF.\n";  
  ::abort();

}



void EvtBToPlnuBKFF::gettensorff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :gettensorff in EvtBToPlnuBKFf.\n";  
  ::abort();

}



void EvtBToPlnuBKFF::getbaryonff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtBToPlnuBKFF.\n";  
  ::abort();

}

void EvtBToPlnuBKFF::getdiracff(EvtId, EvtId, double, double, double*, double*,
				double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtBToPlnuBKFF.\n";
  ::abort();

}

void EvtBToPlnuBKFF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
				 double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtBToPlnuBKFF.\n";
  ::abort();

}

