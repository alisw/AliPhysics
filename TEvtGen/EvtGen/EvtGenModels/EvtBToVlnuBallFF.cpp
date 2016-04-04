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
// Module: EvtHQETFF.cc
//
// Description:   B->Xu l nu with the Ball/Zwicky decay model
//                Xu is a vector (rho, rho0, omega)
//
// Modification history:
//
//    Wells Wulsin      2008 Aug 14         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBToVlnuBallFF.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtPDL.hh"
#include <math.h>
#include <stdlib.h>

EvtBToVlnuBallFF::EvtBToVlnuBallFF(double    r2_A1, 
				   double mfit2_A1,
				   double    r1_A2,
				   double    r2_A2,
				   double mfit2_A2,
				   double    r1_V,
				   double    r2_V,
				   double mfit2_V) {

  _r2_A1    = r2_A1;
  _mfit2_A1 = mfit2_A1;
  _r1_A2    = r1_A2;
  _r2_A2    = r2_A2;
  _mfit2_A2 = mfit2_A2;
  _r1_V     = r1_V;
  _r2_V     = r2_V;
  _mfit2_V  = mfit2_V;
  
  return;
}


void EvtBToVlnuBallFF::getvectorff(EvtId parent,EvtId /*daught*/,
				   double t, double /*mass*/, double *a1f,
				   double *a2f, double *vf, double *a0f ){
  
  // FF calculations taken from the LCSR calculation of 
  // P. Ball, R. Zwicky, Phys.~Rev.~{\bf D71} 014029 (2005), hep-ph/0412079.

  //Define mBstar
  EvtId Bplus = EvtPDL::getId("B+");
  EvtId Bminus = EvtPDL::getId("B-");
  double mBstar = EvtPDL::getMeanMass(EvtPDL::getId("B*0"));
  if (parent==Bplus || parent==Bminus) mBstar = EvtPDL::getMeanMass(EvtPDL::getId("B*+")); 
  
  double q2 = t;
  *a1f = _r2_A1/(1.-q2/_mfit2_A1);
  *a2f = _r1_A2/(1.-q2/_mfit2_A2) + _r2_A2/pow(1.-q2/_mfit2_A2,2.);
  *vf  = _r1_V /(1.-q2/mBstar/mBstar) + _r2_V/(1.-q2/_mfit2_V);
  *a0f = 0.0;
    
  return;
  
  // OLD STUFF from HQETFF

//   double mb=EvtPDL::getMeanMass(parent);
//   double w = ((mb*mb)+(mass*mass)-t)/(2.0*mb*mass);

// Form factors have a general form, with parameters passed in
// from the arguements.

//   double rstar = ( 2.0*sqrt(mb*mass))/(mb+mass);
//   double ha1 = 1-rho2*(w-1);

//   *a1f = (1.0 - (t/((mb+mass)*(mb+mass))))*ha1;
//   *a1f = (*a1f)/rstar;
//   *a2f = (r2/rstar)*ha1;
//   *vf = (r1/rstar)*ha1;

}


void EvtBToVlnuBallFF::getscalarff(EvtId, EvtId, double, double, double*, 
			       double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getvectorff in EvtBToVlnuBallFF.\n";  
  ::abort();

}



void EvtBToVlnuBallFF::gettensorff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :gettensorff in EvtBToVlnuBallFF.\n";  
  ::abort();

}



void EvtBToVlnuBallFF::getbaryonff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtBToVlnuBallFF.\n";  
  ::abort();

}


void EvtBToVlnuBallFF::getdiracff(EvtId, EvtId, double, double, double*, double*,
				  double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtBToVlnuBallFF.\n";
  ::abort();

}

void EvtBToVlnuBallFF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
				   double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtBToVlnuBallFF.\n";
  ::abort();

}
