//--------------------------------------------------------------------------
//
// Module: EvtHQET2FF.cc
//
// Description: form factors for B->D*lnu & B->Dlnu according to HQET 
//              with dispersive FF
//
// Modification history:
//
//    Marco Bomben     March 10, 2003        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtHQET2FF.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtPDL.hh"
#include <math.h>
#include <stdlib.h>


EvtHQET2FF::EvtHQET2FF(double hqetrho2, double hqetha1_1 , double hqetr1_1, double hqetr2_1) {

  rho2 = hqetrho2;
  r1_1 = hqetr1_1;
  r2_1 = hqetr2_1;
  ha1_1 = hqetha1_1;

  return;
}

EvtHQET2FF::EvtHQET2FF(double hqetrho2, double hqetv1_1) {

  rho2 = hqetrho2;
  v1_1 = hqetv1_1;

  return;
}


void EvtHQET2FF::getscalarff(EvtId parent,EvtId,
			    double t, double mass, double *f0p, double *f0m) {


  double mb=EvtPDL::getMeanMass(parent);
  double w = ((mb*mb)+(mass*mass)-t)/(2.0*mb*mass);

// Form factors have a general form, with parameters passed in
// from the arguements.

  // Use disparsion relation parametrization from 
  // I.Caprini, L.Lelluch, M.Neubert, Nucl. Phys. B 530,153(1998)
  const double z = (sqrt(w+1)-sqrt(2.))/(sqrt(w+1)+sqrt(2.));
  double v1 = v1_1*(1.- 8.*rho2*z + (51.*rho2-10.)*z*z - (252.*rho2-84.)*z*z*z)
;

  *f0p=v1;
  *f0m = 0.0;

  return;
 }

void EvtHQET2FF::getvectorff(EvtId parent,EvtId,
			    double t, double mass, double *a1f,
			    double *a2f, double *vf, double *a0f ){


  double mb=EvtPDL::getMeanMass(parent);
  double w = ((mb*mb)+(mass*mass)-t)/(2.0*mb*mass);

// Form factors have a general form, with parameters passed in
// from the arguements.

  double rstar = ( 2.0*sqrt(mb*mass))/(mb+mass);

  // Use disparsion relation parametrization from 
  // I.Caprini, L.Lelluch, M.Neubert, Nucl. Phys. B 530,153(1998)
  const double z = (sqrt(w+1)-sqrt(2.))/(sqrt(w+1)+sqrt(2.));
  double ha1 =ha1_1*(1.- 8.*rho2*z + (53.*rho2-15.)*z*z - (231.*rho2-91.)*z*z*z);
  double r1 = r1_1-0.12*(w-1)+0.05*(w-1)*(w-1);
  double r2 = r2_1+0.11*(w-1)-0.06*(w-1)*(w-1);
;

  *a1f = (1.0 - (t/((mb+mass)*(mb+mass))))*ha1;
  *a1f = (*a1f)/rstar;
  *a2f = (r2/rstar)*ha1;
  *vf = (r1/rstar)*ha1;
  *a0f = 0.0;

  return;
 }

void EvtHQET2FF::gettensorff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :gettensorff in EvtHQET2FF.\n";  
  ::abort();

}



void EvtHQET2FF::getbaryonff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtHQET2FF.\n";  
  ::abort();

}

void EvtHQET2FF::getdiracff(EvtId, EvtId, double, double, double*, double*,
			    double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtHQET2FF.\n";
  ::abort();

}

void EvtHQET2FF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
			     double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtHQET2FF.\n";
  ::abort();

}
