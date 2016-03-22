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
// Module: EvtMelikhovFF.cc
//
// Description: form factors for B->D*lnu according to HQET
//
// Modification history:
//
//    DJL     April 17, 1998        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtMelikhovFF.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include <math.h>
#include "EvtGenBase/EvtPDL.hh"
#include <stdlib.h>

EvtMelikhovFF::EvtMelikhovFF(double arg1) {

  whichfit = int(arg1+0.1);
  return;
}

void EvtMelikhovFF::getvectorff(EvtId parent,EvtId,
                       double t, double mass, double *a1f,
			     double *a2f, double *vf, double *a0f ){

  double ma1(0.0),ra1(0.0),na1(0.0);
  double ma2(0.0),ra2(0.0),na2(0.0);
  double mv(0.0),rv(0.0),nv(0.0);

  if ( whichfit==1) {
     ma1 = 7.07;
     ra1 = 0.27;
     na1 = 2.65;
     ma2 = 6.13;
     ra2 = 0.25;
     na2 = 2.17;
     mv = 6.28;
     rv = 0.30;
     nv = 2.36;
  }
  if ( whichfit==2) {
     ma1 = 6.78;
     ra1 = 0.20;
     na1 = 2.65;
     ma2 = 6.00;
     ra2 = 0.19;
     na2 = 2.34;
     mv = 6.22;
     rv = 0.20;
     nv = 2.46;
  }
  if ( whichfit==3) {
     ma1 = 6.50;
     ra1 = 0.21;
     na1 = 2.70;
     ma2 = 5.90;
     ra2 = 0.20;
     na2 = 2.45;
     mv = 5.90;
     rv = 0.21;
     nv = 2.35;
  }
  if ( whichfit==4) {
     ma1 = 5.68;
     ra1 = 0.29;
     na1 = 1.67;
     ma2 = 5.36;
     ra2 = 0.28;
     na2 = 1.67;
     mv = 5.46;
     rv = 0.29;
     nv = 1.73;
  }

  double mb=EvtPDL::getMeanMass(parent);
  //double w = ((mb*mb)+(mass*mass)-t)/(2.0*mb*mass);

  double melr = mass/mb;
  double mely = t/(mb*mb);
 

  *a1f = ((1.0+melr*melr-mely)/(1+melr))*ra1/pow(1.0-(t/(ma1*ma1)),na1);
  *a2f = (1+melr)*((1.0-melr*melr-mely)/((1+melr)*(1+melr)-mely))
    *ra2/pow(1.0-(t/(ma2*ma2)),na2);
  *vf = (1+melr)*rv/pow(1.0-(t/(mv*mv)),nv);
  *a0f = 0.0;

  return;
 }


void EvtMelikhovFF::getscalarff(EvtId, EvtId, double, double, double*, 
			       double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getvectorff in EvtMelikhovFF.\n";  
  ::abort();

}



void EvtMelikhovFF::gettensorff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :gettensorff in EvtMelikhovFF.\n";  
  ::abort();

}



void EvtMelikhovFF::getbaryonff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtMelikhovFF.\n";  
  ::abort();

}

void EvtMelikhovFF::getdiracff(EvtId, EvtId, double, double, double*, double*,
			       double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtMelikhovFF.\n";
  ::abort();

}

void EvtMelikhovFF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
				double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtMelikhovFF.\n";
  ::abort();

}

