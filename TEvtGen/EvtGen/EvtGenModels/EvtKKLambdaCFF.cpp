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
// Module: EvtSLPoleFF.cc
//
// Description: Routine to implement semileptonic form factors
//              according to the model SLPoles
//
// Modification history:
//
//    DJL       April 17,1998       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtKKLambdaCFF.hh"
#include <string>
#include "EvtGenBase/EvtPDL.hh"
#include <math.h>
#include <stdlib.h>

EvtKKLambdaCFF::EvtKKLambdaCFF(int numarg, double *arglist) {
   _nargs = numarg;
   for (int i=0; i<numarg; i++) {
     _args[i] = arglist[i]; }

   return;
}

void EvtKKLambdaCFF::getbaryonff(EvtId /*parent*/,EvtId /*daught*/,
		 double t, double /*mass*/, double *f1v,
		 double *f1a, double *f2v, double *f2a ) {
  
  *f1v=(_args[0])/(1.0-(t/(_args[1]*_args[1])));

  *f2v=0.;
  *f2a=0.;
  *f1a=-1.0*(*f1v);
  
}

void EvtKKLambdaCFF::getscalarff(EvtId, EvtId,
				 double, double, double*,
				 double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getscalarff in EvtKKLambdaCFF.\n";  
  ::abort();

}

void EvtKKLambdaCFF::getvectorff(EvtId, EvtId,
				 double, double, double*,
				 double*, double*, double* ){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getvectorff in EvtKKLambdaCFF.\n";  
  ::abort();

}

  
void EvtKKLambdaCFF::gettensorff(EvtId, EvtId, double, double, double*, 
			       double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :gettensorff in EvtKKLambdaCFF.\n";  
  ::abort();

}

void EvtKKLambdaCFF::getdiracff(EvtId, EvtId, double, double, double*, double*,
				double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtKKLambdaCFF.\n";
  ::abort();

}

void EvtKKLambdaCFF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
				 double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtKKLambdaCFF.\n";
  ::abort();

}
