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
// Module: EvtSLBKPoleFF.cc
//
// Description: Routine to implement semileptonic form factors
//              according to the model SLBKPoles
//
// Modification history:
//
//    liheng       October 20,2005       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtSLBKPoleFF.hh"//modified
#include <string>
#include "EvtGenBase/EvtPDL.hh"
#include <math.h>
#include <stdlib.h>

EvtSLBKPoleFF::EvtSLBKPoleFF(int numarg, double *arglist) {//modified
   numSLBKPoleargs = numarg;//modified
   for (int i=0; i<numarg; i++) {
      SLBKPoleargs[i] = arglist[i]; }//modified

   return;
}


void EvtSLBKPoleFF::getscalarff(EvtId parent,EvtId daught,
                       double t, double /*mass*/, double *fpf,
			    double *f0f ) {

// Form factors have a general form, with parameters passed in
// from the arguments.

   if ( numSLBKPoleargs != 4 ) {//modified
     report(Severity::Error,"EvtGen") << "Problem in EvtSLBKPoleFF::getscalarff\n";
     report(Severity::Error,"EvtGen") << "wrong number of arguments!\n";
     report(Severity::Error,"EvtGen") << "number args:"<<numSLBKPoleargs<<" (expected 4)\n";
     report(Severity::Error,"EvtGen") << "Parent:"<<EvtPDL::name(parent)<<"\n";
     report(Severity::Error,"EvtGen") << "Daughter:"<<EvtPDL::name(daught)<<"\n";

   }


  double f0,af,powf;

  //double a_0, a_1, a_2, a_3, a_4, a_5, a_6, a_7;

  f0 = SLBKPoleargs[0];//f0
  af = SLBKPoleargs[1];//alpha
  //bf = SLBKPoleargs[2];
  double mass_star2 = SLBKPoleargs[3]*SLBKPoleargs[3];
  powf = 1.0;
  *fpf = f0/(pow( 1.0 - (1.0+af)*(t/mass_star2) + (af*((t/mass_star2)*(t/mass_star2))),powf));//modified

  f0 = SLBKPoleargs[0];//f0
  af = SLBKPoleargs[2];//beta
  //bf = SLBKPoleargs[6];
  powf = 1.0;

  *f0f = f0/(pow( 1.0 - (t/mass_star2/af),powf));//modified

  return;
}

void EvtSLBKPoleFF::getvectorff(EvtId parent,EvtId /*daught*/,
                       double t, double /*mass*/, double *a1f,
			     double *a2f, double *vf, double *a0f ){

   if ( numSLBKPoleargs !=8 ) {//modified
     report(Severity::Error,"EvtGen") << "Problem in EvtSLBKPoleFF::getvectorff\n";//modified
     report(Severity::Error,"EvtGen") << "wrong number of arguements!!!\n";
     report(Severity::Error,"EvtGen") << numSLBKPoleargs<<"\n";//modified
//     printf("\n*********************%d*********************",numSLBKPoleargs);
  }


   report(Severity::Info,"EvtGen")<<"Check the implementation of EvtSLBKPoleFF::getvectorff()!\n";


  double mb=EvtPDL::getMeanMass(parent);
  double mb2 = mb*mb;

//modified-begin
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");  
  static EvtId BP=EvtPDL::getId("B+");
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId BS0=EvtPDL::getId("B_s0");

  static EvtId B0S=EvtPDL::getId("B*0");
  static EvtId BPMS=EvtPDL::getId("B*+");  
  static EvtId BS0S=EvtPDL::getId("B_s*0");
  
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId DM=EvtPDL::getId("D-");
  static EvtId DSP=EvtPDL::getId("D_s+");
  static EvtId DSM=EvtPDL::getId("D_s-");

  static EvtId D0S=EvtPDL::getId("D*0");
  static EvtId DPMS=EvtPDL::getId("D*+");
  static EvtId DSPMS=EvtPDL::getId("D_s*+");

  double mass_star=0.0;
  double mass_star2=0.0;
  if(parent==B0||parent==B0B){
     mass_star=EvtPDL::getMeanMass(B0S);
     mass_star2=mass_star*mass_star;
  }
  if(parent==BP||parent==BM){
     mass_star=EvtPDL::getMeanMass(BPMS);
     mass_star2=mass_star*mass_star;
  }
  if(parent==BS0){
     mass_star=EvtPDL::getMeanMass(BS0S);
     mass_star2=mass_star*mass_star;
  }

  if(parent==D0||parent==D0B){
     mass_star=EvtPDL::getMeanMass(D0S);
     mass_star2=mass_star*mass_star;
  }
  if(parent==DP||parent==DM){
     mass_star=EvtPDL::getMeanMass(DPMS);
     mass_star2=mass_star*mass_star;
  }
  if(parent==DSP||parent==DSM){
     mass_star=EvtPDL::getMeanMass(DSPMS);
     mass_star2=mass_star*mass_star;
  }
//modified-end

  double f0,af,bf,powf;

  f0 = SLBKPoleargs[2];//A1
  af = SLBKPoleargs[6];//b'
  bf = 0;//0
  powf = 1.0;//1.0
  *a1f = f0/(pow( 1.0 - af*t/mass_star2,powf));//modified

  f0 = SLBKPoleargs[3];//A2
  af = SLBKPoleargs[6];//b'
  bf = SLBKPoleargs[7];//b''==0
  powf = 1.0;//1.0

  *a2f = f0/(pow(1.0 - (af+bf)*(t/mass_star2) + (af*bf)*((t/mass_star2)*(t/mass_star2)),powf));//modified

  f0 = SLBKPoleargs[0];//V0
  af = SLBKPoleargs[4];//a
  bf = 0;//0
  powf = 1.0;//1.0

  *vf = f0/(pow( 1.0 - (1.0+af)*(t/mass_star2) + af*(t/mass_star2)*(t/mass_star2),powf));//modified

  f0 = SLBKPoleargs[1];//A0
  af = SLBKPoleargs[5];//a'
  bf = 0;//0
  powf = 1.0;//1.0

  *a0f = f0/(pow( 1.0 - (1.0+af)*(t/mb2) + af*((t/mb2)*(t/mb2)),powf));//modified
  return;
 }



void EvtSLBKPoleFF::gettensorff(EvtId parent,EvtId /*daught*/,
                       double t, double /*mass*/, double *hf,
			     double *kf, double *bpf, double *bmf ){

  if ( numSLBKPoleargs !=16 ) {
     report(Severity::Error,"EvtGen") << "Problem in EvtSLBKPoleFF::gettensorff\n";
     report(Severity::Error,"EvtGen") << "wrong number of arguements!!!\n";
  }

  report(Severity::Info,"EvtGen")<<"Check the implementation of EvtSLBKPoleFF::gettensorff()!\n";

  double mb=EvtPDL::getMeanMass(parent);
  double mb2 = mb*mb;

  double f0,af,bf,powf;

  f0 = SLBKPoleargs[0];
  af = SLBKPoleargs[1];
  bf = SLBKPoleargs[2];
  powf = SLBKPoleargs[3];
  *hf = f0/(pow( 1.0 + (af*t/mb2) + (bf*((t/mb2)*(t/mb2))),powf));

  f0 = SLBKPoleargs[4];
  af = SLBKPoleargs[5];
  bf = SLBKPoleargs[6];
  powf = SLBKPoleargs[7];

  *kf = f0/(pow( 1.0 + (af*t/mb2) + (bf*((t/mb2)*(t/mb2))),powf));

  f0 = SLBKPoleargs[8];
  af = SLBKPoleargs[9];
  bf = SLBKPoleargs[10];
  powf = SLBKPoleargs[11];

  *bpf = f0/(pow( 1.0 + (af*t/mb2) + (bf*((t/mb2)*(t/mb2))),powf));

  f0 = SLBKPoleargs[12];
  af = SLBKPoleargs[13];
  bf = SLBKPoleargs[14];
  powf = SLBKPoleargs[15];

  *bmf = f0/(pow( 1.0 + (af*t/mb2) + (bf*((t/mb2)*(t/mb2))),powf));
  return;
 }

void EvtSLBKPoleFF::getbaryonff(EvtId, EvtId, double, double, double*, 
				double*, double*, double*){
  
  report(Severity::Error,"EvtGen") << "Not implemented :getbaryonff in EvtSLBKPoleFF.\n";  
  ::abort();

}

void EvtSLBKPoleFF::getdiracff(EvtId, EvtId, double, double, double*, double*,
			       double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getdiracff in EvtSLBKPoleFF.\n";
  ::abort();

}

void EvtSLBKPoleFF::getraritaff(EvtId, EvtId, double, double, double*, double*, 
				double*, double*, double*, double*, double*, double*) {
  
  report(Severity::Error,"EvtGen") << "Not implemented :getraritaff in EvtSLBKPoleFF.\n";
  ::abort();

}


