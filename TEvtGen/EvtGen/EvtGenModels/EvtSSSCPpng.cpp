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
// Module: EvtSSSCPpng.cc
//
// Description: Routine to decay B -> 2 scalars taking into account penguin 
//              contributions (assuming single quark dominance for penguins)
//
// Modification history:
//
//    RYD/NK    December 3, 1997          Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtSSSCPpng.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"

EvtSSSCPpng::~EvtSSSCPpng() {}

std::string EvtSSSCPpng::getName(){

  return "SSS_CP_PNG";     

}


EvtDecayBase* EvtSSSCPpng::clone(){

  return new EvtSSSCPpng;

}

void EvtSSSCPpng::init(){

  // check that there are 7 arguments
  checkNArg(7);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);

}

void EvtSSSCPpng::initProbMax(){
 
   setProbMax(getArg(5)*getArg(5)*(1 + getArg(6)*getArg(6)));

}

void EvtSSSCPpng::decay( EvtParticle *p ){

  //added by Lange Jan4,2000
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  double t;
  EvtId other_b;
 
  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtComplex amp;

  EvtComplex A,Abar;
  //EvtComplex ACC, AbarCC;

// assume single (top) quark dominance for the penguin.

// old: a0=alpha, a1=dm, a2=1, a3=1, a4=0, a5=1, a6=0
// new: a0=beta, a1=gamma, a2=delta, a3=dm, a4=1, a5=1=A_{T}, a6=A_{P}/A_{T}

// e.g., for B -> pi pi  
// A_{T} = |V_{ub} V_{ud}| T
// A_{P} = |V_{tb} V_{td}| P
// P and T are purely hadronic matrix elements       
// P/T = 0.055, A_{P}/A_{T} = 0.2 (see Marrocchesi and Paver, hep-ph/9702353)

// A = A_{T}( exp(i(beta+gamma)) + (A_{P}/A_{T}) exp(i(delta))
// A_bar = same, except for the sign of the weak phases 
// here, delta = delta_{p}-delta_{t} (rel. strong phase)

   A=getArg(5)*(EvtComplex(cos(-getArg(0)-getArg(1)),sin(-getArg(0)-getArg(1)))+getArg(6)*EvtComplex(cos(getArg(2)),sin(getArg(2))));
 
   Abar=getArg(5)*(EvtComplex(cos(getArg(0)+getArg(1)),sin(getArg(0)+getArg(1)))+getArg(6)*EvtComplex(cos(getArg(2)),sin(getArg(2))));
 
// get fraction of B0 tags with these amplitudes

  //double xd = 0.65;
  double ratio = 1/(1 + 0.65*0.65);
  
  EvtComplex rf, rbarf;

  rf = EvtComplex(cos(2.0*getArg(0)),sin(2.0*getArg(0)))*Abar/A;
  rbarf = EvtComplex(1.0)/rf;

  double A2 = real(A)*real(A) + imag(A)*imag(A);
  double Abar2 = real(Abar)*real(Abar) + imag(Abar)*imag(Abar);
  
  double rf2 = real(rf)*real(rf) + imag(rf)*imag(rf);    
  double rbarf2 = real(rbarf)*real(rbarf) + imag(rbarf)*imag(rbarf);    

  //fraction of B0 _tags_
  double fract =(Abar2*(1+ rbarf2 + (1 - rbarf2)*ratio))/(Abar2*(1+ rbarf2 + (1 - rbarf2)*ratio) + A2*(1+ rf2 + (1 - rf2)*ratio)); 
  
  EvtCPUtil::getInstance()->OtherB(p,t,other_b,fract);

//this method works just as well -- NK
//randomly generate the tag (B0 or B0B) 

//  double tag = EvtRandom::Flat(0.0,1.0);
//  if (tag < 0.5) {
//
//   EvtCPUtil::OtherB(p,t,other_b,1.0);
//   other_b = B0;
//  }
//  else {
//   
//   EvtCPUtil::OtherB(p,t,other_b,0.0);
//   other_b = B0B;
//  }

//mixing angle = -beta

   if (other_b==B0B){
     amp=A*cos(getArg(3)*t/(2*EvtConst::c))+
       EvtComplex(cos(2.0*getArg(0)),sin(2.0*getArg(0)))*
       getArg(4)*EvtComplex(0.0,1.0)*Abar*sin(getArg(3)*t/(2*EvtConst::c));
   }
   if (other_b==B0){
     amp=A*EvtComplex(cos(-2.0*getArg(0)),sin(-2.0*getArg(0)))*
       EvtComplex(0.0,1.0)*sin(getArg(3)*t/(2*EvtConst::c))+       
       getArg(4)*Abar*cos(getArg(3)*t/(2*EvtConst::c));
   }

   vertex(amp);

  return ;
}


