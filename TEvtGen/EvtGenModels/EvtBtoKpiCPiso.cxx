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
// Module: EvtBtoKpiCPiso.cc
//
// Description: Routine to decay B -> K pi with isospin amplitudes
//
// Modification history:
//
//    RYD/NK     Febuary 7, 1998         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBtoKpiCPiso.hh"
#include <string>

EvtBtoKpiCPiso::~EvtBtoKpiCPiso() {}

std::string EvtBtoKpiCPiso::getName(){

  return "BTOKPI_CP_ISO";     

}


EvtDecayBase* EvtBtoKpiCPiso::clone(){

  return new EvtBtoKpiCPiso;

}


void EvtBtoKpiCPiso::init(){

  // check that there are 15 arguments
  checkNArg(15);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
}


void EvtBtoKpiCPiso::initProbMax() {

//this might need to be revised

  //added by Lange Jan4,2000
  static EvtId PI0=EvtPDL::getId("pi0");
  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PIM=EvtPDL::getId("pi+");
  static EvtId K0=EvtPDL::getId("K0");
  static EvtId KB=EvtPDL::getId("anti-K0");
  static EvtId KP=EvtPDL::getId("K+");
  static EvtId KM=EvtPDL::getId("K-");

if (((getDaug(0)==PI0) && (getDaug(1)==KP)) || ((getDaug(0)==KP) && (getDaug(1)==PI0)))  {
   setProbMax(2.0*(getArg(2)*getArg(2) + getArg(10)*getArg(10)));
 }

if (((getDaug(0)==PI0) && (getDaug(1)==KM)) || ((getDaug(0)==KM) && (getDaug(1)==PI0)))  {
   setProbMax( 2.0*(getArg(4)*getArg(4) + getArg(12)*getArg(12)));
 }

if (((getDaug(0)==PIP) && (getDaug(1)==K0)) || ((getDaug(0)==K0) && (getDaug(1)==PIP))) {
   setProbMax( 4.0*(getArg(6)*getArg(6) + getArg(10)*getArg(10)));
 }

if (((getDaug(0)==PIM) && (getDaug(1)==KB)) || ((getDaug(0)==KB) && (getDaug(1)==PIM))) {
   setProbMax( 4.0*(getArg(8)*getArg(8) + getArg(12)*getArg(12)));
}

if (((getDaug(0)==PI0) && (getDaug(1)==K0)) || ((getDaug(0)==K0) && (getDaug(1)==PI0))) { 
   setProbMax( 2.0*(getArg(2)*getArg(2) + getArg(10)*getArg(10))); 
}

if (((getDaug(0)==PI0) && (getDaug(1)==KB)) || ((getDaug(0)==KB) && (getDaug(1)==PI0))) {
   setProbMax( 2.0*(getArg(4)*getArg(4) + getArg(12)*getArg(12))); 
}

if (((getDaug(0)==PIM) && (getDaug(1)==KP)) || ((getDaug(0)==KP) && (getDaug(1)==PIM))) {
   setProbMax( 4.0*(getArg(6)*getArg(6) + getArg(10)*getArg(10)));
}

if (((getDaug(0)==PIP) && (getDaug(1)==KM)) || ((getDaug(0)==KM) && (getDaug(1)==PIP))) {
   setProbMax( 4.0*(getArg(8)*getArg(8) + getArg(12)*getArg(12)));
}
}

void EvtBtoKpiCPiso::decay( EvtParticle *p ){


  p->initializePhaseSpace(getNDaug(),getDaugs());
  //added by Lange Jan4,2000
  static EvtId PI0=EvtPDL::getId("pi0");
  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PIM=EvtPDL::getId("pi+");
  static EvtId K0=EvtPDL::getId("K0");
  static EvtId KB=EvtPDL::getId("anti-K0");
  static EvtId KP=EvtPDL::getId("K+");
  static EvtId KM=EvtPDL::getId("K-");


   EvtComplex A;
   EvtComplex U, Ubar, V, Vbar, W, Wbar;

   U = EvtComplex(getArg(2)*cos(getArg(3)),getArg(2)*sin(getArg(3)));
   Ubar = EvtComplex(getArg(4)*cos(getArg(5)),getArg(4)*sin(getArg(5)));
   V = EvtComplex(getArg(6)*cos(getArg(7)),getArg(6)*sin(getArg(7)));
   Vbar = EvtComplex(getArg(8)*cos(getArg(9)),getArg(8)*sin(getArg(9)));
   W = EvtComplex(getArg(10)*cos(getArg(11)),getArg(10)*sin(getArg(11)));
   Wbar = EvtComplex(getArg(12)*cos(getArg(13)),getArg(12)*sin(getArg(13)));  

//depending on what combination of K pi we have, there will be different 
//A and Abar (only A given in comments!)

if (((getDaug(0)==PI0) && (getDaug(1)==KP)) || ((getDaug(0)==KP) && (getDaug(1)==PI0)))  {

//pi0 K+, so U - W

   A = U - W;
  }
 
if (((getDaug(0)==PI0) && (getDaug(1)==KM)) || ((getDaug(0)==KM) && (getDaug(1)==PI0)))  {

//pi0 K-, so Ubar - Wbar
  
   A = Ubar - Wbar;
 }

 if (((getDaug(0)==PIP) && (getDaug(1)==K0)) || ((getDaug(0)==K0) && (getDaug(1)==PIP))) {

//pi+ K0, so V + W
  
    A = sqrt(2.0)*(V + W);
 }

 if (((getDaug(0)==PIM) && (getDaug(1)==KB)) || ((getDaug(0)==KB) && (getDaug(1)==PIM))) {

//pi- K0bar, so Vbar + Wbar
   A = sqrt(2.0)*(Vbar + Wbar);
 }

 if (((getDaug(0)==PI0) && (getDaug(1)==K0)) || ((getDaug(0)==K0) && (getDaug(1)==PI0))) { 

//pi0 K0, so U + W

   A= U + W;
  }

 if (((getDaug(0)==PI0) && (getDaug(1)==KB)) || ((getDaug(0)==KB) && (getDaug(1)==PI0))) {

  A = Ubar + Wbar;
  }

 if (((getDaug(0)==PIM) && (getDaug(1)==KP)) || ((getDaug(0)==KP) && (getDaug(1)==PIM))) {

//pi- K+, so V - W
    
    A = sqrt(2.0)*(V - W);
  }

 if (((getDaug(0)==PIP) && (getDaug(1)==KM)) || ((getDaug(0)==KM) && (getDaug(1)==PIP))) {

     A = sqrt(2.0)*(Vbar - Wbar);
  }

   vertex(A);

  return ;
}

