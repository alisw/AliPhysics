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
// Module: EvtVVSPwave.cc
//
// Description: Routine to decay vector-> vector pi pi where the 
//              decay is S-wave dominated.
//
// Modification history:
//
//    RYD       December 11, 1999       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtVVpipi.hh"
#include <string>
using std::endl;

EvtVVpipi::~EvtVVpipi() {}

std::string EvtVVpipi::getName(){

  return "VVPIPI";     

}


EvtDecayBase* EvtVVpipi::clone(){

  return new EvtVVpipi;

}

void EvtVVpipi::init(){

  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PIM=EvtPDL::getId("pi-");
  static EvtId PI0=EvtPDL::getId("pi0");

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(3);

  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);



  if ((!(getDaug(1)==PIP&&getDaug(2)==PIM))&&
      (!(getDaug(1)==PI0&&getDaug(2)==PI0))) {
    report(Severity::Error,"EvtGen") << "EvtVVpipi generator expected "
                           << " pi+ and pi- (or pi0 and pi0) "
			   << "as 2nd and 3rd daughter. "<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }

}

void EvtVVpipi::initProbMax() {

  //Hard coded... should not be hard to calculate...
  setProbMax(0.08);

}      

void EvtVVpipi::decay( EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtParticle *v,*s1,*s2;
  
  v=p->getDaug(0);
  s1=p->getDaug(1);
  s2=p->getDaug(2);

//  Put phase space results into the daughters.
  
  EvtVector4C ep0,ep1,ep2;  
  
  ep0=p->eps(0);
  ep1=p->eps(1);
  ep2=p->eps(2);

  double fac=(s1->getP4()+s2->getP4()).mass2()-4*s1->mass()*s2->mass();

  vertex(0,0,fac*(ep0*v->epsParent(0).conj()));
  vertex(0,1,fac*(ep0*v->epsParent(1).conj()));
  vertex(0,2,fac*(ep0*v->epsParent(2).conj()));
  
  vertex(1,0,fac*(ep1*v->epsParent(0).conj()));
  vertex(1,1,fac*(ep1*v->epsParent(1).conj()));
  vertex(1,2,fac*(ep1*v->epsParent(2).conj()));
  
  vertex(2,0,fac*(ep2*v->epsParent(0).conj()));
  vertex(2,1,fac*(ep2*v->epsParent(1).conj()));
  vertex(2,2,fac*(ep2*v->epsParent(2).conj()));

  return ;

}




