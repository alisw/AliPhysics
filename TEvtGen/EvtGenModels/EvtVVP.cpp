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
// Module: EvtVVP.cc
//
// Description: The decay Vector -> Vector gamma
//              E.g., CHI1->PSI GAMMA
//
// Modification history:
//
//    RYD       September 5, 1997       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <string>
#include "EvtGenBase/EvtVector3C.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenModels/EvtVVP.hh"
#include "EvtGenBase/EvtReport.hh"

EvtVVP::~EvtVVP() {}

std::string EvtVVP::getName(){

  return "VVP";
     
}


EvtDecayBase* EvtVVP::clone(){

  return new EvtVVP;

}

void EvtVVP::init(){

  // check that there are 8 arguments

  checkNArg(8);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::VECTOR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::PHOTON);
}

void EvtVVP::initProbMax(){

  setProbMax(4.0);

}

void EvtVVP::decay(EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtParticle *v,*ph;

  v = p->getDaug(0);
  ph = p->getDaug(1);

  EvtVector3C epsp[3];
  EvtVector3C epsv[3];
  EvtVector3C epsph[2];

  epsp[0]=p->eps(0).vec();
  epsp[1]=p->eps(1).vec();
  epsp[2]=p->eps(2).vec();

  epsv[0]=v->eps(0).vec().conj();
  epsv[1]=v->eps(1).vec().conj();
  epsv[2]=v->eps(2).vec().conj();

  epsph[0]=ph->epsParentPhoton(0).vec().conj();
  epsph[1]=ph->epsParentPhoton(1).vec().conj();

  int i,j,k;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      for(k=0;k<2;k++){
	vertex(i,j,k,epsp[i].cross(epsv[j])*epsph[k]);
      }
    }
  }

  return;

}

