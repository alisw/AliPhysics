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
// Module: EvtTVSPwave.cc
//
// Description: Routine to decay tensor-> vector scalar
//              by specifying the partial waves
//
// Modification history:
//
//    DJL/RYD   August 11, 1997           Module created
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
#include "EvtGenModels/EvtTVSPwave.hh"
#include <string>

EvtTVSPwave::~EvtTVSPwave() {}

std::string EvtTVSPwave::getName(){

  return "TVS_PWAVE";     

}


EvtDecayBase* EvtTVSPwave::clone(){

  return new EvtTVSPwave;

}

void EvtTVSPwave::init(){

  // check that there are 6 arguments
  checkNArg(6);
  checkNDaug(2);
    
  checkSpinParent(EvtSpinType::TENSOR);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
}

void EvtTVSPwave::initProbMax(){

  setProbMax(0.5);


}

void EvtTVSPwave::decay( EvtParticle *p ){

  EvtComplex ap(getArg(0)*cos(getArg(1)),getArg(0)*sin(getArg(1)));
  EvtComplex ad(getArg(2)*cos(getArg(3)),getArg(2)*sin(getArg(3)));
  EvtComplex af(getArg(4)*cos(getArg(5)),getArg(4)*sin(getArg(5)));

  if (ap!=EvtComplex(0.0,0.0)||af!=EvtComplex(0.0,0.0)) {
    report(Severity::Error,"EvtGen") << "dfslkh8945wqh:In EvtTensorToVectorScalar.c\n";
    report(Severity::Error,"EvtGen") << "P or F wave not yet implemented!! (ryd) \n";
  }

  p->initializePhaseSpace(getNDaug(),getDaugs());
    
  EvtParticle *v;
  v = p->getDaug(0);
  EvtVector4R momv = v->getP4();
  double massv = v->mass();

  EvtComplex temp;
  temp = ad;
  double m_parent = p->mass();

  EvtVector4R p_parent;
  
  p_parent.set(m_parent,0.0,0.0,0.0);
  
  EvtVector4C  pep0,pep1,pep2,pep3,pep4;  
  EvtTensor4C pdual;
  
  EvtVector4C epsdual0,epsdual1,epsdual2;

  double norm=massv/(m_parent*momv.get(0)*momv.d3mag()*momv.d3mag());
  pdual=dual(EvtGenFunctions::directProd(norm*p_parent,momv));
  
  epsdual0=pdual.cont1(v->epsParent(0).conj());
  epsdual1=pdual.cont1(v->epsParent(1).conj());
  epsdual2=pdual.cont1(v->epsParent(2).conj());
  
  pep0=p->epsTensor(0).cont1(momv);
  pep1=p->epsTensor(1).cont1(momv);
  pep2=p->epsTensor(2).cont1(momv);
  pep3=p->epsTensor(3).cont1(momv);
  pep4=p->epsTensor(4).cont1(momv);

  vertex(0,0,pep0*epsdual0);
  vertex(1,0,pep1*epsdual0);
  vertex(2,0,pep2*epsdual0);
  vertex(3,0,pep3*epsdual0);
  vertex(4,0,pep4*epsdual0);
  
  vertex(0,1,pep0*epsdual1);
  vertex(1,1,pep1*epsdual1);
  vertex(2,1,pep2*epsdual1);
  vertex(3,1,pep3*epsdual1);
  vertex(4,1,pep4*epsdual1);
  
  vertex(0,2,pep0*epsdual2);
  vertex(1,2,pep1*epsdual2);
  vertex(2,2,pep2*epsdual2);
  vertex(3,2,pep3*epsdual2);
  vertex(4,2,pep4*epsdual2);
  
  return ;
}

std::string EvtTVSPwave::getParamName(int i) {
  switch(i) {
  case 0:
    return "PWave";
  case 1:
    return "PWavePhase";
  case 2:
    return "DWave";
  case 3:
    return "DWavePhase";
  case 4:
    return "FWave";
  case 5:
    return "FWavePhase";
  default:
    return "";
  }
}

std::string EvtTVSPwave::getParamDefault(int i) {
  switch(i) {
  case 0:
    return "0.0";
  case 1:
    return "0.0";
  case 2:
    return "1.0";
  case 3:
    return "0.0";
  case 4:
    return "0.0";
  case 5:
    return "0.0";
  default:
    return "";
  }
}
