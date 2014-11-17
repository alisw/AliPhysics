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
// Module: EvtKstarnunu.cc
//
// Description: B ==> K* nu nubar
//
// Modification history:
//
//    RYD     September 29, 1997        Module created
//
//------------------------------------------------------------------------
//
// Package:     EVT
// Module:      EvtKstarnunu
// 
//              
// Implementation:
//     
//
// Author:      RYD 
// Created:     Sept. 29 1997
//
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <string>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenModels/EvtKstarnunu.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"

EvtKstarnunu::~EvtKstarnunu() {}

std::string EvtKstarnunu::getName(){

  return "KSTARNUNU";     

}


EvtDecayBase* EvtKstarnunu::clone(){

  return new EvtKstarnunu;

}

void EvtKstarnunu::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(3);

  //We expect the parent to be a scalar 
  //and the daughters to be K neutrino netrino

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::NEUTRINO);
  checkSpinDaughter(2,EvtSpinType::NEUTRINO);

}


void EvtKstarnunu::decay(EvtParticle *p){
  
  static EvtId NUE=EvtPDL::getId("nu_e");
  static EvtId NUM=EvtPDL::getId("nu_mu");
  static EvtId NUT=EvtPDL::getId("nu_tau");
  static EvtId NUEB=EvtPDL::getId("anti-nu_e");
  static EvtId NUMB=EvtPDL::getId("anti-nu_mu");
  static EvtId NUTB=EvtPDL::getId("anti-nu_tau");

  p->initializePhaseSpace(getNDaug(),getDaugs());
  
  double m_b = p->mass();

  EvtParticle *meson, *neutrino1, *neutrino2;
  meson = p->getDaug(0);
  neutrino1 = p->getDaug(1);
  neutrino2 = p->getDaug(2);
  EvtVector4R momnu1 = neutrino1->getP4();
  EvtVector4R momnu2 = neutrino2->getP4();
  EvtVector4R momkstar = meson->getP4();
  
  double v0_0, a1_0, a2_0;
  double m2v0, a1_b, a2_b;
  v0_0 = 0.47;
  a1_0 = 0.37;
  a2_0 = 0.40;
  m2v0 = 5.*5.;
  a1_b = -0.023;
  a2_b = 0.034;

  EvtVector4R q = momnu1+momnu2;
  double q2 = q.mass2();   

  double v0, a1, a2;
  v0 = v0_0/(1-q2/m2v0);
  a1 = a1_0*(1+a1_b*q2);
  a2 = a2_0*(1+a2_b*q2);

  EvtVector4R p4b; p4b.set(m_b,0.,0.,0.);  // Do calcs in mother rest frame
  
  double m_k = meson->mass();

  EvtTensor4C tds=(-2*v0/(m_b+m_k))*dual(EvtGenFunctions::directProd(p4b,momkstar))
    - EvtComplex(0.0,1.0)*
    ( (m_b+m_k)*a1*EvtTensor4C::g()
      - (a2/(m_b+m_k))*EvtGenFunctions::directProd(p4b-momkstar,p4b+momkstar)); 
  
  EvtVector4C l;

  if (getDaug(1)==NUE||getDaug(1)==NUM||getDaug(1)==NUT) {
        l=EvtLeptonVACurrent(neutrino1->spParentNeutrino(),
			     neutrino2->spParentNeutrino());
  }
  if (getDaug(1)==NUEB||getDaug(1)==NUMB||getDaug(1)==NUTB) {
        l=EvtLeptonVACurrent(neutrino2->spParentNeutrino(),
			     neutrino1->spParentNeutrino());
  }

  EvtVector4C et0,et1,et2; 
  et0 = tds.cont1( meson->epsParent(0).conj() );
  et1 = tds.cont1( meson->epsParent(1).conj() );
  et2 = tds.cont1( meson->epsParent(2).conj() );

  vertex(0,l*et0);
  vertex(1,l*et1);
  vertex(2,l*et2);

  return;
}

