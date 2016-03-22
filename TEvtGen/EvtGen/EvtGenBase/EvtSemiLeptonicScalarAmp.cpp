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
// Module: EvtSemiLeptonicScalarAmp.cc
//
// Description: Routine to implement semileptonic decays to pseudo-scalar
//              mesons. 
//
// Modification history:
//
//    DJL       April 17,1998       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtSemiLeptonicScalarAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtSemiLeptonicFF.hh"

void EvtSemiLeptonicScalarAmp::CalcAmp( EvtParticle *parent,
					EvtAmp& amp,
					EvtSemiLeptonicFF *FormFactors ) {

  static EvtId EM=EvtPDL::getId("e-");
  static EvtId MUM=EvtPDL::getId("mu-");
  static EvtId TAUM=EvtPDL::getId("tau-");
  static EvtId EP=EvtPDL::getId("e+");
  static EvtId MUP=EvtPDL::getId("mu+");
  static EvtId TAUP=EvtPDL::getId("tau+");
  
  //Add the lepton and neutrino 4 momenta to find q2

  EvtVector4R q = parent->getDaug(1)->getP4() 
                    + parent->getDaug(2)->getP4();
  double q2 = (q.mass2());

  double fpf,f0f;
  double mesonmass = parent->getDaug(0)->mass();
  double parentmass = parent->mass();

  FormFactors->getscalarff(parent->getId(),
                           parent->getDaug(0)->getId(),
                           q2,
                           mesonmass,
                           &fpf, 
                           &f0f);


  EvtVector4R p4b;
  p4b.set(parent->mass(),0.0,0.0,0.0);
  EvtVector4R p4meson = parent->getDaug(0)->getP4();
  double mdiffoverq2;
  mdiffoverq2 = parentmass*parentmass - mesonmass*mesonmass;
  mdiffoverq2 = mdiffoverq2 / q2;

  EvtVector4C l1,l2;

  EvtId l_num = parent->getDaug(1)->getId();
  EvtVector4C tds;

  if (l_num==EM||l_num==MUM||l_num==TAUM){

    tds = EvtVector4C(fpf*(p4b+p4meson - (mdiffoverq2*(p4b-p4meson)))+
		    + f0f*mdiffoverq2*(p4b-p4meson));

    l1=EvtLeptonVACurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParentNeutrino());
    l2=EvtLeptonVACurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParentNeutrino());
  }
  else{
    if (l_num==EP||l_num==MUP||l_num==TAUP){

      tds = EvtVector4C(fpf*(p4b+p4meson - (mdiffoverq2*(p4b-p4meson)))+
		    + f0f*mdiffoverq2*(p4b-p4meson));

      l1=EvtLeptonVACurrent(parent->getDaug(2)->spParentNeutrino(),
			    parent->getDaug(1)->spParent(0));
      l2=EvtLeptonVACurrent(parent->getDaug(2)->spParentNeutrino(),
			    parent->getDaug(1)->spParent(1));
    }
    else{
      report(Severity::Error,"EvtGen") << "dfnb89agngri wrong lepton number\n";
    }
  }

  amp.vertex(0,l1*tds);
  amp.vertex(1,l2*tds);

}

