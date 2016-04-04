//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2000      Caltech, UCSB
//
// Module: EvtbTosllScalarAmp.cc
//
// Description: Routine to implement bTosll decays to pseudo-scalar
//              mesons. 
//
// Modification history:
//
//    Ryd       January 5,2000       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenModels/EvtbTosllScalarAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenModels/EvtbTosllAmp.hh"
#include "EvtGenModels/EvtbTosllFF.hh"

void EvtbTosllScalarAmp::CalcAmp( EvtParticle *parent,
				  EvtAmp& amp,
				  EvtbTosllFF *formFactors ) {

  //Add the lepton and neutrino 4 momenta to find q2

  EvtVector4R q = parent->getDaug(1)->getP4() 
                    + parent->getDaug(2)->getP4();
  double q2 = (q.mass2());

  double fp(0.),f0(0.),ft(0.);
  double mesonmass = parent->getDaug(0)->mass();
  double parentmass = parent->mass();

  formFactors->getScalarFF(parent->getId(),
                           parent->getDaug(0)->getId(),
                           q2,
                           mesonmass,
                           fp,f0,ft);

  EvtId daught = parent->getDaug(0)->getId();
  bool btod = false;
  bool nnlo = true;
  if (daught == EvtPDL::getId(std::string("pi+")) ||
	  daught == EvtPDL::getId(std::string("pi-")) ||
          daught == EvtPDL::getId(std::string("pi0")) ||
          daught == EvtPDL::getId(std::string("eta")) ||
          daught == EvtPDL::getId(std::string("eta'"))
         )  btod = true;


  EvtVector4R p4b;
  p4b.set(parent->mass(),0.0,0.0,0.0);
 
  EvtVector4C l11,l12;
  EvtVector4C l21,l22;

  EvtVector4C a11,a12;
  EvtVector4C a21,a22;

  EvtId l_num = parent->getDaug(1)->getId();

  EvtVector4C T1,T2;

  EvtVector4R phat=p4b/parentmass;
  EvtVector4R qhat=q/parentmass;

  EvtComplex c7eff = EvtbTosllAmp::GetC7Eff(q2,nnlo);
  EvtComplex c9eff = EvtbTosllAmp::GetC9Eff(q2,nnlo,btod);
  EvtComplex c10eff = EvtbTosllAmp::GetC10Eff(q2,nnlo);

  //double mbhat=1;
  double mbhat=4.4/(parentmass);      
  //double mkhat = 0.15;
  double mkhat=mesonmass/(parentmass);
  double shat=q2/(parentmass*parentmass);

  double fm=(f0-fp)*(1-mkhat*mkhat)/shat;

  EvtComplex aprime;
  aprime = c9eff*fp+2.0*mbhat*c7eff*ft/(1+mkhat);
  EvtComplex bprime;
  bprime = c9eff*fm-2*mbhat*c7eff*ft*(1-mkhat)/shat;

  EvtComplex cprime;
  cprime = c10eff*fp;
  EvtComplex dprime; 
  dprime = c10eff*fm;

  static EvtIdSet leptons("e-","mu-","tau-");
  static EvtIdSet antileptons("e+","mu+","tau+");

  if (leptons.contains(l_num)){

    T1=aprime*phat+bprime*qhat;
    T2=cprime*phat+dprime*qhat;

    l11=EvtLeptonVCurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParent(0));
    l21=EvtLeptonVCurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParent(0));
    l12=EvtLeptonVCurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParent(1));
    l22=EvtLeptonVCurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParent(1));
    a11=EvtLeptonACurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParent(0));
    a21=EvtLeptonACurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParent(0));
    a12=EvtLeptonACurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParent(1));
    a22=EvtLeptonACurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParent(1));
  }
  else{
    if (antileptons.contains(l_num)){

    T1=aprime*phat+bprime*qhat;
    T2=cprime*phat+dprime*qhat;

    l11=EvtLeptonVCurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParent(1));
    l21=EvtLeptonVCurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParent(1));
    l12=EvtLeptonVCurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParent(0));
    l22=EvtLeptonVCurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParent(0));

    a11=EvtLeptonACurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParent(1));
    a21=EvtLeptonACurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParent(1));
    a12=EvtLeptonACurrent(parent->getDaug(1)->spParent(1),
                          parent->getDaug(2)->spParent(0));
    a22=EvtLeptonACurrent(parent->getDaug(1)->spParent(0),
                          parent->getDaug(2)->spParent(0));

    }
    else{
      report(Severity::Error,"EvtGen") << "Wrong lepton number\n";
    }
  }

  amp.vertex(0,0,l11*T1+a11*T2);
  amp.vertex(0,1,l12*T1+a12*T2);
  amp.vertex(1,0,l21*T1+a21*T2);
  amp.vertex(1,1,l22*T1+a22*T2);

}















