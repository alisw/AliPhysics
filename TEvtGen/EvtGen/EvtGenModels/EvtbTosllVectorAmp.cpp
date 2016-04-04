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
// Module: EvtbTosllVectorAmp.cc
//
// Description: Routine to implement bTosll decays to vector
//              mesons. 
//
// Modification history:
//
//    Ryd       January 5,2000       Module created
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
#include "EvtGenModels/EvtbTosllVectorAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenModels/EvtbTosllAmp.hh"
#include "EvtGenModels/EvtbTosllFF.hh"

void EvtbTosllVectorAmp::CalcAmp( EvtParticle *parent,
				  EvtAmp& amp,
				  EvtbTosllFF *formFactors ) {

  //Add the lepton and neutrino 4 momenta to find q2

  EvtVector4R q = parent->getDaug(1)->getP4() 
                    + parent->getDaug(2)->getP4();
  double q2 = (q.mass2());


  double a1,a2,a0,v,t1,t2,t3;
  double mesonmass = parent->getDaug(0)->mass();
  double parentmass = parent->mass();

  formFactors->getVectorFF(parent->getId(),
                           parent->getDaug(0)->getId(),
                           q2,
                           mesonmass,
                           a1,a2,a0,v,t1,t2,t3);


  EvtId daught = parent->getDaug(0)->getId();
  bool btod = false;
  bool nnlo = true;
 if 
  (
   daught == EvtPDL::getId(std::string("rho+")) ||
   daught == EvtPDL::getId(std::string("rho-")) ||
   daught == EvtPDL::getId(std::string("rho0")) ||
   daught == EvtPDL::getId(std::string("omega"))
  ) 
   btod = true;

  EvtVector4R p4b;
  p4b.set(parent->mass(),0.0,0.0,0.0);
  EvtVector4R p4meson = parent->getDaug(0)->getP4();
 
  EvtVector4C l11,l12;
  EvtVector4C l21,l22;

  EvtVector4C a11,a12;
  EvtVector4C a21,a22;

  EvtId parentID = parent->getId();

      //EvtId l_num = parent->getDaug(1)->getId();

  EvtVector4R pbhat=p4b/parentmass;
  EvtVector4R qhat=q/parentmass;
  EvtVector4R pkstarhat=p4meson/parentmass;
  EvtVector4R phat=pbhat+pkstarhat;

  EvtComplex c7eff = EvtbTosllAmp::GetC7Eff(q2,nnlo);
  EvtComplex c9eff = EvtbTosllAmp::GetC9Eff(q2,nnlo,btod);
  EvtComplex c10eff = EvtbTosllAmp::GetC10Eff(q2,nnlo);
  EvtComplex uniti(0.0,1.0);

  double mhatb=4.4/(parentmass); 
  double mhatkstar=mesonmass/(parentmass);
  double shat=q2/(parentmass*parentmass);


  EvtComplex a;
  a=c9eff*v*2/(1+mhatkstar)+4*mhatb*c7eff*t1/shat;
  EvtComplex b;
  b=(1+mhatkstar)*(c9eff*a1+2*mhatb*(1-mhatkstar)*c7eff*t2/shat);
  EvtComplex c;
  c=((1-mhatkstar)*c9eff*a2+
	    2*mhatb*c7eff*(t3+(1-mhatkstar*mhatkstar)*t2/shat))/
    (1-mhatkstar*mhatkstar);
  EvtComplex d;
  d=(c9eff*((1+mhatkstar)*a1-(1-mhatkstar)*a2-2*mhatkstar*a0)
	    -2*mhatb*c7eff*t3)/shat;
  EvtComplex e;
  e=2*c10eff*v/(1+mhatkstar);
  EvtComplex f;
  f=(1+mhatkstar)*c10eff*a1;
  EvtComplex g;
  g=c10eff*a2/(1+mhatkstar);
  EvtComplex h;
  h=c10eff*((1+mhatkstar)*a1-(1-mhatkstar)*a2-2*mhatkstar*a0)/shat;
  
  EvtTensor4C T1,T2;
  
  static EvtIdSet bmesons("B-","anti-B0","anti-B_s0");
  static EvtIdSet bbarmesons("B+","B0","B_s0");
  
  EvtParticle* lepPlus=0;
  EvtParticle* lepMinus=0;
  
  int charge1 = EvtPDL::chg3(parent->getDaug(1)->getId());
  int charge2 = EvtPDL::chg3(parent->getDaug(2)->getId());
  
  lepPlus = (charge1 > charge2) ? parent->getDaug(1) : parent->getDaug(2);
  lepMinus = (charge1 < charge2) ? parent->getDaug(1) : parent->getDaug(2);

  if (bmesons.contains(parentID)) {

    T1=a*dual(EvtGenFunctions::directProd(pbhat,pkstarhat))
       -b*uniti*EvtTensor4C::g()
       +c*uniti*EvtGenFunctions::directProd(pbhat,phat)
       +d*uniti*EvtGenFunctions::directProd(pbhat,qhat);
    
    T2=e*dual(EvtGenFunctions::directProd(pbhat,pkstarhat))
       -f*uniti*EvtTensor4C::g()
       +g*uniti*EvtGenFunctions::directProd(pbhat,phat)
       +h*uniti*EvtGenFunctions::directProd(pbhat,qhat);
    
    l11=EvtLeptonVCurrent(lepPlus->spParent(0),
                          lepMinus->spParent(0));
    l21=EvtLeptonVCurrent(lepPlus->spParent(1),
                          lepMinus->spParent(0));
    l12=EvtLeptonVCurrent(lepPlus->spParent(0), 
                          lepMinus->spParent(1));
    l22=EvtLeptonVCurrent(lepPlus->spParent(1),
                          lepMinus->spParent(1));

    a11=EvtLeptonACurrent(lepPlus->spParent(0),
                          lepMinus->spParent(0));
    a21=EvtLeptonACurrent(lepPlus->spParent(1),
                          lepMinus->spParent(0));
    a12=EvtLeptonACurrent(lepPlus->spParent(0),
                          lepMinus->spParent(1));
    a22=EvtLeptonACurrent(lepPlus->spParent(1),
                          lepMinus->spParent(1));

  } else {
    
    if (bbarmesons.contains(parentID)) {

      T1=-a*dual(EvtGenFunctions::directProd(pbhat,pkstarhat))
         -b*uniti*EvtTensor4C::g()
         +c*uniti*EvtGenFunctions::directProd(pbhat,phat)
         +d*uniti*EvtGenFunctions::directProd(pbhat,qhat);
      
      T2=-e*dual(EvtGenFunctions::directProd(pbhat,pkstarhat))
         -f*uniti*EvtTensor4C::g()
         +g*uniti*EvtGenFunctions::directProd(pbhat,phat)
         +h*uniti*EvtGenFunctions::directProd(pbhat,qhat);
      
      l11=EvtLeptonVCurrent(lepPlus->spParent(1),
                            lepMinus->spParent(1));
      l21=EvtLeptonVCurrent(lepPlus->spParent(0),
                            lepMinus->spParent(1));
      l12=EvtLeptonVCurrent(lepPlus->spParent(1),
                            lepMinus->spParent(0));
      l22=EvtLeptonVCurrent(lepPlus->spParent(0),
                            lepMinus->spParent(0));
      
      a11=EvtLeptonACurrent(lepPlus->spParent(1),
                            lepMinus->spParent(1));
      a21=EvtLeptonACurrent(lepPlus->spParent(0),
                            lepMinus->spParent(1));
      a12=EvtLeptonACurrent(lepPlus->spParent(1),
                            lepMinus->spParent(0));
      a22=EvtLeptonACurrent(lepPlus->spParent(0),
                            lepMinus->spParent(0));
      
    }
    else{
      report(Severity::Error,"EvtGen") << "Wrong lepton number\n";
      T1.zero(); T2.zero(); // Set all tensor terms to zero.
    }    
  }


  int i;

  for(i=0;i<3;i++){
    EvtVector4C eps=parent->getDaug(0)->epsParent(i).conj();

    EvtVector4C E1=T1.cont1(eps);
    EvtVector4C E2=T2.cont1(eps);

    amp.vertex(i,0,0,l11*E1+a11*E2);
    amp.vertex(i,0,1,l12*E1+a12*E2);
    amp.vertex(i,1,0,l21*E1+a21*E2);
    amp.vertex(i,1,1,l22*E1+a22*E2);
  } 
}







