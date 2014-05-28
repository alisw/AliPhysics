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
// Module: EvtBTo4piCP.cc
//
// Description: Routine to decay B->pi+ pi- pi+ pi-.
//
// Modification history:
//
//    RYD     March 2, 1997         Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBTo4piCP.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtConst.hh"
#include <string>

EvtBTo4piCP::~EvtBTo4piCP() {}


EvtComplex EvtAmpA2(const EvtVector4R& p4pi1,const EvtVector4R& p4pi2,
		    const EvtVector4R& p4pi3,const EvtVector4R& p4pi4){

  //added by Lange Jan4,2000
  static EvtId A2M=EvtPDL::getId("a_2-");
  static EvtId RHO0=EvtPDL::getId("rho0");

  EvtVector4R p4a2,p4rho,p4b;

  p4rho=p4pi1+p4pi2;

  p4a2=p4rho+p4pi3;

  p4b=p4a2+p4pi4;

  EvtVector4R p4b_a2,p4rho_a2,p4pi1_a2,p4a2_a2;

  p4b_a2=boostTo(p4b,p4a2);
  p4rho_a2=boostTo(p4rho,p4a2);
  p4pi1_a2=boostTo(p4pi1,p4a2);
  p4a2_a2=boostTo(p4a2,p4a2);

  EvtVector4R p4pi1_rho;

  p4pi1_rho=boostTo(p4pi1_a2,p4rho_a2);

  EvtVector4R vb,vrho,vpi,t;

  vb=p4b_a2/p4b_a2.d3mag();
  vrho=p4rho_a2/p4rho_a2.d3mag();
  vpi=p4pi1_rho/p4pi1_rho.d3mag();

  t.set(1.0,0.0,0.0,0.0);

  //  EvtComplex amp_a1,amp_a2;
  EvtComplex amp_a2;
 
  //  double bwm_a1=EvtPDL::getMeanMass(A1M);
  //  double gamma_a1=EvtPDL::getWidth(A1M);
  double bwm_a2=EvtPDL::getMeanMass(A2M);
  double gamma_a2=EvtPDL::getWidth(A2M);
  double bwm_rho=EvtPDL::getMeanMass(RHO0);
  double gamma_rho=EvtPDL::getWidth(RHO0);

  amp_a2=(sqrt(gamma_a2/EvtConst::twoPi)/
    ((p4a2).mass()-bwm_a2-EvtComplex(0.0,0.5*gamma_a2)))*
         (sqrt(gamma_rho/EvtConst::twoPi)/
    ((p4rho).mass()-bwm_rho-EvtComplex(0.0,0.5*gamma_rho)));

  return amp_a2*
    (vb.get(1)*vrho.get(1)+vb.get(2)*vrho.get(2)+vb.get(3)*vrho.get(3))*
    (
     vpi.get(1)*(vb.get(2)*vrho.get(3)-vb.get(3)*vrho.get(2))+
     vpi.get(2)*(vb.get(3)*vrho.get(1)-vb.get(1)*vrho.get(3))+
     vpi.get(3)*(vb.get(1)*vrho.get(2)-vb.get(2)*vrho.get(1))
     );

}

EvtComplex EvtAmpA1(const EvtVector4R& p4pi1,const EvtVector4R& p4pi2,
		    const EvtVector4R& p4pi3,const EvtVector4R& p4pi4){

  //added by Lange Jan4,2000
  static EvtId A1M=EvtPDL::getId("a_1-");
  static EvtId RHO0=EvtPDL::getId("rho0");

  EvtVector4R p4a1,p4rho,p4b;

  p4rho=p4pi1+p4pi2;

  p4a1=p4rho+p4pi3;

  p4b=p4a1+p4pi4;

  EvtVector4R p4b_a1,p4rho_a1,p4pi1_a1,p4a1_a1;

  p4b_a1=boostTo(p4b,p4a1);
  p4rho_a1=boostTo(p4rho,p4a1);
  p4pi1_a1=boostTo(p4pi1,p4a1);
  p4a1_a1=boostTo(p4a1,p4a1);

  EvtVector4R p4pi1_rho;

  p4pi1_rho=boostTo(p4pi1_a1,p4rho_a1);

  EvtVector4R vb,vrho,vpi,t;

  vb=p4b_a1/p4b_a1.d3mag();
  vrho=p4rho_a1/p4rho_a1.d3mag();
  vpi=p4pi1_rho/p4pi1_rho.d3mag();

  t.set(1.0,0.0,0.0,0.0);

  EvtComplex amp_a1;
 
  double bwm_a1=EvtPDL::getMeanMass(A1M);
  double gamma_a1=EvtPDL::getWidth(A1M);
  //  double bwm_a2=EvtPDL::getMeanMass(A2M);
  //  double gamma_a2=EvtPDL::getWidth(A2M);
  double bwm_rho=EvtPDL::getMeanMass(RHO0);
  double gamma_rho=EvtPDL::getWidth(RHO0);

  amp_a1=(sqrt(gamma_a1/EvtConst::twoPi)/
    ((p4a1).mass()-bwm_a1-EvtComplex(0.0,0.5*gamma_a1)))*
         (sqrt(gamma_rho/EvtConst::twoPi)/
    ((p4rho).mass()-bwm_rho-EvtComplex(0.0,0.5*gamma_rho)));

  return amp_a1*
    (vb.get(1)*vpi.get(1)+vb.get(2)*vpi.get(2)+vb.get(3)*vpi.get(3));

}


std::string EvtBTo4piCP::getName(){
 
  return "BTO4PI_CP";     

}


EvtDecayBase* EvtBTo4piCP::clone(){

  return new EvtBTo4piCP;

}

void EvtBTo4piCP::init(){

  // check that there are 18 arguments
  checkNArg(18);
  checkNDaug(4);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
  checkSpinDaughter(2,EvtSpinType::SCALAR);
  checkSpinDaughter(3,EvtSpinType::SCALAR);
}

void EvtBTo4piCP::decay( EvtParticle *p){

  //added by Lange Jan4,2000
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");


  double t;
  EvtId other_b;

  EvtCPUtil::getInstance()->OtherB(p,t,other_b,0.5);
  
  p->initializePhaseSpace(getNDaug(),getDaugs());
  EvtVector4R mom1 = p->getDaug(0)->getP4();
  EvtVector4R mom2 = p->getDaug(1)->getP4();
  EvtVector4R mom3 = p->getDaug(2)->getP4();
  EvtVector4R mom4 = p->getDaug(3)->getP4();

  //  double alpha=getArg(0);
  //double dm=getArg(1);

   EvtComplex amp;

   EvtComplex A,Abar;


   EvtComplex A_a1p,Abar_a1p,A_a2p,Abar_a2p;
   EvtComplex A_a1m,Abar_a1m,A_a2m,Abar_a2m;

   A_a1p=EvtComplex(getArg(2)*cos(getArg(3)),getArg(2)*sin(getArg(3)));
   Abar_a1p=EvtComplex(getArg(4)*cos(getArg(5)),getArg(4)*sin(getArg(5)));

   A_a2p=EvtComplex(getArg(6)*cos(getArg(7)),getArg(6)*sin(getArg(7)));
   Abar_a2p=EvtComplex(getArg(8)*cos(getArg(9)),getArg(8)*sin(getArg(9)));

   A_a1m=EvtComplex(getArg(10)*cos(getArg(11)),getArg(10)*sin(getArg(11)));
   Abar_a1m=EvtComplex(getArg(12)*cos(getArg(13)),getArg(12)*sin(getArg(13)));

   A_a2m=EvtComplex(getArg(14)*cos(getArg(15)),getArg(14)*sin(getArg(15)));
   Abar_a2m=EvtComplex(getArg(16)*cos(getArg(17)),getArg(16)*sin(getArg(17)));

   EvtComplex a2p_amp=EvtAmpA2(mom1,mom2,mom3,mom4)+
                      EvtAmpA2(mom1,mom4,mom3,mom2)+
                      EvtAmpA2(mom3,mom2,mom1,mom4)+
                      EvtAmpA2(mom3,mom4,mom1,mom2);

   EvtComplex a2m_amp=EvtAmpA2(mom2,mom3,mom4,mom1)+
                      EvtAmpA2(mom2,mom1,mom4,mom3)+
                      EvtAmpA2(mom4,mom3,mom2,mom1)+
                      EvtAmpA2(mom4,mom1,mom2,mom3);

   EvtComplex a1p_amp=EvtAmpA1(mom1,mom2,mom3,mom4)+
                      EvtAmpA1(mom1,mom4,mom3,mom2)+
                      EvtAmpA1(mom3,mom2,mom1,mom4)+
                      EvtAmpA1(mom3,mom4,mom1,mom2);

   EvtComplex a1m_amp=EvtAmpA1(mom2,mom3,mom4,mom1)+
                      EvtAmpA1(mom2,mom1,mom4,mom3)+
                      EvtAmpA1(mom4,mom3,mom2,mom1)+
                      EvtAmpA1(mom4,mom1,mom2,mom3);


   A=A_a2p*a2p_amp+A_a1p*a1p_amp+
     A_a2m*a2m_amp+A_a1m*a1m_amp;
   Abar=Abar_a2p*a2p_amp+Abar_a1p*a1p_amp+
        Abar_a2m*a2m_amp+Abar_a1m*a1m_amp;


   if (other_b==B0B){
     amp=A*cos(getArg(1)*t/(2*EvtConst::c))+
       EvtComplex(cos(-2.0*getArg(0)),sin(-2.0*getArg(0)))*
       getArg(2)*EvtComplex(0.0,1.0)*Abar*sin(getArg(1)*t/(2*EvtConst::c));
   }
   if (other_b==B0){
     amp=A*EvtComplex(cos(2.0*getArg(0)),sin(2.0*getArg(0)))*
       EvtComplex(0.0,1.0)*sin(getArg(1)*t/(2*EvtConst::c))+       
       getArg(2)*Abar*cos(getArg(1)*t/(2*EvtConst::c));
   }
 
   vertex(amp);

   return ;
}

