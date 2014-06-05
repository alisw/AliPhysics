//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtSVPCP.cc
//
// Description: Routine to decay scalar -> vectors+photon
//              including CP violation effects
//
// Modification history:
//
//    Maurizio pierini   Nov 11, 2003       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtTensor3C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenModels/EvtSVPCP.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"

EvtSVPCP::~EvtSVPCP() {}

std::string EvtSVPCP::getName(){

  return "SVP_CP";     

}


EvtDecayBase* EvtSVPCP::clone(){

  return new EvtSVPCP;

}

void EvtSVPCP::initProbMax(){

  setProbMax(2*(getArg(3)*getArg(3)+getArg(5)*getArg(5)));

}


void EvtSVPCP::init(){

  // check that there are 7 arguments
  checkNArg(7);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::PHOTON);

}

void EvtSVPCP::decay( EvtParticle *p ){

  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  double t;
  EvtId other_b;

  EvtCPUtil::getInstance()->OtherB(p,t,other_b,0.5);

  EvtComplex G1P,G1M, G1_T_even, G1_T_odd;

  double norm = getArg(3)*getArg(3)+getArg(5)*getArg(5);
  
  G1P=EvtComplex(getArg(3)*cos(getArg(4))/norm,getArg(3)*sin(getArg(4))/norm);
  G1M=EvtComplex(getArg(5)*cos(getArg(6))/norm,getArg(5)*sin(getArg(6))/norm);

  G1_T_even = (G1P+G1M)/sqrt(2.0);
  G1_T_odd  = (G1P-G1M)/sqrt(2.0);
  
  EvtComplex lambda_km =EvtComplex(cos(-2*getArg(0)),sin(-2*getArg(0)));
  
  double cdmt=cos(getArg(1)*t/(2*EvtConst::c));
  double sdmt=sin(getArg(1)*t/(2*EvtConst::c));

  EvtComplex cG1_T_even,cG1_T_odd;
  
  if (other_b==B0B){
    cG1_T_even = G1_T_even*(cdmt+lambda_km*EvtComplex(0.0,getArg(2)*sdmt));
    cG1_T_odd  = G1_T_odd*(cdmt-lambda_km*EvtComplex(0.0,getArg(2)*sdmt));
  }
  if (other_b==B0){
    cG1_T_even = G1_T_even*(cdmt+(1.0/lambda_km)*EvtComplex(0.0,getArg(2)*sdmt));
    cG1_T_odd  =-G1_T_odd*(cdmt-(1.0/lambda_km)*EvtComplex(0.0,getArg(2)*sdmt));
  }
  
  EvtComplex hp, hm, h0;

  // This part is adopted from EvtSVVHel and since there is
  // a photon that can not have helicity 0 this is put in by 
  // setting the h0 amplitude to 0.
  hm=(cG1_T_even-cG1_T_odd)/sqrt(2.0); 
  hp=(cG1_T_even+cG1_T_odd)/sqrt(2.0); 
  h0=EvtComplex(0.0,0.0);

  EvtParticle *v1,*ph;

  p->initializePhaseSpace(getNDaug(),getDaugs());
  v1 = p->getDaug(0);
  ph = p->getDaug(1);
  EvtVector4R momv1 = v1->getP4();
  EvtVector4R momph = ph->getP4();

  EvtTensor4C d,g;

  g.setdiag(1.0,-1.0,-1.0,-1.0);

  EvtVector4R v,vp;

  v=momv1/momv1.d3mag();
  vp=(momv1+momph)/(momv1+momph).mass();   

  d=((1.0/sqrt(3.0))*(h0-(hp+hm))*(-1.0/sqrt(3.0)))*g+
    ((1.0/sqrt(2.0))*(hp-hm)*EvtComplex(0.0,1.0)*(sqrt(1.0/2.0)))*dual(EvtGenFunctions::directProd(v,vp))+
    (sqrt(2.0/3.0)*(h0+0.5*(hp+hm))*sqrt(3.0/2.0))*(EvtGenFunctions::directProd(v,v)+(1.0/3.0)*g);

  EvtVector4C ep0,ep1,ep2;  
  
  ep0=d.cont1(v1->eps(0).conj());
  ep1=d.cont1(v1->eps(1).conj());
  ep2=d.cont1(v1->eps(2).conj());

  EvtVector4C ep20,ep21,ep22;

  ep20=ph->epsParentPhoton(0).conj();  
  ep21=ph->epsParentPhoton(1).conj();  

  vertex(0,0,ep0*ep20);
  vertex(0,1,ep0*ep21);
  
  vertex(1,0,ep1*ep20);
  vertex(1,1,ep1*ep21);
   
  vertex(2,0,ep2*ep20);
  vertex(2,1,ep2*ep21);

				 
  return ;

}

std::string EvtSVPCP::getParamName(int i) {
  switch(i) {
  case 0:
    return "weakPhase";
  case 1:
    return "deltaM";
  case 2:
    return "finalStateCP";
  case 3:
    return "Af";
  case 4:
    return "AfPhase";
  case 5:
    return "Abarf";
  case 6:
    return "AbarfPhase";
  default:
    return "";
  }
}

std::string EvtSVPCP::getParamDefault(int i) {
  switch(i) {
  case 3:
    return "1.0";
  case 4:
    return "0.0";
  case 5:
    return "1.0";
  case 6:
    return "0.0";
  default:
    return "";
  }
}
