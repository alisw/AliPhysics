//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1999      Caltech, UCSB
//
// Module: EvtSVVCPLH.cc
//
// Description: The decay of a scalar to two vector particles are 
//              performed with CP violation and different widths for
//              the cpe even and od states. E.g. Bs->J/psi phi.
//
// Modification history:
//
//    RYD       November 5, 1999       Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtSVVHelAmp.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtSVVCPLH.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"
using std::endl;

EvtSVVCPLH::~EvtSVVCPLH() {}

std::string EvtSVVCPLH::getName(){

  return "SVV_CPLH";     

}


EvtDecayBase* EvtSVVCPLH::clone(){

  return new EvtSVVCPLH;

}

void EvtSVVCPLH::init(){

  // check that there are 9 arguments
  checkNArg(9);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::VECTOR);

}

void EvtSVVCPLH::initProbMax(){

  //This is probably not quite right, but it should do as a start...
  //Anders

  setProbMax(2*(getArg(3)*getArg(3)+getArg(5)*getArg(5)+getArg(7)*getArg(7)));

}

void EvtSVVCPLH::decay( EvtParticle *p){

  //added by Lange Jan4,2000
  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");

  double t;
  EvtId other_b;

  EvtCPUtil::getInstance()->OtherB(p,t,other_b);

  EvtComplex G0P,G1P,G1M;
  
  G1P=EvtComplex(getArg(3)*cos(getArg(4)),getArg(3)*sin(getArg(4)));
  G0P=EvtComplex(getArg(5)*cos(getArg(6)),getArg(5)*sin(getArg(6)));
  G1M=EvtComplex(getArg(7)*cos(getArg(8)),getArg(7)*sin(getArg(8)));

  EvtComplex lambda_km=EvtComplex(cos(2*getArg(0)),sin(2*getArg(0)));

  double cdmt=cos(getArg(1)*t/(2*EvtConst::c));
  double sdmt=sin(getArg(1)*t/(2*EvtConst::c));

  EvtComplex cG0P,cG1P,cG1M;

  static double ctauL=EvtPDL::getctau(EvtPDL::getId("B_s0L"));
  static double ctauH=EvtPDL::getctau(EvtPDL::getId("B_s0H"));

  //I'm not sure if the fabs() is right when t can be
  //negative as in the case of Bs produced coherently.
  double pt=1;
  double mt=exp(-fabs(t*(ctauL-ctauH)/(ctauL*ctauH)));

  if (other_b==BSB){
    cG0P=pt*G0P*(cdmt+lambda_km*EvtComplex(0.0,getArg(2)*sdmt));
    cG1P=pt*G1P*(cdmt+lambda_km*EvtComplex(0.0,getArg(2)*sdmt));
    cG1M=mt*G1M*(cdmt-lambda_km*EvtComplex(0.0,getArg(2)*sdmt));
  }
  else if (other_b==BS0){
    cG0P=pt*G0P*(cdmt+(1.0/lambda_km)*EvtComplex(0.0,getArg(2)*sdmt));
    cG1P=pt*G1P*(cdmt+(1.0/lambda_km)*EvtComplex(0.0,getArg(2)*sdmt));
    cG1M=-mt*G1M*(cdmt-(1.0/lambda_km)*EvtComplex(0.0,getArg(2)*sdmt));
  }
  else{
    report(Severity::Error,"EvtGen") << "other_b was not BSB or BS0!"<<endl;
    ::abort();
  }

   EvtComplex A0,AP,AM;

   A0=cG0P/sqrt(2.0);
   AP=(cG1P+cG1M)/sqrt(2.0); 
   AM=(cG1P-cG1M)/sqrt(2.0); 

   EvtSVVHelAmp::SVVHel(p,_amp2,getDaug(0),getDaug(1),AP,A0,AM);

  return ;
}



