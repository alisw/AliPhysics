//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2001      Caltech, UCSB
//
// Module: EvtSVSCPLH.cc
//
// Description: The decay of a scalar to a scalar and a vector particle are 
//              performed with CP violation and different widths for
//              the cp even and odd states. E.g. B->J/psi K_S.
//
// Modification history:
//
//    Ryd       March 29, 2001       Module created
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
#include "EvtGenModels/EvtSVSCPLH.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtVector4C.hh"
using std::endl;

EvtSVSCPLH::~EvtSVSCPLH() {}

std::string EvtSVSCPLH::getName(){

  return "SVS_CPLH";     

}


EvtDecayBase* EvtSVSCPLH::clone(){

  return new EvtSVSCPLH;

}

void EvtSVSCPLH::init(){

  // check that there are 8 arguments
  checkNArg(8);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);

  static double ctau=EvtPDL::getctau(EvtPDL::getId("B0"));

  // hbar/s
  _dm=getArg(0);
  _dgamma=EvtConst::c*getArg(1)/ctau;

  _qop=getArg(2)*exp(EvtComplex(0.0,getArg(3)));

  _poq=1.0/_qop;

  _Af=getArg(4)*exp(EvtComplex(0.0,getArg(5)));
  _Abarf=getArg(6)*exp(EvtComplex(0.0,getArg(7)));
  
  if (verbose()){
    report(Severity::Info,"EvtGen")<<":EvtSVSCPLH:dm="<<_dm<<endl;
    report(Severity::Info,"EvtGen")<<":EvtSVSCPLH:dGamma="<<_dgamma<<endl;
    report(Severity::Info,"EvtGen")<<":EvtSVSCPLH:q/p="<<_qop<<endl;
    report(Severity::Info,"EvtGen")<<":EvtSVSCPLH:Af="<<_Af<<endl;
    report(Severity::Info,"EvtGen")<<":EvtSVSCPLH:Abarf="<<_Abarf<<endl;
  }


}

void EvtSVSCPLH::initProbMax(){

  //This is probably not quite right, but it should do as a start...
  //Anders

  setProbMax(4.0*(getArg(4)*getArg(4)+getArg(6)*getArg(6)));

}

void EvtSVSCPLH::decay( EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  double t;
  EvtId other_b;

  EvtCPUtil::getInstance()->OtherB(p,t,other_b,0.5);

  //convert time from mm to seconds
  t/=EvtConst::c;

  //sign convention is dm=Mheavy-Mlight
  //                   dGamma=Gammalight-Gammaheavy
  //such that in the standard model both of these are positive.
  EvtComplex gp=0.5*(exp(EvtComplex(0.25*t*_dgamma,-0.5*t*_dm))+exp(EvtComplex(-0.25*t*_dgamma,0.5*t*_dm)));
  EvtComplex gm=0.5*(exp(EvtComplex(0.25*t*_dgamma,-0.5*t*_dm))-exp(EvtComplex(-0.25*t*_dgamma,0.5*t*_dm)));

  EvtComplex amp;

  if (other_b==B0B){
    amp=gp*_Af+_qop*gm*_Abarf;
  }
  else if (other_b==B0){
    amp=gp*_Abarf+_poq*gm*_Af;
  }
  else{
    report(Severity::Error,"EvtGen") << "other_b was not B0 or B0B!"<<endl;
    ::abort();
  }

  EvtVector4R p4_parent=p->getP4Restframe();;
  
  double norm=p->getDaug(0)->mass()/(p->getDaug(0)->getP4().d3mag()*p4_parent.mass());

  EvtParticle* v=p->getDaug(0);

  vertex(0,amp*norm*(p4_parent*(v->epsParent(0))));
  vertex(1,amp*norm*(p4_parent*(v->epsParent(1))));
  vertex(2,amp*norm*(p4_parent*(v->epsParent(2))));
  

  return ;
}



