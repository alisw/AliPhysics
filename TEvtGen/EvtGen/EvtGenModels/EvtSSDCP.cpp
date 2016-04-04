//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2001      Caltech
//
// Module: EvtSSDCP.cc
//
// Description: See EvtSSDCP.hh
//
// Modification history:
//
//    RYD       August 12, 2001        Module created
//    F. Sandrelli, Fernando M-V  March 1, 2002     Debugged and added z parameter (CPT violation)
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenModels/EvtSSDCP.hh"
#include "EvtGenBase/EvtIncoherentMixing.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"
using std::endl;

EvtSSDCP::~EvtSSDCP() {}

std::string EvtSSDCP::getName(){

  return "SSD_CP";     

}


EvtDecayBase* EvtSSDCP::clone(){

  return new EvtSSDCP;

}

void EvtSSDCP::init(){

  // check that there are 8 or 12 or 14 arguments

  checkNArg(14,12,8);
  checkNDaug(2);

  EvtSpinType::spintype d1type=EvtPDL::getSpinType(getDaug(0));
  EvtSpinType::spintype d2type=EvtPDL::getSpinType(getDaug(1));

  // Check it is a B0 or B0s
  if ( ( getParentId() != EvtPDL::getId( "B0" ) )
       && ( getParentId() != EvtPDL::getId( "anti-B0" ) ) 
       && ( getParentId() != EvtPDL::getId( "B_s0" ) )
       && ( getParentId() != EvtPDL::getId( "anti-B_s0" ) ) ) {
    report( Severity::Error , "EvtGen" ) << "EvtSSDCP only decays B0 and B0s" 
                               << std::endl ;
    ::abort() ;
  }

  
  if ( (!(d1type == EvtSpinType::SCALAR || d2type == EvtSpinType::SCALAR))||
       (!((d2type==EvtSpinType::SCALAR)||(d2type==EvtSpinType::VECTOR)||(d2type==EvtSpinType::TENSOR)))||
       (!((d1type==EvtSpinType::SCALAR)||(d1type==EvtSpinType::VECTOR)||(d1type==EvtSpinType::TENSOR)))
       ) {
    report(Severity::Error,"EvtGen") << "EvtSSDCP generator expected "
                           << "one of the daugters to be a scalar, the other either scalar, vector, or tensor, found:"
			   << EvtPDL::name(getDaug(0)).c_str()<<" and "<<EvtPDL::name(getDaug(1)).c_str()<<endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<endl;
    ::abort();
  }

  _dm=getArg(0)/EvtConst::c; //units of 1/mm

  _dgog=getArg(1);

  
  _qoverp=getArg(2)*EvtComplex(cos(getArg(3)),sin(getArg(3)));
  _poverq=1.0/_qoverp;

  _A_f=getArg(4)*EvtComplex(cos(getArg(5)),sin(getArg(5)));

  _Abar_f=getArg(6)*EvtComplex(cos(getArg(7)),sin(getArg(7)));
  
  if (getNArg()>=12){
    _eigenstate=false;
    _A_fbar=getArg(8)*EvtComplex(cos(getArg(9)),sin(getArg(9)));
    _Abar_fbar=getArg(10)*EvtComplex(cos(getArg(11)),sin(getArg(11)));
  }
  else{
    //I'm somewhat confused about this. For a CP eigenstate set the
    //amplitudes to the same. For a non CP eigenstate CPT invariance
    //is enforced. (ryd)
    if (
	(getDaug(0)==EvtPDL::chargeConj(getDaug(0))&&
	 getDaug(1)==EvtPDL::chargeConj(getDaug(1)))||
	(getDaug(0)==EvtPDL::chargeConj(getDaug(1))&&
	 getDaug(1)==EvtPDL::chargeConj(getDaug(0)))){
      _eigenstate=true;
    }else{
      _eigenstate=false;
      _A_fbar=conj(_Abar_f);
      _Abar_fbar=conj(_A_f);
    }
  }

  //FS: new check for z 
  if (getNArg()==14){ //FS Set _z parameter if provided else set it 0
    _z=EvtComplex(getArg(12),getArg(13));
  }
  else{
    _z=EvtComplex(0.0,0.0);
  }

  // FS substituted next 2 lines...


  //
  //  _gamma=EvtPDL::getctau(EvtPDL::getId("B0"));  //units of 1/mm
  //_dgamma=_gamma*0.5*_dgog;
  //
  // ...with:

  if ( ( getParentId() == EvtPDL::getId("B0") ) || 
       ( getParentId() == EvtPDL::getId("anti-B0") ) ) {
    _gamma=1./EvtPDL::getctau(EvtPDL::getId("B0")); //gamma/c (1/mm)
  }
  else{
    _gamma=1./EvtPDL::getctau(EvtPDL::getId("B_s0")) ;
  }

  _dgamma=_gamma*_dgog;  //dgamma/c (1/mm) 

  if (verbose()){
    report(Severity::Info,"EvtGen") << "SSD_CP will generate CP/CPT violation:"
			  << endl << endl
			  << "    " << EvtPDL::name(getParentId()).c_str() << " --> "
			  << EvtPDL::name(getDaug(0)).c_str() << " + "
			  << EvtPDL::name(getDaug(1)).c_str() << endl << endl
			  << "using parameters:" << endl << endl
			  << "  delta(m)  = " << _dm << " hbar/ps" << endl
			  << "dGamma      = "  << _dgamma <<" ps-1" <<endl
			  << "       q/p  = " << _qoverp << endl  
			  << "        z  = " << _z << endl  
			  << "       tau  = " << 1./_gamma << " ps" << endl;
  }

}

void EvtSSDCP::initProbMax() {
  double theProbMax = 
    abs(_A_f) * abs(_A_f) +
    abs(_Abar_f) * abs(_Abar_f) +
    abs(_A_fbar) * abs(_A_fbar) +
    abs(_Abar_fbar) * abs(_Abar_fbar);

  if (_eigenstate) theProbMax*=2;

  EvtSpinType::spintype d2type=EvtPDL::getSpinType(getDaug(1));
  EvtSpinType::spintype d1type=EvtPDL::getSpinType(getDaug(0));
  if (d1type==EvtSpinType::TENSOR||d2type==EvtSpinType::TENSOR) theProbMax*=10;

  setProbMax(theProbMax);
}

void EvtSSDCP::decay( EvtParticle *p){

  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");
  
  static EvtId B0s = EvtPDL::getId("B_s0");
  static EvtId B0Bs = EvtPDL::getId("anti-B_s0");

  double t;
  EvtId other_b;
  EvtId daugs[2];

  int flip=0;
  if (!_eigenstate){
    if (EvtRandom::Flat(0.0,1.0)<0.5) flip=1;
  }

  if (!flip) {
    daugs[0]=getDaug(0);
    daugs[1]=getDaug(1);
  }
  else{
    daugs[0]=EvtPDL::chargeConj(getDaug(0));
    daugs[1]=EvtPDL::chargeConj(getDaug(1));
  }

  EvtParticle *d;
  p->initializePhaseSpace(2, daugs);

  EvtComplex amp;


  EvtCPUtil::getInstance()->OtherB(p,t,other_b,0.5); // t is c*Dt (mm)
//  EvtIncoherentMixing::OtherB( p , t , other_b , 0.5 ) ;
  
  //if (flip) t=-t;

  //FS We assume DGamma=GammaLow-GammaHeavy and Dm=mHeavy-mLow
  EvtComplex expH=exp(-EvtComplex(-0.25*_dgamma*t,0.5*_dm*t));
  EvtComplex expL=exp(EvtComplex(-0.25*_dgamma*t,0.5*_dm*t));
  //FS Definition of gp and gm
  EvtComplex gp=0.5*(expL+expH);
  EvtComplex gm=0.5*(expL-expH);
  //FS Calculation os sqrt(1-z^2) 
  EvtComplex sqz=sqrt(abs(1-_z*_z))*exp(EvtComplex(0,arg(1-_z*_z)/2));
  
  //EvtComplex BB=0.5*(expL+expH);                  // <B0|B0(t)>
  //EvtComplex barBB=_qoverp*0.5*(expL-expH);       // <B0bar|B0(t)>
  //EvtComplex BbarB=_poverq*0.5*(expL-expH);       // <B0|B0bar(t)>
  //EvtComplex barBbarB=BB;                         // <B0bar|B0bar(t)>
  //  FS redefinition of these guys... (See BAD #188 eq.35 for ref.)
  //  q/p is taken as in the BaBar Phys. Book (opposite sign wrt ref.)
  EvtComplex BB=gp+_z*gm;                 // <B0|B0(t)>
  EvtComplex barBB=sqz*_qoverp*gm;       // <B0bar|B0(t)>
  EvtComplex BbarB=sqz*_poverq*gm;       // <B0|B0bar(t)>
  EvtComplex barBbarB=gp-_z*gm;           // <B0bar|B0bar(t)>

  if (!flip){
    if (other_b==B0B||other_b==B0Bs){
      //at t=0 we have a B0
      //report(Severity::Info,"EvtGen") << "B0B"<<endl;
      amp=BB*_A_f+barBB*_Abar_f;
      //std::cout << "noflip B0B tag:"<<amp<<std::endl;
      //amp=0.0;
    }
    if (other_b==B0||other_b==B0s){
      //report(Severity::Info,"EvtGen") << "B0"<<endl;
      amp=BbarB*_A_f+barBbarB*_Abar_f;
    }
  }else{
    if (other_b==B0||other_b==B0s){
      amp=BbarB*_A_fbar+barBbarB*_Abar_fbar;
      //std::cout << "flip B0 tag:"<<amp<<std::endl;
      //amp=0.0;
    }
    if (other_b==B0B||other_b==B0Bs){
      amp=BB*_A_fbar+barBB*_Abar_fbar;
    }
  }


  EvtVector4R p4_parent=p->getP4Restframe();
  double m_parent=p4_parent.mass();

  EvtSpinType::spintype d2type=EvtPDL::getSpinType(getDaug(1));

  EvtVector4R momv;
  EvtVector4R moms;

  if (d2type==EvtSpinType::SCALAR){
    d2type=EvtPDL::getSpinType(getDaug(0));
    d= p->getDaug(0);
    momv = d->getP4();
    moms = p->getDaug(1)->getP4();
  }
  else{
    d= p->getDaug(1);
    momv = d->getP4();
    moms = p->getDaug(0)->getP4();
  }



  if (d2type==EvtSpinType::SCALAR) {
    vertex(amp);
  }
  
  if (d2type==EvtSpinType::VECTOR) {
    
    double norm=momv.mass()/(momv.d3mag()*p->mass());

    //std::cout << amp << " " << norm << " " << p4_parent << d->getP4()<< std::endl;
    //    std::cout << EvtPDL::name(d->getId()) << " " << EvtPDL::name(p->getDaug(0)->getId()) << 
    //  " 1and2 " << EvtPDL::name(p->getDaug(1)->getId()) << std::endl;
    //std::cout << d->eps(0) << std::endl;
    //std::cout << d->epsParent(0) << std::endl;
    vertex(0,amp*norm*p4_parent*(d->epsParent(0)));
    vertex(1,amp*norm*p4_parent*(d->epsParent(1)));
    vertex(2,amp*norm*p4_parent*(d->epsParent(2)));
  
  }

  if (d2type==EvtSpinType::TENSOR) {

    double norm=d->mass()*d->mass()/(m_parent*d->getP4().d3mag()*d->getP4().d3mag());
 
   
   vertex(0,amp*norm*d->epsTensorParent(0).cont1(p4_parent)*p4_parent);
   vertex(1,amp*norm*d->epsTensorParent(1).cont1(p4_parent)*p4_parent);
   vertex(2,amp*norm*d->epsTensorParent(2).cont1(p4_parent)*p4_parent);
   vertex(3,amp*norm*d->epsTensorParent(3).cont1(p4_parent)*p4_parent);
   vertex(4,amp*norm*d->epsTensorParent(4).cont1(p4_parent)*p4_parent);
  
  }


  return ;
}

std::string EvtSSDCP::getParamName(int i) {
  switch(i) {
  case 0:
    return "deltaM";
  case 1:
    return "deltaGammaOverGamma";
  case 2:
    return "qOverP";
  case 3:
    return "qOverPPhase";
  case 4:
    return "Af";
  case 5:
    return "AfPhase";
  case 6:
    return "Abarf";
  case 7:
    return "AbarfPhase";
  case 8:
    return "Afbar";
  case 9:
    return "AfbarPhase";
  case 10:
    return "Abarfbar";
  case 11:
    return "AbarfbarPhase";
  case 12:
    return "Z";
  case 13:
    return "ZPhase";
  default:
    return "";
  }
}

std::string EvtSSDCP::getParamDefault(int i) {
  switch(i) {
  case 3:
    return "0.0";
  case 4:
    return "1.0";
  case 5:
    return "0.0";
  case 6:
    return "1.0";
  case 7:
    return "0.0";
  default:
    return "";
  }
}
