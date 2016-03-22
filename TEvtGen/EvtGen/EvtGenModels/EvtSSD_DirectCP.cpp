// $Id: EvtSSD_DirectCP.cpp,v 1.2 2009-03-16 16:24:05 robbep Exp $
// Generation of direct CP violation in hadronic environment
// Patrick Robbe, LHCb,  08 Nov 2006
// 
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenModels/EvtSSD_DirectCP.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"

EvtSSD_DirectCP::~EvtSSD_DirectCP() {}

std::string EvtSSD_DirectCP::getName( ){

  return "SSD_DirectCP" ;

}


EvtDecayBase* EvtSSD_DirectCP::clone(){

  return new EvtSSD_DirectCP;

}

void EvtSSD_DirectCP::init(){

  // check that there is 1 argument and 2-body decay

  checkNArg(1);
  checkNDaug(2);

  EvtSpinType::spintype d1type=EvtPDL::getSpinType(getDaug(0));
  EvtSpinType::spintype d2type=EvtPDL::getSpinType(getDaug(1));
  
  if ( (!(d1type == EvtSpinType::SCALAR || d2type == EvtSpinType::SCALAR))||
       (!((d2type==EvtSpinType::SCALAR)||(d2type==EvtSpinType::VECTOR)||
          (d2type==EvtSpinType::TENSOR)))||
       (!((d1type==EvtSpinType::SCALAR)||(d1type==EvtSpinType::VECTOR)||
          (d1type==EvtSpinType::TENSOR)))
       ) {
    report(Severity::Error,"EvtGen") << "EvtSSD_DirectCP generator expected "
                           << "one of the daugters to be a scalar, "
                           << "the other either scalar, vector, or tensor, "
                           << "found:"
                           << EvtPDL::name(getDaug(0)).c_str()
                           <<" and "
                           <<EvtPDL::name(getDaug(1)).c_str()<<std::endl;
    report(Severity::Error,"EvtGen") << "Will terminate execution!"<<std::endl;
    ::abort();
  }
  
  _acp = getArg( 0 ) ; // A_CP defined as A_CP = (BR(fbar)-BR(f))/(BR(fbar)+BR(f))

}

void EvtSSD_DirectCP::initProbMax() {
  double theProbMax = 1. ;

  EvtSpinType::spintype d2type=EvtPDL::getSpinType(getDaug(1));
  EvtSpinType::spintype d1type=EvtPDL::getSpinType(getDaug(0));
  if (d1type==EvtSpinType::TENSOR||d2type==EvtSpinType::TENSOR) theProbMax*=10;

  setProbMax(theProbMax);
}

void EvtSSD_DirectCP::decay( EvtParticle *p) {

  bool flip = false ;
  EvtId daugs[2];
  
  // decide it is B or Bbar:
  if ( EvtRandom::Flat(0.,1.) < ( ( 1. - _acp ) / 2. ) ) {
    // it is a B
    if ( EvtPDL::getStdHep( getParentId() ) < 0 ) flip = true ;
  } else {
    // it is a Bbar
    if ( EvtPDL::getStdHep( getParentId() ) > 0 ) flip = true ;
  }
  
  if ( flip ) {
    if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
      p->getParent()
        ->setId( EvtPDL::chargeConj( p->getParent()->getId() ) ) ;
      p->setId( EvtPDL::chargeConj( p->getId() ) ) ;
    }
    else {
      p->setId( EvtPDL::chargeConj( p->getId() ) ) ;
    }
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
    vertex(1.);
  }
  
  if (d2type==EvtSpinType::VECTOR) {
    
    double norm=momv.mass()/(momv.d3mag()*p->mass());
    
    vertex(0,norm*p4_parent*(d->epsParent(0)));
    vertex(1,norm*p4_parent*(d->epsParent(1)));
    vertex(2,norm*p4_parent*(d->epsParent(2)));
  
  }

  if (d2type==EvtSpinType::TENSOR) {

    double norm=
      d->mass()*d->mass()/(m_parent*d->getP4().d3mag()*d->getP4().d3mag());
 
   
   vertex(0,norm*d->epsTensorParent(0).cont1(p4_parent)*p4_parent);
   vertex(1,norm*d->epsTensorParent(1).cont1(p4_parent)*p4_parent);
   vertex(2,norm*d->epsTensorParent(2).cont1(p4_parent)*p4_parent);
   vertex(3,norm*d->epsTensorParent(3).cont1(p4_parent)*p4_parent);
   vertex(4,norm*d->epsTensorParent(4).cont1(p4_parent)*p4_parent);  
  }
}

bool EvtSSD_DirectCP::isB0Mixed ( EvtParticle * p ) {
  if ( ! ( p->getParent() ) ) return false ;

  static EvtId B0 =EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  if ( ( p->getId() != B0 ) && ( p->getId() != B0B ) ) return false ;

  if ( ( p->getParent()->getId() == B0 ) ||
       ( p->getParent()->getId() == B0B ) ) return true ;

  return false ;
}

bool EvtSSD_DirectCP::isBsMixed ( EvtParticle * p ) {
  if ( ! ( p->getParent() ) ) return false ;

  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");

  if ( ( p->getId() != BS0 ) && ( p->getId() != BSB ) ) return false ;

  if ( ( p->getParent()->getId() == BS0 ) ||
       ( p->getParent()->getId() == BSB ) ) return true ;

  return false ;
}

std::string EvtSSD_DirectCP::getParamName(int i) {
  switch(i) {
  case 0:
    return "ACP";
  default:
    return "";
  }
}
