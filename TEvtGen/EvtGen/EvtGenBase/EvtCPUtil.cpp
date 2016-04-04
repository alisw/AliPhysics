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
// Module: EvtCPUtil.cc
//
// Description: Utilities needed for generation of CP violating
//              decays.
//
// Modification history:
//
//    RYD     March 24, 1998         Module created
//
//    COWAN   June 10, 2009	     Added methods for getting dGamma(s)
//				     and dm(s) using B(s)0H and B(s)0L.
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSymTable.hh"
#include "EvtGenBase/EvtConst.hh"
#include <stdio.h>
#include <stdlib.h>

#include <assert.h>
using std::endl;

EvtCPUtil::EvtCPUtil(int mixingType) {
  _enableFlip = false;
  _mixingType = mixingType;
}

EvtCPUtil::~EvtCPUtil() {
}

EvtCPUtil* EvtCPUtil::getInstance() {

  static EvtCPUtil* theCPUtil = 0;

  if (theCPUtil == 0) {
    theCPUtil = new EvtCPUtil(1);
  }

  return theCPUtil;

}

//added two functions for finding the fraction of B0 tags for decays into 
//both CP eigenstates and non-CP eigenstates -- NK, Jan. 27th, 1998

void EvtCPUtil::fractB0CP(EvtComplex Af, EvtComplex Abarf, 
			  double /*deltam*/, double beta, double &fract) {

  //This function returns the number of B0 tags for decays into CP-eigenstates
  //(the "probB0" in the new EvtOtherB)

  //double gamma_B = EvtPDL::getWidth(B0);   
  //double xd = deltam/gamma_B;
  //double xd = 0.65;
  double ratio = 1/(1 + 0.65*0.65);
  
  EvtComplex rf, rbarf;

  rf = EvtComplex(cos(2.0*beta),sin(2.0*beta))*Abarf/Af;
  rbarf = EvtComplex(1.0)/rf;

  double A2 = real(Af)*real(Af) + imag(Af)*imag(Af);
  double Abar2 = real(Abarf)*real(Abarf) + imag(Abarf)*imag(Abarf);
  
  double rf2 = real(rf)*real(rf) + imag(rf)*imag(rf);    
  double rbarf2 = real(rbarf)*real(rbarf) + imag(rbarf)*imag(rbarf);    

  fract = (Abar2*(1+ rbarf2 + (1 - rbarf2)*ratio))/(Abar2*(1+ rbarf2 + (1 - rbarf2)*ratio) + A2*(1+ rf2 + (1 - rf2)*ratio));  
  return; 

}

void EvtCPUtil::fractB0nonCP(EvtComplex Af, EvtComplex Abarf, 
			     EvtComplex Afbar, EvtComplex Abarfbar, 
			     double deltam, double beta, 
			     int flip, double &fract) {

  //this function returns the number of B0 tags for decays into non-CP eigenstates
  //(the "probB0" in the new EvtOtherB)
  //this needs more thought... 

  //double gamma_B = EvtPDL::getWidth(B0);
  //report(Severity::Info,"EvtGen") << "gamma " << gamma_B<< endl;
  //double xd = deltam/gamma_B;

  //why is the width of B0 0 in PDL??

  double xd = 0.65;
  double gamma_B = deltam/xd;
  double IAf, IAfbar, IAbarf, IAbarfbar;
  EvtComplex rf, rfbar, rbarf, rbarfbar;
  double rf2, rfbar2, rbarf2, rbarfbar2;
  double Af2, Afbar2, Abarf2, Abarfbar2;

  rf = EvtComplex(cos(2.0*beta),sin(2.0*beta))*Abarf/Af;
  rfbar = EvtComplex(cos(2.0*beta),sin(2.0*beta))*Abarfbar/Afbar; 
  rbarf = EvtComplex(cos(-2.0*beta),sin(-2.0*beta))*Af/Abarf;
  rbarfbar = EvtComplex(cos(-2.0*beta),sin(-2.0*beta))*Afbar/Abarfbar;
  
  
  rf2 = real(rf)*real(rf) + imag(rf)*imag(rf);
  rfbar2 = real(rfbar)*real(rfbar) + imag(rfbar)*imag(rfbar);
  rbarf2 = real(rbarf)*real(rbarf) + imag(rbarf)*imag(rbarf);
  rbarfbar2 = real(rbarfbar)*real(rbarfbar) + imag(rbarfbar)*imag(rbarfbar);

  Af2 = real(Af)*real(Af) + imag(Af)*imag(Af);
  Afbar2 = real(Afbar)*real(Afbar) + imag(Afbar)*imag(Afbar); 
  Abarf2 = real(Abarf)*real(Abarf) + imag(Abarf)*imag(Abarf);
  Abarfbar2 = real(Abarfbar)*real(Abarfbar) + imag(Abarfbar)*imag(Abarfbar);


  //
  //IAf = integral(gamma(B0->f)), etc.
  //

  IAf = (Af2/(2*gamma_B))*(1+rf2+(1-rf2)/(1+xd*xd));
  IAfbar = (Afbar2/(2*gamma_B))*(1+rfbar2+(1-rfbar2)/(1+xd*xd));
  IAbarf = (Abarf2/(2*gamma_B))*(1+rbarf2+(1-rbarf2)/(1+xd*xd));
  IAbarfbar = (Abarfbar2/(2*gamma_B))*(1+rbarfbar2+(1-rbarfbar2)/(1+xd*xd));
  
  //flip specifies the relative fraction of fbar events
 
  fract = IAbarf/(IAbarf+IAf) + flip*IAbarfbar/(IAfbar+IAbarfbar);


  return;  
} 

void EvtCPUtil::OtherB( EvtParticle *p,double &t, EvtId &otherb, double probB0){

  if (_mixingType == EvtCPUtil::Coherent) {

    OtherCoherentB(p, t, otherb, probB0);

  } else if (_mixingType == EvtCPUtil::Incoherent) {

    OtherIncoherentB(p, t, otherb, probB0);

  }

}

void EvtCPUtil::OtherCoherentB( EvtParticle *p,double &t, EvtId &otherb, double probB0){

  //Can not call this recursively!!!
  static int entryCount=0;
  entryCount++;

  //added by Lange Jan4,2000
  static EvtId B0B=EvtPDL::getId("anti-B0");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BS=EvtPDL::getId("B_s0");

  static EvtId UPS4S=EvtPDL::getId("Upsilon(4S)");

  int isB0=EvtRandom::Flat(0.0,1.0)<probB0;
  
  int idaug;
  
  p->setLifetime();
  
  // now get the time between the decay of this B and the other B!
  
  EvtParticle *parent=p->getParent();
  
  EvtParticle *other;

  bool incoherentmix=false;
    
  if ((parent!=0)&&(parent->getId()==B0||
		    parent->getId()==B0B||
		    parent->getId()==BS||
		    parent->getId()==BSB)) {
    incoherentmix=true;
  }

  if (incoherentmix) parent=parent->getParent();
  
  if (parent==0||parent->getId()!=UPS4S) {
    //Need to make this more general, but for now
    //assume no parent. If we have parent of B we
    //need to charge conj. full decay tree.
    
    
    if (parent!=0) {
      report(Severity::Info,"EvtGen") << "p="<<EvtPDL::name(p->getId())
			    << " parent="<<EvtPDL::name(parent->getId())
			    << endl;
    }
    assert(parent==0);
    p->setLifetime();
    t=p->getLifetime();
    bool needToChargeConj=false;
    if (p->getId()==B0B&&isB0) needToChargeConj=true;
    if (p->getId()==B0&&!isB0) needToChargeConj=true;
    if (p->getId()==BSB&&isB0) needToChargeConj=true;
    if (p->getId()==BS&&!isB0) needToChargeConj=true;

    if (needToChargeConj) {
      p->setId( EvtPDL::chargeConj(p->getId()));
      if (incoherentmix) {
	p->getDaug(0)->setId(EvtPDL::chargeConj(p->getDaug(0)->getId()));
      }
    }
    otherb=EvtPDL::chargeConj(p->getId());
    
    entryCount--;
    return;
  }
  else{
    if (parent->getDaug(0)!=p){
      other=parent->getDaug(0);
      idaug=0;
    }
    else{
      other=parent->getDaug(1);
      idaug=1;
    }
  }
  
  if (parent != 0 ) {

    //if (entryCount>1){
    //  report(Severity::Info,"EvtGen") << "Double CP decay:"<<entryCount<<endl;
    //}

    //kludge!! Lange Mar21, 2003 	 
    // if the other B is an alias... don't change the flavor.. 	 
    if ( other->getId().isAlias() ) { 	 
      OtherB(p,t,otherb);
      entryCount--;
      return; 	 
      
    }
    
    if (entryCount==1){
    
      EvtVector4R p_init=other->getP4();
      //int decayed=other->getNDaug()>0;
      bool decayed = other->isDecayed();

      other->deleteTree();
    
      EvtScalarParticle* scalar_part;
      
      scalar_part=new EvtScalarParticle;
      if (isB0) {
	scalar_part->init(B0,p_init);
      }
      else{
	scalar_part->init(B0B,p_init);
      }
      other=(EvtParticle *)scalar_part;
      //    other->set_type(EvtSpinType::SCALAR);
      other->setDiagonalSpinDensity();      
    
      parent->insertDaugPtr(idaug,other);
    
      if (decayed){
	//report(Severity::Info,"EvtGen") << "In CP Util calling decay \n";
	other->decay();
      }

    }

    otherb=other->getId();

    other->setLifetime();
    t=p->getLifetime()-other->getLifetime();
    
    otherb = other->getId();
    
  }
  else {
    report(Severity::Info,"EvtGen") << "We have an error here!!!!"<<endl;
    otherb = EvtId(-1,-1); 
  }
  
  entryCount--;
  return ;
}

// ========================================================================
bool EvtCPUtil::isBsMixed ( EvtParticle * p )
{
  if ( ! ( p->getParent() ) ) return false ;

  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");

  if ( ( p->getId() != BS0 ) && ( p->getId() != BSB ) ) return false ;

  if ( ( p->getParent()->getId() == BS0 ) ||
       ( p->getParent()->getId() == BSB ) ) return true ;

  return false ;
}

// ========================================================================
bool EvtCPUtil::isB0Mixed ( EvtParticle * p )
{
  if ( ! ( p->getParent() ) ) return false ;

  static EvtId B0 =EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  if ( ( p->getId() != B0 ) && ( p->getId() != B0B ) ) return false ;

  if ( ( p->getParent()->getId() == B0 ) ||
       ( p->getParent()->getId() == B0B ) ) return true ;

  return false ;
}
//============================================================================
// Return the tag of the event (ie the anti-flavour of the produced 
// B meson). Flip the flavour of the event with probB probability
//============================================================================
void EvtCPUtil::OtherIncoherentB( EvtParticle * p ,
                                  double & t ,
                                  EvtId & otherb ,
                                  double probB )
{
  //std::cout<<"New routine running"<<endl;
  //if(p->getId() == B0 || p->getId() == B0B) 
  //added by liming Zhang
  enableFlip();
  if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
    p->getParent()->setLifetime() ;
    t = p->getParent()->getLifetime() ;
  }
  else {
    p->setLifetime() ;
    t = p->getLifetime() ;
  }

  if ( flipIsEnabled() ) {
    //std::cout << " liming << flipIsEnabled " << std::endl;
    // Flip the flavour of the particle with probability probB
    bool isFlipped = ( EvtRandom::Flat( 0. , 1. ) < probB ) ;

    if ( isFlipped ) {
      if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
        p->getParent()
          ->setId( EvtPDL::chargeConj( p->getParent()->getId() ) ) ;
        p->setId( EvtPDL::chargeConj( p->getId() ) ) ;
      }
      else {
        p->setId( EvtPDL::chargeConj( p->getId() ) ) ;
      }
    }
  }

  if ( ( isB0Mixed( p ) ) || ( isBsMixed( p ) ) ) {
    // if B has mixed, tag flavour is charge conjugate of parent of B-meson
    otherb = EvtPDL::chargeConj( p->getParent()->getId() ) ;
  }
  else {
    // else it is opposite flavour than this B hadron
    otherb = EvtPDL::chargeConj( p->getId() ) ;
  }

  return ;
}
//============================================================================
void EvtCPUtil::OtherB( EvtParticle *p,double &t, EvtId &otherb){

  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId B0B=EvtPDL::getId("anti-B0");
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId UPS4=EvtPDL::getId("Upsilon(4S)");

  if (p->getId()==BS0||p->getId()==BSB){
    static double ctauL=EvtPDL::getctau(EvtPDL::getId("B_s0L"));
    static double ctauH=EvtPDL::getctau(EvtPDL::getId("B_s0H"));
    static double ctau=ctauL<ctauH?ctauH:ctauL;
    t=-log(EvtRandom::Flat())*ctau;
    EvtParticle* parent=p->getParent();
    if (parent!=0&&(parent->getId()==BS0||parent->getId()==BSB)){
      if (parent->getId()==BS0) otherb=BSB;
      if (parent->getId()==BSB) otherb=BS0;
      parent->setLifetime(t);
      return;
    }
    if (p->getId()==BS0) otherb=BSB;
    if (p->getId()==BSB) otherb=BS0;
    p->setLifetime(t);
    return;
  }

  if (p->getId()==D0||p->getId()==D0B){
    static double ctauL=EvtPDL::getctau(EvtPDL::getId("D0L"));
    static double ctauH=EvtPDL::getctau(EvtPDL::getId("D0H"));
    static double ctau=ctauL<ctauH?ctauH:ctauL;
    t=-log(EvtRandom::Flat())*ctau;
    EvtParticle* parent=p->getParent();
    if (parent!=0&&(parent->getId()==D0||parent->getId()==D0B)){
      if (parent->getId()==D0) otherb=D0B;
      if (parent->getId()==D0B) otherb=D0;
      parent->setLifetime(t);
      return;
    }
    if (p->getId()==D0) otherb=D0B;
    if (p->getId()==D0B) otherb=D0;
    p->setLifetime(t);
    return;
  }

  p->setLifetime();

  // now get the time between the decay of this B and the other B!
  
  EvtParticle *parent=p->getParent();

  if (parent==0||parent->getId()!=UPS4) {
    //report(Severity::Error,"EvtGen") << 
    //  "Warning CP violation with B having no parent!"<<endl;
    t=p->getLifetime();
    if (p->getId()==B0) otherb=B0B;
    if (p->getId()==B0B) otherb=B0;
    if (p->getId()==BS0) otherb=BSB;
    if (p->getId()==BSB) otherb=BS0;
    return;
  }
  else{
    if (parent->getDaug(0)!=p){
      otherb=parent->getDaug(0)->getId();
      parent->getDaug(0)->setLifetime();
      t=p->getLifetime()-parent->getDaug(0)->getLifetime();
    }
    else{
     otherb=parent->getDaug(1)->getId();
      parent->getDaug(1)->setLifetime();
      t=p->getLifetime()-parent->getDaug(1)->getLifetime();
   }
  }


  return ;
}

// No CP violation is assumed
void EvtCPUtil::incoherentMix(const EvtId id, double &t, int &mix){

  int stdHepNum=EvtPDL::getStdHep(id);
  stdHepNum=abs(stdHepNum);
  
  EvtId partId=EvtPDL::evtIdFromStdHep(stdHepNum);

  std::string partName=EvtPDL::name(partId);
  std::string hname=partName+std::string("H");
  std::string lname=partName+std::string("L");

  EvtId lId=EvtPDL::getId(lname);
  EvtId hId=EvtPDL::getId(hname);

  double ctauL=EvtPDL::getctau(lId);
  double ctauH=EvtPDL::getctau(hId);

  // Bug Fixed: Corrected the average as gamma is the relevent parameter
  double ctau=2.0*(ctauL*ctauH)/(ctauL+ctauH);
  //double ctau=0.5*(ctauL+ctauH);

  // Bug Fixed: ctau definition changed above
  //double y=(ctauH-ctauL)/(2*ctau);
  double y=(ctauH-ctauL)/(ctauH+ctauL);

  //deltam and qoverp defined in DECAY.DEC

  std::string qoverpParmName=std::string("qoverp_incohMix_")+partName;
  std::string mdParmName=std::string("dm_incohMix_")+partName;
  int ierr;
  double qoverp=atof(EvtSymTable::get(qoverpParmName,ierr).c_str());
  double x=atof(EvtSymTable::get(mdParmName,ierr).c_str())*ctau/EvtConst::c;
  double fac;

  if(id==partId){
    fac=1.0/(qoverp*qoverp);
  }
  else{
    fac=qoverp*qoverp;
  }

  double mixprob=(x*x+y*y)/(x*x+y*y+fac*(2.0+x*x-y*y));

  int mixsign;

  mixsign=(mixprob>EvtRandom::Flat(0.0,1.0))?-1:1;

  double prob;

  // Find the longest of the two lifetimes
  double ctaulong = ctauL<=ctauH?ctauH:ctauL;

  // Bug fixed: Ensure cosine argument is dimensionless so /ctau
  do{
    t=-log(EvtRandom::Flat())*ctaulong;
    prob=1.0+exp(-2.0*fabs(y)*t/ctau)+mixsign*2.0*exp(-fabs(y)*t/ctau)*cos(x*t/ctau);
  }while(prob<4.0*EvtRandom::Flat());

  mix=0;

  if (mixsign==-1) mix=1;
   
  return;
}


double EvtCPUtil::getDeltaGamma(const EvtId id){

  int stdHepNum = EvtPDL::getStdHep(id);
  stdHepNum = abs(stdHepNum);
  EvtId partId = EvtPDL::evtIdFromStdHep(stdHepNum);

  std::string partName = EvtPDL::name(partId);
  std::string hname = partName + std::string("H");
  std::string lname = partName + std::string("L");
  
  EvtId lId = EvtPDL::getId(lname);
  EvtId hId = EvtPDL::getId(hname);

  double ctauL = EvtPDL::getctau(lId);  
  double ctauH = EvtPDL::getctau(hId);
  
  double dGamma = (1/ctauL - 1/ctauH)*EvtConst::c;
  return dGamma;
}

double EvtCPUtil::getDeltaM(const EvtId id){

  int stdHepNum = EvtPDL::getStdHep(id);  
  stdHepNum = abs(stdHepNum); 
  EvtId partId = EvtPDL::evtIdFromStdHep(stdHepNum);
  
  std::string partName = EvtPDL::name(partId);
  std::string parmName = std::string("dm_incohMix_") + partName;

  int ierr;  
  double dM = atof(EvtSymTable::get(parmName,ierr).c_str());
  return dM;
}

bool EvtCPUtil::flipIsEnabled() { return _enableFlip ; }
void EvtCPUtil::enableFlip() { _enableFlip = true ; }
void EvtCPUtil::disableFlip() { _enableFlip = false ; }

