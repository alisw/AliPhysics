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
// Module: EvtParticle.cc
//
// Description: Class to describe all particles
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtNeutrinoParticle.hh"
#include "EvtGenBase/EvtRaritaSchwingerParticle.hh"
#include "EvtGenBase/EvtStringParticle.hh"
#include "EvtGenBase/EvtStdHep.hh"
#include "EvtGenBase/EvtSecondary.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtIdSet.hh"
#include "EvtGenBase/EvtStatus.hh"

using std::endl;



EvtParticle::~EvtParticle() {
  delete _decayProb;
}

EvtParticle::EvtParticle() {
   _ndaug=0;
   _parent=0;
   _channel=-10;
   _t=0.0;
   _genlifetime=1;
   _first=1;
   _isInit=false;
   _validP4=false;
   _isDecayed=false;
   _decayProb=0;
   //   _mix=false;
}

void EvtParticle::setFirstOrNot() {
  _first=0;
}
void EvtParticle::resetFirstOrNot() {
  _first=1;
}

void EvtParticle::setChannel( int i ) { 
  _channel=i;
}

EvtParticle *EvtParticle::getDaug(int i) { return _daug[i]; }

EvtParticle *EvtParticle::getParent() const { return _parent;}

void EvtParticle::setLifetime(double tau){
  _t=tau;
}

void EvtParticle::setLifetime(){
  if (_genlifetime){
    _t=-log(EvtRandom::Flat())*EvtPDL::getctau(getId());
  }
}

double EvtParticle::getLifetime(){

  return _t;
}

void EvtParticle::addDaug(EvtParticle *node) {
  node->_daug[node->_ndaug++]=this;
  _ndaug=0;
  _parent=node; 
}


int EvtParticle::firstornot() const { return _first;}

EvtId EvtParticle::getId() const { return _id;}

int EvtParticle::getPDGId() const {return EvtPDL::getStdHep(_id);}

EvtSpinType::spintype EvtParticle::getSpinType() const 
      { return EvtPDL::getSpinType(_id);}

int EvtParticle::getSpinStates() const 
  { return EvtSpinType::getSpinStates(EvtPDL::getSpinType(_id));}

const EvtVector4R& EvtParticle::getP4() const { return _p;}

int EvtParticle::getChannel() const { return _channel;}

size_t EvtParticle::getNDaug() const { return _ndaug;}

double EvtParticle::mass() const {

     return _p.mass();
}


void EvtParticle::setDiagonalSpinDensity(){

  _rhoForward.setDiag(getSpinStates());
}

void EvtParticle::setVectorSpinDensity(){

  if (getSpinStates()!=3) {
    report(Severity::Error,"EvtGen")<<"Error in EvtParticle::setVectorSpinDensity"<<endl;
    report(Severity::Error,"EvtGen")<<"spin_states:"<<getSpinStates()<<endl;
    report(Severity::Error,"EvtGen")<<"particle:"<<EvtPDL::name(_id).c_str()<<endl;
    ::abort();
  }

  EvtSpinDensity rho;

  //Set helicity +1 and -1 to 1.
  rho.setDiag(getSpinStates());
  rho.set(1,1,EvtComplex(0.0,0.0));

  setSpinDensityForwardHelicityBasis(rho);

}


void EvtParticle::setSpinDensityForwardHelicityBasis(const EvtSpinDensity& rho){

  EvtSpinDensity R=rotateToHelicityBasis();

  assert(R.getDim()==rho.getDim());

  int n=rho.getDim();

  _rhoForward.setDim(n);

  int i,j,k,l;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      EvtComplex tmp=0.0;
      for(k=0;k<n;k++){
	for(l=0;l<n;l++){
	  tmp+=R.get(l,i)*rho.get(l,k)*conj(R.get(k,j));
	}
      }
      _rhoForward.set(i,j,tmp);
    }
  }

}

void EvtParticle::setSpinDensityForwardHelicityBasis(const EvtSpinDensity& rho,
						     double alpha,
						     double beta,
						     double gamma){

  EvtSpinDensity R=rotateToHelicityBasis(alpha,beta,gamma);

  assert(R.getDim()==rho.getDim());

  int n=rho.getDim();

  _rhoForward.setDim(n);

  int i,j,k,l;

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      EvtComplex tmp=0.0;
      for(k=0;k<n;k++){
	for(l=0;l<n;l++){
	  tmp+=R.get(l,i)*rho.get(l,k)*conj(R.get(k,j));
	}
      }
      _rhoForward.set(i,j,tmp);
    }
  }

}

void EvtParticle::initDecay(bool useMinMass) {

  EvtParticle* p=this;
  // carefull - the parent mass might be fixed in stone..
  EvtParticle* par=p->getParent();
  double parMass=-1.;
  if ( par != 0 ) {
    if ( par->hasValidP4() ) parMass=par->mass();
    for (size_t i=0;i<par->getNDaug();i++) {
      EvtParticle *tDaug=par->getDaug(i);
      if ( p != tDaug )
        parMass-=EvtPDL::getMinMass(tDaug->getId());
    }
  }
  
  if ( _isInit ) {
    //we have already been here - just reroll the masses!
    if ( _ndaug>0) {
      for(size_t ii=0;ii<_ndaug;ii++){
        if ( _ndaug==1 ||  EvtPDL::getWidth(p->getDaug(ii)->getId()) > 0.0000001)
          p->getDaug(ii)->initDecay(useMinMass);
        else p->getDaug(ii)->setMass(EvtPDL::getMeanMass(p->getDaug(ii)->getId()));
      }
    }
    
    EvtId *dauId=0;
    double *dauMasses=0;
    if ( _ndaug > 0) {
      dauId=new EvtId[_ndaug];
      dauMasses=new double[_ndaug];
      for (size_t j=0;j<_ndaug;j++) { 
        dauId[j]=p->getDaug(j)->getId();
        dauMasses[j]=p->getDaug(j)->mass();
      }
    }
    EvtId *parId=0;
    EvtId *othDauId=0;
    EvtParticle *tempPar=p->getParent();
    if (tempPar) {
      parId=new EvtId(tempPar->getId());
      if ( tempPar->getNDaug()==2 ) {
        if ( tempPar->getDaug(0) == this ) othDauId=new EvtId(tempPar->getDaug(1)->getId());
        else othDauId=new EvtId(tempPar->getDaug(0)->getId());
      }
    }
    if ( p->getParent() && _validP4==false ) {
      if ( !useMinMass ) {
        p->setMass(EvtPDL::getRandMass(p->getId(),parId,_ndaug,dauId,othDauId,parMass,dauMasses));
      }
      else p->setMass(EvtPDL::getMinMass(p->getId()));
    }
    if ( parId) delete parId;
    if ( othDauId) delete othDauId;
    if ( dauId) delete [] dauId;
    if ( dauMasses) delete [] dauMasses;
    return;
  }
  
  
  //Will include effects of mixing here
  //added by Lange Jan4,2000
  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BD0=EvtPDL::getId("B0");
  static EvtId BDB=EvtPDL::getId("anti-B0");
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId U4S=EvtPDL::getId("Upsilon(4S)");
  static EvtIdSet borUps(BS0,BSB,BD0,BDB,U4S);
  
  //only makes sense if there is no parent particle which is a B or an Upsilon
  bool hasBorUps=false;
  if ( getParent() && borUps.contains(getParent()->getId()) ) hasBorUps=true;
  //    if ( (getNDaug()==0)&&(getParent()==0) && (getId()==BS0||getId()==BSB||getId()==BD0||getId()==BDB)){
  EvtId thisId=getId();
  // remove D0 mixing for now.
  //  if ( (getNDaug()==0 && !hasBorUps) && (thisId==BS0||thisId==BSB||thisId==BD0||thisId==BDB||thisId==D0||thisId==D0B)){
  if ( (getNDaug()==0 && !hasBorUps) && (thisId==BS0||thisId==BSB||thisId==BD0||thisId==BDB)){
    double t;
    int mix;
    EvtCPUtil::getInstance()->incoherentMix(getId(), t, mix);
    setLifetime(t);
    
    if (mix) {

      EvtScalarParticle* scalar_part;
    
      scalar_part=new EvtScalarParticle;
      if (getId()==BS0) {
        EvtVector4R p_init(EvtPDL::getMass(BSB),0.0,0.0,0.0);
        scalar_part->init(EvtPDL::chargeConj(getId()),p_init);
      }
      else if (getId()==BSB) {
        EvtVector4R p_init(EvtPDL::getMass(BS0),0.0,0.0,0.0);
        scalar_part->init(EvtPDL::chargeConj(getId()),p_init);
      }
      else if (getId()==BD0) {
        EvtVector4R p_init(EvtPDL::getMass(BDB),0.0,0.0,0.0);
        scalar_part->init(EvtPDL::chargeConj(getId()),p_init);
      }
      else if (getId()==BDB) {
        EvtVector4R p_init(EvtPDL::getMass(BD0),0.0,0.0,0.0);
        scalar_part->init(EvtPDL::chargeConj(getId()),p_init);
      }
      else if (getId()==D0) {
        EvtVector4R p_init(EvtPDL::getMass(D0B),0.0,0.0,0.0);
        scalar_part->init(EvtPDL::chargeConj(getId()),p_init);
      }
      else if (getId()==D0B) {
        EvtVector4R p_init(EvtPDL::getMass(D0),0.0,0.0,0.0);
        scalar_part->init(EvtPDL::chargeConj(getId()),p_init);
      }
      
      scalar_part->setLifetime(0);
      scalar_part->setDiagonalSpinDensity();      
      
      insertDaugPtr(0,scalar_part);

      _ndaug=1;
      _isInit=true;
      p=scalar_part;
      p->initDecay(useMinMass);
      return;


    }
  }

  EvtDecayBase *decayer;
  decayer = EvtDecayTable::getInstance()->getDecayFunc(p);

  if ( decayer ) {
    p->makeDaughters(decayer->nRealDaughters(),decayer->getDaugs());
    //then loop over the daughters and init their decay
    for(size_t i=0;i<p->getNDaug();i++){
      //      std::cout << EvtPDL::name(p->getDaug(i)->getId()) << " " << i << " " << p->getDaug(i)->getSpinType() << " " << EvtPDL::name(p->getId()) << std::endl;
      if ( EvtPDL::getWidth(p->getDaug(i)->getId()) > 0.0000001)
        p->getDaug(i)->initDecay(useMinMass);
      else p->getDaug(i)->setMass(EvtPDL::getMeanMass(p->getDaug(i)->getId()));
    }
  }
  
  int j;
  EvtId *dauId=0;
  double *dauMasses=0;
  int nDaugT=p->getNDaug();
  if ( nDaugT > 0) {
    dauId=new EvtId[nDaugT];
    dauMasses=new double[nDaugT];
    for (j=0;j<nDaugT;j++) { 
      dauId[j]=p->getDaug(j)->getId();
      dauMasses[j]=p->getDaug(j)->mass();
    }
  }

  EvtId *parId=0;
  EvtId *othDauId=0;
  EvtParticle *tempPar=p->getParent();
  if (tempPar) {
    parId=new EvtId(tempPar->getId());
    if ( tempPar->getNDaug()==2 ) {
      if ( tempPar->getDaug(0) == this ) othDauId=new EvtId(tempPar->getDaug(1)->getId());
      else othDauId=new EvtId(tempPar->getDaug(0)->getId());
    }
  }
  if ( p->getParent() && p->hasValidP4()==false ) {
    if ( !useMinMass ) {
      p->setMass(EvtPDL::getRandMass(p->getId(),parId,p->getNDaug(),dauId,othDauId,parMass,dauMasses));
    }
    else {
      p->setMass(EvtPDL::getMinMass(p->getId()));
    }
  }
  if ( parId) delete parId;
  if ( othDauId) delete othDauId;
  if ( dauId) delete [] dauId;
  if ( dauMasses) delete [] dauMasses;
  _isInit=true;
}


void EvtParticle::decay(){
  //P is particle to decay, typically 'this' but sometime
  //modified by mixing 
  EvtParticle* p=this;
  //Did it mix?
  //if ( p->getMixed() ) {
    //should take C(p) - this should only
    //happen the first time we call decay for this
    //particle
  //p->takeCConj();
  // p->setUnMixed();
  //}

  EvtDecayBase *decayer;
  decayer = EvtDecayTable::getInstance()->getDecayFunc(p);
  //  if ( decayer ) {
  //    report(Severity::Info,"EvtGen") << "calling decay for " << EvtPDL::name(p->getId()) << " " << p->mass() << " " << p->getP4() << " " << p->getNDaug() << " " << p << endl;
  //    report(Severity::Info,"EvtGen") << "NDaug= " << decayer->getNDaug() << endl;
  //    int ti;
  //    for ( ti=0; ti<decayer->getNDaug(); ti++) 
  //      report(Severity::Info,"EvtGen") << "Daug " << ti << " " << EvtPDL::name(decayer->getDaug(ti)) << endl;
  //  }
  //if (p->_ndaug>0) {
  //      report(Severity::Info,"EvtGen") <<"Is decaying particle with daughters!!!!!"<<endl;
  //     ::abort();
    //return;
    //call initdecay first - April 29,2002 - Lange
  //}

  //if there are already daughters, then this step is already done!
  // figure out the masses
  bool massTreeOK(true);
  if ( _ndaug == 0 ) {
    massTreeOK = generateMassTree();
  }

  if (massTreeOK == false) {
    report(Severity::Info,"EvtGen")<<"Could not decay "<<EvtPDL::name(p->getId())
			 <<" with mass "<<p->mass()
			 <<" to decay channel number "<<_channel<<endl;
    _isDecayed = false;
    return;
  }
  
  static EvtId BS0=EvtPDL::getId("B_s0");
  static EvtId BSB=EvtPDL::getId("anti-B_s0");
  static EvtId BD0=EvtPDL::getId("B0");
  static EvtId BDB=EvtPDL::getId("anti-B0"); 
  // static EvtId D0=EvtPDL::getId("D0");
  // static EvtId D0B=EvtPDL::getId("anti-D0");

  EvtId thisId=getId();
  // remove D0 mixing for now..
  //  if ( _ndaug==1 &&  (thisId==BS0||thisId==BSB||thisId==BD0||thisId==BDB||thisId==D0||thisId==D0B) ) {
  if ( _ndaug==1 &&  (thisId==BS0||thisId==BSB||thisId==BD0||thisId==BDB) ) {
    p=p->getDaug(0);
    decayer = EvtDecayTable::getInstance()->getDecayFunc(p);
  }
  //now we have accepted a set of masses - time
  if ( decayer != 0) {
    decayer->makeDecay(p);
  }
  else{
    p->_rhoBackward.setDiag(p->getSpinStates());
  }

  _isDecayed=true;
  return;  
}

bool EvtParticle::generateMassTree() {

  bool isOK(true);

  double massProb=1.;
  double ranNum=2.;
  int counter=0;
  EvtParticle *p=this;
  while (massProb<ranNum) {
    //check it out the first time.
    p->initDecay();
    massProb=p->compMassProb();
    ranNum=EvtRandom::Flat();
    counter++;

    if ( counter > 10000 ) {
      if ( counter == 10001 ) {
	report(Severity::Info,"EvtGen") << "Too many iterations to determine the mass tree. Parent mass= "<< p->mass() << " " << massProb <<endl;
	p->printTree();
	report(Severity::Info,"EvtGen") << "will take next combo with non-zero likelihood\n"; 
      }
      if ( massProb>0. ) massProb=2.0;
      if ( counter > 20000 ) {
	// one last try - take the minimum masses
	p->initDecay(true);
	p->printTree();
	massProb=p->compMassProb();
	if ( massProb>0. ) {
	  massProb=2.0;
	  report(Severity::Info,"EvtGen") << "Taking the minimum mass of all particles in the chain\n";
	}
	else {
	  report(Severity::Info,"EvtGen") << "Sorry, no luck finding a valid set of masses.  This may be a pathological combo\n";
	  isOK = false;
	  break;
	}
      }
    }
  }

  return isOK;

}

double EvtParticle::compMassProb() {

  EvtParticle *p=this;
  double mass=p->mass();
  double parMass=0.;
  if ( p->getParent()) { 
    parMass=p->getParent()->mass();
  }

  int nDaug=p->getNDaug();
  double *dMasses=0;

  int i;
  if ( nDaug>0 ) {
    dMasses=new double[nDaug];
    for (i=0; i<nDaug; i++) dMasses[i]=p->getDaug(i)->mass();
  }

  double temp=1.0;
  temp=EvtPDL::getMassProb(p->getId(), mass, parMass, nDaug, dMasses);

  //If the particle already has a mass, we dont need to include
  //it in the probability calculation
  if ( (!p->getParent() || _validP4 ) && temp>0.0 ) temp=1.; 

  delete [] dMasses;
  for (i=0; i<nDaug; i++) {
    temp*=p->getDaug(i)->compMassProb();
  }
  return temp;
}

void EvtParticle::deleteDaughters(bool keepChannel){

  for(size_t i=0;i<_ndaug;i++){
    _daug[i]->deleteTree();
  }
  
  _ndaug=0;
  if ( !keepChannel) _channel=-10;
  _first=1;
  _isInit=false;
}

void EvtParticle::deleteTree(){

  this->deleteDaughters();
  
  delete this;
  
}

EvtVector4C EvtParticle::epsParent(int i) const {
  EvtVector4C temp;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th polarization vector."
			 <<" I.e. you thought it was a"
			 <<" vector particle!" << endl;
  ::abort();
  return temp;
}

EvtVector4C EvtParticle::eps(int i) const {
  EvtVector4C temp;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th polarization vector."
			 <<" I.e. you thought it was a"
			 <<" vector particle!" << endl;
  ::abort();
  return temp;
}

EvtVector4C EvtParticle::epsParentPhoton(int i){
  EvtVector4C temp;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th polarization vector of photon."
			 <<" I.e. you thought it was a"
			 <<" photon particle!" << endl;
  ::abort();
  return temp;
}

EvtVector4C EvtParticle::epsPhoton(int i){
  EvtVector4C temp;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th polarization vector of a photon."
			 <<" I.e. you thought it was a"
			 <<" photon particle!" << endl;
  ::abort();
  return temp;
}

EvtDiracSpinor EvtParticle::spParent(int i) const {
  EvtDiracSpinor tempD;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th dirac spinor."
			 <<" I.e. you thought it was a"
			 <<" Dirac particle!" << endl;
  ::abort();
  return tempD;
}

EvtDiracSpinor EvtParticle::sp(int i) const {
  EvtDiracSpinor tempD;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th dirac spinor."
			 <<" I.e. you thought it was a"
			 <<" Dirac particle!" << endl;
  ::abort();
  return tempD;
}

EvtDiracSpinor EvtParticle::spParentNeutrino() const {
  EvtDiracSpinor tempD;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the "
			 <<"dirac spinor."
			 <<" I.e. you thought it was a"
			 <<" neutrino particle!" << endl;
  ::abort();
  return tempD;
}

EvtDiracSpinor EvtParticle::spNeutrino() const {
  EvtDiracSpinor tempD;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the "
			 <<"dirac spinor."
			 <<" I.e. you thought it was a"
			 <<" neutrino particle!" << endl;
  ::abort();
  return tempD;
}

EvtTensor4C EvtParticle::epsTensorParent(int i) const {
  EvtTensor4C tempC; 
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th tensor."
			 <<" I.e. you thought it was a"
			 <<" Tensor particle!" << endl;
  ::abort();
  return tempC;
}

EvtTensor4C EvtParticle::epsTensor(int i) const {
  EvtTensor4C tempC; 
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th tensor."
			 <<" I.e. you thought it was a"
			 <<" Tensor particle!" << endl;
  ::abort();
  return tempC;
}


EvtRaritaSchwinger EvtParticle::spRSParent(int i) const {
  EvtRaritaSchwinger tempD;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th Rarita-Schwinger spinor."
			 <<" I.e. you thought it was a"
			 <<" RaritaSchwinger particle!" << std::endl;
  ::abort();
  return tempD;
}

EvtRaritaSchwinger EvtParticle::spRS(int i) const {
  EvtRaritaSchwinger tempD;
  printParticle();
  report(Severity::Error,"EvtGen") << "and you have asked for the:"<<i
			 <<"th Rarita-Schwinger spinor."
			 <<" I.e. you thought it was a"
			 <<" RaritaSchwinger particle!" << std::endl;
  ::abort();
  return tempD;
}



EvtVector4R EvtParticle::getP4Lab() const {
  EvtVector4R temp,mom;
  const EvtParticle *ptemp;
  
  temp=this->getP4();
  ptemp=this;
  
  while (ptemp->getParent()!=0) {
    ptemp=ptemp->getParent();
    mom=ptemp->getP4();
    temp=boostTo(temp,mom);   
  } 
  return temp;
}

EvtVector4R EvtParticle::getP4LabBeforeFSR() {
  EvtVector4R temp,mom;
  EvtParticle *ptemp;

  temp=this->_pBeforeFSR;
  ptemp=this;

  while (ptemp->getParent()!=0) {
    ptemp=ptemp->getParent();
    mom=ptemp->getP4();
    temp=boostTo(temp,mom);
  }
  return temp;
}



EvtVector4R EvtParticle::getP4Restframe() const {

  return EvtVector4R(mass(),0.0,0.0,0.0);

}

EvtVector4R EvtParticle::get4Pos() const {

  EvtVector4R temp,mom;
  EvtParticle *ptemp;
  
  temp.set(0.0,0.0,0.0,0.0);
  ptemp=getParent();

  if (ptemp==0) return temp;

  temp=(ptemp->_t/ptemp->mass())*(ptemp->getP4());

  while (ptemp->getParent()!=0) {
    ptemp=ptemp->getParent();
    mom=ptemp->getP4();
    temp=boostTo(temp,mom);
    temp=temp+(ptemp->_t/ptemp->mass())*(ptemp->getP4());
  } 
  
  return temp;
}


EvtParticle * EvtParticle::nextIter(EvtParticle *rootOfTree) {

  EvtParticle *bpart;
  EvtParticle *current;

  current=this;
  size_t i;

  if (_ndaug!=0) return _daug[0];

  do{
    bpart=current->_parent;
    if (bpart==0) return 0;
    i=0;
    while (bpart->_daug[i]!=current) {i++;}

    if ( bpart==rootOfTree ) {
      if ( i+1 == bpart->_ndaug ) return 0;
    }

    i++;
    current=bpart;

  }while(i>=bpart->_ndaug);

  return bpart->_daug[i];

}


void EvtParticle::makeStdHep(EvtStdHep& stdhep,EvtSecondary& secondary,
			     EvtId *list_of_stable){

  //first add particle to the stdhep list;
  stdhep.createParticle(getP4Lab(),get4Pos(),-1,-1,
			EvtPDL::getStdHep(getId()));

  int ii=0;

  //lets see if this is a longlived particle and terminate the 
  //list building!
  
  while (list_of_stable[ii]!=EvtId(-1,-1)) {
    if (getId()==list_of_stable[ii]){
      secondary.createSecondary(0,this);
      return;
    }
    ii++;
  }



  for(size_t i=0;i<_ndaug;i++){
    stdhep.createParticle(_daug[i]->getP4Lab(),_daug[i]->get4Pos(),0,0,
			  EvtPDL::getStdHep(_daug[i]->getId()));
  }

  for(size_t i=0;i<_ndaug;i++){
    _daug[i]->makeStdHepRec(1+i,1+i,stdhep,secondary,list_of_stable);
  }
  return;

}

void EvtParticle::makeStdHep(EvtStdHep& stdhep){

  //first add particle to the stdhep list;
  stdhep.createParticle(getP4Lab(),get4Pos(),-1,-1,
			EvtPDL::getStdHep(getId()));

  for(size_t i=0;i<_ndaug;i++){
    stdhep.createParticle(_daug[i]->getP4Lab(),_daug[i]->get4Pos(),0,0,
			  EvtPDL::getStdHep(_daug[i]->getId()));
  }

  for(size_t i=0;i<_ndaug;i++){
    _daug[i]->makeStdHepRec(1+i,1+i,stdhep);
  }
  return;

}


void EvtParticle::makeStdHepRec(int firstparent,int lastparent,
				EvtStdHep& stdhep,
				EvtSecondary& secondary,
				EvtId *list_of_stable){


  int ii=0;

  //lets see if this is a longlived particle and terminate the 
  //list building!
  
  while (list_of_stable[ii]!=EvtId(-1,-1)) {
    if (getId()==list_of_stable[ii]){
      secondary.createSecondary(firstparent,this);
      return;
    }
    ii++;
  }



  int parent_num=stdhep.getNPart();
  for(size_t i=0;i<_ndaug;i++){
    stdhep.createParticle(_daug[i]->getP4Lab(),_daug[i]->get4Pos(),
			  firstparent,lastparent,
			  EvtPDL::getStdHep(_daug[i]->getId()));
  }

  for(size_t i=0;i<_ndaug;i++){
    _daug[i]->makeStdHepRec(parent_num+i,parent_num+i,stdhep,
			   secondary,list_of_stable);
  }
  return;

}

void EvtParticle::makeStdHepRec(int firstparent,int lastparent,
				EvtStdHep& stdhep){

  int parent_num=stdhep.getNPart();
  for(size_t i=0;i<_ndaug;i++){
    stdhep.createParticle(_daug[i]->getP4Lab(),_daug[i]->get4Pos(),
			  firstparent,lastparent,
			  EvtPDL::getStdHep(_daug[i]->getId()));
  }

  for(size_t i=0;i<_ndaug;i++){
    _daug[i]->makeStdHepRec(parent_num+i,parent_num+i,stdhep);
  }
  return;

}

void EvtParticle::printTreeRec(unsigned int level) const {

  size_t newlevel,i;
  newlevel = level +1;

  
  if (_ndaug!=0) {
    if ( level > 0 ) {
      for (i=0;i<(5*level);i++) {
	report(Severity::Info,"") <<" ";
      }
    }
    report(Severity::Info,"") << EvtPDL::name(_id).c_str();  
    report(Severity::Info,"") << " -> ";
    for(i=0;i<_ndaug;i++){
      report(Severity::Info,"") << EvtPDL::name(_daug[i]->getId()).c_str()<<" ";
    }
    for(i=0;i<_ndaug;i++){
      report(Severity::Info,"") << _daug[i]->mass()<< " " << _daug[i]->getP4() << " " <<_daug[i]->getSpinStates() << "; ";
    }
    report(Severity::Info,"")<<endl;
    for(i=0;i<_ndaug;i++){
      _daug[i]->printTreeRec(newlevel);
    }
  }
}

void EvtParticle::printTree() const {
  
  report(Severity::Info,"EvtGen") << "This is the current decay chain"<<endl;
  report(Severity::Info,"") << "This top particle is "<<
    EvtPDL::name(_id).c_str()<<" " << this->mass() << " " << this->getP4() << endl;  
  
  this->printTreeRec(0);
  report(Severity::Info,"EvtGen") << "End of decay chain."<<endl;

}

std::string EvtParticle::treeStrRec(unsigned int level) const {

  size_t newlevel,i;
  newlevel = level +1;

  std::string retval="";

  for(i=0;i<_ndaug;i++){
    retval+=EvtPDL::name(_daug[i]->getId());
    if ( _daug[i]->getNDaug() > 0 ) {
      retval+= " (";
      retval+= _daug[i]->treeStrRec(newlevel);
      retval+= ") ";
    }
    else{
      if ( i+1 !=_ndaug) retval+=" ";
    }
  }

  return retval;
}


std::string EvtParticle::treeStr() const {

  std::string retval=EvtPDL::name(_id);
  retval+=" -> ";

  retval+=treeStrRec(0);

  return retval;
}

void EvtParticle::printParticle() const {

  switch (EvtPDL::getSpinType(_id)){ 
  case EvtSpinType::SCALAR:
    report(Severity::Info,"EvtGen") << "This is a scalar particle:"<<EvtPDL::name(_id).c_str()<<"\n";
    break;     
  case EvtSpinType::VECTOR:
    report(Severity::Info,"EvtGen") << "This is a vector particle:"<<EvtPDL::name(_id).c_str()<<"\n";
    break;     
  case EvtSpinType::TENSOR:
    report(Severity::Info,"EvtGen") << "This is a tensor particle:"<<EvtPDL::name(_id).c_str()<<"\n";
    break;
  case EvtSpinType::DIRAC:
    report(Severity::Info,"EvtGen") << "This is a dirac particle:"<<EvtPDL::name(_id).c_str()<<"\n";
    break;
  case EvtSpinType::PHOTON:
    report(Severity::Info,"EvtGen") << "This is a photon:"<<EvtPDL::name(_id).c_str()<<"\n";
    break;
  case EvtSpinType::NEUTRINO:
    report(Severity::Info,"EvtGen") << "This is a neutrino:"<<EvtPDL::name(_id).c_str()<<"\n";
    break;
  case EvtSpinType::STRING:
    report(Severity::Info,"EvtGen") << "This is a string:"<<EvtPDL::name(_id).c_str()<<"\n";
    break;
  default:
    report(Severity::Info,"EvtGen") <<"Unknown particle type in EvtParticle::printParticle()"<<endl;
    break;
  }
  report(Severity::Info,"EvtGen") << "Number of daughters:"<<_ndaug<<"\n";


}



void init_vector( EvtParticle **part ){
  *part = (EvtParticle *) new EvtVectorParticle;
} 


void init_scalar( EvtParticle **part ){
  *part = (EvtParticle *) new EvtScalarParticle;
} 

void init_tensor( EvtParticle **part ){
  *part = (EvtParticle *) new EvtTensorParticle;
} 

void init_dirac( EvtParticle **part ){
  *part = (EvtParticle *) new EvtDiracParticle;
} 

void init_photon( EvtParticle **part ){
  *part = (EvtParticle *) new EvtPhotonParticle;
} 

void init_neutrino( EvtParticle **part ){
  *part = (EvtParticle *) new EvtNeutrinoParticle;
} 

void init_string( EvtParticle **part ){
  *part = (EvtParticle *) new EvtStringParticle;
} 

double EvtParticle::initializePhaseSpace(
                   unsigned int numdaughter,EvtId *daughters, 
		   bool forceDaugMassReset, double poleSize,
		   int whichTwo1, int whichTwo2) {

  double m_b;
  unsigned int i;
  //lange
  //  this->makeDaughters(numdaughter,daughters);

  static EvtVector4R p4[100];
  static double mass[100];

  m_b = this->mass();

  //lange - Jan2,2002 - Need to check to see if the daughters of the parent
  // have changed. If so, delete them and start over.
  //report(Severity::Info,"EvtGen") << "the parent is\n";
  //if ( this->getParent() ) {
  //  if ( this->getParent()->getParent() ) this->getParent()->getParent()->printTree();
    //    this->getParent()->printTree();
  //}
  //report(Severity::Info,"EvtGen") << "and this is\n";
  //if ( this) this->printTree();
  bool resetDaughters=false;
  
  if ( numdaughter != this->getNDaug() && this->getNDaug() > 0 ) resetDaughters=true;
  if ( numdaughter == this->getNDaug() ) 
    for (i=0; i<numdaughter;i++) {
      if ( this->getDaug(i)->getId() != daughters[i] ) resetDaughters=true;
      //report(Severity::Info,"EvtGen") << EvtPDL::name(this->getDaug(i)->getId()) 
      //		    << " " << EvtPDL::name(daughters[i]) << endl;
    }

  if ( resetDaughters || forceDaugMassReset) {
    bool t1=true;
    //but keep the decay channel of the parent.
    this->deleteDaughters(t1);
    this->makeDaughters(numdaughter,daughters);
    bool massTreeOK = this->generateMassTree();
    if (massTreeOK == false) {return 0.0;}
  }

  double weight=0.;
  for (i=0; i<numdaughter;i++) {
    mass[i]=this->getDaug(i)->mass();
  }

  if ( poleSize<-0.1) {
    //special case to enforce 4-momentum conservation in 1->1 decays
    if (numdaughter==1) {
      this->getDaug(0)->init(daughters[0],EvtVector4R(m_b,0.0,0.0,0.0));
    }
    else{
      EvtGenKine::PhaseSpace( numdaughter, mass, p4, m_b );
      for(i=0;i<numdaughter;i++){
	this->getDaug(i)->init(daughters[i],p4[i]);
      }
    }
  }
  else  {
    if ( numdaughter != 3 ) {
      report(Severity::Error,"EvtGen") << "Only can generate pole phase space "
			     << "distributions for 3 body final states"
			     << endl<<"Will terminate."<<endl;
      ::abort();
    }
    bool ok=false;
    if ( (whichTwo1 == 1 && whichTwo2 == 0 ) ||
	 (whichTwo1 == 0 && whichTwo2 == 1 ) ) {
      weight=EvtGenKine::PhaseSpacePole( m_b, mass[0], mass[1], mass[2], 
					  poleSize, p4);
      this->getDaug(0)->init(daughters[0],p4[0]);
      this->getDaug(1)->init(daughters[1],p4[1]);
      this->getDaug(2)->init(daughters[2],p4[2]);
      ok=true;
    }
    if ( (whichTwo1 == 1 && whichTwo2 == 2 ) ||
	 (whichTwo1 == 2 && whichTwo2 == 1 ) ) {
      weight=EvtGenKine::PhaseSpacePole( m_b, mass[2], mass[1], mass[0], 
					  poleSize, p4);
      this->getDaug(0)->init(daughters[0],p4[2]);
      this->getDaug(1)->init(daughters[1],p4[1]);
      this->getDaug(2)->init(daughters[2],p4[0]);
      ok=true;
    }
    if ( (whichTwo1 == 0 && whichTwo2 == 2 ) ||
	 (whichTwo1 == 2 && whichTwo2 == 0 ) ) {
      weight=EvtGenKine::PhaseSpacePole( m_b, mass[1], mass[0], mass[2], 
					  poleSize, p4);
      this->getDaug(0)->init(daughters[0],p4[1]);
      this->getDaug(1)->init(daughters[1],p4[0]);
      this->getDaug(2)->init(daughters[2],p4[2]);
      ok=true;
    }
    if ( !ok) {
      report(Severity::Error,"EvtGen") << "Invalid pair of particle to generate a pole dist "
			     << whichTwo1 << " " << whichTwo2
			     << endl<<"Will terminate."<<endl;
      ::abort();
    }
  }

  return weight;
}

void EvtParticle::makeDaughters(unsigned int ndaugstore, std::vector<EvtId> idVector) {

  // Convert the STL vector method to use the array method for now, since the
  // array method pervades most of the EvtGen code...

  unsigned int nVector = idVector.size();
  if (nVector < ndaugstore) {
    report(Severity::Error,"EvtGen") << "Asking to make "<<ndaugstore<<" daughters when there "
			   << "are only "<<nVector<<" EvtId values available"<<endl;
    return;
  }

  EvtId *idArray=new EvtId[ndaugstore];
  unsigned int i;
  for (i = 0; i < ndaugstore; i++) {
    idArray[i] = idVector[i];
  }

  this->makeDaughters(ndaugstore, idArray);

  delete[] idArray;
}

void EvtParticle::makeDaughters( unsigned int ndaugstore, EvtId *id){

  unsigned int i;
  if ( _channel < 0 ) {
    setChannel(0);
  }
  EvtParticle* pdaug;  
  if (_ndaug!=0 ){
    if (_ndaug!=ndaugstore){
      report(Severity::Error,"EvtGen") << "Asking to make a different number of "
			     << "daughters than what was previously created."<<endl;
      report(Severity::Error,"EvtGen") << "Original parent:"<<EvtPDL::name(_id)<<endl;
      for (size_t i=0;i<_ndaug;i++){
	  report(Severity::Error,"EvtGen") << "Original daugther:"<<EvtPDL::name(getDaug(i)->getId())<<endl;
      }
      for (size_t i=0;i<ndaugstore;i++){
	  report(Severity::Error,"EvtGen") << "New Daug:"<<EvtPDL::name(id[i])<<endl;
      }
      report(Severity::Error,"EvtGen") << "Will terminate."<<endl;
      ::abort();
    }
  } 
  else{
    for(i=0;i<ndaugstore;i++){
      pdaug=EvtParticleFactory::particleFactory(EvtPDL::getSpinType(id[i]));
      pdaug->setId(id[i]);
      pdaug->addDaug(this);	
    }

  } //else
} //makeDaughters


void EvtParticle::setDecayProb(double prob) {

  if ( _decayProb == 0 )  _decayProb=new double;
  *_decayProb=prob;
}

std::string EvtParticle::getName() {
  
  std::string theName = _id.getName();
  return theName;

}
