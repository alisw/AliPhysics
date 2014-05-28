#include "EvtGenBase/EvtPatches.hh"
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
// Module: EvtPartProp.cc
//
// Description: Store particle properties for one particle.
//
// Modification history:
//
//    RYD     April 4, 1997        Module created
//
//------------------------------------------------------------------------
//
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include "EvtGenBase/EvtPartProp.hh"
#include "EvtGenBase/EvtAbsLineShape.hh"
#include "EvtGenBase/EvtFlatLineShape.hh"
#include "EvtGenBase/EvtManyDeltaFuncLineShape.hh"
#include "EvtGenBase/EvtRelBreitWignerBarrierFact.hh"
#include <string>
using std::fstream;

EvtPartProp::EvtPartProp():
  _id(-1,-1)
  ,_idchgconj(-1,-1)
  ,_chg3(0)
  ,_stdhep(0)
  ,_lundkc(0)
{
  _lineShape=0;
  _ctau=0.0;
  _name="*******";
  _spintype=EvtSpinType::SCALAR;
}

EvtPartProp::EvtPartProp(const EvtPartProp& x){

  if (0!=x._lineShape){
    _lineShape=x._lineShape->clone();
  }
  else{
    _lineShape=0;
  }
  _ctau=x._ctau;
  _name=x._name;
  _spintype=x._spintype;
  _id=x._id;
  _idchgconj=x._idchgconj;
  _chg3=x._chg3;
  _stdhep=x._stdhep;
  _lundkc=x._lundkc;

}

EvtPartProp::~EvtPartProp() {
  if ( _lineShape ) delete _lineShape;
  _lineShape=0;
}


void EvtPartProp::setName(std::string pname) {

  _name=pname;

}


EvtPartProp& EvtPartProp::operator=(const EvtPartProp& x){

  _lineShape=x._lineShape->clone();

  _ctau=x._ctau;
  _name=x._name;
  _chg3=x._chg3;
  _spintype=x._spintype;
  return *this;
}

void EvtPartProp::initLineShape(double mass, double width, double maxRange){

  _lineShape=new EvtRelBreitWignerBarrierFact(mass,width,maxRange,_spintype);

}

void EvtPartProp::newLineShape(std::string type){

  double m=_lineShape->getMass();
  double w=_lineShape->getWidth();
  double mR=_lineShape->getMaxRange();
  EvtSpinType::spintype  st=_lineShape->getSpinType();
  delete _lineShape;
  if ( type == "RELBW" ) {
    _lineShape=new EvtRelBreitWignerBarrierFact(m,w,mR,st);
  }
  if ( type == "NONRELBW" ) {
    _lineShape = new EvtAbsLineShape(m,w,mR,st);
  }
  if ( type == "FLAT" ) {
    _lineShape = new EvtFlatLineShape(m,w,mR,st);
  }
  if ( type == "MANYDELTAFUNC" ) {
    _lineShape = new EvtManyDeltaFuncLineShape(m,w,mR,st);
  }
}


void EvtPartProp::reSetMass(double mass) {
  if (!_lineShape) ::abort();
  _lineShape->reSetMass(mass);
}
void EvtPartProp::reSetWidth(double width){
  if (!_lineShape) ::abort();
  _lineShape->reSetWidth(width);
}

void EvtPartProp::setPWForDecay( int spin, EvtId d1, EvtId d2) { 
  if (!_lineShape) ::abort();
  _lineShape->setPWForDecay(spin,d1,d2);
}

void EvtPartProp::setPWForBirthL( int spin, EvtId par, EvtId othD) { 
  if (!_lineShape) ::abort();
  _lineShape->setPWForBirthL(spin,par,othD);
}


void EvtPartProp::reSetMassMin(double mass){
  if (!_lineShape) ::abort();
  _lineShape->reSetMassMin(mass);
}
void EvtPartProp::reSetMassMax(double mass){
  if (!_lineShape) ::abort();
  _lineShape->reSetMassMax(mass);
}
void EvtPartProp::reSetBlatt(double blatt){
  if (!_lineShape) ::abort();
  _lineShape->reSetBlatt(blatt);
}
void EvtPartProp::reSetBlattBirth(double blatt){
  if (!_lineShape) ::abort();
  _lineShape->reSetBlattBirth(blatt);
}
void EvtPartProp::includeBirthFactor(bool yesno){
  if (!_lineShape) ::abort();
  _lineShape->includeBirthFactor(yesno);
}
void EvtPartProp::includeDecayFactor(bool yesno){
  if (!_lineShape) ::abort();
  _lineShape->includeDecayFactor(yesno);
}






