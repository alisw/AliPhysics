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
// Module: EvtGen/EvtPartProp.hh
//
// Description: Class to keep the particle properties for
//              one particle
//
// Modification history:
//
//    RYD     April 4, 1997         Module created
//
//------------------------------------------------------------------------

#ifndef EVTPARTPROP_HH
#define EVTPARTPROP_HH

#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtAbsLineShape.hh"


class EvtPartProp {

public:

  EvtPartProp(); 
  EvtPartProp(const EvtPartProp& x); 

  ~EvtPartProp(); 

  double getMass() {return _lineShape->getMass();} 
  double getMassMin() {return _lineShape->getMassMin();} 
  double getMassMax() {return _lineShape->getMassMax();} 
  double getMaxRange() {return _lineShape->getMaxRange();} 
  double getWidth() {return _lineShape->getWidth();} 

  double getRandMass(EvtId *parId, int nDaug, EvtId *dauId,EvtId *othDauId,double maxMass, double *dauMasses) {return _lineShape->getRandMass(parId,nDaug,dauId,othDauId,maxMass,dauMasses);}
  double getMassProb(double mass, double massPar, int nDaug, double *massDau) { return _lineShape->getMassProb(mass,massPar,nDaug,massDau);}

  double getctau() {return _ctau; } 
  void   setctau(double tau) { _ctau=tau; }

  int    getChg3() {return _chg3; } 
  void   setChg3(int c3) { _chg3=c3; }

  EvtSpinType::spintype  getSpinType() {return _spintype; }
  void   setSpinType(EvtSpinType::spintype stype ) { _spintype=stype; }

  const std::string&  getName() {return _name;}
  void   setName(std::string pname);

  EvtId  getId() {return _id;}
  void   setId(EvtId id) {_id=id;}

  EvtId  getIdChgConj() {return _idchgconj;}
  void   setIdChgConj(EvtId idchgconj) {_idchgconj=idchgconj;}

  int  getStdHep() {return _stdhep;}
  void   setStdHep(int stdhep) {_stdhep=stdhep;}

  int  getLundKC() {return _lundkc;}
  void   setLundKC(int lundkc) {_lundkc=lundkc;}

  EvtAbsLineShape* getLineShape() {return _lineShape;}
  void initLineShape(double mass, double width, double maxRange);
  //  void initLineShape(double mass, double width, double maxRange, double mDaug1, double mDaug2, int l);

  // setLineShape takes ownership of l
  void setLineShape(EvtAbsLineShape *l) { _lineShape=l;}
  double rollMass(){return _lineShape->rollMass();}

  EvtPartProp& operator=(const EvtPartProp& x);

  void reSetMass(double mass);
  void reSetWidth(double width);

  void reSetMassMin(double mass);
  void reSetMassMax(double mass);
  void reSetBlatt(double blatt);
  void reSetBlattBirth(double blatt);
  void includeBirthFactor(bool yesno);
  void includeDecayFactor(bool yesno);
  void newLineShape(std::string type);
  void setPWForDecay( int spin, EvtId d1, EvtId d2);
  void setPWForBirthL( int spin, EvtId par, EvtId othD);

private:

  EvtAbsLineShape *_lineShape;

  double _ctau;
  EvtId  _id;
  EvtId  _idchgconj;
  EvtSpinType::spintype  _spintype;
  int _chg3;
  int _stdhep;
  int _lundkc;
  std::string _name;
  
}; 

#endif

