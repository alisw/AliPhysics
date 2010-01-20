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
// Module: EvtLineShape.cc
//
// Description: Store particle properties for one particle.
//
// Modification history:
//
//    Lange      March 10, 2001        Module created
//    Dvoretskii June  03, 2002        Reimplemented rollMass()
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtManyDeltaFuncLineShape.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"
#include "EvtGenBase/EvtBlattWeisskopf.hh"
#include "EvtGenBase/EvtPropBreitWignerRel.hh"
#include "EvtGenBase/EvtPropBreitWigner.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSpinType.hh"

EvtManyDeltaFuncLineShape::EvtManyDeltaFuncLineShape() {

}

EvtManyDeltaFuncLineShape::~EvtManyDeltaFuncLineShape() {
}

EvtManyDeltaFuncLineShape::EvtManyDeltaFuncLineShape(double mass, double width, double maxRange, EvtSpinType::spintype sp) { 

  _mass=mass;
  _width=width;
  _spin=sp;
  _maxRange=maxRange;

  double maxdelta = width;

  _massMax=mass+maxdelta;
  _massMin=mass-maxdelta;

  if ( _massMin< 0. ) _massMin=0.;

}

EvtManyDeltaFuncLineShape::EvtManyDeltaFuncLineShape(const EvtManyDeltaFuncLineShape& x):
EvtAbsLineShape( x ) {
  _mass=x._mass;
  _width=x._width;
  _spin=x._spin;
  _massMax=x._massMax;
  _massMin=x._massMin;
  _maxRange=x._maxRange;

}

EvtManyDeltaFuncLineShape& EvtManyDeltaFuncLineShape::operator=(const EvtManyDeltaFuncLineShape& x){
  _mass=x._mass;
  _massMax=x._massMax;
  _massMin=x._massMin;
  _width=x._width;
  _maxRange=x._maxRange;
  _spin=x._spin;
  return *this;

}

EvtAbsLineShape* EvtManyDeltaFuncLineShape::clone() {

  return new EvtManyDeltaFuncLineShape(*this);
}


double EvtManyDeltaFuncLineShape::getMassProb(double mass, double massPar,int nDaug, double *massDau) {

  
  double dTotMass=0.;

  int i;
  for (i=0; i<nDaug; i++) {
    dTotMass+=massDau[i];
  }
  if ( (mass<dTotMass) ) return 0.;

  if ( massPar>0.0000000001 ) {
    if ( mass > massPar) return 0.;
  }

  return 1.;
}

double EvtManyDeltaFuncLineShape::getRandMass(EvtId*,int, EvtId*, EvtId*, double, double *) {

  int nDelta = int((_massMax - _massMin)/_width);
  nDelta++;
  double rand=EvtRandom::Flat(0.,float(nDelta));
  int randI=int(rand);
  return _massMin+randI*_width;
}




