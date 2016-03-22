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
//    Lange     March 10, 2001        Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "EvtGenBase/EvtAbsLineShape.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"
#include "EvtGenBase/EvtPropBreitWigner.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"

using namespace std;

EvtAbsLineShape::EvtAbsLineShape() {
}

EvtAbsLineShape::~EvtAbsLineShape() {
}

EvtAbsLineShape::EvtAbsLineShape(double mass, double width, double maxRange, EvtSpinType::spintype sp) {

  _includeDecayFact = false;
  _includeBirthFact = false;
  _mass=mass;
  _width=width;
  _spin=sp;
  _maxRange=maxRange;
  double maxdelta=15.0*width;
  //if ( width>0.001 ) {
  //  if ( 5.0*width < 0.6 ) maxdelta = 0.6;
  //}
  if ( maxRange > 0.00001 ) {
    _massMax=mass+maxdelta;
    _massMin=mass-maxRange;
  }
  else{
    _massMax=mass+maxdelta;
    _massMin=mass-15.0*width;
  }
  if ( _massMin< 0. ) _massMin=0.;
  _massMax=mass+maxdelta;
}

EvtAbsLineShape::EvtAbsLineShape(const EvtAbsLineShape& x){

  _includeDecayFact = x._includeDecayFact;
  _includeBirthFact = x._includeBirthFact;
  _mass=x._mass;
  _massMax=x._massMax;
  _massMin=x._massMin;
  _width=x._width;
  _spin=x._spin;
  _maxRange=x._maxRange;
}

EvtAbsLineShape& EvtAbsLineShape::operator=(const EvtAbsLineShape& x){

  _includeDecayFact = x._includeDecayFact;
  _includeBirthFact = x._includeBirthFact;
  _mass=x._mass;
  _massMax=x._massMax;
  _massMin=x._massMin;
  _width=x._width;
  _spin=x._spin;
  _maxRange=x._maxRange;
  return *this; 
}

EvtAbsLineShape* EvtAbsLineShape::clone() {

  return new EvtAbsLineShape(*this);
}


double EvtAbsLineShape::rollMass() {

  double ymin, ymax;
  double temp;

  if ( _width < 0.0001 ) {
    return _mass;
  }
  else{
    ymin = atan( 2.0*(_massMin-_mass)/_width);
    ymax = atan( 2.0*(_massMax-_mass)/_width);

    temp= ( _mass + ((_width/2.0)*tan(EvtRandom::Flat(ymin,ymax))));

    return temp;
  }
}
double EvtAbsLineShape::getRandMass(EvtId *parId, int /* nDaug */, EvtId * /*dauId*/, EvtId */*othDaugId*/, double maxMass, double */*dauMasses*/) {

  if ( _width< 0.0001) return _mass;
  //its not flat - but generated according to a BW

  if (maxMass>0&&maxMass<_massMin) {
    report(Severity::Debug,"EvtGen") << "In EvtAbsLineShape::getRandMass:"<<endl;
    report(Severity::Debug,"EvtGen") << "Cannot create a particle with a minimal mass of "
                           << _massMin << " from a "<<EvtPDL::name(*parId)
                           << " decay with available left-over mass-energy " << maxMass
                           << ". Returning 0.0 mass. The rest of this decay chain will probably fail..." << endl;
    return 0.0;
  }

  double mMin=_massMin;
  double mMax=_massMax;
  if ( maxMass>-0.5 && maxMass< mMax) mMax=maxMass;
  double ymin = atan( 2.0*(mMin-_mass)/_width);
  double ymax = atan( 2.0*(mMax-_mass)/_width);
  
  return ( _mass + ((_width/2.0)*tan(EvtRandom::Flat(ymin,ymax))));
  //  return EvtRandom::Flat(_massMin,_massMax);
};

double EvtAbsLineShape::getMassProb(double mass, double massPar, int nDaug, double *massDau) {

  double dTotMass=0.;
  if ( nDaug>1) {
    int i;
    for (i=0; i<nDaug; i++) {
      dTotMass+=massDau[i];
    }
    //report(Severity::Info,"EvtGen") << mass << " " << massPar << " " << dTotMass << " "<< endl;
    //    if ( (mass-dTotMass)<0.0001 ) return 0.;
    if ( (mass<dTotMass) ) return 0.;
  }
  if ( _width< 0.0001) return 1.;

  // no parent - lets not panic
  if ( massPar>0.0000000001 ) {
    if ( mass > massPar) return 0.;
  }
  //Otherwise return the right value.
  //Fortunately we have generated events according to a non-rel BW, so 
  //just return..
  //EvtPropBreitWigner bw(_mass,_width);
  //EvtPropFactor<EvtTwoBodyVertex> f(bw);
  //EvtComplex fm=f.eval(mass);
  //EvtComplex fm0=f.eval(_mass);
  //return (abs(fm)*abs(fm))/(abs(fm0)*abs(fm0));
  return 1.0;
}


