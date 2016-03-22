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

#include "EvtGenBase/EvtPredGen.hh"

#include "EvtGenBase/EvtRelBreitWignerBarrierFact.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"
#include "EvtGenBase/EvtPropBreitWignerRel.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtAmpPdf.hh"
#include "EvtGenBase/EvtMassAmp.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtIntervalFlatPdf.hh"
#include <algorithm>
EvtRelBreitWignerBarrierFact::EvtRelBreitWignerBarrierFact() {

}

EvtRelBreitWignerBarrierFact::~EvtRelBreitWignerBarrierFact() {
}

EvtRelBreitWignerBarrierFact::EvtRelBreitWignerBarrierFact(double mass, double width, double maxRange, EvtSpinType::spintype sp) :
  EvtAbsLineShape(mass,width,maxRange,sp)
{ // double mDaug1, double mDaug2, int l) {

  _includeDecayFact=true;
  _includeBirthFact=true;
  _mass=mass;
  _width=width;
  _spin=sp;
  _blattDecay=3.0;
  _blattBirth=1.0;
  _maxRange=maxRange;
  _errorCond=false;

  double maxdelta = 15.0*width;

  if ( maxRange > 0.00001 ) {
    _massMax=mass+maxdelta;
    _massMin=mass-maxRange;
  }
  else{
    _massMax=mass+maxdelta;
    _massMin=mass-15.0*width;
  }
 
  _massMax=mass+maxdelta;
  if ( _massMin< 0. ) _massMin=0.;
}

EvtRelBreitWignerBarrierFact::EvtRelBreitWignerBarrierFact(const EvtRelBreitWignerBarrierFact& x) :
  EvtAbsLineShape(x)
{
  _massMax=x._massMax;
  _massMin=x._massMin;
  _blattDecay=x._blattDecay;
  _blattBirth=x._blattBirth;
  _maxRange=x._maxRange;
  _includeDecayFact=x._includeDecayFact;
  _includeBirthFact=x._includeBirthFact;
  _errorCond=x._errorCond;

}

EvtRelBreitWignerBarrierFact& EvtRelBreitWignerBarrierFact::operator=(const EvtRelBreitWignerBarrierFact& x) {
  _mass=x._mass;
  _width=x._width;
  _spin=x._spin;
  _massMax=x._massMax;
  _massMin=x._massMin;
  _blattDecay=x._blattDecay;
  _blattBirth=x._blattBirth;
  _maxRange=x._maxRange;
  _includeDecayFact=x._includeDecayFact;
  _includeBirthFact=x._includeBirthFact;
  _errorCond=x._errorCond;
  
  return *this;
}

EvtAbsLineShape* EvtRelBreitWignerBarrierFact::clone() {

  return new EvtRelBreitWignerBarrierFact(*this);
}


double EvtRelBreitWignerBarrierFact::getMassProb(double mass, double massPar,int nDaug, double *massDau) {

  _errorCond=false;
  //return EvtAbsLineShape::getMassProb(mass,massPar,nDaug,massDau);
  if (nDaug!=2) return EvtAbsLineShape::getMassProb(mass,massPar,nDaug,massDau);

  double dTotMass=0.;

  int i;
  for (i=0; i<nDaug; i++) {
    dTotMass+=massDau[i];
  }
  //report(Severity::Info,"EvtGen") << mass << " " << massPar << " " << dTotMass << " "<< endl;
  //    if ( (mass-dTotMass)<0.0001 ) return 0.;
  //report(Severity::Info,"EvtGen") << mass << " " << dTotMass << endl;
  if ( (mass<dTotMass) ) return 0.;

  if ( _width< 0.0001) return 1.;

  if ( massPar>0.0000000001 ) {
    if ( mass > massPar) return 0.;
  }

  if ( _errorCond ) return 0.;

  // we did all the work in getRandMass
  return 1.;
}

double EvtRelBreitWignerBarrierFact::getRandMass(EvtId *parId,int nDaug, EvtId *dauId, EvtId *othDaugId, double maxMass, double *dauMasses) {
  if ( nDaug!=2) return EvtAbsLineShape::getRandMass(parId,nDaug,dauId,othDaugId,maxMass,dauMasses);

  if ( _width< 0.00001) return _mass;

  //first figure out L - take the lowest allowed.

  EvtSpinType::spintype spinD1=EvtPDL::getSpinType(dauId[0]);
  EvtSpinType::spintype spinD2=EvtPDL::getSpinType(dauId[1]);

  int t1=EvtSpinType::getSpin2(spinD1);
  int t2=EvtSpinType::getSpin2(spinD2);
  int t3=EvtSpinType::getSpin2(_spin);

  int Lmin=-10;


  // the user has overridden the partial wave to use.
  for (unsigned int vC=0; vC<_userSetPW.size(); vC++) {
    if ( dauId[0]==_userSetPWD1[vC] &&  dauId[1]==_userSetPWD2[vC] ) {
      Lmin=2*_userSetPW[vC]; 
    }
    if ( dauId[0]==_userSetPWD2[vC] &&  dauId[1]==_userSetPWD1[vC] ) {
      Lmin=2*_userSetPW[vC]; 
    }
  }
  
  // allow for special cases.
  if (Lmin<-1 ) {
    
    //There are some things I don't know how to deal with
    if ( t3>4) return EvtAbsLineShape::getRandMass(parId,nDaug,dauId,othDaugId,maxMass,dauMasses);
    if ( t1>4) return EvtAbsLineShape::getRandMass(parId,nDaug,dauId,othDaugId,maxMass,dauMasses);
    if ( t2>4) return EvtAbsLineShape::getRandMass(parId,nDaug,dauId,othDaugId,maxMass,dauMasses);
    
    //figure the min and max allowwed "spins" for the daughters state
    Lmin=std::max(t3-t2-t1,std::max(t2-t3-t1,t1-t3-t2));
    if (Lmin<0) Lmin=0;
    assert(Lmin==0||Lmin==2||Lmin==4);
  }

  //double massD1=EvtPDL::getMeanMass(dauId[0]);
  //double massD2=EvtPDL::getMeanMass(dauId[1]);
  double massD1=dauMasses[0];
  double massD2=dauMasses[1];

  // I'm not sure how to define the vertex factor here - so retreat to nonRel code.
  if ( (massD1+massD2)> _mass ) return  EvtAbsLineShape::getRandMass(parId,nDaug,dauId,othDaugId,maxMass,dauMasses);

  //parent vertex factor not yet implemented
  double massOthD=-10.;
  double massParent=-10.;
  int birthl=-10;
  if ( othDaugId) {
    EvtSpinType::spintype spinOth=EvtPDL::getSpinType(*othDaugId);
    EvtSpinType::spintype spinPar=EvtPDL::getSpinType(*parId);
    
    int tt1=EvtSpinType::getSpin2(spinOth);
    int tt2=EvtSpinType::getSpin2(spinPar);
    int tt3=EvtSpinType::getSpin2(_spin);
    
    
    //figure the min and max allowwed "spins" for the daughters state
    if ( (tt1<=4) && ( tt2<=4) ) {
      birthl=std::max(tt3-tt2-tt1,std::max(tt2-tt3-tt1,tt1-tt3-tt2));
      if (birthl<0) birthl=0;
    
      massOthD=EvtPDL::getMeanMass(*othDaugId);
      massParent=EvtPDL::getMeanMass(*parId);
    
    }

    // allow user to override
    for (size_t vC=0; vC<_userSetBirthPW.size(); vC++) {
      if ( *othDaugId==_userSetBirthOthD[vC] && *parId==_userSetBirthPar[vC]){
	birthl=2*_userSetBirthPW[vC];
      } 
    }

  }
  double massM=_massMax;
  if ( (maxMass > -0.5) && (maxMass < massM) ) massM=maxMass;

  //special case... if the parent mass is _fixed_ we can do a little better
  //and only for a two body decay as that seems to be where we have problems

  // Define relativistic propagator amplitude

  EvtTwoBodyVertex vd(massD1,massD2,_mass,Lmin/2);
  vd.set_f(_blattDecay);
  EvtPropBreitWignerRel bw(_mass,_width);
  EvtMassAmp amp(bw,vd);


  if ( _includeDecayFact) {
    amp.addDeathFact();
    amp.addDeathFactFF();
  }
  if ( massParent>-1.) {
    if ( _includeBirthFact ) {

      EvtTwoBodyVertex vb(_mass,massOthD,massParent,birthl/2);
      vb.set_f(_blattBirth);
      amp.setBirthVtx(vb);
      amp.addBirthFact();
      amp.addBirthFactFF();
    }
  }


  EvtAmpPdf<EvtPoint1D> pdf(amp);

  // Estimate maximum and create predicate for accept reject


  double tempMaxLoc=_mass;
  if ( maxMass>-0.5 && maxMass<_mass) tempMaxLoc=maxMass;
  double tempMax=_massMax;
  if ( maxMass>-0.5 && maxMass<_massMax) tempMax=maxMass;
  double tempMinMass=_massMin;
  if ( massD1+massD2 > _massMin) tempMinMass=massD1+massD2;

  //redo sanity check - is there a solution to our problem.
  //if not return an error condition that is caught by the
  //mass prob calculation above.
  if ( tempMinMass > tempMax ) {
    _errorCond=true;
    return tempMinMass;
  }
  
  if ( tempMaxLoc < tempMinMass) tempMaxLoc=tempMinMass;

  double safetyFactor=1.2;

  EvtPdfMax<EvtPoint1D> max(safetyFactor*pdf.evaluate(EvtPoint1D(tempMinMass,tempMax,tempMaxLoc)));

  EvtPdfPred<EvtPoint1D> pred(pdf);
  pred.setMax(max);

  EvtIntervalFlatPdf flat(tempMinMass,tempMax);
  EvtPdfGen<EvtPoint1D> gen(flat);
  EvtPredGen<EvtPdfGen<EvtPoint1D>,EvtPdfPred<EvtPoint1D> > predgen(gen,pred);

  EvtPoint1D point = predgen();
  return point.value();

};









