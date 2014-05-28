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
// Module: EvtGen/EvtAbsLineShape.hh
//
// Description: Class to keep the particle properties for
//              one particle
//
// Modification history:
//
//    Lange     March 10, 2001         Module created
//
//------------------------------------------------------------------------

#ifndef EVTRELBREITWIGNERBARRIERFACT_HH
#define EVTRELBREITWIGNERBARRIERFACT_HH

#include "EvtGenBase/EvtAbsLineShape.hh"

class EvtRelBreitWignerBarrierFact :public EvtAbsLineShape {

public:

  EvtRelBreitWignerBarrierFact(); 
  EvtRelBreitWignerBarrierFact(double mass, double width, double maxRange, EvtSpinType::spintype sp);
    //figure the m1 and l on the fly
    //			       double mDaug1, double mDaug2, int l); 
  ~EvtRelBreitWignerBarrierFact();
  EvtRelBreitWignerBarrierFact& operator=(const EvtRelBreitWignerBarrierFact& x);
  EvtRelBreitWignerBarrierFact(const EvtRelBreitWignerBarrierFact& x); 

  EvtAbsLineShape* clone();

  double getMassProb(double mass, double massPar, int nDaug, double *massDau);
  // othDaugId is the other daughter of the parent in the case of a two body decay (only!)
  // ie B->rho K rho->pipi, othDaugId = K
   double getRandMass(EvtId *parId, int nDaug, EvtId *dauId, EvtId *othDaugId, double maxMass, double *dauMasses);

  virtual void reSetBlatt(double blatt) { _blattDecay = blatt; }
  virtual void reSetBlattBirth(double blatt) { _blattBirth = blatt; }

protected:

  double _blattDecay;
  double _blattBirth;
  bool _errorCond;

}; 

#endif

