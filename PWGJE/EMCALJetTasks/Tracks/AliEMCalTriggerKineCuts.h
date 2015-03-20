#ifndef ALIEMCALTRIGGERKINECUTS_H
#define ALIEMCALTRIGGERKINECUTS_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TObject.h>
#include "AliCutValueRange.h"

class AliVParticle;

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerKineCuts: public TObject {
public:
  AliEMCalTriggerKineCuts();
  virtual ~AliEMCalTriggerKineCuts() {}

  void SetPtRange(double ptmin, double ptmax) { fPtCut.SetLimits(ptmin, ptmax); }
  void SetEtaRange(double etamin, double etamax) { fEtaCut.SetLimits(etamin, etamax); }
  void SetPhiRange(double phimin, double phimax) { fPhiCut.SetLimits(phimin, phimax); }

  bool IsSelected(const AliVParticle *const track) const;

protected:
  AliCutValueRange<double>     fPtCut;
  AliCutValueRange<double>     fEtaCut;
  AliCutValueRange<double>     fPhiCut;

  ClassDef(AliEMCalTriggerKineCuts, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERKINECUTS_H */
