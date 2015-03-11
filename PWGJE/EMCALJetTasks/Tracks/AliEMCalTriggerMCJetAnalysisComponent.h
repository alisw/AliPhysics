#ifndef ALIEMCALTRIGGERMCJETANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERMCJETANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliEMCalTriggerTracksAnalysisComponent.h"

class TString;
class AliVParticle;
class AliEmcalJet;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;

class AliEMCalTriggerMCJetAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerMCJetAnalysisComponent();
  AliEMCalTriggerMCJetAnalysisComponent(const char * name);
  virtual ~AliEMCalTriggerMCJetAnalysisComponent() {}

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  void SetUsePatches(Bool_t doUse = kTRUE) { fUsePatches = doUse; }
  void SetMinimumJetPt(Double_t minpt) { fMinimumJetPt = minpt; }

protected:
  void FillHistogram(const TString &histname, const AliVParticle *track, const AliEmcalJet *jet, double vz, double weight);
  void FillJetHistogram(const TString &histname, const AliEmcalJet *recjet, double vz, double weight);

  Double_t                fMinimumJetPt;                      // Min. pt request for the jet
  Bool_t                  fUsePatches;                        // Use patches for trigger decision

  ClassDef(AliEMCalTriggerMCJetAnalysisComponent, 1);         // Analysis component for MC Jets
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERMCJETANALYSISCOMPONENT_H */
