#ifndef ALIJETTRIGGERSELECTIONTASK_H
#define ALIJETTRIGGERSELECTIONTASK_H

// $Id$

class AliEmcalJet;

#include "AliAnalysisTaskEmcalJet.h"

class AliJetTriggerSelectionTask : public AliAnalysisTaskEmcalJet {
 public:

  AliJetTriggerSelectionTask();
  AliJetTriggerSelectionTask(const char *name);
  virtual ~AliJetTriggerSelectionTask() {;}

  void                        SetMaxDistance(Double_t d) { fMaxDistance2    = d*d ; }
  void                        SetEnergyThreshold(TF1 *f) { fEnergyThreshold = f   ; }
  void                        SetTriggerBits(UInt_t d)   { fTriggerBits     = d   ; }

 protected:
  Bool_t                      Run();
  void                        ExecOnce();
  Bool_t                      RetrieveEventObjects();
  void                        FindTriggers();
  void                        SelectJets();
  Bool_t                      IsTriggerJet(AliEmcalJet *jet);

  TF1                        *fEnergyThreshold;                // energy threshold vs. centrality
  Double_t                    fMaxDistance2;                   // max distance square between trigger patch and jet
  UInt_t                      fTriggerBits;                    // trigger bit to be set

  Bool_t                      fTaskSettingsOk;                 //!if false, don't execute task  
  Int_t                       fNTriggers;                      //!number of triggers in the current event
  Double_t                    fTrigPos[999][2];                //!(eta,phi) trigger positions in the current event
  AliVVZERO                  *fVZERO;                          //!Event V0 object
  Double_t                    fV0ATotMult;                     //!Event V0A total multiplicity
  Double_t                    fV0CTotMult;                     //!Event V0C total multiplicity
 
 private:
  AliJetTriggerSelectionTask(const AliJetTriggerSelectionTask&);            // not implemented
  AliJetTriggerSelectionTask &operator=(const AliJetTriggerSelectionTask&); // not implemented

  ClassDef(AliJetTriggerSelectionTask, 1) // jet trigger selection task
};
#endif
