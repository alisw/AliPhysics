#ifndef ALIANALYSISTASKEMCALTRIGGERBASE_H
#define ALIANALYSISTASKEMCALTRIGGERBASE_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>
#include <vector>

#include "AliAnalysisTaskEmcal.h"
#include "AliCutValueRange.h"

class TClonesArray;
class THistManager;
class AliOADBContainer;

namespace EMCalTriggerPtAnalysis {

class AliEmcalTriggerOfflineSelection;

class AliAnalysisTaskEmcalTriggerBase : public AliAnalysisTaskEmcal{
public:
  enum OnlineTrigger_t{
    kCPREG1 = 0,
    kCPREG2 = 1,
    kCPREJ1 = 2,
    kCPREJ2 = 3,
    kCPREL0 = 4
  };
  AliAnalysisTaskEmcalTriggerBase();
  AliAnalysisTaskEmcalTriggerBase(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerBase();

  void SetDownscaleOADB(TString oadbname) { fNameDownscaleOADB = oadbname; }
  void SetExcludeNoiseEvents(Bool_t doExclude = true) { fRejectNoiseEvents = doExclude; }
  void SetSelectNoiseEvents(Bool_t doSelect = true) { fSelectNoiseEvents = doSelect; }
  void AddMaskedFastor(int fastorID) { fMaskedFastors.push_back(fastorID); }
  void SetMaskedFastorOADB(TString oadbname) { fNameMaskedFastorOADB = oadbname; }
  void SetOnlineSelectionMethodV1(bool doselect) { fUseOnlineTriggerSelectorV1 = true; }
  void SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }
  void SetVertexCut(double zmin, double zmax) { fVertexCut.SetLimits(zmin, zmax); }
  void UseTriggerPatches(Bool_t doUse) { fTriggerStringFromPatches = doUse; }

protected:
  virtual void UserCreateOutputObjects();
  virtual bool IsEventSelected();
  virtual bool IsUserEventSelected() { return true; }
  virtual void CreateUserHistos() = 0;
  virtual void CreateUserObjects() = 0;
  virtual void ExecOnce();
  virtual void RunChanged(Int_t runnuber);

  virtual void UserFillHistosBeforeEventSelection() { }
  virtual void UserFillHistosAfterEventSelection() { }

  std::vector<TString> GetSupportedTriggers();

  Double_t GetTriggerWeight(const TString &triggerclass) const;
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;
  void TriggerSelection();

  bool SelectOnlineTrigger(OnlineTrigger_t trigger);
  bool SelectOnlineTriggerV1(OnlineTrigger_t trigger);


  THistManager                    *fHistos;                   ///< Task Histogram container

  AliEmcalTriggerOfflineSelection *fTriggerSelection;         ///< Offline trigger selection
  Bool_t                          fTriggerStringFromPatches;  ///< Do rebuild the trigger string from trigger patches
  std::vector<TString>            fSelectedTriggers;          //!<! Triggers selected for given event

  AliCutValueRange<double>        fVertexCut;                 ///< Cut on the z-position of the primary vertex

  TString                         fNameDownscaleOADB;         ///< Name of the downscale OADB container
  AliOADBContainer                *fDownscaleOADB;            //!<! Container with downscale factors for different triggers
  TObjArray                       *fDownscaleFactors;         //!<! Downscalfactors for given run
  TString                         fNameMaskedFastorOADB;      ///< Name of the masked fastor OADB container
  AliOADBContainer                *fMaskedFastorOADB;         //!<! Container with masked fastors
  std::vector<int>                fMaskedFastors;             ///< List of masked fastors

  Bool_t                          fSelectNoiseEvents;         ///< Explicitly select events triggered only by noisy fastors
  Bool_t                          fRejectNoiseEvents;         ///< Reject events triggered by noisy fastors
  Bool_t                          fUseOnlineTriggerSelectorV1;  ///< Switch separating between different trigger selection methods


private:
  AliAnalysisTaskEmcalTriggerBase(const AliAnalysisTaskEmcalTriggerBase &);
  AliAnalysisTaskEmcalTriggerBase &operator=(const AliAnalysisTaskEmcalTriggerBase &);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALTRIGGERBASE_H */
