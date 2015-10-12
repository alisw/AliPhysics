#ifndef ALIANALYSISTASKEVENTSELECTIONREF_H
#define ALIANALYSISTASKEVENTSELECTIONREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TArrayD;
class TClonesArray;
class TObjArray;
class TString;

class AliAnalysisUtils;
class AliAODTrack;
class AliESDtrack;
class AliESDtrackCuts;
class AliEmcalTriggerPatchInfo;
class AliEMCALGeometry;
class AliVCluster;
class AliVTrack;

namespace EMCalTriggerPtAnalysis {

class AliEMCalHistoContainer;

class AliAnalysisTaskEventSelectionRef : public AliAnalysisTaskSE {
public:
  enum EmcalTriggerClass{
    kCPREL0 = 0,
    kCPREG1,
    kCPREG2,
    kCPREJ1,
    kCPREJ2,
    kCPRntrig
  };
  AliAnalysisTaskEventSelectionRef();
  AliAnalysisTaskEventSelectionRef(const char *name);
  virtual ~AliAnalysisTaskEventSelectionRef();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);

  void SetOfflineEnergyThreshold(EmcalTriggerClass trgcls, double threshold) { fOfflineEnergyThreshold[trgcls] = threshold; }
  void SetClusterContainer(TString name) { fClusterContainerName = name; }

protected:
  void FillEventCounterHists(const char *triggerclass, double vtxz, bool isSelected, bool isOfflineSelected);
  Bool_t IsOfflineSelected(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const;

  void ProcessTrack(const char *triggerclass, const AliVTrack * track, bool isOfflineSelected);
  void ProcessCluster(const char *triggerclass, const AliVCluster *clust, bool isOfflineSelected);
  void ProcessOfflinePatch(const char * triggerclass, const AliEmcalTriggerPatchInfo * patch, bool isOfflineSelected);

  Bool_t TrackSelectionESD(AliESDtrack* track) ;
  Bool_t TrackSelectionAOD(AliAODTrack* track);

  void CreatePtBinning(TArrayD& binning) const;
  void CreateEnergyBinning(TArrayD& binning) const;

  TString                       fClusterContainerName;
  AliAnalysisUtils              *fAnalysisUtils;
  AliESDtrackCuts               *fTrackCuts;
  AliEMCalHistoContainer        *fHistos;                   //!
  AliEMCALGeometry              *fGeometry;                 //!
  TClonesArray                  *fTriggerPatchContainer;    //!
  TClonesArray                  *fClusterContainer;         //!
  TObjArray                     *fTrackContainer;           //!
  Double_t                      fOfflineEnergyThreshold[kCPRntrig];

  ClassDef(AliAnalysisTaskEventSelectionRef, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif
