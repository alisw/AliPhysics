#ifndef ALIANALYSISTASKEVENTSELECTIONREF_H
#define ALIANALYSISTASKEVENTSELECTIONREF_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class TArrayD;
class TClonesArray;
class THistManager;
class TObjArray;
class TString;

class AliAnalysisUtils;
class AliAODTrack;
class AliESDtrack;
class AliESDtrackCuts;
class AliEMCALTriggerPatchInfo;
class AliEMCALGeometry;
class AliVCluster;
class AliVTrack;

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEmcalTriggerOfflineSelection;

class AliAnalysisTaskEventSelectionRef : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEventSelectionRef();
  AliAnalysisTaskEventSelectionRef(const char *name);
  virtual ~AliAnalysisTaskEventSelectionRef();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);

  void SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }
  void SetClusterContainer(TString name) { fClusterContainerName = name; }

protected:
  void FillEventCounterHists(const char *triggerclass, double vtxz, bool isSelected, bool isOfflineSelected);

  void ProcessTrack(const char *triggerclass, const AliVTrack * track, bool isOfflineSelected);
  void ProcessCluster(const char *triggerclass, const AliVCluster *clust, bool isOfflineSelected);
  void ProcessOfflinePatch(const char * triggerclass, const AliEMCALTriggerPatchInfo * patch, bool isOfflineSelected);

  Bool_t TrackSelectionESD(AliESDtrack* track) ;
  Bool_t TrackSelectionAOD(AliAODTrack* track);

  void CreatePtBinning(TArrayD& binning) const;
  void CreateEnergyBinning(TArrayD& binning) const;

  TString                       fClusterContainerName;
  AliAnalysisUtils              *fAnalysisUtils;
  AliEmcalTriggerOfflineSelection *fTriggerSelection;
  AliESDtrackCuts               *fTrackCuts;
  THistManager                  *fHistos;                   //!
  AliEMCALGeometry              *fGeometry;                 //!
  TClonesArray                  *fTriggerPatchContainer;    //!
  TClonesArray                  *fClusterContainer;         //!
  TObjArray                     *fTrackContainer;           //!

private:

  AliAnalysisTaskEventSelectionRef(const AliAnalysisTaskEventSelectionRef &);
  AliAnalysisTaskEventSelectionRef &operator=(const AliAnalysisTaskEventSelectionRef &);

  ClassDef(AliAnalysisTaskEventSelectionRef, 1);
};

} /* namespace EMCALTriggerJets */

} /* namespace PWGJE */

#endif
