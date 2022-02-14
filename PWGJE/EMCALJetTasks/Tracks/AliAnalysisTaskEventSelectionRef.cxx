/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <TArrayD.h>
#include <TClonesArray.h>
#include <THashList.h>
#include <THistManager.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TString.h>

#include "AliAnalysisUtils.h"
#include "AliAODTrack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliInputEventHandler.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

#include "AliAnalysisTaskEventSelectionRef.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEventSelectionRef)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEventSelectionRef::AliAnalysisTaskEventSelectionRef():
  AliAnalysisTaskSE(),
  fClusterContainerName(""),
  fAnalysisUtils(nullptr),
  fTriggerSelection(nullptr),
  fTrackCuts(nullptr),
  fHistos(nullptr),
  fGeometry(nullptr),
  fTriggerPatchContainer(nullptr),
  fClusterContainer(nullptr),
  fTrackContainer(nullptr)
{
}

AliAnalysisTaskEventSelectionRef::AliAnalysisTaskEventSelectionRef(const char *name):
  AliAnalysisTaskSE(name),
  fClusterContainerName(""),
  fAnalysisUtils(nullptr),
  fTriggerSelection(nullptr),
  fTrackCuts(nullptr),
  fHistos(nullptr),
  fGeometry(nullptr),
  fTriggerPatchContainer(nullptr),
  fClusterContainer(nullptr),
  fTrackContainer(nullptr)
{
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEventSelectionRef::~AliAnalysisTaskEventSelectionRef() {
  if(fAnalysisUtils) delete fAnalysisUtils;
  if(fTriggerSelection) delete fTriggerSelection;
  if(fTrackCuts) delete fTrackCuts;
  if(fTrackContainer) delete fTrackContainer;
}

void AliAnalysisTaskEventSelectionRef::UserCreateOutputObjects(){
  fAnalysisUtils = new AliAnalysisUtils;

  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(true, 1);
  fTrackCuts->SetName("Standard Track cuts");
  fTrackCuts->SetMinNCrossedRowsTPC(120);
  fTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

  fHistos = new THistManager("Ref");

  TArrayD ptbinning, energybinning;
  CreatePtBinning(ptbinning);
  CreateEnergyBinning(energybinning);

  TString triggers[6] = {"MB", "EMC7", "EJ1", "EJ2", "EG1", "EG2"};
  TString selectionStatus[3] = {"All","Accepted","Rejected"};
  for(TString *trgit = triggers; trgit < triggers + sizeof(triggers)/sizeof(TString); ++trgit){
    fHistos->CreateTH1(Form("hEventCount%sBeforeEventSelection", trgit->Data()), Form("Event count for trigger %s before event selection", trgit->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCount%sBeforeOfflineTrigger", trgit->Data()), Form("Event count for trigger %s before offline selection", trgit->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hEventCount%sAfterOfflineTrigger", trgit->Data()), Form("Event count for trigger %s before offline selection", trgit->Data()), 1, 0.5, 1.5);
    fHistos->CreateTH1(Form("hVertexTrigger%sBeforeEventSelection", trgit->Data()), Form("Vertex Distribution for trigger %s before event selection", trgit->Data()), 400, -40, 40);
    fHistos->CreateTH1(Form("hVertexTrigger%sBeforeOfflineTrigger", trgit->Data()), Form("Vertex Distribution for trigger %s before offline trigger", trgit->Data()), 400, -40, 40);
    fHistos->CreateTH1(Form("hVertexTrigger%sAfterOfflineTrigger", trgit->Data()), Form("Vertex Distribution for trigger %s after offline trigger", trgit->Data()), 400, -40, 40);
    for(TString *statit = selectionStatus; statit < selectionStatus + sizeof(selectionStatus)/sizeof(TString); ++statit){
      fHistos->CreateTH1(Form("hTrackPt%s%s", statit->Data(), trgit->Data()), Form("Pt spectrum of tracks in %s events of trigger %s", statit->Data(), trgit->Data()), ptbinning);
      fHistos->CreateTH1(Form("hClusterEnergy%s%s", statit->Data(), trgit->Data()), Form("Cluster energy spectrum in %s events of trigger %s", statit->Data(), trgit->Data()), energybinning);
      if(trgit->CompareTo("MB"))
        fHistos->CreateTH1(Form("hPatchEnergy%s%s", statit->Data(), trgit->Data()), Form("Patch energy spectrum in %s events of trigger %s", statit->Data(), trgit->Data()), energybinning);
    }
    // Make histos for patch energy spectra for Min. Bias for different triggers
    fHistos->CreateTH1(Form("hPatchEnergy%sMinBias", trgit->Data()), Form("Patch energy spectrum for %s patches found in MB events", trgit->Data()), energybinning);
  }


  fTrackContainer = new TObjArray;
  fTrackContainer->SetOwner(kFALSE);

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEventSelectionRef::UserExec(Option_t *){
  // Select event
  if(!fGeometry){
    fGeometry = AliEMCALGeometry::GetInstanceFromRunNumber(fInputEvent->GetRunNumber());
  }
  fTriggerPatchContainer = static_cast<TClonesArray *>(fInputEvent->FindListObject("EmcalTriggers"));
  fClusterContainer = static_cast<TClonesArray *>(fInputEvent->FindListObject(fClusterContainerName.Data()));
  TString triggerstring = fInputEvent->GetFiredTriggerClasses();
  UInt_t selectionstatus = fInputHandler->IsEventSelected();
  Bool_t isMinBias = selectionstatus & AliVEvent::kINT7,
      isEJ1 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ1"),
      isEJ2 = (selectionstatus & AliVEvent::kEMCEJE) && triggerstring.Contains("EJ2"),
      isEG1 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG1"),
      isEG2 = (selectionstatus & AliVEvent::kEMCEGA) && triggerstring.Contains("EG2"),
      isEMC7 = (selectionstatus & AliVEvent::kEMC7) && triggerstring.Contains("CEMC7");
  if(!(isMinBias || isEMC7 || isEG1 || isEG2 || isEJ1 || isEJ2)) return;
  const AliVVertex *vtx = fInputEvent->GetPrimaryVertex();
  //if(!fInputEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return;         // reject pileup event
  if(vtx->GetNContributors() < 1) return;
  if(fInputEvent->IsA() == AliESDEvent::Class() && fAnalysisUtils->IsFirstEventInChunk(fInputEvent)) return;
  bool isSelected  = kTRUE;
  if(!fAnalysisUtils->IsVertexSelected2013pA(fInputEvent)) isSelected = kFALSE;       // Apply new vertex cut
  if(fAnalysisUtils->IsPileUpEvent(fInputEvent)) isSelected = kFALSE;       // Apply new vertex cut
  // Apply vertex z cut
  if(vtx->GetZ() < -10. || vtx->GetZ() > 10.) isSelected = kFALSE;

  // prepare tracks: Apply track selection and extrapolate to EMCAL surface
  fTrackContainer->Clear();
  AliVTrack *trk(nullptr); Double_t etaEMCAL(0.), phiEMCAL(0.); Int_t supermoduleID(-1);
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); itrk++){
    trk = static_cast<AliESDtrack *>(fInputEvent->GetTrack(itrk));
    if(TMath::Abs(trk->Eta()) > 0.6) continue;
    if(trk->IsA() == AliESDtrack::Class()){
      AliESDtrack *origtrack = static_cast<AliESDtrack *>(trk);
      if(!TrackSelectionESD(origtrack)) continue;
      AliESDtrack copytrack(*origtrack);
      AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&copytrack);
      etaEMCAL = copytrack.GetTrackEtaOnEMCal();
      phiEMCAL = copytrack.GetTrackPhiOnEMCal();
    } else {
      AliAODTrack *origtrack = static_cast<AliAODTrack *>(trk);
      if(!TrackSelectionAOD(origtrack)) continue;
      AliAODTrack copytrack(*origtrack);
      AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&copytrack);
      etaEMCAL = copytrack.GetTrackEtaOnEMCal();
      phiEMCAL = copytrack.GetTrackPhiOnEMCal();
    }
    if(!fGeometry->SuperModuleNumberFromEtaPhi(etaEMCAL, phiEMCAL, supermoduleID)) continue;
    if(supermoduleID < 0 || supermoduleID >= 10) continue;
    fTrackContainer->Add(trk);
  }

  // Fill Event counter and reference vertex distributions for the different trigger classes
  if(isMinBias){
    FillEventCounterHists("MB", vtx->GetZ(), isSelected, true);
  }
  if(isEMC7){
    FillEventCounterHists("EMC7", vtx->GetZ(), isSelected, fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEL0, fInputEvent));
  }
  if(isEJ2){
    FillEventCounterHists("EJ2", vtx->GetZ(), isSelected, fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ2, fInputEvent));
  }
  if(isEJ1){
    FillEventCounterHists("EJ1", vtx->GetZ(), isSelected, fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEJ1, fInputEvent));
  }
  if(isEG2){
    FillEventCounterHists("EG2", vtx->GetZ(), isSelected, fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG2, fInputEvent));
  }
  if(isEG1){
    FillEventCounterHists("EG1", vtx->GetZ(), isSelected, fTriggerSelection->IsOfflineSelected(AliEmcalTriggerOfflineSelection::kTrgEG1, fInputEvent));
  }

  PostData(1, fHistos->GetListOfHistograms());
}

void AliAnalysisTaskEventSelectionRef::FillEventCounterHists(
    const char *triggerclass,
    double vtxz,
    bool isSelected,
    bool isOfflineSelected
)
{
  // Fill reference distribution for the primary vertex before any z-cut
  fHistos->FillTH1(Form("hVertexTrigger%sBeforeEventSelection", triggerclass), vtxz);
  fHistos->FillTH1(Form("hEventCount%sBeforeEventSelection", triggerclass), 1.);
  if(isSelected){
    // Fill Event counter and reference vertex distributions after event selection
    fHistos->FillTH1(Form("hEventCount%sBeforeOfflineTrigger", triggerclass), 1);
    fHistos->FillTH1(Form("hVertexTrigger%sBeforeOfflineTrigger", triggerclass), vtxz);
    if(isOfflineSelected){
      fHistos->FillTH1(Form("hEventCount%sAfterOfflineTrigger", triggerclass), 1);
      fHistos->FillTH1(Form("hVertexTrigger%sAfterOfflineTrigger", triggerclass), vtxz);
    }

    // Now make some distributions of quantities in selected and rejected event
    //... tracks
    for(TIter trackIter = TIter(fTrackContainer).Begin(); trackIter != TIter::End(); ++trackIter){
      ProcessTrack(triggerclass, static_cast<AliVTrack *>(*trackIter), isOfflineSelected);
    }

    // ... clusters
    if(fClusterContainer){
      for(TIter clsit = TIter(fClusterContainer).Begin(); clsit != TIter::End(); ++clsit){
        ProcessCluster(triggerclass, static_cast<AliVCluster *>(*clsit), isOfflineSelected);
      }
    }

    // ... patches
    for(TIter patchiter = TIter(fTriggerPatchContainer).Begin(); patchiter != TIter::End(); ++patchiter){
      ProcessOfflinePatch(triggerclass, static_cast<AliEMCALTriggerPatchInfo *>(*patchiter), isOfflineSelected);
    }
  }
}

void AliAnalysisTaskEventSelectionRef::ProcessTrack(
    const char *triggerclass,
    const AliVTrack * track,
    bool isOfflineSelected
)
{
  fHistos->FillTH1(Form("hTrackPtAll%s", triggerclass), TMath::Abs(track->Pt()));
  if(isOfflineSelected){
    fHistos->FillTH1(Form("hTrackPtAccepted%s", triggerclass), TMath::Abs(track->Pt()));
  } else {
    fHistos->FillTH1(Form("hTrackPtRejected%s", triggerclass), TMath::Abs(track->Pt()));
  }
}

void AliAnalysisTaskEventSelectionRef::ProcessCluster(
    const char *triggerclass,
    const AliVCluster *clust,
    bool isOfflineSelected
)
{
  if(!clust->IsEMCAL()) return;
  fHistos->FillTH1(Form("hClusterEnergyAll%s", triggerclass), clust->E());
  if(isOfflineSelected){
    fHistos->FillTH1(Form("hClusterEnergyAccepted%s", triggerclass), clust->E());
  } else {
    fHistos->FillTH1(Form("hClusterEnergyRejected%s", triggerclass), clust->E());
  }
}

void AliAnalysisTaskEventSelectionRef::ProcessOfflinePatch(
    const char * triggerclass,
    const AliEMCALTriggerPatchInfo * patch,
    bool isOfflineSelected
)
{
  bool isSingleShower = patch->IsGammaLowSimple();
  bool isJetPatch = patch->IsJetLowSimple();
  if(!(isSingleShower || isJetPatch)) return;
  TString triggerstring(triggerclass);
  if(!triggerstring.CompareTo("MB")){
    if(isSingleShower){
      fHistos->FillTH1("hPatchEnergyEMC7MinBias", patch->GetPatchE());
      fHistos->FillTH1("hPatchEnergyEG1MinBias", patch->GetPatchE());
      fHistos->FillTH1("hPatchEnergyEG2MinBias", patch->GetPatchE());
    } else {
      fHistos->FillTH1("hPatchEnergyEJ1MinBias", patch->GetPatchE());
      fHistos->FillTH1("hPatchEnergyEJ2MinBias", patch->GetPatchE());
    }
  } else {
    bool singleShowerTrigger = !triggerstring.CompareTo("EMC7") || !triggerstring.CompareTo("EG1") || !triggerstring.CompareTo("EG2");
    if((isSingleShower && singleShowerTrigger) || (isJetPatch && !singleShowerTrigger)){
      fHistos->FillTH1(Form("hPatchEnergyAll%s", triggerclass), patch->GetPatchE());
      if(isOfflineSelected)
        fHistos->FillTH1(Form("hPatchEnergyAccepted%s", triggerclass), patch->GetPatchE());
      else
        fHistos->FillTH1(Form("hPatchEnergyRejected%s", triggerclass), patch->GetPatchE());
    }
  }
}

/**
 * Create new energy binning
 * @param binning
 */
void AliAnalysisTaskEventSelectionRef::CreateEnergyBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(32, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Create new Pt binning
 * @param binning
 */
void AliAnalysisTaskEventSelectionRef::CreatePtBinning(TArrayD& binning) const {
  std::vector<double> mybinning;
  std::map<double,double> definitions;
  definitions.insert(std::pair<double, double>(1, 0.05));
  definitions.insert(std::pair<double, double>(2, 0.1));
  definitions.insert(std::pair<double, double>(4, 0.2));
  definitions.insert(std::pair<double, double>(7, 0.5));
  definitions.insert(std::pair<double, double>(16, 1));
  definitions.insert(std::pair<double, double>(36, 2));
  definitions.insert(std::pair<double, double>(40, 4));
  definitions.insert(std::pair<double, double>(50, 5));
  definitions.insert(std::pair<double, double>(100, 10));
  definitions.insert(std::pair<double, double>(200, 20));
  double currentval = 0.;
  mybinning.push_back(currentval);
  for(std::map<double,double>::iterator id = definitions.begin(); id != definitions.end(); ++id){
    double limit = id->first, binwidth = id->second;
    while(currentval < limit){
      currentval += binwidth;
      mybinning.push_back(currentval);
    }
  }
  binning.Set(mybinning.size());
  int ib = 0;
  for(std::vector<double>::iterator it = mybinning.begin(); it != mybinning.end(); ++it)
    binning[ib++] = *it;
}

/**
 * Run track selection for ESD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskEventSelectionRef::TrackSelectionESD(AliESDtrack* track) {
  return fTrackCuts->AcceptTrack(track);
}

/**
 * Run track selection for AOD tracks
 * @param track The track to check
 * @return True if the track is selected, false otherwise
 */
Bool_t AliAnalysisTaskEventSelectionRef::TrackSelectionAOD(AliAODTrack* track) {
  if(!track->TestFilterBit(AliAODTrack::kTrkGlobal)) return false;
  if(track->GetTPCNCrossedRows() < 120) return false;
  return true;
}
