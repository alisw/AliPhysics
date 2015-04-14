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
#include <map>
#include <vector>

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TString.h>
#include <TTree.h>

#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliEmcalTriggerPatchInfo.h"

#include "AliEMCalPtTaskVTrackSelection.h"
#include "AliReducedEmcalCluster.h"
#include "AliReducedGeneratedParticle.h"
#include "AliReducedHighPtEvent.h"
#include "AliReducedMCHeader.h"
#include "AliReducedPatchContainer.h"
#include "AliReducedPatchInfo.h"
#include "AliReducedReconstructedTrack.h"
#include "AliReducedHighPtEventCreator.h"

/// \cond CLASSIMP
ClassImp(HighPtTracks::AliReducedTrackSelectionContainer)
ClassImp(HighPtTracks::AliReducedHighPtEventCreator)
/// \endcond

namespace HighPtTracks {

/**
 * Dummy constructor, initializing default values
 */
AliReducedHighPtEventCreator::AliReducedHighPtEventCreator():
  AliAnalysisTaskEmcal(),
  fOutputTree(NULL),
  fOutputEvent(NULL),
  fTrackSelections(NULL),
  fSwapTriggerThresholds(kFALSE),
  fMinClusterE(-1),
  fMaxClusterE(1000),
  fMinPt(-1),
  fMaxPt(1000),
  fMinEta(-1000),
  fMaxEta(1000)
{

}

/**
 * Named constructor, linking also output container and allocated internal containers
 * \param name Name of the task
 */
AliReducedHighPtEventCreator::AliReducedHighPtEventCreator(const char* name):
  AliAnalysisTaskEmcal(name, kTRUE),
  fOutputTree(NULL),
  fOutputEvent(NULL),
  fTrackSelections(NULL),
  fSwapTriggerThresholds(kFALSE),
  fMinClusterE(-1),
  fMaxClusterE(1000),
  fMinPt(-1),
  fMaxPt(1000),
  fMinEta(-1000),
  fMaxEta(1000)
{
  DefineOutput(2, TTree::Class());

  fTrackSelections = new TObjArray;
  fTrackSelections->SetOwner(kTRUE);
}

/**
 * Destructor, deleting cut list
 */
AliReducedHighPtEventCreator::~AliReducedHighPtEventCreator() {
  if(fTrackSelections) delete fTrackSelections;
}

/**
 * Creates the output tree and the output event
 */
void AliReducedHighPtEventCreator::UserCreateOutputObjects() {
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  OpenFile(2);

  fOutputTree = new TTree("ReducedEvent", "A simple reduced event");
  fOutputEvent = new AliReducedHighPtEvent();
  fOutputTree->Branch("ReducedEvent", "AliReducedHighPtEvent", fOutputEvent, 128000, 0);

  PostData(1, fOutputTree);
}

/**
 * Run filter:
 * -# Do event selection
 * -# Create output event
 * -# Fill basic event information
 * -# If Monte-Carlo is available, add MC header and filter true particles
 * -# Filter clusters
 * -# Filter reconstructed tracks
 * \return True if the event was selected, false otherwise
 */
Bool_t AliReducedHighPtEventCreator::Run() {
  if(!fCaloClusters){
    AliError("Cluster container missing");
    return kFALSE;
  }
  if(!fTracks){
    AliError("Track container missing");
    return kFALSE;
  }
  if(!SelectEvent(fInputEvent)) return kFALSE;
  new(fOutputEvent) AliReducedHighPtEvent(kTRUE);

  // Write event specific information
  fOutputEvent->SetVertexZ(fInputEvent->GetPrimaryVertex()->GetZ());
  AliCentrality *centralityHandler = fInputEvent->GetCentrality();
  if(centralityHandler) fOutputEvent->SetCentralityPercentile(centralityHandler->GetCentralityPercentile("V0A"));
  ConvertTriggerPatches(fTriggerPatchInfo, fOutputEvent->GetPatchContainer());
  TString triggerString(fInputEvent->GetFiredTriggerClasses());
  fOutputEvent->SetDecisionFromTriggerString(triggerString.Contains("EG1"), triggerString.Contains("EG2"), triggerString.Contains("EJ1"), triggerString.Contains("EJ2"));
  fOutputEvent->SetMinBiasEvent(fInputHandler->IsEventSelected() & AliVEvent::kINT7);

  std::map<int, int> mcindexmap;
  if(MCEvent()){
    // Generate Monte-Carlo header
    AliReducedMCHeader *mcheader = new AliReducedMCHeader();
    if(fIsPythia){
      mcheader->SetCrossSection(fXsection);
      mcheader->SetNumberOfTrials(fNTrials);
      mcheader->SetPtHard(fPtHard);
    }
    fOutputEvent->SetMonteCarloHeader(mcheader);

    // Convert Monte Carlo particles, fill index map
    AliMCEvent *mcev = MCEvent();
    int npart(0);
    for(Int_t ipart = 0; ipart < mcev->GetNumberOfTracks(); ipart++){
      AliVParticle *part = mcev->GetTrack(ipart);
      Double_t pt(TMath::Abs(part->Pt())), eta(part->Eta());
      if(pt < fMinPt || pt > fMaxPt) return 0;
      if(eta < fMinEta || eta > fMaxEta) return 0;
      if(part->IsA() == AliAODMCParticle::Class()){
        AliAODMCParticle *aodpart = static_cast<AliAODMCParticle *>(part);
        if(!aodpart->IsPhysicalPrimary()) continue;
      } else {
        // ESD - get information from stack
        if(!mcev->Stack()->IsPhysicalPrimary(ipart)) continue;
      }

      AliReducedGeneratedParticle *reducedgen = new AliReducedGeneratedParticle(npart, part->PdgCode(), part->Px(), part->Py(), part->Pz(), part->E());
      fOutputEvent->AddReducedGeneratedParticle(reducedgen);
      mcindexmap.insert(std::pair<int,int>(ipart, npart));
      npart++;
    }
  }

  // Filter Clusters
  std::map<int, int> clusterindexmap;
  int ncluster(0);
  Double_t vtxpos[3];fInputEvent->GetPrimaryVertex()->GetXYZ(vtxpos);
  for(TIter cliter = TIter(fCaloClusters).Begin(); cliter != TIter::End(); ++cliter){
    AliVCluster *incluster = static_cast<AliVCluster *>(*cliter);
    if(!SelectCluster(incluster)) continue;
    TLorentzVector clustervec;
    incluster->GetMomentum(clustervec, vtxpos);
    AliReducedEmcalCluster *redcluster = new AliReducedEmcalCluster(ncluster, incluster->E(), clustervec.Eta(), clustervec.Phi(), incluster->GetM02(), incluster->GetM20());
    fOutputEvent->AddReducedCluster(redcluster);
    clusterindexmap.insert(std::pair<int,int>(incluster->GetID(), ncluster));
    ncluster++;
  }

  // Filter tracks
  for(TIter trkiter = TIter(fTracks).Begin(); trkiter != TIter::End(); ++trkiter){
    AliVTrack *trackToCheck = static_cast<AliVTrack *>(*trkiter);
    TArrayI cutIndices;
    if(!SelectTrack(trackToCheck,cutIndices)) continue;
    AliReducedReconstructedTrack *rectrack = new AliReducedReconstructedTrack;
    double pvec[3];
    trackToCheck->GetPxPyPz(pvec);
    rectrack->SetMomentumVector(pvec[0], pvec[1], pvec[2]);
    for(Int_t icut = 0; icut < cutIndices.GetSize(); icut++) rectrack->SetTrackCuts(cutIndices[icut]);
    if(MCEvent()){
      // handle label
      Int_t label = TMath::Abs(trackToCheck->GetLabel());
      if(label < fMCEvent->GetNumberOfTracks()){
        std::map<int,int>::iterator found = mcindexmap.find(label);
        if(found != mcindexmap.end()){
          rectrack->SetMatchedParticleIndex(found->second);
        }
        if(label >= 0) rectrack->SetGoodTrackLabel(kTRUE);
      }
    }
    // handle clusters
    Int_t clusterIndex = trackToCheck->GetEMCALcluster();
    if(clusterIndex >= 0 && clusterIndex < fCaloClusters->GetEntries()){
      std::map<int,int>::iterator found = clusterindexmap.find(clusterIndex);
      if(found != clusterindexmap.end()){
        rectrack->SetMatchedClusterIndex(found->second);
      }
    }
  }

  fOutputTree->Fill();
  PostData(1, fOutputTree);
  return kTRUE;
}

/**
 * Add new virtual track selection to the list of track selections
 * \param sel The track selection
 * \param cutindex The index of the track selection
 */
void AliReducedHighPtEventCreator::AddVirtualTrackSelection(
    EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection* sel, int cutindex) {
  fTrackSelections->Add(new AliReducedTrackSelectionContainer(cutindex, sel));
}

/**
 * Apply standard event selection
 * \param event The event to check
 * \return True if the event was selected, false otherwise
 */
Bool_t AliReducedHighPtEventCreator::SelectEvent(AliVEvent* event) const {
  const AliVVertex *primvtx = event->GetPrimaryVertex();
  if(!primvtx || !primvtx->GetNContributors()) return kFALSE;
  if(TMath::Abs(primvtx->GetZ()) > 10.) return kFALSE;
  if(event->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return kFALSE;
  AliAnalysisUtils eventSelUtil;
  if(!eventSelUtil.IsVertexSelected2013pA(event)) return kTRUE;
  return kTRUE;
}

/**
 * Select cluster based on energy and whether it is in EMCAL
 * \param clust
 * \return True if the cluster was selected, false otherwise
 */
Bool_t AliReducedHighPtEventCreator::SelectCluster(AliVCluster* clust) const {
    if(!clust->IsEMCAL()) return kFALSE;
    if(clust->E() < fMinClusterE || clust->E() > fMaxClusterE) return kFALSE;
    return kTRUE;
}

/**
 * Run track selection: First apply kinematics selection, then test different cut combinations
 * \param track The track to check
 * \param cutindices Indicies of cut combinations by which the track was selected
 * \return Number of successfull cut combinations
 */
Int_t AliReducedHighPtEventCreator::SelectTrack(AliVTrack* track, TArrayI& cutindices) const {
  std::vector<int> selected;

  double pt = TMath::Abs(track->Pt()), eta = track->Eta();
  if(pt < fMinPt || pt > fMaxPt) return 0;
  if(eta < fMinEta || eta > fMaxEta) return 0;

  for(TIter cutiter = TIter(fTrackSelections).Begin(); cutiter != TIter::End(); ++cutiter){
    AliReducedTrackSelectionContainer *mycut = static_cast<AliReducedTrackSelectionContainer *>(*cutiter);
    if(!mycut->GetTrackSelection()->IsTrackAccepted(track)) continue;
    selected.push_back(mycut->GetIndex());
  }
  if(!selected.size()) return 0;
  cutindices.Set(selected.size());
  int counter = 0;
  for(std::vector<int>::iterator inditer = selected.begin(); inditer != selected.end(); ++inditer){
    cutindices[counter++] = *inditer;
  }
  return selected.size();
}

/**
 * Convert trigger patches to the reduced format
 * \param patches Input patches
 * \param cont Output container
 */
void AliReducedHighPtEventCreator::ConvertTriggerPatches(TClonesArray* patches,
    AliReducedPatchContainer* cont) {
  for(TIter patchIter = TIter(patches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEmcalTriggerPatchInfo *mypatch = static_cast<AliEmcalTriggerPatchInfo *>(*patchIter);
    if(!mypatch->IsOfflineSimple() && mypatch->IsLevel0()) continue;
    AliReducedPatchContainer::PatchType_t triggertype;
    if(mypatch->IsOfflineSimple()){
      if(mypatch->IsGammaHighSimple()) triggertype = AliReducedPatchContainer::kEMCGammaHigh;
      if(mypatch->IsGammaLowSimple()) triggertype = AliReducedPatchContainer::kEMCGammaLow;
      if(mypatch->IsJetHighSimple()) triggertype = AliReducedPatchContainer::kEMCJetHigh;
      if(mypatch->IsJetLowSimple()) triggertype = AliReducedPatchContainer::kEMCJetLow;
      cont->AddTriggerPatch(kTRUE, triggertype, mypatch->GetPatchE(), mypatch->GetADCAmp(), mypatch->GetEtaGeo(), mypatch->GetPhiGeo());
    } else {
      if(mypatch->IsGammaHigh()) triggertype = fSwapTriggerThresholds ? AliReducedPatchContainer::kEMCGammaLow : AliReducedPatchContainer::kEMCGammaHigh;
      if(mypatch->IsGammaLow()) triggertype = fSwapTriggerThresholds ? AliReducedPatchContainer::kEMCGammaHigh : AliReducedPatchContainer::kEMCGammaLow;
      if(mypatch->IsJetHigh()) triggertype = fSwapTriggerThresholds ? AliReducedPatchContainer::kEMCJetLow : AliReducedPatchContainer::kEMCJetHigh;
      if(mypatch->IsJetLow()) triggertype = fSwapTriggerThresholds ? AliReducedPatchContainer::kEMCJetHigh : AliReducedPatchContainer::kEMCJetLow;
      cont->AddTriggerPatch(kTRUE, triggertype, mypatch->GetPatchE(), mypatch->GetADCOfflineAmp(), mypatch->GetEtaGeo(), mypatch->GetPhiGeo());
    }
  }
}

/**
 * Dummy constructor
 */
AliReducedTrackSelectionContainer::AliReducedTrackSelectionContainer():
  fIndex(-1),
  fTrackSelection(NULL)
{
}

/**
 * Constuctor fully defining the object
 * \param index Cut index
 * \param sel Track selection object
 */
AliReducedTrackSelectionContainer::AliReducedTrackSelectionContainer(
    Int_t index, EMCalTriggerPtAnalysis::AliEMCalPtTaskVTrackSelection* sel):
  fIndex(index),
  fTrackSelection(sel)
{
}

/**
 * Destructor, cleaning memory
 */
AliReducedTrackSelectionContainer::~AliReducedTrackSelectionContainer() {
  if(fTrackSelection) delete fTrackSelection;
}


} /* namespace HighPtTracks */

