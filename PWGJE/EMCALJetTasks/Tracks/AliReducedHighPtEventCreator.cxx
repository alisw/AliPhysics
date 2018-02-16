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
#include <cfloat>
#include <map>
#include <vector>

#include <TArrayD.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TString.h>
#include <TTree.h>

#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliEmcalTrackSelection.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliPicoTrack.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliEMCALTriggerPatchInfo.h"

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
  fEventSelectionBits(AliVEvent::kAny),
  fSwapTriggerThresholds(kFALSE),
  fMinClusterE(-1),
  fMaxClusterE(1000),
  fMinPt(-1),
  fMaxPt(1000),
  fMinEta(-1000),
  fMaxEta(1000),
  fKeepFractionEvents(1.),
  fApplyCentralitySelection(kFALSE),
  fCentralityMethod("V0A"),
  fTriggerSetup("4classes"),
  fMinBiasSelection(AliVEvent::kINT7)
{
  fSelectCentralityRange[0] = 0.;
  fSelectCentralityRange[1] = 100.;
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
  fEventSelectionBits(AliVEvent::kAny),
  fSwapTriggerThresholds(kFALSE),
  fMinClusterE(-1),
  fMaxClusterE(1000),
  fMinPt(-1),
  fMaxPt(1000),
  fMinEta(-1000),
  fMaxEta(1000),
  fKeepFractionEvents(1.),
  fApplyCentralitySelection(kFALSE),
  fCentralityMethod("V0A"),
  fTriggerSetup("4classes"),
  fMinBiasSelection(AliVEvent::kINT7)
{
  DefineOutput(2, TTree::Class());

  fSelectCentralityRange[0] = 0.;
  fSelectCentralityRange[1] = 100.;

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

  PostData(1, fOutput);
  PostData(2, fOutputTree);
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
  AliDebug(1, "Starting creation of the reconstructed event");
  if(!fCaloClusters){
    AliError("Cluster container missing");
    return kFALSE;
  }
  if(!fTracks){
    AliError("Track container missing");
    return kFALSE;
  }
  if(!(this->fInputHandler->IsEventSelected() & fEventSelectionBits)){
    AliDebug(1, "Trigger not selected");
    return kFALSE;
  } else {
    AliDebug(1, "Trigger selected");
  }
  if(!SelectEvent(fInputEvent)) return kFALSE;
  new(fOutputEvent) AliReducedHighPtEvent(kTRUE);

  // Write event specific information
  fOutputEvent->SetVertexZ(fInputEvent->GetPrimaryVertex()->GetZ());
  AliCentrality *centralityHandler = fInputEvent->GetCentrality();
  if(centralityHandler) fOutputEvent->SetCentralityPercentile(centralityHandler->GetCentralityPercentile(fCentralityMethod.Data()));
  ConvertTriggerPatches(fTriggerPatchInfo, fOutputEvent->GetPatchContainer());
  TString triggerString(fInputEvent->GetFiredTriggerClasses());
  if(!fTriggerSetup.CompareTo("4classes"))
    fOutputEvent->SetDecisionFromTriggerString(triggerString.Contains("EG1"), triggerString.Contains("EG2"), triggerString.Contains("EJ1"), triggerString.Contains("EJ2"));
  else
    fOutputEvent->SetDecisionFromTriggerString(triggerString.Contains("EGA"), kFALSE, triggerString.Contains("EJE"), kFALSE);
  fOutputEvent->SetMinBiasEvent(fInputHandler->IsEventSelected() & fMinBiasSelection);
  fOutputEvent->SetRunNumber(fInputEvent->GetRunNumber());

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
      if(pt < fMinPt || pt > fMaxPt) continue;
      if(eta < fMinEta || eta > fMaxEta) continue;
      if(!part->Charge()) continue;
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
    // Get leading 3 cell energies
    TArrayD cellEnergies;
    GetCellEnergies(incluster, cellEnergies);
    TArrayI indices(cellEnergies.GetSize());
    TMath::Sort(cellEnergies.GetSize(), cellEnergies.GetArray(), indices.GetArray(), kTRUE);
    redcluster->SetLeadingCellEnergies(
        cellEnergies[indices[0]],
        cellEnergies.GetSize() > 1 ? cellEnergies[indices[1]] : 0,
        cellEnergies.GetSize() > 2 ? cellEnergies[indices[2]] : 0
    );
    // Assing MC particles
    if(MCEvent()){
      for(Int_t ilab = 0; ilab < incluster->GetNLabels(); ilab++){
        AliVParticle *assigned = MCEvent()->GetTrack(TMath::Abs(incluster->GetLabels()[ilab]));
        if(!assigned) continue;
        redcluster->AddTrueContributor(assigned->PdgCode(), assigned->Px(), assigned->Py(), assigned->Pz(), assigned->E());
      }
    }
    fOutputEvent->AddReducedCluster(redcluster);
    clusterindexmap.insert(std::pair<int,int>(incluster->GetID(), ncluster));
    ncluster++;
  }

  // Filter tracks
  int nsel(0);
  AliDebug(1, Form("Number of tracks in container: %d", fTracks->GetEntries()));
  for(TIter trkiter = TIter(fTracks).Begin(); trkiter != TIter::End(); ++trkiter){
    AliVTrack *trackRaw = static_cast<AliVTrack *>(*trkiter), *trackToCheck(NULL);
    // handlinf for pico tracks as input
    if(trackRaw->IsA() == AliPicoTrack::Class()){
      AliPicoTrack *picotrack = static_cast<AliPicoTrack *>(trackRaw);
      trackToCheck = picotrack->GetTrack();
    } else {
      trackToCheck = trackRaw;
    }
    FixTrackInputEvent(trackToCheck);
    TArrayI cutIndices;
    if(!SelectTrack(trackToCheck,cutIndices)) continue;
    AliReducedReconstructedTrack *rectrack = new AliReducedReconstructedTrack;
    double pvec[3];
    trackToCheck->GetPxPyPz(pvec);
    rectrack->SetMomentumVector(pvec[0], pvec[1], pvec[2]);
    rectrack->SetCharge(static_cast<Char_t>(trackToCheck->Charge()));
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
        rectrack->SetTPCClusters(trackToCheck->GetTPCNcls());
        rectrack->SetTPCCrossedRows(GetTPCCrossedRows(trackToCheck));
        rectrack->SetTPCSharedClusters(trackToCheck->GetTPCSharedMapPtr()->CountBits());
        rectrack->SetTPCFindableClusters(trackToCheck->GetTPCNclsF());
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
    //std::cout << "Track selected" << std::endl;
    fOutputEvent->AddReducedReconstructedParticle(rectrack);
    nsel++;
  }
  AliDebug(1, Form("Number of selected tracks :%d ", nsel ));

  fOutputTree->Fill();
  PostData(2, fOutputTree);
  return kTRUE;
}

/**
 * Add new virtual track selection to the list of track selections
 * \param sel The track selection
 * \param cutindex The index of the track selection
 */
void AliReducedHighPtEventCreator::AddVirtualTrackSelection(AliEmcalTrackSelection* sel, int cutindex) {
  fTrackSelections->Add(new AliReducedTrackSelectionContainer(cutindex, sel));
}

/**
 * Apply standard event selection. Includes also downscaling if requested.
 * Note: Downscaling always applied after full event selection
 * \param event The event to check
 * \return True if the event was selected, false otherwise
 */
Bool_t AliReducedHighPtEventCreator::SelectEvent(AliVEvent* event) const {
  if(fApplyCentralitySelection && InputEvent()->GetCentrality()){
    AliDebug(1, Form("Centrality selection applied in range %f - %f\n", fSelectCentralityRange[0], fSelectCentralityRange[1]));
    Float_t centrality = InputEvent()->GetCentrality()->GetCentralityPercentile(fCentralityMethod.Data());
    AliDebug(1, Form("Centrality before: %f, estimator %s\n", centrality, fCentralityMethod.Data()));
    if(centrality < fSelectCentralityRange[0] || centrality > fSelectCentralityRange[1]) return kFALSE;
    AliDebug(1, "Centrality selected");
  }
  const AliVVertex *primvtx = event->GetPrimaryVertex();
  if(!primvtx || !primvtx->GetNContributors()) return kFALSE;
  if(TMath::Abs(primvtx->GetZ()) > 10.) return kFALSE;
  if(event->IsPileupFromSPD(3, 0.8, 3., 2., 5.)) return kFALSE;
  AliAnalysisUtils eventSelUtil;
  if(!eventSelUtil.IsVertexSelected2013pA(event)) return kFALSE;
  if(fKeepFractionEvents < (1. - FLT_EPSILON)){
    // Apply downscaling
    Float_t choice = gRandom->Uniform(0., 1.);
    if(choice > fKeepFractionEvents) return kFALSE;     // event rejected by downscaling
  }
  return kTRUE;
}

/**
 * Select cluster based on energy and whether it is in EMCAL
 * \param clust
 * \return True if the cluster was selected, false otherwise
 */
Bool_t AliReducedHighPtEventCreator::SelectCluster(const AliVCluster* clust) const {
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
 * Access number of TPC crossed rows in a transparent way both for
 * \param trk Track to check
 * \return Number of TPC crossed rows
 */
Int_t AliReducedHighPtEventCreator::GetTPCCrossedRows(const AliVTrack* trk) const {
  if(trk->IsA() == AliESDtrack::Class()){
    return (static_cast<const AliESDtrack *>(trk))->GetTPCCrossedRows();
  } else if(trk->IsA() == AliAODTrack::Class()){
    return (static_cast<const AliAODTrack *>(trk))->GetTPCNCrossedRows();
  }
  return 0;
}

/**
 * Get the cluster cell energies
 * \param emccluster EMCAL cluster to check
 * \param energies Array storing the cell energies
 */
void AliReducedHighPtEventCreator::GetCellEnergies(AliVCluster* emccluster, TArrayD& energies) const {
  if(!fInputEvent->GetEMCALCells()) {
    AliError("No EMCAL cells");
    return;
  }
  AliDebug(2, Form("Number of cells: %d, array: %p", emccluster->GetNCells(), emccluster->GetCellsAbsId()));
  energies.Set(emccluster->GetNCells());
  for(int icell = 0; icell < emccluster->GetNCells(); icell++){
    // printf("Cell ID: %d\n", emccluster->GetCellsAbsId()[icell]);
    energies[icell] = fInputEvent->GetEMCALCells()->GetCellAmplitude(emccluster->GetCellsAbsId()[icell]);
  }
}

/**
 * Convert trigger patches to the reduced format
 * \param patches Input patches
 * \param cont Output container
 */
void AliReducedHighPtEventCreator::ConvertTriggerPatches(TClonesArray* patches,
    AliReducedPatchContainer* cont) {
  if(!patches){
    AliError("Trigger patch container not found\n");
    return;
  }
  for(TIter patchIter = TIter(patches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEMCALTriggerPatchInfo *mypatch = static_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(!mypatch->IsOfflineSimple() && mypatch->IsLevel0()) continue;
    AliReducedPatchContainer::PatchType_t triggertype;
    Bool_t isDefined(kFALSE);
    if(mypatch->IsOfflineSimple()){
      if(mypatch->IsGammaHighSimple()) { triggertype = AliReducedPatchContainer::kEMCGammaHigh; isDefined = kTRUE; }
      if(mypatch->IsGammaLowSimple()) { triggertype = AliReducedPatchContainer::kEMCGammaLow; isDefined = kTRUE; }
      if(mypatch->IsJetHighSimple()) { triggertype = AliReducedPatchContainer::kEMCJetHigh; isDefined = kTRUE; }
      if(mypatch->IsJetLowSimple()) { triggertype = AliReducedPatchContainer::kEMCJetLow; isDefined = kTRUE; }
      if(!isDefined){
        AliDebug(2, "Unknown offline patch type");
        continue;
      }
      AliDebug(2, Form("Adding offline patch of type %d", int(triggertype)));
      cont->AddTriggerPatch(kTRUE, triggertype, mypatch->GetPatchE(), mypatch->GetADCAmp(), mypatch->GetEtaGeo(), mypatch->GetPhiGeo());
    } else {
      if(mypatch->IsGammaHigh()) { triggertype = fSwapTriggerThresholds ? AliReducedPatchContainer::kEMCGammaLow : AliReducedPatchContainer::kEMCGammaHigh; isDefined=kTRUE; }
      if(mypatch->IsGammaLow()) { triggertype = fSwapTriggerThresholds ? AliReducedPatchContainer::kEMCGammaHigh : AliReducedPatchContainer::kEMCGammaLow; isDefined=kTRUE; }
      if(mypatch->IsJetHigh()) { triggertype = fSwapTriggerThresholds ? AliReducedPatchContainer::kEMCJetLow : AliReducedPatchContainer::kEMCJetHigh; isDefined=kTRUE; }
      if(mypatch->IsJetLow()) { triggertype = fSwapTriggerThresholds ? AliReducedPatchContainer::kEMCJetHigh : AliReducedPatchContainer::kEMCJetLow; isDefined=kTRUE; }
      if(!isDefined){
        AliDebug(2, "Unknown online patch type");
        continue;
      }
      AliDebug(2, Form("Adding online patch of type %d", int(triggertype)));
      cont->AddTriggerPatch(kFALSE, triggertype, mypatch->GetPatchE(), mypatch->GetADCAmp(), mypatch->GetEtaGeo(), mypatch->GetPhiGeo());
    }
  }
}

/**
 * Set the corresponding pointer to the original event to the track in a transparent way for
 * ESD, AOD and pico tracks
 * \param trk The track to handle
 */
void AliReducedHighPtEventCreator::FixTrackInputEvent(AliVTrack* trk) {
  if(!trk->GetEvent()){
    if(trk->IsA() == AliESDtrack::Class())
      (static_cast<AliESDtrack *>(trk))->SetESDEvent(static_cast<AliESDEvent *>(fInputEvent));
    else if(trk->IsA() == AliAODTrack::Class())
      (static_cast<AliAODTrack *>(trk))->SetAODEvent(static_cast<AliAODEvent *>(fInputEvent));
    else if(trk->IsA() == AliPicoTrack::Class()){
      AliPicoTrack *mytrk = static_cast<AliPicoTrack *>(trk);
      if(!mytrk->GetEvent()){
        if(mytrk->GetTrack()->IsA() == AliESDtrack::Class())
          (static_cast<AliESDtrack *>(mytrk->GetTrack()))->SetESDEvent(static_cast<AliESDEvent *>(fInputEvent));
        else
          (static_cast<AliAODTrack *>(mytrk->GetTrack()))->SetAODEvent(static_cast<AliAODEvent *>(fInputEvent));
      }
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
AliReducedTrackSelectionContainer::AliReducedTrackSelectionContainer(Int_t index, AliEmcalTrackSelection* sel):
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

