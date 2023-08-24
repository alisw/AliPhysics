/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
#include <string>
#include <vector>

#include <TAxis.h>
#include <TClonesArray.h>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>

#include "AliAODMCParticle.h"
#include "AliLog.h"
#include "AliEmcalTrackSelection.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliMCEvent.h"
#include "AliPicoTrack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliEMCalTriggerAnaClassManager.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerKineCuts.h"
#include "AliEMCalTriggerRecTrackAnalysisComponent.h"
#include "AliEMCalTriggerWeightHandler.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerRecTrackAnalysisComponent)

using namespace PWGJE::EMCALJetTasks;

/**
 * \brief Dummy constructor
 *
 * Dummy constructor. For ROOT I/O. Not intended to be used by the users.
 */
AliEMCalTriggerRecTrackAnalysisComponent::AliEMCalTriggerRecTrackAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent(),
  fTrackSelection(NULL),
  fSwapEta(kFALSE),
  fRequestMCtrue(kFALSE),
  fDoMatchPatches(kFALSE)
{
}

/**
 * \brief Main constructor
 *
 * Main constructor, initialises component with a component name. To be used by the users to
 * create the component object.
 */
AliEMCalTriggerRecTrackAnalysisComponent::AliEMCalTriggerRecTrackAnalysisComponent(const char *name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fTrackSelection(NULL),
  fSwapEta(kFALSE),
  fRequestMCtrue(kFALSE),
  fDoMatchPatches(kFALSE)
{
}

/**
 * \brief Destructor
 *
 * Destructor, taking care of the track selection.
 */
AliEMCalTriggerRecTrackAnalysisComponent::~AliEMCalTriggerRecTrackAnalysisComponent() {
  if(fTrackSelection) delete fTrackSelection;
}

/**
 * \brief Creating histograms for the analysis component
 *
 * Create histograms of the track analysis component and add it to the list of
 * histograms. For each trigger class we have
 * - tracks with esd information
 * - tracks with MC information
 * - tracks with clusters and esd information
 * - tracks with clusters and MC information
 * In addition, a correlation matrix (THnSparse) for the correlation beween generated
 * and reconstructed \f$ p_{t} \f$ is created.
 *
 * This function is the implementation of the abstract method CreateHistos
 * declared in AliEMCalTriggerTracksAnalysisComponent.
 */
void AliEMCalTriggerRecTrackAnalysisComponent::CreateHistos() {
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  // Create trigger definitions
  std::map<std::string, std::string> triggerCombinations;
  GetAllTriggerNamesAndTitles(triggerCombinations);

  // Create axis definitions
  const TBinning *ptbinning = fBinning->GetBinning("pt"),
      *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi"),
      *vertexbinning = fBinning->GetBinning("zvertex");

  const TAxis *trackaxes[5] = {
      DefineAxis("pt", *ptbinning),
      DefineAxis("eta", *etabinning),
      DefineAxis("phi", *phibinning),
      DefineAxis("zvertex", *vertexbinning),
      DefineAxis("mbtrigger", TLinearBinning(2, -0.5, 1.5))
  };

  // Build histograms
  for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
    const std::string name = it->first, &title = it->second;
    fHistos->CreateTHnSparse(Form("hTrackHist%s", name.c_str()), Form("Track-based data for %s events", title.c_str()), 5, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hTrackInAcceptanceHist%s", name.c_str()), Form("Track-based data for %s events for tracks matched to EMCal clusters", title.c_str()), 5, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hMCTrackHist%s", name.c_str()), Form("Track-based data for %s events with MC kinematics", title.c_str()), 5, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hMCTrackInAcceptanceHist%s", name.c_str()), Form("Track-based data for %s events for tracks matched to EMCal clusters with MC kinematics", title.c_str()), 5, trackaxes, "s");
    if(fDoMatchPatches){
      fHistos->CreateTHnSparse(Form("hTrackHistPatchMatch%s", name.c_str()), Form("Track-based data for %s events for tracks matched to EMCal patches", title.c_str()), 5, trackaxes, "s");
      fHistos->CreateTHnSparse(Form("hMCTrackHistPatchMatch%s", name.c_str()), Form("Track-based data for %s events for tracks matched to EMCal patches with MC kinematics", title.c_str()), 5, trackaxes, "s");
    }
  }

  // Correlation Matrix
  // 1. gen pt
  // 2. rec pt
  // 3. rec eta
  // 4. rec phi
  const TAxis *corraxes[4] = {
	  DefineAxis("ptgen", *ptbinning),
	  DefineAxis("ptrec", *ptbinning),
	  DefineAxis("eta", *etabinning),
	  DefineAxis("phi", *phibinning)
  };
  fHistos->CreateTHnSparse("hTrackPtCorrelation", "Correlation matrix for track pt", 4, corraxes);

  for(int iaxis = 0; iaxis < 5; iaxis++) delete trackaxes[iaxis];
}

/**
 * \brief Run track loop on list of matching tracks
 *
 * Analyses all particles at reconstruction level in an event, and fill a track based
 * histogram for all triggers which selected the event. Following steps are performed:
 *  -# Get a list of all triggers which selected the event
 *  -# Iterate over pre-selected tracks
 *      - Check kinematics and track selection cuts
 *      - Check whether the track is a true signal track (if MC information is available). If yes,
 *        fill correlation matrix
 *      - Fill histogram at track level (and if an associated particle is available also at MC-truth level)
 *      - Find cluster matched to track. If cluster is available, fill histograms for tracks with clusters
 *      - Find trigger patches mathced to track. If trigger patch is available, fill histogram for tracks with patches
 * This function is the implementation of the abstract method Process declared
 * in AliEMCalTriggerTracksAnalysisComponent.
 *
 * \param data the event data
 */
void AliEMCalTriggerRecTrackAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  AliDebug(1, Form("Number of matched tracks: %d", data->GetMatchedTrackContainer()->GetEntries()));
  if(fRequestMCtrue && !data->GetMCEvent()) return;

  std::vector<std::string> triggernames;
  this->GetMachingTriggerNames(triggernames);

  AliVTrack *track(NULL);
  const AliVParticle *assocMC(NULL);
  if(!data->GetMatchedTrackContainer()){
    AliError("No container for matched tracks");
    return;
  }
  TIter trackIter(data->GetMatchedTrackContainer());

  double weight = 1.;
  if(fWeightHandler && data->GetMCEvent()){
    weight = fWeightHandler->GetEventWeight(data->GetMCEvent());
  }
  while((track = dynamic_cast<AliVTrack *>(trackIter()))){
    // Apply track selection
    assocMC = NULL;
    if(fKineCuts && !fKineCuts->IsSelected(track)) continue;
    if(fTrackSelection && !fTrackSelection->IsTrackAccepted(track)){
      AliDebug(2, "Track not accepted");
      continue;
    }

    if(fRequestMCtrue){
    	if(!(assocMC = IsMCTrueTrack(track, data->GetMCEvent()))) continue;	// Not a true track
    	this->FillCorrelation(assocMC, track);
    }
    // Try to match the cluster
    Bool_t hasCluster = kFALSE;
    AliVCluster *clust(NULL);
    AliVTrack *testtrack = track;
    AliPicoTrack *pictrack = dynamic_cast<AliPicoTrack *>(track);
    if(pictrack) testtrack = pictrack->GetTrack();
    if(testtrack->GetEMCALcluster() >= 0 && (clust = dynamic_cast<AliVCluster *>(data->GetClusterContainer()->At(testtrack->GetEMCALcluster()))))
      hasCluster = kTRUE;

    // Try to match to EMCAL patches (only in case the event is triggered)
    TList patches;
    if(fDoMatchPatches && triggernames.size()) MatchTriggerPatches(track, data->GetTriggerPatchContainer(), patches);

    // Fill histograms
    for(std::vector<std::string>::iterator name = triggernames.begin(); name != triggernames.end(); name++){
      FillHistogram(Form("hTrackHist%s", name->c_str()), track, NULL, data->GetRecEvent(), kFALSE, weight);
      if(hasCluster) FillHistogram(Form("hTrackInAcceptanceHist%s", name->c_str()), track, NULL, data->GetRecEvent(), kFALSE, weight);
      Bool_t hasPatch = HasMatchedPatchOfType(*name, patches);
      if(fDoMatchPatches && hasPatch) FillHistogram(Form("hTrackHistPatchMatch%s", name->c_str()), track, NULL, data->GetRecEvent(), kFALSE, weight);
      if(assocMC){
        FillHistogram(Form("hMCTrackHist%s", name->c_str()), track, assocMC, data->GetRecEvent(), kTRUE, weight);
        if(hasCluster) FillHistogram(Form("hMCTrackInAcceptanceHist%s", name->c_str()), track, assocMC, data->GetRecEvent(), kTRUE, weight);
        if(fDoMatchPatches && hasPatch) FillHistogram(Form("hMCTrackHistPatchMatch%s", name->c_str()), track, NULL, data->GetRecEvent(), kFALSE, weight);
      }
    }
  }
}

/**
 * \brief Check whether track is a true signal track
 *
 * Check according to the associated MC information whether the track is a MC true track,
 * and whether it is physical primary
 *
 * \param trk track to check
 * \param evnt MC event information necessary for the check
 * \return the associated MC particle (NULL if not MC true)
 */
const AliVParticle * AliEMCalTriggerRecTrackAnalysisComponent::IsMCTrueTrack(
    const AliVTrack* const trk, const AliMCEvent* evnt) const {
  int label = TMath::Abs(trk->GetLabel());
  const AliVParticle *mcpart = evnt->GetTrack(label);
  if(!mcpart) return NULL;
  const AliAODMCParticle *aodpart = dynamic_cast<const AliAODMCParticle *>(mcpart);
  if(aodpart){
    if(!aodpart->IsPhysicalPrimary()) return NULL;
  } else {
    if(!evnt->IsPhysicalPrimary(mcpart->GetLabel())) return NULL;
  }
  return mcpart;
}

/**
 * \brief Fill main track-based histogram
 *
 * Fill main track-based histogram defined by its name with
 *  -# \f$ p_{t} \f$
 *  -# \f$ \eta \f$
 *  -# \f$ \phi \f$
 *  -# z-position of the primary vertex
 *  -# status flag signalizing that the event was also a minimum bias event
 * If useMCkine is set to true, the kinematic quantities will be obtained from the
 * associated MC particle.
 *
 * \param histname Name of the THnSparse to fill
 * \param trk Reconstructed track
 * \param assocMC The associated MC track
 * \param recev Reconstructed event
 * \param useMCkine If true we fill histogram with MC truth information
 * \param weight Event weight (optional)
 */
void AliEMCalTriggerRecTrackAnalysisComponent::FillHistogram(
    const TString& histname, const AliVTrack* const trk,
    const AliVParticle* assocMC, const AliVEvent* const recev,
    Bool_t useMCkine, Double_t weight) {
  if(useMCkine && !assocMC) return;
  double data[5];
  data[0] = useMCkine ? TMath::Abs(assocMC->Pt()) : TMath::Abs(trk->Pt());
  data[1] = (fSwapEta ? -1. : 1.) * (useMCkine ? assocMC->Eta() : trk->Eta());
  data[2] = useMCkine ? assocMC->Phi() : trk->Phi();
  data[3] = recev->GetPrimaryVertex()->GetZ();
  data[4] = fTriggerClassManager->HasMinBiasTrigger();
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}

/**
 * \brief Fill correlation matrix between generated and reconstructed particle information
 *
 * Fills the correlation matrix between \f$ p_{t} \f$ at generation and at reconstruction level
 * using information from the reconstructed track and the associated generated track. In
 * addition the postion of the track in \f$ \eta \f$ and \$ \phi \f$ is saved as well.
 *
 * \param genparticle Particle at generation level
 * \param recparticle Particle at reconstruction level
 * \param weight Event weight (optional)
 */
void AliEMCalTriggerRecTrackAnalysisComponent::FillCorrelation(
		const AliVParticle* const genparticle,
		const AliVParticle* const recparticle, double weight) {
	double data[4] = {TMath::Abs(genparticle->Pt()), TMath::Abs(recparticle->Pt()), recparticle->Eta(), recparticle->Phi()};
	fHistos->FillTHnSparse("hTrackPtCorrelation", data, weight);
}

/**
 * Match track to trigger patches. A patch is matched to the track if the track \f$ \eta \f$ and
 * \f$ phi \f$ are within the maximum and minimum of the EMCAL trigger patches.
 * \param rectrack Reconstructed track to be matched to EMCAL patches
 * \param inputpatches Container with reconstructed patches in the event
 * \param outputpatches Output list where matched patches are stored
 */
void AliEMCalTriggerRecTrackAnalysisComponent::MatchTriggerPatches(
    const AliVTrack* rectrack, const TClonesArray* inputpatches,
    TList& outputpatches) const {
  outputpatches.Clear();
  Double_t etaEMCAL = rectrack->GetTrackEtaOnEMCal(),
           phiEMCAL = rectrack->GetTrackPhiOnEMCal();
  for(TIter patchiter = TIter(inputpatches).Begin(); patchiter != TIter::End(); ++patchiter){
    AliEMCALTriggerPatchInfo *recpatch = static_cast<AliEMCALTriggerPatchInfo *>(*patchiter);
    double etamin = TMath::Min(recpatch->GetEtaMin(), recpatch->GetEtaMax()),
        etamax = TMath::Max(recpatch->GetEtaMin(), recpatch->GetEtaMax()),
        phimin = TMath::Min(recpatch->GetPhiMin(), recpatch->GetPhiMax()),
        phimax = TMath::Max(recpatch->GetPhiMin(), recpatch->GetPhiMax());
    if(etaEMCAL >= etamin && etaEMCAL <= etamax && phiEMCAL >= phimin && phiEMCAL <= phimax){
      // Patch accepted
      outputpatches.Add(recpatch);
    }
  }
}

/**
 * Check if the track is matched to any patch that fired the given tigger class
 * \param triggertype Name of the triggertype
 * \param patches List of matched patches
 * \return True if a matched patch of the given type is found
 */
Bool_t AliEMCalTriggerRecTrackAnalysisComponent::HasMatchedPatchOfType(
    TString triggertype, const TList& patches) const {
  Bool_t hasGammaLow = kFALSE,
      hasGammaHigh = kFALSE,
      hasJetLow = kFALSE,
      hasJetHigh = kFALSE;
  for(TIter patchIter = TIter(&patches).Begin(); patchIter != TIter::End(); ++patchIter){
    AliEMCALTriggerPatchInfo *matchedpad = static_cast<AliEMCALTriggerPatchInfo *>(*patchIter);
    if(matchedpad->IsOfflineSimple()) continue;
    if(matchedpad->IsJetHigh()) hasJetHigh = kTRUE;
    if(matchedpad->IsJetLow()) hasJetLow = kTRUE;
    if(matchedpad->IsGammaHigh()) hasGammaHigh = kTRUE;
    if(matchedpad->IsGammaLow()) hasGammaLow = kTRUE;
  }
  // Evaluate status
  if((triggertype.Contains("EMCJHigh") || triggertype.Contains("EMCHighJetOnly")) && hasJetHigh) return kTRUE;
  if((triggertype.Contains("EMCJLow") || triggertype.Contains("EMCLowJetOnly")) && hasJetLow) return kTRUE;
  if((triggertype.Contains("EMCGHigh") || triggertype.Contains("EMCHighGammaOnly")) && hasGammaHigh) return kTRUE;
  if((triggertype.Contains("EMCGLow") || triggertype.Contains("EMCLowGammaOnly")) && hasGammaLow) return kTRUE;
  if(triggertype.Contains("EMCHighBoth") && (hasJetHigh || hasGammaHigh)) return kTRUE;
  if(triggertype.Contains("EMCLowBoth") && (hasJetLow || hasGammaLow)) return kTRUE;
  return kFALSE;
}
