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
/*
 * Track analysis component: Loops over tracks from the EMCal track container and
 * counts the tracks in histograms
 *
 *   Author: Markus Fasel
 */
#include <map>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TString.h>

#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerKineCuts.h"
#include "AliEMCalPtTaskVTrackSelection.h"
#include "AliEMCalTriggerRecTrackAnalysisComponent.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerRecTrackAnalysisComponent)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerRecTrackAnalysisComponent::AliEMCalTriggerRecTrackAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent(),
  fTrackSelection(NULL),
  fSwapEta(kFALSE),
  fUsePatches(kFALSE),
  fRequestMCtrue(kFALSE)
{
  /*
   * Dummy (I/O) constructor
   */
}

//______________________________________________________________________________
AliEMCalTriggerRecTrackAnalysisComponent::AliEMCalTriggerRecTrackAnalysisComponent(const char *name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fTrackSelection(NULL),
  fSwapEta(kFALSE),
  fUsePatches(kFALSE),
  fRequestMCtrue(kFALSE)
{
  /*
   * Main constructor
   */
}

//______________________________________________________________________________
AliEMCalTriggerRecTrackAnalysisComponent::~AliEMCalTriggerRecTrackAnalysisComponent() {
  /*
   * Destructor, taking care of the track selection
   */
  if(fTrackSelection) delete fTrackSelection;
}

//______________________________________________________________________________
void AliEMCalTriggerRecTrackAnalysisComponent::CreateHistos() {
  /*
   * Create histograms of the track analysis component. For each trigger class we have
   * - tracks with esd information
   * - tracks with MC information
   * - tracks with clusters and esd information
   * - tracks with clusters and MC information
   */
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  // Create trigger definitions
  std::map<std::string, std::string> triggerCombinations;
  const char *triggernames[11] = {"MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh",
      "EMCGLow", "EMCHighBoth", "EMCHighGammaOnly", "EMCHighJetOnly",
      "EMCLowBoth", "EMCLowGammaOnly", "EMCLowJetOnly"};
  // Define names and titles for different triggers in the histogram container
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[0], "min. bias events"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[1], "jet-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[2], "jet-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[3], "gamma-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[4], "gamma-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[5], "jet and gamma triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[6], "exclusively gamma-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[7], "exclusively jet-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[8], "jet and gamma triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[9], "exclusively gamma-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[10], "exclusively-triggered events (low threshold)"));

  // Create axis definitions
  const AliEMCalTriggerBinningDimension *ptbinning = fBinning->GetBinning("pt"),
      *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi"),
      *vertexbinning = fBinning->GetBinning("zvertex");

  const TAxis *trackaxes[5] = {
      DefineAxis("pt", ptbinning),
      DefineAxis("eta", etabinning),
      DefineAxis("phi", phibinning),
      DefineAxis("zvertex", vertexbinning),
      DefineAxis("mbtrigger", 2, -0.5, 1.5)
  };

  // Build histograms
  for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
    const std::string name = it->first, &title = it->second;
    fHistos->CreateTHnSparse(Form("hTrackHist%s", name.c_str()), Form("Track-based data for %s events", title.c_str()), 5, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hTrackInAcceptanceHist%s", name.c_str()), Form("Track-based data for %s events  for tracks matched to EMCal clusters", title.c_str()), 5, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hMCTrackHist%s", name.c_str()), Form("Track-based data for %s events with MC kinematics", title.c_str()), 5, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hMCTrackInAcceptanceHist%s", name.c_str()), Form("Track-based data for %s events with MC kinematics for tracks matched to EMCal clusters", title.c_str()), 5, trackaxes, "s");
  }

  for(int iaxis = 0; iaxis < 5; iaxis++) delete trackaxes[iaxis];
}

//______________________________________________________________________________
void AliEMCalTriggerRecTrackAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  /*
   * Run track loop on list of matching tracks
   *
   * @param data: the event data
   */
  AliDebug(1, Form("Number of matched tracks: %d", data->GetMatchedTrackContainer()->GetEntries()));

  std::vector<std::string> triggernames;
  this->GetMachingTriggerNames(triggernames, fUsePatches);

  AliVTrack *track(NULL);
  AliVParticle *assocMC(NULL);
  TIter trackIter(data->GetMatchedTrackContainer());
  while((track = dynamic_cast<AliVTrack *>(trackIter()))){
    // Apply track selection
    assocMC = NULL;
    if(fKineCuts && !fKineCuts->IsSelected(track)) continue;
    if(fTrackSelection && !fTrackSelection->IsTrackAccepted(track)) continue;

    if(fRequestMCtrue && data->GetMCEvent() && !(assocMC = IsMCTrueTrack(track, data->GetMCEvent()))) continue;
    // Try to match the cluster
    Bool_t hasCluster = kFALSE;
    AliVCluster *clust(NULL);
    if(track->GetEMCALcluster() >= 0 && (clust = dynamic_cast<AliVCluster *>(data->GetClusterContainer()->At(track->GetEMCALcluster()))))
      hasCluster = kTRUE;

    // Fill histograms
    for(std::vector<std::string>::iterator name = triggernames.begin(); name != triggernames.end(); name++){
      FillHistogram(Form("hTrackHist%s", name->c_str()), track, NULL, data->GetRecEvent(), kFALSE);
      if(hasCluster) FillHistogram(Form("hTrackInAcceptanceHist%s", name->c_str()), track, NULL, data->GetRecEvent(), kFALSE);
      if(assocMC){
        FillHistogram(Form("hMCTrackHist%s", name->c_str()), track, NULL, data->GetRecEvent(), kTRUE);
        if(hasCluster) FillHistogram(Form("hMCTrackInAcceptanceHist%s", name->c_str()), track, NULL, data->GetRecEvent(), kTRUE);
      }
    }
  }
}

//______________________________________________________________________________
AliVParticle * AliEMCalTriggerRecTrackAnalysisComponent::IsMCTrueTrack(
    const AliVTrack* const trk, const AliMCEvent* evnt) const {
  /*
   * Check according to the associated MC information whether the track is a MC true track,
   * and whether it is physical primary
   *
   * @param trk: track to check
   * @param evnt: MC event information necessary for the check
   *
   * @return: the associated MC particle (NULL if not MC true)
   */
  int label = TMath::Abs(trk->GetLabel());
  AliVParticle *mcpart = evnt->GetTrack(label);
  if(!mcpart) return NULL;
  if(!evnt->IsPhysicalPrimary(label)) return NULL;
  return mcpart;
}

//______________________________________________________________________________
void AliEMCalTriggerRecTrackAnalysisComponent::FillHistogram(
    const TString& histname, const AliVTrack* const trk,
    const AliVParticle* assocMC, const AliVEvent* const recev,
    Bool_t useMCkine) {
  /*
   *
   */
  if(useMCkine && !assocMC) return;
  double data[5];
  data[0] = useMCkine ? TMath::Abs(assocMC->Pt()) : TMath::Abs(trk->Pt());
  data[1] = (fSwapEta ? -1. : 1.) * (useMCkine ? assocMC->Eta() : trk->Eta());
  data[2] = useMCkine ? assocMC->Phi() : trk->Phi();
  data[3] = recev->GetPrimaryVertex()->GetZ();
  data[4] = fTriggerDecision->IsMinBias();
  fHistos->FillTHnSparse(histname.Data(), data);
}

} /* namespace EMCalTriggerPtAnalysis */
