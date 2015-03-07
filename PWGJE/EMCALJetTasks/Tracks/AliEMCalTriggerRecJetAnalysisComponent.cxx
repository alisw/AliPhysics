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
 * Analysis component inspecting tracks with a minimum given jet pt
 *
 *   Author: Markus Fasel
 */
#include <iostream>
#include <string>
#include <vector>

#include <TMath.h>
#include <TString.h>

#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#include "AliEMCalHistoContainer.h"
#include "AliEMCalPtTaskVTrackSelection.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerKineCuts.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerRecJetAnalysisComponent.h"
#include "AliEMCalTriggerWeightHandler.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerRecJetAnalysisComponent)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerRecJetAnalysisComponent::AliEMCalTriggerRecJetAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent(),
  fTrackSelection(NULL),
  fMinimumJetPt(20.),
  fRequestMCtrue(kFALSE),
  fSwapEta(kFALSE),
  fUsePatches(kFALSE)
{
  /*
   * Dummy (I/O) constructor, not to be used
   */
}

//______________________________________________________________________________
AliEMCalTriggerRecJetAnalysisComponent::AliEMCalTriggerRecJetAnalysisComponent(const char* name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fTrackSelection(NULL),
  fMinimumJetPt(20.),
  fRequestMCtrue(kFALSE),
  fSwapEta(kFALSE),
  fUsePatches(kFALSE)
{
  /*
   * Main constructor for the users
   */
}

//______________________________________________________________________________
AliEMCalTriggerRecJetAnalysisComponent::~AliEMCalTriggerRecJetAnalysisComponent() {
  /*
   * Destructor
   */
  if(fTrackSelection) delete fTrackSelection;
}

//______________________________________________________________________________
void AliEMCalTriggerRecJetAnalysisComponent::CreateHistos() {
  /*
   * Create histrogram for the jet pt analysis
   */
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  TString jetptstring = Form("jetPt%03d", int(fMinimumJetPt));
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
      *jetptbinning = fBinning->GetBinning("jetpt"),
      *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi"),
      *vertexbinning = fBinning->GetBinning("zvertex"),
	  *centralitybinning = fBinning->GetBinning("centrality");

  const TAxis *trackaxes[6] = {
      DefineAxis("trackpt", ptbinning),
      DefineAxis("jettpt", jetptbinning ? jetptbinning : ptbinning),
      DefineAxis("eta", etabinning),
      DefineAxis("phi", phibinning),
      DefineAxis("zvertex", vertexbinning),
      DefineAxis("mbtrigger", 2, -0.5, 1.5)
  };

  const TAxis *trackaxes1[5] = {
		  DefineAxis("trackpt", ptbinning),
		  DefineAxis("jetpt", jetptbinning),
		  DefineAxis("eta", etabinning),
		  DefineAxis("centrality", centralitybinning),
		  DefineAxis("dR", 20, 0., 0.5),
  };

  const TAxis *jetaxes[4] = {
      DefineAxis("jetpt", jetptbinning ? jetptbinning : ptbinning),
      DefineAxis("jeteta", etabinning),
      DefineAxis("jetphi", phibinning),
      DefineAxis("zvertex", vertexbinning)
  };

  // Build histograms
  for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
    const std::string name = it->first, &title = it->second;
    fHistos->CreateTHnSparse(Form("hTrackJetHist%s%s", jetptstring.Data(), name.c_str()), Form("Track-based data for tracks in jets in %s events", title.c_str()), 6, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hTrackJetCentralityHist%s%s", jetptstring.Data(), name.c_str()), Form("Track-based histogram for tracks in jets in %s events and centrality", title.c_str()), 5, trackaxes1, "s");
    fHistos->CreateTHnSparse(Form("hMCTrackJetHist%s%s", jetptstring.Data(), name.c_str()), Form("Track-based data for tracks in jets in %s events with MC kinematics", title.c_str()), 6, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hRecJetHist%s%s", jetptstring.Data(), name.c_str()), Form("Reconstructed jets in %s-triggered events", name.c_str()), 4, jetaxes);
  }

  for(int iaxis = 0; iaxis < 6; iaxis++) delete trackaxes[iaxis];

}

//______________________________________________________________________________
void AliEMCalTriggerRecJetAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  /*
   * Analyse tracks from jets with a given minimum pt
   */
  std::vector<std::string> triggernames;
  this->GetMachingTriggerNames(triggernames, fUsePatches);
  TString jetptstring = Form("jetPt%03d", int(fMinimumJetPt));

  // Debugging:
  if(fComponentDebugLevel > 1){
    PrintTriggerNames(triggernames, "RecJets");
    fTriggerDecision->Print();
    if(!fTriggerDecision->CheckConsistency())
      std::cout << "Decision from patches and from strings do not match" << std::endl;
  }

  double weight = 1.;
  if(fWeightHandler && data->GetMCEvent()){
    weight = fWeightHandler->GetEventWeight(data->GetMCEvent());
  }

  AliJetContainer *cont = data->GetJetContainerData();
  AliEmcalJet *reconstructedJet = cont->GetNextAcceptJet(0);
  AliVTrack *foundtrack(NULL);
  const AliVParticle *assocMC(NULL);
  AliCentrality *centralityHandler = data->GetRecEvent()->GetCentrality();
  while(reconstructedJet){
    if(TMath::Abs(reconstructedJet->Pt()) > fMinimumJetPt){
      for(std::vector<std::string>::iterator name = triggernames.begin(); name != triggernames.end(); ++name)
        FillJetHistogram(Form("hRecJetHist%s%s", jetptstring.Data(), name->c_str()), reconstructedJet, data->GetRecEvent()->GetPrimaryVertex()->GetZ(), weight);
      // Jet selected, loop over particles
      for(int ipart = 0; ipart < reconstructedJet->GetNumberOfTracks(); ipart++){
        foundtrack = dynamic_cast<AliVTrack *>(reconstructedJet->TrackAt(ipart, cont->GetParticleContainer()->GetArray()));
        if(!fKineCuts->IsSelected(foundtrack)) continue;
        if(fRequestMCtrue && data->GetMCEvent() && (assocMC = IsMCTrueTrack(foundtrack, data->GetMCEvent()))) continue;
        if(fTrackSelection && !fTrackSelection->IsTrackAccepted(foundtrack)) continue;
        // track selected, fill histogram
        for(std::vector<std::string>::iterator name = triggernames.begin(); name != triggernames.end(); ++name){
          FillHistogram(Form("hTrackJetHist%s%s", jetptstring.Data(), name->c_str()), foundtrack,  reconstructedJet, data->GetRecEvent()->GetPrimaryVertex()->GetZ(), weight);
          FillTrackHistogramCentrality(Form("hTrackJetCentralityHist%s%s", jetptstring.Data(), name->c_str()), foundtrack, reconstructedJet, centralityHandler->GetCentralityPercentile("V0A"), weight);
          if(assocMC){
            FillHistogram(Form("hMCTrackJetHist%s%s", jetptstring.Data(), name->c_str()), assocMC, reconstructedJet, data->GetRecEvent()->GetPrimaryVertex()->GetZ(), weight);
          }
        }
      }
    }
    reconstructedJet = cont->GetNextAcceptJet();
  }
}

//______________________________________________________________________________
const AliVParticle * AliEMCalTriggerRecJetAnalysisComponent::IsMCTrueTrack(
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
  const AliVParticle *mcpart = evnt->GetTrack(label);
  if(!mcpart) return NULL;
  const AliAODMCParticle *aodpart = dynamic_cast<const AliAODMCParticle *>(mcpart);
  if(aodpart) {
    if(!aodpart->IsPhysicalPrimary()) return NULL;
  } else {
    if(!evnt->IsPhysicalPrimary(label)) return NULL;
  }
  return mcpart;
}

//______________________________________________________________________________
void AliEMCalTriggerRecJetAnalysisComponent::FillHistogram(
    const TString& histname, const AliVParticle* track, const AliEmcalJet* jet,
    double vz, double weight) {
  /*
   * Fill Histogram with relevant information
   */
  if(!fTriggerDecision) return;
  double data[6] = {TMath::Abs(track->Pt()), TMath::Abs(jet->Pt()), (fSwapEta ? -1. : 1.) * track->Eta(), track->Phi(), vz, fTriggerDecision->IsMinBias() ? 1. : 0.};
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}

//______________________________________________________________________________
void AliEMCalTriggerRecJetAnalysisComponent::FillJetHistogram(
    const TString& histname, const AliEmcalJet* recjet, double vz, double weight) {
  /*
   * Fill histogram for reconstructed jets with the relevant information
   */
  double data[4] = {TMath::Abs(recjet->Pt()), (fSwapEta ? -1. : 1.) * recjet->Eta(), recjet->Phi(), vz};
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}

//______________________________________________________________________________
void AliEMCalTriggerRecJetAnalysisComponent::FillTrackHistogramCentrality(
		const TString& histname, const AliVTrack* const trk, const AliEmcalJet* jet, double centpercent, double weight) {
	/*
	 * Fill Histogram for tracks with:
	 * - Track pt
	 * - Jet pt
	 * - Track eta
	 * - distance to the main jet axis
	 * - centrality percentile
	 */
	double data[5] = { TMath::Abs(trk->Pt()), TMath::Abs(jet->Pt()), (fSwapEta ? -1. : 1.) * trk->Eta(), centpercent, jet->DeltaR(trk)};
	fHistos->FillTHnSparse(histname.Data(), data, weight);
}

} /* namespace EMCalTriggerPtAnalysis */
