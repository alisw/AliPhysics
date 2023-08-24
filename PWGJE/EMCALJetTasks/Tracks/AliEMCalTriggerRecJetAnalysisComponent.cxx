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
#include <iostream>
#include <string>
#include <vector>

#include <THistManager.h>
#include <TLinearBinning.h>
#include <TMath.h>
#include <TString.h>

#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliEmcalTrackSelection.h"
#include "AliJetContainer.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#include "AliEMCalTriggerAnaClassManager.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerKineCuts.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerRecJetAnalysisComponent.h"
#include "AliEMCalTriggerWeightHandler.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerRecJetAnalysisComponent)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy (I/O) constructor, not to be used
 */
AliEMCalTriggerRecJetAnalysisComponent::AliEMCalTriggerRecJetAnalysisComponent() :
  AliEMCalTriggerTracksAnalysisComponent(),
  fTrackSelection(NULL),
  fMinimumJetPt(20.),
  fRequestMCtrue(kFALSE),
  fSwapEta(kFALSE)
{
}

/**
 * Main constructor for the users. Initializes all fields with default values.
 * \param name Name of the analysis component
 */
AliEMCalTriggerRecJetAnalysisComponent::AliEMCalTriggerRecJetAnalysisComponent(const char* name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fTrackSelection(NULL),
  fMinimumJetPt(20.),
  fRequestMCtrue(kFALSE),
  fSwapEta(kFALSE)
{
}

/**
 * Destructor, deletes elements the component has ownership over.
 */
AliEMCalTriggerRecJetAnalysisComponent::~AliEMCalTriggerRecJetAnalysisComponent() {
  if(fTrackSelection) delete fTrackSelection;
}

/**
 * Create histrogram for the jet pt analysis
 */
void AliEMCalTriggerRecJetAnalysisComponent::CreateHistos() {
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  TString jetptstring = Form("jetPt%03d", int(fMinimumJetPt));
  // Create trigger definitions
  std::map<std::string, std::string> triggerCombinations;
  GetAllTriggerNamesAndTitles(triggerCombinations);
  // Create axis definitions
  const TBinning *ptbinning = fBinning->GetBinning("pt"),
      *jetptbinning = fBinning->GetBinning("jetpt"),
      *etabinning = fBinning->GetBinning("eta"),
      *phibinning = fBinning->GetBinning("phi"),
      *vertexbinning = fBinning->GetBinning("zvertex"),
	  *centralitybinning = fBinning->GetBinning("centrality");

  const TAxis *trackaxes[6] = {
      DefineAxis("trackpt", *ptbinning),
      DefineAxis("jettpt", jetptbinning ? *jetptbinning : *ptbinning),
      DefineAxis("eta", *etabinning),
      DefineAxis("phi", *phibinning),
      DefineAxis("zvertex", *vertexbinning),
      DefineAxis("mbtrigger", TLinearBinning(2, -0.5, 1.5))
  };

  const TAxis *trackaxes1[5] = {
		  DefineAxis("trackpt", *ptbinning),
		  DefineAxis("jetpt", jetptbinning ? *jetptbinning : *ptbinning),
		  DefineAxis("eta", *etabinning),
		  DefineAxis("centrality", *centralitybinning),
		  DefineAxis("dR", TLinearBinning(20, 0., 0.5)),
  };

  const TAxis *jetaxes[4] = {
      DefineAxis("jetpt", jetptbinning ? *jetptbinning : *ptbinning),
      DefineAxis("jeteta", *etabinning),
      DefineAxis("jetphi", *phibinning),
      DefineAxis("zvertex", *vertexbinning)
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

/**
 * Run the event loop
 * -# Select reconstructed jets
 * -# Fill jet based histogram
 * -# Associated particles to the jet and make a track selection
 * -# Fill track based histogram with reconstructed information
 * -# If we have Monte-Carlo information an the track has an associated true particle, fill track based histogram with true information.
 * \param data All data for the given event
 */
void AliEMCalTriggerRecJetAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  /*
   * Analyse tracks from jets with a given minimum pt
   */
  std::vector<std::string> triggernames;
  this->GetMachingTriggerNames(triggernames);
  TString jetptstring = Form("jetPt%03d", int(fMinimumJetPt));

  // Debugging:
  if(fComponentDebugLevel > 1){
    PrintTriggerNames(triggernames, "RecJets");
  }

  double weight = 1.;
  if(fWeightHandler && data->GetMCEvent()){
    weight = fWeightHandler->GetEventWeight(data->GetMCEvent());
  }

  AliJetContainer *cont = data->GetJetContainerData();
  cont->ResetCurrentID();
  AliEmcalJet *reconstructedJet = cont->GetNextAcceptJet();
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

/**
 * Check according to the associated MC information whether the track is a MC true track,
 * and whether it is physical primary
 * \param trk track to check
 * \param evnt MC event information necessary for the check
 * \return the associated MC particle (NULL if not MC true)
 */
const AliVParticle * AliEMCalTriggerRecJetAnalysisComponent::IsMCTrueTrack(
    const AliVTrack* const trk, const AliMCEvent* evnt) const {
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

/**
 * Fill histogram for tracks with relevant information
 * @param histname Name of the histogram to be fille
 * @param track Reconstructed selected particle
 * @param jet Jet in which the particle was found
 * @param vz z-position of the primary vertex
 * @param weight Event weight
 */
void AliEMCalTriggerRecJetAnalysisComponent::FillHistogram(
    const TString& histname, const AliVParticle* track, const AliEmcalJet* jet,
    double vz, double weight) {

  if(!fTriggerClassManager) return;
  double data[6] = {TMath::Abs(track->Pt()), TMath::Abs(jet->Pt()), (fSwapEta ? -1. : 1.) * track->Eta(), track->Phi(), vz, fTriggerClassManager->HasMinBiasTrigger() ? 1. : 0.};
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}

/**
 * Fill histogram for reconstructed jets with the relevant information
 * \param histname Name of the histogram to be filled
 * \param recjet Reconstructed jet
 * \param vz z-position of the primary vertex
 * \param weight Event weight
 */
void AliEMCalTriggerRecJetAnalysisComponent::FillJetHistogram(
    const TString& histname, const AliEmcalJet* recjet, double vz, double weight) {
  double data[4] = {TMath::Abs(recjet->Pt()), (fSwapEta ? -1. : 1.) * recjet->Eta(), recjet->Phi(), vz};
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}

/**
 * Fill Histogram for tracks with:
 * - Track pt
 * - Jet pt
 * - Track eta
 * - distance to the main jet axis
 * - centrality percentile
 * \param histname Name of the histogram to be filled
 * \param trk Reconstructed and selected tracl
 * \param jet Reconstricted jet
 * \param centpercent Centrality percentile
 * \param weight Event weight
 */
void AliEMCalTriggerRecJetAnalysisComponent::FillTrackHistogramCentrality(
		const TString& histname, const AliVTrack* const trk, const AliEmcalJet* jet, double centpercent, double weight) {
	/*
	 */
	double data[5] = { TMath::Abs(trk->Pt()), TMath::Abs(jet->Pt()), (fSwapEta ? -1. : 1.) * trk->Eta(), centpercent, jet->DeltaR(trk)};
	fHistos->FillTHnSparse(histname.Data(), data, weight);
}
