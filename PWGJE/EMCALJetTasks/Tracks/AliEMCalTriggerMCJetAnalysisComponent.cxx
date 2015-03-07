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
 * Analysis component for tracks in jets of MC particles where the jet has a given
 * minimum pt
 *
 *   Author: Markus Fasel
 */
#include <string>
#include <vector>

#include <TMath.h>
#include <TString.h>

#include "AliAODMCParticle.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVVertex.h"

#include "AliEMCalHistoContainer.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerKineCuts.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerMCJetAnalysisComponent.h"
#include "AliEMCalTriggerWeightHandler.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerMCJetAnalysisComponent)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerMCJetAnalysisComponent::AliEMCalTriggerMCJetAnalysisComponent():
  AliEMCalTriggerTracksAnalysisComponent(),
  fMinimumJetPt(20.),
  fUsePatches(kFALSE)
{
  /*
   * Dummy (I/O) constructor, not to be used
   */
}

//______________________________________________________________________________
AliEMCalTriggerMCJetAnalysisComponent::AliEMCalTriggerMCJetAnalysisComponent(const char* name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fMinimumJetPt(20.),
  fUsePatches(kFALSE)
{
  /*
   * Main constructor, to be used
   */
}

//______________________________________________________________________________
void AliEMCalTriggerMCJetAnalysisComponent::CreateHistos() {
  /*
   * Create histograms for the MC jet analysis
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
      *vertexbinning = fBinning->GetBinning("zvertex");

  const TAxis *trackaxes[6] = {
      DefineAxis("trackpt", ptbinning),
      DefineAxis("jettpt", jetptbinning ? jetptbinning : ptbinning),
      DefineAxis("eta", etabinning),
      DefineAxis("phi", phibinning),
      DefineAxis("zvertex", vertexbinning),
      DefineAxis("mbtrigger", 2, -0.5, 1.5)
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
    fHistos->CreateTHnSparse(Form("hParticleJetHist%s%s", jetptstring.Data(), name.c_str()), Form("Track-based data for tracks in jets in %s events", title.c_str()), 6, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hMCJetHist%s%s", jetptstring.Data(), name.c_str()), Form("Reconstructed jets in %s-triggered events", name.c_str()), 4, jetaxes);
  }

  for(int iaxis = 0; iaxis < 6; iaxis++) delete trackaxes[iaxis];
}

//______________________________________________________________________________
void AliEMCalTriggerMCJetAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  /*
   * Analyse particles in a jet with a given minimum jet pt
   */
  if(!data->GetJetContainerMC()){
    AliError("No Jet container for MC found");
    return;
  }
  std::vector<std::string> triggernames;
  this->GetMachingTriggerNames(triggernames, fUsePatches);
  TString jetptstring = Form("jetPt%03d", int(fMinimumJetPt));

  double weight = 1.;
  if(fWeightHandler && data->GetMCEvent()){
    weight = fWeightHandler->GetEventWeight(data->GetMCEvent());
  }

  AliJetContainer *cont = data->GetJetContainerMC();
  AliEmcalJet *reconstructedJet = cont->GetNextAcceptJet(0);
  AliAODMCParticle *foundtrack(NULL);
  while(reconstructedJet){
    if(TMath::Abs(reconstructedJet->Pt()) > fMinimumJetPt){
      for(std::vector<std::string>::iterator name = triggernames.begin(); name != triggernames.end(); ++name)
        FillJetHistogram(Form("hMCJetHist%s%s", jetptstring.Data(), name->c_str()), reconstructedJet, data->GetRecEvent()->GetPrimaryVertex()->GetZ(), weight);
      // Jet selected, loop over particles
      for(int ipart = 0; ipart < reconstructedJet->GetNumberOfTracks(); ipart++){
        foundtrack = dynamic_cast<AliAODMCParticle *>(reconstructedJet->TrackAt(ipart, cont->GetParticleContainer()->GetArray()));
        if(!fKineCuts->IsSelected(foundtrack)) continue;
        if(!foundtrack->Charge()) continue;
        if(!foundtrack->IsPhysicalPrimary()) continue;
        // track selected, fill histogram
        for(std::vector<std::string>::iterator name = triggernames.begin(); name != triggernames.end(); ++name){
          FillHistogram(Form("hParticleJetHist%s%s", jetptstring.Data(), name->c_str()), foundtrack,  reconstructedJet, data->GetRecEvent()->GetPrimaryVertex()->GetZ(), weight);
        }
      }
    }
    reconstructedJet = cont->GetNextAcceptJet();
  }
}

//______________________________________________________________________________
void AliEMCalTriggerMCJetAnalysisComponent::FillHistogram(
    const TString& histname, const AliVParticle* track, const AliEmcalJet* jet,
    double vz, double weight) {
  /*
   * Fill Histogram with relevant information
   */
  if(!fTriggerDecision) return;
  double data[6] = {TMath::Abs(track->Pt()), TMath::Abs(jet->Pt()), track->Eta(), track->Phi(), vz, fTriggerDecision->IsMinBias() ? 1. : 0.};
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}

//______________________________________________________________________________
void AliEMCalTriggerMCJetAnalysisComponent::FillJetHistogram(
    const TString& histname, const AliEmcalJet* recjet, double vz, double weight) {
  /*
   * Fill histogram for reconstructed jets with the relevant information
   */
  double data[4] = {TMath::Abs(recjet->Pt()), recjet->Eta(), recjet->Phi(), vz};
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}
} /* namespace EMCalTriggerPtAnalysis */
