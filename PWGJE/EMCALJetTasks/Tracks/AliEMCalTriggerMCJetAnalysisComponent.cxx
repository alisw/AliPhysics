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
#include <string>
#include <vector>

#include <THistManager.h>
#include <TLinearBinning.h>
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

#include "AliEMCalTriggerAnaClassManager.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerKineCuts.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerMCJetAnalysisComponent.h"
#include "AliEMCalTriggerWeightHandler.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerMCJetAnalysisComponent)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy (I/O) constructor, not to be used
 */
AliEMCalTriggerMCJetAnalysisComponent::AliEMCalTriggerMCJetAnalysisComponent():
  AliEMCalTriggerTracksAnalysisComponent(),
  fMinimumJetPt(20.)
{
}

/**
 * Main constructor, initialising component with default values.
 * \param name Name of the component
 */
AliEMCalTriggerMCJetAnalysisComponent::AliEMCalTriggerMCJetAnalysisComponent(const char* name) :
  AliEMCalTriggerTracksAnalysisComponent(name),
  fMinimumJetPt(20.)
{
}

/**
 * Create histograms for the MC jet analysis
 */
void AliEMCalTriggerMCJetAnalysisComponent::CreateHistos() {
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
      *vertexbinning = fBinning->GetBinning("zvertex");

  const TAxis *trackaxes[6] = {
      DefineAxis("trackpt", *ptbinning),
      DefineAxis("jettpt", jetptbinning ? *jetptbinning : *ptbinning),
      DefineAxis("eta", *etabinning),
      DefineAxis("phi", *phibinning),
      DefineAxis("zvertex", *vertexbinning),
      DefineAxis("mbtrigger", TLinearBinning(2, -0.5, 1.5))
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
    fHistos->CreateTHnSparse(Form("hParticleJetHist%s%s", jetptstring.Data(), name.c_str()), Form("Track-based data for tracks in jets in %s events", title.c_str()), 6, trackaxes, "s");
    fHistos->CreateTHnSparse(Form("hMCJetHist%s%s", jetptstring.Data(), name.c_str()), Form("Reconstructed jets in %s-triggered events", name.c_str()), 4, jetaxes);
  }

  for(int iaxis = 0; iaxis < 6; iaxis++) delete trackaxes[iaxis];
}

/**
 * Analyse particles in a jet with a given minimum jet \f$ p_{t} \f$
 *  -# Select jets with a minimum \f$ p_{t} \f$
 *  -# Fill track-based histograms only for particles found in the given jet
 * \param data Event data
 */
void AliEMCalTriggerMCJetAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  if(!data->GetJetContainerMC()){
    AliError("No Jet container for MC found");
    return;
  }
  std::vector<std::string> triggernames;
  this->GetMachingTriggerNames(triggernames);
  TString jetptstring = Form("jetPt%03d", int(fMinimumJetPt));

  double weight = 1.;
  if(fWeightHandler && data->GetMCEvent()){
    weight = fWeightHandler->GetEventWeight(data->GetMCEvent());
  }

  AliJetContainer *cont = data->GetJetContainerMC();
  cont->ResetCurrentID();
  AliEmcalJet *reconstructedJet = cont->GetNextAcceptJet();
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

/**
 * Fill track-based Histogram with relevant information
 * \param histname Name of the histogram
 * \param track Particle associated to the jet
 * \param jet Reconstructed jet
 * \param vz z-position of the primary vertex
 * \param weight Event weight
 */
void AliEMCalTriggerMCJetAnalysisComponent::FillHistogram(
    const TString& histname, const AliVParticle* track, const AliEmcalJet* jet,
    double vz, double weight) {
  if(!fTriggerClassManager) return;
  double data[6] = {TMath::Abs(track->Pt()), TMath::Abs(jet->Pt()), track->Eta(), track->Phi(), vz, fTriggerClassManager->HasMinBiasTrigger() ? 1. : 0.};
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}

/**
 * Fill jet-based histogram for reconstructed jets with the relevant information
 * \param histname Name of the histogram
 * \param recjet Reconstructed jet
 * \param vz z-position of the primary vertex
 * \param weight Event weight
 */
void AliEMCalTriggerMCJetAnalysisComponent::FillJetHistogram(
    const TString& histname, const AliEmcalJet* recjet, double vz, double weight) {
  double data[4] = {TMath::Abs(recjet->Pt()), recjet->Eta(), recjet->Phi(), vz};
  fHistos->FillTHnSparse(histname.Data(), data, weight);
}
