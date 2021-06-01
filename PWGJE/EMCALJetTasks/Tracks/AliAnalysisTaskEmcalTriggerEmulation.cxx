 /**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
#include <sstream>

#include <TArrayD.h>
#include <THashList.h>
#include <THistManager.h>
#include <TMath.h>

#include "AliAODMCParticle.h"
#include "AliClusterContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#include "AliAnalysisTaskEmcalTriggerEmulation.h"
#include "AliEmcalTriggerEmulation.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerBinningFactory.h"
#include "AliEMCalTriggerWeightHandler.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalTriggerEmulation)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalTriggerEmulation::AliAnalysisTaskEmcalTriggerEmulation() :
    AliAnalysisTaskEmcal(),
    fkWeightHandler(nullptr),
    fkTriggerEmulation(nullptr),
    fHistManager(nullptr),
    fNameMCParticles(),
    fNameClusters(),
    fNameTracks(),
    fIsMC(false),
    fEtaRange(-0.8, 0.8),
    fPhiRange(0, 2*TMath::Pi())
{
  /*
   * See header file for details
   */
}

AliAnalysisTaskEmcalTriggerEmulation::AliAnalysisTaskEmcalTriggerEmulation(const char *name) :
    AliAnalysisTaskEmcal(name, true),
    fkWeightHandler(nullptr),
    fkTriggerEmulation(nullptr),
    fHistManager(nullptr),
    fNameMCParticles(),
    fNameClusters(),
    fNameTracks(),
    fIsMC(false),
    fEtaRange(-0.8, 0.8),
    fPhiRange(0, 2*TMath::Pi())
{
  /*
   * See header file for details
   */
}

AliAnalysisTaskEmcalTriggerEmulation::~AliAnalysisTaskEmcalTriggerEmulation() {
  /*
   * See header file for details
   */
  if(fkWeightHandler) delete fkWeightHandler;
  if(fkTriggerEmulation) delete fkTriggerEmulation;
  if(fHistManager) delete fHistManager;
}

void AliAnalysisTaskEmcalTriggerEmulation::UserCreateOutputObjects(){
  /*
   * See header file for details
   */
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  AliEMCalTriggerBinningComponent binhandler;
  AliEMCalTriggerBinningFactory binfact;
  binfact.Create(&binhandler);

  std::stringstream kinestring;
  kinestring << "eta " << fEtaRange << ", phi " << fPhiRange;

  fHistManager = new THistManager("histos");
  fHistManager->CreateTH1("hEvents", "Event counter", 1, 0.5, 1.5);
  fHistManager->CreateTH1("hVertexZ", "z-component of the primary vertex", 100, -40, 40);
  fHistManager->CreateTH1("hParticles", Form("Particles (%s)", kinestring.str().c_str()), *binhandler.GetBinning("pt"));
  fHistManager->CreateTH1("hClusters", "Clusters", *binhandler.GetBinning("pt"));
  fHistManager->CreateTH1("hTracks", Form("Tracks (%s)", kinestring.str().c_str()), *binhandler.GetBinning("pt"));

  for(auto histo : *fHistManager) fOutput->Add(histo);
}

bool AliAnalysisTaskEmcalTriggerEmulation::Run(){
  /*
   * See header file for details
   */
  double weight = 1;
  if(fIsMC){
    if(fkWeightHandler) weight = fkWeightHandler->GetEventWeight(this->fPythiaHeader);
    const AliAODMCParticle *mcpart(nullptr);
    // Loop MC particles
    for(auto mcen : GetMCParticleContainer(fNameMCParticles.Data())->accepted()){
      mcpart =static_cast<const AliAODMCParticle *>(mcen);
      if(!mcpart->IsPhysicalPrimary()) continue;
      if(mcpart->Charge() == 0) continue;
      if(!fEtaRange.IsInRange(mcpart->Eta())) continue;
      if(!fPhiRange.IsInRange(mcpart->Phi())) continue;
      fHistManager->FillTH1("hParticles", TMath::Abs(mcpart->Pt()), weight);
    }
    if(fkTriggerEmulation){
      if(!fkTriggerEmulation->SelectEvent(GetClusterContainer(fNameClusters.Data()))) return false;
    }
  }

  fHistManager->FillTH1("hEvents", 1, weight);
  fHistManager->FillTH1("hVertexZ", InputEvent()->GetPrimaryVertex()->GetZ(), weight);

  // Loop clusters
  const AliVCluster *clust(nullptr);
  for(auto clusten : GetClusterContainer(fNameClusters.Data())->accepted()){
    clust = static_cast<const AliVCluster *>(clusten);
    fHistManager->FillTH1("hClusters", TMath::Abs(clust->E()), weight);
  }

  // Loop tracks
  // Assume: Track quality cuts done by the track container
  const AliVTrack *track(nullptr);
  const AliAODMCParticle *assocMC(nullptr), *tmp(nullptr);
  for(auto tracken : GetTrackContainer(fNameTracks.Data())->accepted()){
    track = static_cast<const AliVTrack *>(tracken);
    double pt = TMath::Abs(track->Pt());
    if(fIsMC){
      for(auto mcen : GetMCParticleContainer(fNameMCParticles.Data())->accepted()){
        tmp = static_cast<const AliAODMCParticle *>(mcen);
        if(TMath::Abs(tmp->GetLabel()) == TMath::Abs(track->GetLabel())){
          assocMC = tmp;
          break;
        }
      }
      if(!assocMC) continue;
      if(!assocMC->IsPhysicalPrimary()) continue;
      pt = TMath::Abs(assocMC->Pt());
    }
    if(!fEtaRange.IsInRange(track->Eta())) continue;
    if(!fPhiRange.IsInRange(track->Phi())) continue;
    fHistManager->FillTH1("hTracks", pt, weight);
  }

  return true;
}
