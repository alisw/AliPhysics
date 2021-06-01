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
#include <TClonesArray.h>
#include <TList.h>
#include <TMath.h>
#include <TNtuple.h>

#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliEmcalJet.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetContainer.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliParticleContainer.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliEmcalTrackSelection.h"

#include "AliAnalysisTaskPtEfficiencyJets.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskPtEfficiencyJets)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskPtEfficiencyJets::AliAnalysisTaskPtEfficiencyJets() :
    AliAnalysisTaskEmcalJet(),
    fAnalysisUtils(NULL),
    fMCJetContainer(""),
    fTrackCuts(NULL),
    fTrackNtuple(NULL)
{
}

AliAnalysisTaskPtEfficiencyJets::AliAnalysisTaskPtEfficiencyJets(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fAnalysisUtils(NULL),
    fMCJetContainer(""),
    fTrackCuts(NULL),
    fTrackNtuple(NULL)
{
  fAnalysisUtils = new AliAnalysisUtils();
  SetMakeGeneralHistograms(kTRUE);
}

AliAnalysisTaskPtEfficiencyJets::~AliAnalysisTaskPtEfficiencyJets() {
  if(fAnalysisUtils) delete fAnalysisUtils;
  if(fTrackCuts) delete fTrackCuts;
  if(fTrackNtuple) delete fTrackNtuple;
}

void AliAnalysisTaskPtEfficiencyJets::UserCreateOutputObjects() {
  OpenFile(1);
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  fTrackNtuple = new TNtuple("tracktuple", "track tuple", "partpt:trackpt:jetpt:pthard:parteta:partphi:tracketa:trackphi:vertexz");
  fOutput->Add(fTrackNtuple);
  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskPtEfficiencyJets::Run() {
  // Apply event selection first
  if(!fAnalysisUtils->IsVertexSelected2013pA(fInputEvent)) return kFALSE;
  const AliVVertex *primvertex = fInputEvent->GetPrimaryVertex();

  Float_t pthard = 0;
  AliGenPythiaEventHeader *evheader = dynamic_cast<AliGenPythiaEventHeader *>(fMCEvent->GenEventHeader());
  if(evheader) pthard = evheader->GetPtHard();
  // get the Monte-Carlo jet container
  AliJetContainer *jcmc = dynamic_cast<AliJetContainer *>(fJetCollArray.FindObject(fMCJetContainer.Data()));
  if(!jcmc) return kFALSE;
  TClonesArray *particleArray = jcmc->GetParticleContainer()->GetArray();

  for(TIter partIter = TIter(particleArray).Begin(); partIter != TIter::End(); ++partIter){
    AliVParticle *part = static_cast<AliVParticle *>(*partIter);
    if(!SelectTrueParticle(part)) continue;
    AliVTrack *reconstructed = FindAssociatedTrack(part);
    AliEmcalJet *jet = FindAssociatedJet(part, jcmc);

    // Fill ntuple with
    // particle pt
    // track pt
    // jet pt
    // particle eta
    // particle phi
    // track eta
    // track phi
    // vertex z
    Float_t data[] = {(Float_t)TMath::Abs(part->Pt()), reconstructed ? static_cast<Float_t>(TMath::Abs(reconstructed->Pt())) : 0,
        jet ? static_cast<Float_t>(TMath::Abs(jet->Pt())) : 0, pthard, (Float_t)part->Eta(), (Float_t)part->Phi(), static_cast<Float_t>(reconstructed ? reconstructed->Eta() : -1000.),
            static_cast<Float_t>(reconstructed ? reconstructed->Phi() : -1000.), static_cast<Float_t>(primvertex ? primvertex->GetZ() : -1000.)
    };
    fTrackNtuple->Fill(data);
  }
  PostData(1, fOutput);
  return kTRUE;
}

AliVTrack* AliAnalysisTaskPtEfficiencyJets::FindAssociatedTrack(AliVParticle* trueParticle) {
  AliVTrack *result = NULL;
  for(int itrk = 0; itrk < fInputEvent->GetNumberOfTracks(); itrk++){
    AliVTrack *trk = static_cast<AliVTrack *>(fInputEvent->GetTrack(itrk));
    if(TMath::Abs(trk->GetLabel()) != TMath::Abs(trueParticle->GetLabel())) continue;
    if(!fTrackCuts->IsTrackAccepted(trk)) continue;
    result = trk;
    break;
  }
  return result;
}

AliEmcalJet* AliAnalysisTaskPtEfficiencyJets::FindAssociatedJet(AliVParticle* trueParticle, AliJetContainer* jets) {
  AliEmcalJet *result = NULL;
  for(int ijet = 0; ijet < jets->GetNJets(); ijet++){
      AliEmcalJet *nextjet = jets->GetJet(ijet);
      if(nextjet->ContainsTrack(trueParticle, jets->GetParticleContainer()->GetArray())){
        result = nextjet;
        break;
      }
  }
  return result;
}

bool AliAnalysisTaskPtEfficiencyJets::SelectTrueParticle(AliVParticle* part) {
  if(part->IsA() == AliAODMCParticle::Class()){
    AliAODMCParticle *aodpart = static_cast<AliAODMCParticle *>(part);
    if(!aodpart->IsPhysicalPrimary()) return kFALSE;
    if(!aodpart->Charge()) return kFALSE;
    if(TMath::Abs(aodpart->Eta()) > 0.8) return kFALSE;
    if(TMath::Abs(aodpart->Pt()) < 5.) return kFALSE;
    return kTRUE;
  } else {
    AliMCParticle *esdpart = static_cast<AliMCParticle *>(part);
    if(!fMCEvent->IsPhysicalPrimary(esdpart->GetLabel())) return kFALSE;
    if(!esdpart->Charge()) return kFALSE;
    if(TMath::Abs(esdpart->Eta()) > 0.8) return kFALSE;
    if(TMath::Abs(esdpart->Pt()) < 5.) return kFALSE;
    return kTRUE;
  }
}
