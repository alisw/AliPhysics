// $Id$
//
// Emcal jet response matrix maker task.
//
// Author: S. Aiola

#include "AliJetResponseMaker.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliFJWrapper.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliMCEvent.h"

ClassImp(AliJetResponseMaker)

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker() : 
  AliAnalysisTaskEmcal("AliJetResponseMaker"),
  fMCTracksName("MCParticles"),
  fMCJetsName("MCJets"),
  fMaxDistance(0.2),
  fMCTracks(0),
  fMCJets(0),
  fHistMCJetPhiEta(0),
  fHistMCJetsPt(0),
  fHistMCJetsPtNonBias(0),
  fHistMCJetsNEFvsPt(0),
  fHistMCJetsZvsPt(0),
  fHistJetPhiEta(0),
  fHistJetsPt(0),
  fHistJetsPtNonBias(0),
  fHistJetsNEFvsPt(0),
  fHistJetsZvsPt(0),
  fHistClosestDistance(0),
  fHistClosestDeltaPhi(0),
  fHistClosestDeltaEta(0),
  fHistClosestDeltaPt(0),
  fHistPartvsDetecPt(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker(const char *name) : 
  AliAnalysisTaskEmcal(name),
  fMCTracksName("MCParticles"),
  fMCJetsName("MCJets"),
  fMaxDistance(0.2),
  fMCTracks(0),
  fMCJets(0),
  fHistMCJetPhiEta(0),
  fHistMCJetsPt(0),
  fHistMCJetsPtNonBias(0),
  fHistMCJetsNEFvsPt(0),
  fHistMCJetsZvsPt(0),
  fHistJetPhiEta(0),
  fHistJetsPt(0),
  fHistJetsPtNonBias(0),
  fHistJetsNEFvsPt(0),
  fHistJetsZvsPt(0),
  fHistClosestDistance(0),
  fHistClosestDeltaPhi(0),
  fHistClosestDeltaEta(0),
  fHistClosestDeltaPt(0),
  fHistPartvsDetecPt(0)
{
  // Standard constructor.

}

//________________________________________________________________________
AliJetResponseMaker::~AliJetResponseMaker()
{
  // Destructor
}

//________________________________________________________________________
void AliJetResponseMaker::UserCreateOutputObjects()
{
  // Create user objects.

  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

  fHistJetPhiEta = new TH2F("fHistJetPhiEta", "fHistJetPhiEta", 20, -2, 2, 32, 0, 6.4);
  fHistJetPhiEta->GetXaxis()->SetTitle("#eta");
  fHistJetPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJetPhiEta);
  
  fHistJetsPt = new TH1F("fHistJetsPt", "fHistJetsPt", fNbins, fMinPt, fMaxPt);
  fHistJetsPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistJetsPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPt);
  
  fHistJetsPtNonBias = new TH1F("fHistJetsPtNonBias", "fHistJetsPtNonBias", fNbins, fMinPt, fMaxPt);
  fHistJetsPtNonBias->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistJetsPtNonBias->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPtNonBias);
  
  fHistJetsZvsPt = new TH2F("fHistJetsZvsPt", "fHistJetsZvsPt", fNbins, 0, 1.2, fNbins, fMinPt, fMaxPt);
  fHistJetsZvsPt->GetXaxis()->SetTitle("Z");
  fHistJetsZvsPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  fOutput->Add(fHistJetsZvsPt);
  
  if (fAnaType == kEMCAL) {
    fHistJetsNEFvsPt = new TH2F("fHistJetsNEFvsPt", "fHistJetsNEFvsPt", fNbins, 0, 1.2, fNbins, fMinPt, fMaxPt);
    fHistJetsNEFvsPt->GetXaxis()->SetTitle("NEF");
    fHistJetsNEFvsPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutput->Add(fHistJetsNEFvsPt);
  }

  fHistMCJetPhiEta = new TH2F("fHistMCJetPhiEta", "fHistMCJetPhiEta", 20, -2, 2, 32, 0, 6.4);
  fHistMCJetPhiEta->GetXaxis()->SetTitle("#eta");
  fHistMCJetPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistMCJetPhiEta);
  
  fHistMCJetsPt = new TH1F("fHistMCJetsPt", "fHistMCJetsPt", fNbins, fMinPt, fMaxPt);
  fHistMCJetsPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistMCJetsPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistMCJetsPt);
  
  fHistMCJetsPtNonBias = new TH1F("fHistMCJetsPtNonBias", "fHistMCJetsPtNonBias", fNbins, fMinPt, fMaxPt);
  fHistMCJetsPtNonBias->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistMCJetsPtNonBias->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistMCJetsPtNonBias);
  
  fHistMCJetsZvsPt = new TH2F("fHistMCJetsZvsPt", "fHistMCJetsZvsPt", fNbins, 0, 1.2, fNbins, fMinPt, fMaxPt);
  fHistMCJetsZvsPt->GetXaxis()->SetTitle("Z");
  fHistMCJetsZvsPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  fOutput->Add(fHistMCJetsZvsPt);
  
  if (fAnaType == kEMCAL) {
    fHistMCJetsNEFvsPt = new TH2F("fHistMCJetsNEFvsPt", "fHistMCJetsNEFvsPt", fNbins, 0, 1.2, fNbins, fMinPt, fMaxPt);
    fHistMCJetsNEFvsPt->GetXaxis()->SetTitle("NEF");
    fHistMCJetsNEFvsPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutput->Add(fHistMCJetsNEFvsPt);
  }

  fHistClosestDistance = new TH1F("fHistClosestDistance", "fHistClosestDistance", 50, 0, fMaxDistance * 1.2);
  fHistClosestDistance->GetXaxis()->SetTitle("d");
  fHistClosestDistance->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClosestDistance);

  fHistClosestDeltaPhi = new TH1F("fHistClosestDeltaPhi", "fHistClosestDeltaPhi", 64, -1.6, 4.8);
  fHistClosestDeltaPhi->GetXaxis()->SetTitle("#Delta#phi");
  fHistClosestDeltaPhi->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClosestDeltaPhi);

  fHistClosestDeltaEta = new TH1F("fHistClosestDeltaEta", "fHistClosestDeltaEta", TMath::CeilNint(fMaxEta - fMinEta) * 20, fMinEta * 2, fMaxEta * 2);
  fHistClosestDeltaEta->GetXaxis()->SetTitle("#Delta#eta");
  fHistClosestDeltaEta->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClosestDeltaEta);

  fHistClosestDeltaPt = new TH1F("fHistClosestDeltaPt", "fHistClosestDeltaPt", fNbins, -fMaxPt / 2, fMaxPt / 2);
  fHistClosestDeltaPt->GetXaxis()->SetTitle("#Delta p_{T}");
  fHistClosestDeltaPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClosestDeltaPt);

  fHistPartvsDetecPt = new TH2F("fHistPartvsDetecPt", "fHistPartvsDetecPt", fNbins, fMinPt, fMaxPt, fNbins, fMinPt, fMaxPt);
  fHistPartvsDetecPt->GetXaxis()->SetTitle("p_{T}^{det}");
  fHistPartvsDetecPt->GetYaxis()->SetTitle("p_{T}^{rec}");
  fOutput->Add(fHistPartvsDetecPt);

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliJetResponseMaker::FillHistograms()
{
  // Fill histograms.

  // Find the closest jets
  DoJetLoop(fJets, fMCJets);
  DoJetLoop(fMCJets, fJets);

  Int_t nJets = fJets->GetEntriesFast();

  for (Int_t i = 0; i < nJets; i++) {

    AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fJets->At(i));

    if (!jet) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    fHistJetsPtNonBias->Fill(jet->Pt());

    if (jet->MaxTrackPt() < fPtBiasJetTrack && (fAnaType == kTPC || jet->MaxClusterPt() < fPtBiasJetClus))
	continue;

    if (jet->ClosestJet() && jet->ClosestJet()->ClosestJet() == jet && jet->ClosestJetDistance() < fMaxDistance) {    // Matched jet found
      jet->SetMatchedToClosest();
      jet->ClosestJet()->SetMatchedToClosest();
      fHistClosestDistance->Fill(jet->ClosestJetDistance());
      Double_t deta = jet->Eta() - jet->MatchedJet()->Eta();
      fHistClosestDeltaEta->Fill(deta);
      Double_t dphi = jet->Phi() - jet->MatchedJet()->Phi();
      fHistClosestDeltaPhi->Fill(dphi);
      Double_t dpt = jet->Pt() - jet->MatchedJet()->Pt();
      fHistClosestDeltaPt->Fill(dpt);
      fHistPartvsDetecPt->Fill(jet->Pt(), jet->MatchedJet()->Pt());
    }

    fHistJetsPt->Fill(jet->Pt());

    fHistJetPhiEta->Fill(jet->Eta(), jet->Phi());

    if (fAnaType == kEMCAL)
      fHistJetsNEFvsPt->Fill(jet->NEF(), jet->Pt());

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, fTracks);
      if (track)
	fHistJetsZvsPt->Fill(track->Pt() / jet->Pt(), jet->Pt());
    }

    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);
      if (cluster) {
	TLorentzVector nP;
	cluster->GetMomentum(nP, fVertex);
	fHistJetsZvsPt->Fill(nP.Pt() / jet->Pt(), jet->Pt());
      }
    }
  }

 Int_t nMCJets = fMCJets->GetEntriesFast();

  for (Int_t i = 0; i < nMCJets; i++) {

    AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fMCJets->At(i));

    if (!jet) {
      AliError(Form("Could not receive mc jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    fHistMCJetsPtNonBias->Fill(jet->Pt());

    if (jet->MaxTrackPt() < fPtBiasJetTrack)
	continue;

    fHistMCJetsPt->Fill(jet->Pt());

    fHistMCJetPhiEta->Fill(jet->Eta(), jet->Phi());

    if (fAnaType == kEMCAL)
      fHistMCJetsNEFvsPt->Fill(jet->NEF(), jet->Pt());

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, fMCTracks);
      if (track)
	fHistMCJetsZvsPt->Fill(track->Pt() / jet->Pt(), jet->Pt());
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop(TClonesArray *jets1, TClonesArray *jets2)
{
  // Do the jet loop.

  Int_t nJets1 = jets1->GetEntriesFast();
  Int_t nJets2 = jets2->GetEntriesFast();

  for (Int_t i = 0; i < nJets1; i++) {

    AliEmcalJet* jet1 = dynamic_cast<AliEmcalJet*>(jets1->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet1))
      continue;
    
    if (jet1->MaxTrackPt() < fPtBiasJetTrack && (fAnaType == kTPC || jet1->IsMC() || jet1->MaxClusterPt() < fPtBiasJetClus))
	continue;

    for (Int_t j = 0; j < nJets2; j++) {
      
      AliEmcalJet* jet2 = dynamic_cast<AliEmcalJet*>(jets2->At(j));
      
      if (!jet2) {
	AliError(Form("Could not receive jet %d", j));
	continue;
      }  
      
      if (!AcceptJet(jet2))
	continue;
      
      if (jet2->MaxTrackPt() < fPtBiasJetTrack && (fAnaType == kTPC || jet2->IsMC() || jet2->MaxClusterPt() < fPtBiasJetClus))
	continue;
      
      Double_t deta = jet2->Eta() - jet1->Eta();
      Double_t dphi = jet2->Phi() - jet1->Phi();
      Double_t d = TMath::Sqrt(deta * deta + dphi * dphi);

      if (d < jet1->ClosestJetDistance()) {
	jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
	jet1->SetClosestJet(jet2, d);
      }
      else if (d < jet1->SecondClosestJetDistance()) {
	jet1->SetSecondClosestJet(jet2, d);
      }
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::RetrieveEventObjects()
{
  // Retrieve event objects.

  AliAnalysisTaskEmcal::RetrieveEventObjects();
  
  if (!fMCJetsName.IsNull()) {
    fMCJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCJetsName));
    if (!fMCJets) {
      AliWarning(Form("Could not retrieve MC jets %s!", fMCJetsName.Data())); 
    }
  }

  if (!fMCTracksName.IsNull()) {
    fMCTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCTracksName));
    if (!fMCJets) {
      AliWarning(Form("Could not retrieve MC tracks %s!", fMCTracksName.Data())); 
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
