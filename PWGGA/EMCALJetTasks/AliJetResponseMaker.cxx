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
#include "AliGenPythiaEventHeader.h"
#include "AliMCEvent.h"

ClassImp(AliJetResponseMaker)

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker() : 
  AliAnalysisTaskEmcalJet("AliJetResponseMaker", kTRUE),
  fMCTracksName("MCParticles"),
  fMCJetsName("MCJets"),
  fMaxDistance(0.25),
  fDoWeighting(kFALSE),
  fVertexCut(10),
  fPythiaHeader(0),
  fEventWeight(0),
  fPtHardBin(0),
  fNTrials(0),
  fMCTracks(0),
  fMCJets(0),
  fHistNTrials(0),
  fHistAcceptedEvents(0),
  fHistEvents(0),
  fHistMCJetPhiEta(0),
  fHistMCJetsPt(0),
  fHistMCJetsNEFvsPt(0),
  fHistMCJetsZvsPt(0),
  fHistJetPhiEta(0),
  fHistJetsPt(0),
  fHistJetsNEFvsPt(0),
  fHistJetsZvsPt(0),
  fHistClosestDistance(0),
  fHistClosestDeltaPhi(0),
  fHistClosestDeltaEta(0),
  fHistClosestDeltaPt(0),
  fHistNonMatchedMCJetPt(0),
  fHistNonMatchedJetPt(0),
  fHistPartvsDetecPt(0),
  fHistMissedMCJets(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fMCTracksName("MCParticles"),
  fMCJetsName("MCJets"),
  fMaxDistance(0.25),
  fDoWeighting(kFALSE),
  fVertexCut(10),
  fPythiaHeader(0),
  fEventWeight(0),
  fPtHardBin(0),
  fNTrials(0),
  fMCTracks(0),
  fMCJets(0),
  fHistNTrials(0),
  fHistAcceptedEvents(0),
  fHistEvents(0),
  fHistMCJetPhiEta(0),
  fHistMCJetsPt(0),
  fHistMCJetsNEFvsPt(0),
  fHistMCJetsZvsPt(0),
  fHistJetPhiEta(0),
  fHistJetsPt(0),
  fHistJetsNEFvsPt(0),
  fHistJetsZvsPt(0),
  fHistClosestDistance(0),
  fHistClosestDeltaPhi(0),
  fHistClosestDeltaEta(0),
  fHistClosestDeltaPt(0),
  fHistNonMatchedMCJetPt(0),
  fHistNonMatchedJetPt(0),
  fHistPartvsDetecPt(0),
  fHistMissedMCJets(0)
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

  const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  fHistNTrials = new TH1F("fHistNTrials", "fHistNTrials", 11, 0, 11);
  fHistNTrials->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistNTrials->GetYaxis()->SetTitle("trials");
  fOutput->Add(fHistNTrials);

  fHistAcceptedEvents = new TH1F("fHistAcceptedEvents", "fHistAcceptedEvents", 11, 0, 11);
  fHistAcceptedEvents->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistAcceptedEvents->GetYaxis()->SetTitle("accepted events");
  fOutput->Add(fHistAcceptedEvents);

  fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0, 11);
  fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistEvents->GetYaxis()->SetTitle("total events");
  fOutput->Add(fHistEvents);

  for (Int_t i = 1; i < 12; i++) {
    fHistNTrials->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
    fHistAcceptedEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
    fHistEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  }

  fHistJetPhiEta = new TH2F("fHistJetPhiEta", "fHistJetPhiEta", 20, -2, 2, 32, 0, 6.4);
  fHistJetPhiEta->GetXaxis()->SetTitle("#eta");
  fHistJetPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJetPhiEta);
  
  fHistJetsPt = new TH1F("fHistJetsPt", "fHistJetsPt", fNbins, fMinBinPt, fMaxBinPt);
  fHistJetsPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistJetsPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPt);
  
  fHistJetsZvsPt = new TH2F("fHistJetsZvsPt", "fHistJetsZvsPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJetsZvsPt->GetXaxis()->SetTitle("Z");
  fHistJetsZvsPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  fOutput->Add(fHistJetsZvsPt);
  
  if (fAnaType == kEMCAL) {
    fHistJetsNEFvsPt = new TH2F("fHistJetsNEFvsPt", "fHistJetsNEFvsPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetsNEFvsPt->GetXaxis()->SetTitle("NEF");
    fHistJetsNEFvsPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutput->Add(fHistJetsNEFvsPt);
  }

  fHistMCJetPhiEta = new TH2F("fHistMCJetPhiEta", "fHistMCJetPhiEta", 20, -2, 2, 32, 0, 6.4);
  fHistMCJetPhiEta->GetXaxis()->SetTitle("#eta");
  fHistMCJetPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistMCJetPhiEta);
  
  fHistMCJetsPt = new TH1F("fHistMCJetsPt", "fHistMCJetsPt", fNbins, fMinBinPt, fMaxBinPt);
  fHistMCJetsPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistMCJetsPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistMCJetsPt);
  
  fHistMCJetsZvsPt = new TH2F("fHistMCJetsZvsPt", "fHistMCJetsZvsPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistMCJetsZvsPt->GetXaxis()->SetTitle("Z");
  fHistMCJetsZvsPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  fOutput->Add(fHistMCJetsZvsPt);
  
  if (fAnaType == kEMCAL) {
    fHistMCJetsNEFvsPt = new TH2F("fHistMCJetsNEFvsPt", "fHistMCJetsNEFvsPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
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

  fHistClosestDeltaPt = new TH1F("fHistClosestDeltaPt", "fHistClosestDeltaPt", fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
  fHistClosestDeltaPt->GetXaxis()->SetTitle("#Delta p_{T}");
  fHistClosestDeltaPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistClosestDeltaPt);

  fHistNonMatchedMCJetPt = new TH1F("fHistNonMatchedMCJetPt", "fHistNonMatchedMCJetPt", fNbins, fMinBinPt, fMaxBinPt);
  fHistNonMatchedMCJetPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistNonMatchedMCJetPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedMCJetPt);

  fHistNonMatchedJetPt = new TH1F("fHistNonMatchedJetPt", "fHistNonMatchedJetPt", fNbins, fMinBinPt, fMaxBinPt);
  fHistNonMatchedJetPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistNonMatchedJetPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedJetPt);

  fHistPartvsDetecPt = new TH2F("fHistPartvsDetecPt", "fHistPartvsDetecPt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistPartvsDetecPt->GetXaxis()->SetTitle("p_{T}^{rec}");
  fHistPartvsDetecPt->GetYaxis()->SetTitle("p_{T}^{gen}");
  fOutput->Add(fHistPartvsDetecPt);

  fHistMissedMCJets = new TH1F("fHistMissedMCJets", "fHistMissedMCJets", fNbins, fMinBinPt, fMaxBinPt);
  fHistMissedMCJets->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistMissedMCJets->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistMissedMCJets);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::Run()
{
  // Find the closest jets

  fHistEvents->SetBinContent(fPtHardBin, fHistEvents->GetBinContent(fPtHardBin) + 1);

  if (TMath::Abs(fVertex[2]) > fVertexCut)
    return kFALSE;
  
  fHistAcceptedEvents->SetBinContent(fPtHardBin, fHistEvents->GetBinContent(fPtHardBin) + 1);
  fHistNTrials->SetBinContent(fPtHardBin, fHistEvents->GetBinContent(fPtHardBin) + fNTrials);

  DoJetLoop(fJets, fMCJets, kFALSE);
  DoJetLoop(fMCJets, fJets, kTRUE);

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliJetResponseMaker::FillHistograms()
{
  // Fill histograms.

  const Int_t nMCJets = fMCJets->GetEntriesFast();

  for (Int_t i = 0; i < nMCJets; i++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fMCJets->At(i));

    if (!jet) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet, kTRUE, kFALSE))
      continue;

    if (jet->Pt() > fMaxBinPt)
      continue;

    if (jet->ClosestJet() && jet->ClosestJet()->ClosestJet() == jet && 
        jet->ClosestJetDistance() < fMaxDistance) {    // Matched jet found
      jet->SetMatchedToClosest();
      jet->ClosestJet()->SetMatchedToClosest();
      if (jet->MatchedJet()->Pt() > fMaxBinPt) {
	fHistMissedMCJets->Fill(jet->Pt(), fEventWeight);
      }
      else {
	fHistClosestDistance->Fill(jet->ClosestJetDistance(), fEventWeight);

	Double_t deta = jet->MatchedJet()->Eta() - jet->Eta();
	fHistClosestDeltaEta->Fill(deta, fEventWeight);

	Double_t dphi = jet->MatchedJet()->Phi() - jet->Phi();
	fHistClosestDeltaPhi->Fill(dphi, fEventWeight);

	Double_t dpt = jet->MatchedJet()->Pt() - jet->Pt();
	fHistClosestDeltaPt->Fill(dpt, fEventWeight);

	fHistPartvsDetecPt->Fill(jet->MatchedJet()->Pt(), jet->Pt(), fEventWeight);
      }
    }
    else {
      fHistNonMatchedMCJetPt->Fill(jet->Pt(), fEventWeight);
      fHistMissedMCJets->Fill(jet->Pt(), fEventWeight);
    }

    fHistMCJetsPt->Fill(jet->Pt(), fEventWeight);

    fHistMCJetPhiEta->Fill(jet->Eta(), jet->Phi(), fEventWeight);

    if (fAnaType == kEMCAL)
      fHistMCJetsNEFvsPt->Fill(jet->NEF(), jet->Pt(), fEventWeight);

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, fMCTracks);
      if (track)
	fHistMCJetsZvsPt->Fill(track->Pt() / jet->Pt(), jet->Pt(), fEventWeight);
    }
  }

  const Int_t nJets = fJets->GetEntriesFast();

  for (Int_t i = 0; i < nJets; i++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));

    if (!jet) {
      AliError(Form("Could not receive mc jet %d", i));
      continue;
    }  
    
    if (!AcceptJet(jet))
      continue;

    if (!jet->MatchedJet()) {
      fHistNonMatchedJetPt->Fill(jet->Pt(), fEventWeight);
    }

    fHistJetsPt->Fill(jet->Pt(), fEventWeight);

    fHistJetPhiEta->Fill(jet->Eta(), jet->Phi(), fEventWeight);

    if (fAnaType == kEMCAL)
      fHistJetsNEFvsPt->Fill(jet->NEF(), jet->Pt(), fEventWeight);

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, fTracks);
      if (track)
	fHistJetsZvsPt->Fill(track->Pt() / jet->Pt(), jet->Pt(), fEventWeight);
    }

    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);
      if (cluster) {
	TLorentzVector nP;
	cluster->GetMomentum(nP, fVertex);
	fHistJetsZvsPt->Fill(nP.Pt() / jet->Pt(), jet->Pt(), fEventWeight);
      }
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop(TClonesArray *jets1, TClonesArray *jets2, Bool_t mc)
{
  // Do the jet loop.

  Int_t nJets1 = jets1->GetEntriesFast();
  Int_t nJets2 = jets2->GetEntriesFast();

  for (Int_t i = 0; i < nJets1; i++) {

    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(jets1->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet1, kTRUE, mc))
      continue;

    for (Int_t j = 0; j < nJets2; j++) {
      
      AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(jets2->At(j));
      
      if (!jet2) {
	AliError(Form("Could not receive jet %d", j));
	continue;
      }  
      
      if (!AcceptJet(jet2, kTRUE, !mc))
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
Bool_t AliJetResponseMaker::RetrieveEventObjects()
{
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;
  
  fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());

  if (!fPythiaHeader)
    return kFALSE;

  if (fDoWeighting)
    fEventWeight = fPythiaHeader->EventWeight();
  else
    fEventWeight = 1;

  const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  Double_t pthard = fPythiaHeader->GetPtHard();
  
  for (fPtHardBin = 0; fPtHardBin < 11; fPtHardBin++) {
    if (ptHardLo[fPtHardBin] < pthard && ptHardHi[fPtHardBin] > pthard)
      break;
  }

  fNTrials = fPythiaHeader->Trials();

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::ExecOnce()
{
  // Retrieve event objects.

  if (!fMCJetsName.IsNull() && !fMCJets) {
    fMCJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCJetsName));
    if (!fMCJets) {
      AliError(Form("%s: Could not retrieve mc jets %s!", GetName(), fMCJetsName.Data()));
      return;
    }
    else if (!fMCJets->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fMCJetsName.Data())); 
      fMCJets = 0;
      return;
    }
  }

  if (!fMCTracksName.IsNull() && !fMCTracks) {
    fMCTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCTracksName));
    if (!fMCTracks) {
      AliError(Form("%s: Could not retrieve mc tracks %s!", GetName(), fMCTracksName.Data())); 
      return;
    }
    else {
      TClass *cl = fMCTracks->GetClass();
      if (!cl->GetBaseClass("AliVParticle") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fMCTracksName.Data())); 
	fMCTracks = 0;
	return;
      }
    }
  }

  AliAnalysisTaskEmcalJet::ExecOnce();
}

//________________________________________________________________________
void AliJetResponseMaker::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
