// $Id: AliJetResponseMaker.cxx  $
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
  fMCPartName("Tracks"),
  fMCJetsName("mcJets"),
  fMCParts(0),
  fMCJets(0),
  fHistMCJetPhiEta(0),
  fHistMCJetsPt(0),
  fHistMCJetsPtTrack(0),
  fHistMCJetsPtClus(0),
  fHistMCJetsPtNonBias(0),
  fHistMCLeadingJetPt(0),
  fHistMCJetsNEFvsPt(0),
  fHistMCJetsZvsPt(0),
  fHistJetPhiEta(0),
  fHistJetsPt(0),
  fHistJetsPtTrack(0),
  fHistJetsPtClus(0),
  fHistJetsPtNonBias(0),
  fHistLeadingJetPt(0),
  fHistJetsNEFvsPt(0),
  fHistJetsZvsPt(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker(const char *name) : 
  AliAnalysisTaskEmcal(name),
  fMCPartName("Tracks"),
  fMCJetsName("mcJets"),
  fMCParts(0),
  fMCJets(0),
  fHistMCJetPhiEta(0),
  fHistMCJetsPt(0),
  fHistMCJetsPtTrack(0),
  fHistMCJetsPtClus(0),
  fHistMCJetsPtNonBias(0),
  fHistMCLeadingJetPt(0),
  fHistMCJetsNEFvsPt(0),
  fHistMCJetsZvsPt(0),
  fHistJetPhiEta(0),
  fHistJetsPt(0),
  fHistJetsPtTrack(0),
  fHistJetsPtClus(0),
  fHistJetsPtNonBias(0),
  fHistLeadingJetPt(0),
  fHistJetsNEFvsPt(0),
  fHistJetsZvsPt(0)
{
  // Standard constructor.

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
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

  fHistJetsPt = new TH1F("fHistJetsPt", "fHistJetsPt", fNbins, fMinPt, fMaxPt);
  fHistJetsPt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistJetsPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPt);
  
  if (fAnaType == kEMCAL) {
    fHistJetsPtClus = new TH1F("fHistJetsPtClus", "fHistJetsPtClus", fNbins, fMinPt, fMaxPt);
    fHistJetsPtClus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistJetsPtClus->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsPtClus);
  }
  
  fHistJetsPtTrack = new TH1F("fHistJetsPtTrack", "fHistJetsPtTrack", fNbins, fMinPt, fMaxPt);
  fHistJetsPtTrack->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistJetsPtTrack->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPtTrack);
  
  fHistJetsPtNonBias = new TH1F("fHistJetsPtNonBias", "fHistJetsPtNonBias", fNbins, fMinPt, fMaxPt);
  fHistJetsPtNonBias->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  fHistJetsPtNonBias->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPtNonBias);
  
  fHistLeadingJetPt = new TH1F("fHistLeadingJetPt", "fHistLeadingJetPt", fNbins, fMinPt, fMaxPt);
  fHistLeadingJetPt->GetXaxis()->SetTitle("p_{T} [GeV]");
  fHistLeadingJetPt->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistLeadingJetPt);
  
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
}

//________________________________________________________________________
void AliJetResponseMaker::FillHistograms()
{
  // Fill histograms.

  Int_t maxJetIndex  = -1;
  Int_t maxMCJetIndex  = -1;

  DoJetLoop(maxJetIndex, fJets, fTracks, fCaloClusters); 
  if (maxJetIndex < 0) 
    return;
  AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fJets->At(maxJetIndex));
  if (!jet) 
    return;

  DoJetLoop(maxMCJetIndex, fMCJets, fMCParts);  
  if (maxMCJetIndex < 0) 
    return;
  AliEmcalJet* mcjet = dynamic_cast<AliEmcalJet*>(fMCJets->At(maxMCJetIndex));
  if (!mcjet) 
    return;

  fHistLeadingJetPt->Fill(jet->Pt());
  fHistMCLeadingJetPt->Fill(mcjet->Pt());
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop(Int_t &maxJetIndex, TClonesArray *jets, TClonesArray *tracks, TClonesArray *clusters)
{
  // Do the jet loop.

  if (!fJets)
    return;

  Int_t njets = jets->GetEntriesFast();

  Float_t maxJetPt = 0;
  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(jets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    fHistJetsPtNonBias->Fill(jet->Pt());

    if (jet->MaxTrackPt() > fPtBiasJetTrack)
      fHistJetsPtTrack->Fill(jet->Pt());
    
    if (fAnaType == kEMCAL && jet->MaxClusterPt() > fPtBiasJetClus) 
      fHistJetsPtClus->Fill(jet->Pt());
    
    if (jet->MaxTrackPt() < fPtBiasJetTrack && (fAnaType == kTPC || jet->MaxClusterPt() < fPtBiasJetClus))
	continue;

    if (maxJetPt < jet->Pt()) {
      maxJetPt = jet->Pt();
      maxJetIndex = ij;
    }

    fHistJetsPt->Fill(jet->Pt());

    fHistJetPhiEta->Fill(jet->Eta(), jet->Phi());

    if (fAnaType == kEMCAL)
      fHistJetsNEFvsPt->Fill(jet->NEF(), jet->Pt());

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, tracks);
      if (track)
	fHistJetsZvsPt->Fill(track->Pt() / jet->Pt(), jet->Pt());
    }

    if (fAnaType != kEMCAL || !clusters)
      continue;

    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      AliVCluster *cluster = jet->ClusterAt(ic, clusters);

      if (cluster) {
	TLorentzVector nPart;
	cluster->GetMomentum(nPart, fVertex);
	fHistJetsZvsPt->Fill(nPart.Et() / jet->Pt(), jet->Pt());
      }
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::RetrieveEventObjects()
{
  // Retrieve event objects.

  AliAnalysisTaskEmcal::RetrieveEventObjects();
  
  if (strcmp(fMCJetsName,"")) {
    fMCJets = dynamic_cast<TClonesArray*>(MCEvent()->FindListObject(fMCJetsName));
    if (!fMCJets) {
      AliWarning(Form("Could not retrieve MC jets %s!", fMCJetsName.Data())); 
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
