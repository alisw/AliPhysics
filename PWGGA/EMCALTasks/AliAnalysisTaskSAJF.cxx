// $Id$
//
// Jet finder analysis task (S.Aiola).
//
//

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TRandom3.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliLog.h"

#include "AliAnalysisTaskSAJF.h"

ClassImp(AliAnalysisTaskSAJF)

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF() : 
  AliAnalysisTaskSE("AliAnalysisTaskSAJF"),
  fAnaType(kEMCAL),
  fMinEta(-1),
  fMaxEta(1),
  fMinPhi(0),
  fMaxPhi(2 * TMath::Pi()),
  fJetRadius(0.4),
  fOutput(0),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fKtJetsName("KtJets"),
  fEmbJetsName("EmbJets"),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fKtJets(0),
  fEmbJets(0),
  fTrgClusters(0),
  fCent(0),
  fCentBin(-1),
  fHistCentrality(0),
  fHistJetPhiEta(0),
  fHistRhoPartVSleadJetPt(0),
  fHistMedKtVSRhoPart(0),
  fHistRCPtVSRhoPart(0),
  fNbins(500),
  fMinPt(0),
  fMaxPt(250)
{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistJetsPt[i] = 0;
    fHistJetsNEF[i] = 0;
    fHistJetsZ[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHist2LeadingJetPt[i] = 0;
    fHistTracksPtLJ[i] = 0;
    fHistClusEtLJ[i] = 0;
    fHistTracksPtBkg[i] = 0;
    fHistClusEtBkg[i] = 0;
    fHistBkgClusPhiCorr[i] = 0;
    fHistBkgTracksPhiCorr[i] = 0;
    fHistBkgClusMeanRho[i] = 0;
    fHistBkgTracksMeanRho[i] = 0;
    fHistBkgLJetPhiCorr[i] = 0;
    fHistMedianPtKtJet[i] = 0;
    fHistDeltaPtRC[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistEmbJets[i] = 0;
    fHistEmbPart[i] = 0;
    fHistDeltaPtMedKtEmb[i] = 0;
  }

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF(const char *name) : 
  AliAnalysisTaskSE(name),
  fAnaType(kEMCAL),
  fMinEta(-1),
  fMaxEta(1),
  fMinPhi(0),
  fMaxPhi(2 * TMath::Pi()),
  fJetRadius(0.4),
  fOutput(0),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fJetsName("Jets"),
  fKtJetsName("KtJets"),
  fEmbJetsName("EmbJets"),
  fTrgClusName("ClustersL1GAMMAFEE"),
  fTracks(0),
  fCaloClusters(0),
  fJets(0),
  fKtJets(0),
  fEmbJets(0),
  fTrgClusters(0),
  fCent(0),
  fCentBin(-1),
  fHistCentrality(0),
  fHistJetPhiEta(0),
  fHistRhoPartVSleadJetPt(0),
  fHistMedKtVSRhoPart(0),
  fHistRCPtVSRhoPart(0),
  fNbins(500),
  fMinPt(0),
  fMaxPt(250)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistJetsPt[i] = 0;
    fHistJetsNEF[i] = 0;
    fHistJetsZ[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHist2LeadingJetPt[i] = 0;
    fHistTracksPtLJ[i] = 0;
    fHistClusEtLJ[i] = 0;
    fHistTracksPtBkg[i] = 0;
    fHistClusEtBkg[i] = 0;
    fHistBkgClusPhiCorr[i] = 0;
    fHistBkgTracksPhiCorr[i] = 0;
    fHistBkgClusMeanRho[i] = 0;
    fHistBkgTracksMeanRho[i] = 0;
    fHistBkgLJetPhiCorr[i] = 0;
    fHistMedianPtKtJet[i] = 0;
    fHistDeltaPtRC[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistEmbJets[i] = 0;
    fHistEmbPart[i] = 0;
    fHistDeltaPtMedKtEmb[i] = 0;
  }

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class()); 
}

//________________________________________________________________________
AliAnalysisTaskSAJF::~AliAnalysisTaskSAJF()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserCreateOutputObjects()
{
  // Create histograms
  
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
  
  fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", fNbins, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistJetPhiEta = new TH2F("fHistJetPhiEta","Phi-Eta distribution of jets", 20, -2, 2, 32, 0, 6.4);
  fHistJetPhiEta->GetXaxis()->SetTitle("Eta");
  fHistJetPhiEta->GetYaxis()->SetTitle("Phi");
  fOutput->Add(fHistJetPhiEta);

  fHistRhoPartVSleadJetPt = new TH2F("fHistRhoPartVSleadJetPt","fHistRhoPartVSleadJetPt", fNbins, fMinPt, fMaxPt, fNbins, fMinPt, fMaxPt);
  fHistRhoPartVSleadJetPt->GetXaxis()->SetTitle("#rho [GeV]");
  fHistRhoPartVSleadJetPt->GetYaxis()->SetTitle("Leading jet energy [GeV]");
  fOutput->Add(fHistRhoPartVSleadJetPt);

  fHistMedKtVSRhoPart = new TH2F("fHistMedKtVSRhoPart","fHistMedKtVSRhoPart", fNbins, fMinPt, fMaxPt, fNbins, fMinPt, fMaxPt);
  fHistMedKtVSRhoPart->GetXaxis()->SetTitle("median kt P_{T} [GeV]");
  fHistMedKtVSRhoPart->GetYaxis()->SetTitle("#rho [GeV]");
  fOutput->Add(fHistMedKtVSRhoPart);
  
  fHistRCPtVSRhoPart = new TH2F("fHistRCPtVSRhoPart","fHistRCPtVSRhoPart", fNbins, fMinPt, fMaxPt, fNbins, fMinPt, fMaxPt);
  fHistRCPtVSRhoPart->GetXaxis()->SetTitle("rigid cone P_{T} [GeV]");
  fHistRCPtVSRhoPart->GetYaxis()->SetTitle("#rho [GeV]");
  fOutput->Add(fHistRCPtVSRhoPart);

  TString histname;

  for (Int_t i = 0; i < 4; i++) {
    histname = "fHistJetsPt_";
    histname += i;
    fHistJetsPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistJetsPt[i]->GetXaxis()->SetTitle("E [GeV]");
    fHistJetsPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsPt[i]);
    
    histname = "fHistJetsNEF_";
    histname += i;
    fHistJetsNEF[i] = new TH1F(histname.Data(), histname.Data(), fNbins, 0, 1.2);
    fHistJetsNEF[i]->GetXaxis()->SetTitle("NEF");
    fHistJetsNEF[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsNEF[i]);

    histname = "fHistJetsZ_";
    histname += i;
    fHistJetsZ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, 0, 1.2);
    fHistJetsZ[i]->GetXaxis()->SetTitle("Z");
    fHistJetsZ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsZ[i]);

    histname = "fHistLeadingJetPt_";
    histname += i;
    fHistLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistLeadingJetPt[i]->GetXaxis()->SetTitle("P_{T} [GeV]");
    fHistLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJetPt[i]);
    
    histname = "fHistTracksPtLJ_";
    histname += i;
    fHistTracksPtLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistTracksPtLJ[i]->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    fHistTracksPtLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksPtLJ[i]);
    
    histname = "fHistClusEtLJ_";
    histname += i;
    fHistClusEtLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistClusEtLJ[i]->GetXaxis()->SetTitle("E_{T} [GeV]");
    fHistClusEtLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistClusEtLJ[i]);

    histname = "fHistTracksPtBkg_";
    histname += i;
    fHistTracksPtBkg[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistTracksPtBkg[i]->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    fHistTracksPtBkg[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksPtBkg[i]);
    
    histname = "fHistClusEtBkg_";
    histname += i;
    fHistClusEtBkg[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistClusEtBkg[i]->GetXaxis()->SetTitle("E_{T} [GeV]");
    fHistClusEtBkg[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistClusEtBkg[i]);

    histname = "fHistBkgClusPhiCorr_";
    histname += i;
    fHistBkgClusPhiCorr[i] = new TH1F(histname.Data(), histname.Data(), 128, -1.6, 4.8);
    fHistBkgClusPhiCorr[i]->GetXaxis()->SetTitle("#Delta#phi");
    fHistBkgClusPhiCorr[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgClusPhiCorr[i]);

    histname = "fHistBkgTracksPhiCorr_";
    histname += i;
    fHistBkgTracksPhiCorr[i] = new TH1F(histname.Data(), histname.Data(), 128, -1.6, 4.8);
    fHistBkgTracksPhiCorr[i]->GetXaxis()->SetTitle("#Delta#phi");
    fHistBkgTracksPhiCorr[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgTracksPhiCorr[i]);

    histname = "fHistBkgClusMeanRho_";
    histname += i;
    fHistBkgClusMeanRho[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistBkgClusMeanRho[i]->GetXaxis()->SetTitle("#rho [GeV]");
    fHistBkgClusMeanRho[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgClusMeanRho[i]);

    histname = "fHistBkgTracksMeanRho_";
    histname += i;
    fHistBkgTracksMeanRho[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistBkgTracksMeanRho[i]->GetXaxis()->SetTitle("#rho [GeV/c]");
    fHistBkgTracksMeanRho[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgTracksMeanRho[i]);

    histname = "fHistBkgLJetPhiCorr_";
    histname += i;
    fHistBkgLJetPhiCorr[i] = new TH1F(histname.Data(), histname.Data(), 128, -1.6, 4.8);
    fHistBkgLJetPhiCorr[i]->GetXaxis()->SetTitle("#Delta#phi");
    fHistBkgLJetPhiCorr[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistBkgLJetPhiCorr[i]);

    histname = "fHistMedianPtKtJet_";
    histname += i;
    fHistMedianPtKtJet[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistMedianPtKtJet[i]->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    fHistMedianPtKtJet[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistMedianPtKtJet[i]);

    histname = "fHistDeltaPtRC_";
    histname += i;
    fHistDeltaPtRC[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt - fMaxPt / 2, fMaxPt / 2);
    fHistDeltaPtRC[i]->GetXaxis()->SetTitle("#deltaP_{T} [GeV/c]");
    fHistDeltaPtRC[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRC[i]);

    histname = "fHistDeltaPtRCExLJ_";
    histname += i;
    fHistDeltaPtRCExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt - fMaxPt / 2, fMaxPt / 2);
    fHistDeltaPtRCExLJ[i]->GetXaxis()->SetTitle("#deltaP_{T} [GeV/c]");
    fHistDeltaPtRCExLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRCExLJ[i]);

    histname = "fHistRCPt_";
    histname += i;
    fHistRCPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistRCPt[i]->GetXaxis()->SetTitle("rigid cone P_{T} [GeV/c]");
    fHistRCPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPt[i]);

    histname = "fHistRCPtExLJ_";
    histname += i;
    fHistRCPtExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistRCPtExLJ[i]->GetXaxis()->SetTitle("rigid cone P_{T} [GeV/c]");
    fHistRCPtExLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPtExLJ[i]);

    histname = "fHistEmbJets_";
    histname += i;
    fHistEmbJets[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistEmbJets[i]->GetXaxis()->SetTitle("embedded jet P_{T} [GeV/c]");
    fHistEmbJets[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEmbJets[i]);

    histname = "fHistEmbPart_";
    histname += i;
    fHistEmbPart[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt, fMaxPt);
    fHistEmbPart[i]->GetXaxis()->SetTitle("embedded particle P_{T} [GeV/c]");
    fHistEmbPart[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEmbPart[i]);

    histname = "fHistDeltaPtMedKtEmb_";
    histname += i;
    fHistDeltaPtMedKtEmb[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinPt - fMaxPt / 2, fMaxPt / 2);
    fHistDeltaPtMedKtEmb[i]->GetXaxis()->SetTitle("#deltaP_{T} [GeV/c]");
    fHistDeltaPtMedKtEmb[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtMedKtEmb[i]);
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::RetrieveEventObjects()
{
  fCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
  if (!fCaloClusters) {
    AliWarning(Form("Could not retrieve clusters!")); 
  }

  fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
  if (!fTracks) {
    AliWarning(Form("Could not retrieve tracks!")); 
  }

  fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
  if (!fJets) {
    AliWarning(Form("Could not retrieve jets!")); 
  }

  if (strcmp(fKtJetsName,"")) {
    fKtJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fKtJetsName));
    if (!fKtJets) {
      AliWarning(Form("Could not retrieve Kt jets!")); 
    }
  }

  if (strcmp(fEmbJetsName,"")) {
    fEmbJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fEmbJetsName));
    if (!fEmbJets) {
      AliWarning(Form("Could not retrieve emb jets!")); 
    }
  }
    
  if (strcmp(fTrgClusName,"")) {
    fTrgClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrgClusName));
    if (!fTrgClusters) {
      AliWarning(Form("Could not retrieve trigger clusters!")); 
    }
  }

  fCent = InputEvent()->GetCentrality();
  if (fCent) {
    Float_t cent = fCent->GetCentralityPercentile("V0M");
    if      (cent >=  0 && cent <   10) fCentBin = 0;
    else if (cent >= 10 && cent <   30) fCentBin = 1;
    else if (cent >= 30 && cent <   50) fCentBin = 2;
    else if (cent >= 50 && cent <= 100) fCentBin = 3; 
    else {
      AliWarning(Form("Negative centrality: %f. Assuming 99", cent));
      fCentBin = 3;
    }
  }
  else {
    AliWarning(Form("Could not retrieve centrality information! Assuming 99"));
    fCentBin = 3;
  }
}

//________________________________________________________________________
AliVTrack* AliAnalysisTaskSAJF::GetTrack(const Int_t i) const
{
  if (fTracks)
    return dynamic_cast<AliVTrack*>(fTracks->At(i));
  else
    return 0;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAJF::GetNumberOfTracks() const
{
  if (fTracks)
    return fTracks->GetEntriesFast();
  else
    return 0;
}

//________________________________________________________________________
AliVCluster* AliAnalysisTaskSAJF::GetCaloCluster(const Int_t i) const
{ 
  if (fCaloClusters)
    return dynamic_cast<AliVCluster*>(fCaloClusters->At(i));
  else
    return 0;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAJF::GetNumberOfCaloClusters() const
{ 
  if (fCaloClusters)
    return fCaloClusters->GetEntriesFast();
  else
    return 0;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskSAJF::GetJet(const Int_t i) const
{
  if (fJets)
    return dynamic_cast<AliEmcalJet*>(fJets->At(i));
  else
    return 0;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAJF::GetNumberOfJets() const
{
  if (fJets)
    return fJets->GetEntriesFast();
  else
    return 0;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskSAJF::GetKtJet(const Int_t i) const
{
  if (fKtJets)
    return dynamic_cast<AliEmcalJet*>(fKtJets->At(i));
  else
    return 0;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAJF::GetNumberOfKtJets() const
{
  if (fKtJets)
    return fKtJets->GetEntriesFast();
  else
    return 0;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskSAJF::GetEmbJet(const Int_t i) const
{
  if (fEmbJets)
    return dynamic_cast<AliEmcalJet*>(fEmbJets->At(i));
  else
    return 0;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAJF::GetNumberOfEmbJets() const
{
  if (fEmbJets)
    return fEmbJets->GetEntriesFast();
  else
    return 0;
}

//________________________________________________________________________
AliVCluster* AliAnalysisTaskSAJF::GetTrgCluster(const Int_t i) const
{
  if (fTrgClusters)
    return dynamic_cast<AliVCluster*>(fTrgClusters->At(i));
  else
    return 0;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAJF::GetNumberOfTrgClusters() const
{
  if (fTrgClusters)
    return fTrgClusters->GetEntriesFast();
  else
    return 0;
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::FillHistograms()
{
  Float_t A = fJetRadius * fJetRadius * TMath::Pi();

  Float_t cent = 100;

  if (fCent)
    cent = fCent->GetCentralityPercentile("V0M");

  fHistCentrality->Fill(cent);

  Int_t maxJetIndex  = -1;
  Int_t max2JetIndex = -1;

  DoJetLoop(maxJetIndex, max2JetIndex);
  
  if (maxJetIndex < 0) 
    return;

  AliEmcalJet* jet = GetJet(maxJetIndex);
  if (!jet) 
    return;

  fHistLeadingJetPt[fCentBin]->Fill(jet->Pt());
  jet->SortConstituents();
  
  AliEmcalJet* jet2 = 0;
  if (max2JetIndex >= 0)
    jet2 = GetJet(max2JetIndex);

  if (jet2) {
    //fHistLeadingJetPt[fCentBin]->Fill(jet2->Pt());
    jet2->SortConstituents();
  }

  Float_t rhoKt = DoKtJetLoop();
  fHistMedianPtKtJet[fCentBin]->Fill(rhoKt);
  
  Float_t rhoTracks = DoTrackLoop(maxJetIndex, max2JetIndex);
  fHistBkgTracksMeanRho[fCentBin]->Fill(rhoTracks);

  Float_t rhoClus = 0;
  if (fAnaType == kEMCAL || fAnaType == kEMCALFiducial) 
    rhoClus = DoClusterLoop(maxJetIndex, max2JetIndex);
  
  fHistBkgClusMeanRho[fCentBin]->Fill(rhoClus);

  fHistRhoPartVSleadJetPt->Fill(A * (rhoClus + rhoTracks), jet->Pt());

  fHistMedKtVSRhoPart->Fill(rhoKt, rhoClus + rhoTracks);
  
  Int_t nRCs = GetArea() / A - 1;

  for (Int_t i = 0; i < nRCs; i++) {
    Float_t RCpt = GetRigidConePt(0);
    Float_t RCptExLJ = GetRigidConePt(jet);
    
    fHistDeltaPtRC[fCentBin]->Fill(RCpt - A * rhoKt);
    fHistDeltaPtRCExLJ[fCentBin]->Fill(RCptExLJ - A * rhoKt);

    fHistRCPt[fCentBin]->Fill(RCpt / A);
    fHistRCPtExLJ[fCentBin]->Fill(RCptExLJ / A);
    fHistRCPtVSRhoPart->Fill(RCptExLJ / A, rhoClus + rhoTracks);
  }

  Float_t maxEmbJetPt = 0;
  Float_t maxEmbPartPt = 0;

  Bool_t embOK = DoEmbJetLoop(maxEmbJetPt, maxEmbPartPt);

  if (embOK) {
    fHistEmbJets[fCentBin]->Fill(maxEmbJetPt);
    fHistEmbPart[fCentBin]->Fill(maxEmbPartPt);
    fHistDeltaPtMedKtEmb[fCentBin]->Fill(maxEmbJetPt - A * rhoKt - maxEmbPartPt);
  }
  else {
    if (maxEmbPartPt != 0)
      AliWarning("Embedded particle is not the leading particle of the leading jet!");
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::DoJetLoop(Int_t &maxJetIndex, Int_t &max2JetIndex)
{
  Double_t vertex[3] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  Int_t njets = GetNumberOfJets();

  Float_t maxJetPt = 0;
  Float_t max2JetPt = 0;
  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = GetJet(ij);

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  
    
    if (jet->Pt() <= 0)
      continue;

    if (!AcceptJet(jet))
      continue;

    fHistJetPhiEta->Fill(jet->Eta(), jet->Phi());
    fHistJetsPt[fCentBin]->Fill(jet->Pt());
    fHistJetsNEF[fCentBin]->Fill(jet->NEF());

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVTrack *track = jet->TrackAt(it, fTracks);
      if (track)
	fHistJetsZ[fCentBin]->Fill(track->Pt() / jet->Pt());
    }

    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);

      if (cluster) {
	TLorentzVector nPart;
	cluster->GetMomentum(nPart, vertex);
	fHistJetsZ[fCentBin]->Fill(nPart.Et() / jet->Pt());
      }
    }

    if (maxJetPt < jet->Pt()) {
      max2JetPt = maxJetPt;
      max2JetIndex = maxJetIndex;
      maxJetPt = jet->Pt();
      maxJetIndex = ij;
    }
    else if (max2JetPt < jet->Pt()) {
      max2JetPt = jet->Pt();
      max2JetIndex = ij;
    }
  } //jet loop 
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAJF::DoKtJetLoop(Int_t nLJs)
{
  Float_t ktJetsMedian = 0;
  Int_t nktjets =  GetNumberOfKtJets();

  Int_t NoOfZeroJets = 0;
  if (nktjets > 0) {
    
    TArrayF ktJets(nktjets);
    for (Int_t ij = 0; ij < nktjets; ij++) {
      
      AliEmcalJet* jet = GetKtJet(ij);
      
      if (!jet) {
	AliError(Form("Could not receive jet %d", ij));
	continue;
      } 
      
      if (jet->Pt() <= 0) {
	NoOfZeroJets++;
	continue;
      }
      
      if (!AcceptJet(jet)) {
	NoOfZeroJets++;
	continue;
      }

      Float_t rho = jet->Pt() / jet->Area();
      Int_t i = nktjets - 1;
      while (rho < ktJets[i] && i > 0)
	i--;
      memmove(ktJets.GetArray() + nktjets - ij - 1, ktJets.GetArray() + nktjets - ij, (ij + i - nktjets + 1) * sizeof(Float_t));
      ktJets[i] = rho;
    } //kt jet loop 

    nktjets -= NoOfZeroJets;
    
    if (nktjets < 1) return 0;

    memmove(ktJets.GetArray(), ktJets.GetArray() + NoOfZeroJets, nktjets * sizeof(Float_t));

    nktjets -= nLJs;

    if (nktjets < 1) return 0;

    if (nktjets % 2)
      ktJetsMedian = ktJets[nktjets / 2];
    else
      ktJetsMedian = (ktJets[nktjets / 2] + ktJets[nktjets / 2 - 1]) / 2;
  }

  return ktJetsMedian;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::DoEmbJetLoop(Float_t &maxJetPt, Float_t &maxPartPt)
{
  Double_t vertex[3] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  Int_t nembjets = GetNumberOfEmbJets();

  Int_t maxJetIdx = -1;
  Int_t maxTrackIdx = -1;
  Int_t maxClusIdx = -1;

  Float_t maxTrackPt = 0;
  Float_t maxClusEt = 0;

  for (Int_t ij = 0; ij < nembjets; ij++) {
      
    AliEmcalJet* jet = GetEmbJet(ij);
      
    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    } 
      
    if (jet->Pt() <= 0)
	continue;
 
    if (!AcceptJet(jet))
      continue;
    
    if (jet->Pt() > maxJetPt) {
      maxJetPt = jet->Pt();
      maxJetIdx = ij;
    }
  }

  if (maxJetPt <= 0)
    return kFALSE;

  AliEmcalJet* jet = GetEmbJet(maxJetIdx);

  for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
    AliVTrack *track = jet->TrackAt(it, fTracks);

    if (!track) continue;
     
    if (track->Pt() > maxTrackPt) {
      maxTrackPt = track->Pt();
      maxTrackIdx = it;
    }
  }

  for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
    AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);
    
    if (!cluster) continue;
    
    TLorentzVector nPart;
    cluster->GetMomentum(nPart, vertex);

    if (nPart.Et() > maxClusEt) {
      maxClusEt = nPart.Et();
      maxClusIdx = ic;
    } 
  }

  if (maxClusEt > maxTrackPt) {
    AliVCluster *cluster = jet->ClusterAt(maxClusIdx, fCaloClusters);
    maxPartPt = maxClusEt;
    return (Bool_t)(cluster->Chi2() == 100);
  }
  else {
    AliVTrack *track = jet->TrackAt(maxTrackIdx, fTracks);
    maxPartPt = maxTrackPt;
    return (Bool_t)(track->GetLabel() == 100);
  }
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAJF::DoTrackLoop(Int_t maxJetIndex, Int_t max2JetIndex)
{ 
  AliEmcalJet* jet = 0;
  if (max2JetIndex >= 0)
    jet = GetJet(maxJetIndex);

  AliEmcalJet* jet2 = 0;
  if (max2JetIndex >= 0)
    jet2 = GetJet(max2JetIndex);

  Float_t rhoTracks = 0;
  Int_t ntracks = GetNumberOfTracks();

  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliVTrack* track = GetTrack(iTracks);         
    if(!track) {
      AliError(Form("Could not retrieve track %d",iTracks)); 
      continue; 
    }
    
    if (!AcceptTrack(track)) continue;
    
    if (jet && IsJetTrack(jet, iTracks)) {
      fHistTracksPtLJ[fCentBin]->Fill(track->Pt());
    }
    else if (!jet2 || !IsJetTrack(jet2, iTracks)) {
      fHistTracksPtBkg[fCentBin]->Fill(track->Pt());
      rhoTracks += track->Pt();
      
      if (jet) {
	Float_t dphijet = jet->Phi() - track->Phi();
	if (dphijet < -1.6) dphijet += TMath::Pi() * 2;
	if (dphijet > 4.8) dphijet -= TMath::Pi() * 2;
	fHistBkgLJetPhiCorr[fCentBin]->Fill(dphijet);
      }

      for(Int_t it2 = iTracks+1; it2 < ntracks; it2++) {
	AliVTrack* track2 = GetTrack(it2);         
	if(!track2) {
	  AliError(Form("Could not retrieve track %d", it2)); 
	  continue; 
	}
	
	if (!AcceptTrack(track2)) continue;
	
	if (jet && IsJetTrack(jet, it2)) continue;

	if (jet2 && IsJetTrack(jet2, it2)) continue;
	
	Float_t dphi = track->Phi() - track2->Phi();
	if (dphi < -1.6) dphi += TMath::Pi() * 2;
	if (dphi > 4.8) dphi -= TMath::Pi() * 2;
	fHistBkgTracksPhiCorr[fCentBin]->Fill(dphi);
      } // second track loop
    }
  }
  rhoTracks /= GetArea();
  return rhoTracks;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAJF::DoClusterLoop(Int_t maxJetIndex, Int_t max2JetIndex)
{
  Double_t vertex[3] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  AliEmcalJet* jet = 0;
  if (max2JetIndex >= 0)
    jet = GetJet(maxJetIndex);

  AliEmcalJet* jet2 = 0;
  if (max2JetIndex >= 0)
    jet2 = GetJet(max2JetIndex);

  Float_t rhoClus = 0;
  Int_t nclusters =  GetNumberOfCaloClusters();
  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = GetCaloCluster(iClusters);
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;
    
    if (!AcceptCluster(cluster)) continue;

    TLorentzVector nPart;
    cluster->GetMomentum(nPart, vertex);

    if (jet && IsJetCluster(jet, iClusters)) {
      fHistClusEtLJ[fCentBin]->Fill(nPart.Et());
    }
    else if (!jet2 || !IsJetCluster(jet2, iClusters)) {
      fHistClusEtBkg[fCentBin]->Fill(nPart.Et());
      rhoClus += nPart.Et();

      Float_t pos1[3];
      cluster->GetPosition(pos1);
      TVector3 clusVec1(pos1);

      if (jet) {
	Float_t dphijet = jet->Phi() - clusVec1.Phi();
	if (dphijet < -1.6) dphijet += TMath::Pi() * 2;
	if (dphijet > 4.8) dphijet -= TMath::Pi() * 2;
	fHistBkgLJetPhiCorr[fCentBin]->Fill(dphijet);
      }

      for(Int_t ic2 = iClusters+1; ic2 < nclusters; ic2++) {
	AliVCluster* cluster2 = GetCaloCluster(ic2);
	if (!cluster2) {
	  AliError(Form("Could not receive cluster %d", ic2));
	  continue;
	}  
	
	if (!(cluster2->IsEMCAL())) continue;

	if (!AcceptCluster(cluster)) continue;
	
	if (jet && IsJetCluster(jet, ic2)) continue;
	
	if (jet2 && IsJetCluster(jet2, ic2)) continue;

	Float_t pos2[3];
	cluster2->GetPosition(pos2);
	TVector3 clusVec2(pos2);

	Float_t dphi = clusVec1.Phi() - clusVec2.Phi();
	if (dphi < -1.6) dphi += TMath::Pi() * 2;
	if (dphi > 4.8) dphi -= TMath::Pi() * 2;
	fHistBkgClusPhiCorr[fCentBin]->Fill(dphi);
      }
    }
  } //cluster loop 
  rhoClus /= GetArea();
  return rhoClus;
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::Init()
{
  if (fAnaType == kFullAcceptance) {
    fMinEta = -1;
    fMaxEta = 1;
    fMinPhi = 0;
    fMaxPhi = 2 * TMath::Pi();
  }
  else if (fAnaType == kEMCAL || fAnaType == kEMCALFiducial) {
    fMinEta = -0.7;
    fMaxEta = 0.7;
    fMinPhi = 80 * TMath::DegToRad();
    fMaxPhi = 180 * TMath::DegToRad();
  }
  else {
    AliWarning("Analysis type not recognized! Assuming kFullAcceptance...");
    fMinEta = -1;
    fMaxEta = 1;
    fMinPhi = 0;
    fMaxPhi = 2 * TMath::Pi();
  }
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAJF::GetRigidConePt(AliEmcalJet *jet, Float_t minD)
{
  static TRandom3 random;

  Double_t vertex[3] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  Float_t eta = 0;
  Float_t phi = 0;

  Float_t LJeta = 999;
  Float_t LJphi = 999;

  if (jet) {
    LJeta = jet->Eta();
    LJphi = jet->Phi();
  }

  Float_t dLJ = 0;
  do {
    eta = random.Rndm() * (fMaxEta - fMinEta) + fMinEta;
    phi = random.Rndm() * (fMaxPhi - fMinPhi) + fMinPhi;
    dLJ = TMath::Sqrt((LJeta - eta) * (LJeta - eta) + (LJphi - phi) * (LJphi - phi));
  } while (dLJ < minD && !AcceptJet(eta, phi));
  
  Float_t pt = 0;

  Int_t nclusters =  GetNumberOfCaloClusters();
  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = GetCaloCluster(iClusters);
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;

    if (!AcceptCluster(cluster)) continue;
    
    TLorentzVector nPart;
    cluster->GetMomentum(nPart, vertex);

    Float_t pos[3];
    cluster->GetPosition(pos);
    TVector3 clusVec(pos);

    Float_t d = TMath::Sqrt((clusVec.Eta() - eta) * (clusVec.Eta() - eta) + (clusVec.Phi() - phi) * (clusVec.Phi() - phi));

    if (d <= fJetRadius)
      pt += nPart.Pt();
  }

  Int_t ntracks = GetNumberOfTracks();
  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliVTrack* track = GetTrack(iTracks);         
    if(!track) {
      AliError(Form("Could not retrieve track %d",iTracks)); 
      continue; 
    }
    
    if (!AcceptTrack(track)) continue;

    Float_t tracketa = track->Eta();
    Float_t trackphi = track->Phi();
    
    if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi + 2 * TMath::Pi()))
      trackphi += 2 * TMath::Pi();
    if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi - 2 * TMath::Pi()))
      trackphi -= 2 * TMath::Pi();

    Float_t d = TMath::Sqrt((tracketa - eta) * (tracketa - eta) + (trackphi - phi) * (trackphi - phi));

    if (d <= fJetRadius)
      pt += track->Pt();
  }

  return pt;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSAJF::GetArea() const
{
  return (fMaxEta - fMinEta) * (fMaxPhi - fMinPhi);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted) const
{
  for (Int_t i = 0; i < jet->GetNumberOfTracks(); i++) {
    Int_t ijettrack = jet->TrackAt(i);
    if (sorted && ijettrack > itrack)
      return kFALSE;
    if (ijettrack == itrack)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted) const
{
  for (Int_t i = 0; i < jet->GetNumberOfClusters(); i++) {
    Int_t ijetclus = jet->ClusterAt(i);
    if (sorted && ijetclus > iclus)
      return kFALSE;
    if (ijetclus == iclus)
      return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::AcceptJet(Float_t eta, Float_t phi) const
{
  if (fAnaType == kFullAcceptance) {
    return kTRUE;
  }
  else if (fAnaType == kEMCAL) {
    return (Bool_t)(eta > fMinEta && eta < fMaxEta && phi > fMinPhi && phi < fMaxPhi);
  }
  else if (fAnaType == kEMCALFiducial) {
    return (Bool_t)(eta > fMinEta + fJetRadius && eta < fMaxEta - fJetRadius && phi > fMinPhi + fJetRadius && phi < fMaxPhi - fJetRadius);
  }
  else {
    AliWarning("Analysis type not recognized! Assuming kFullAcceptance...");
    return kTRUE;
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::AcceptJet(AliEmcalJet* jet) const
{
  return AcceptJet(jet->Eta(), jet->Phi());
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::AcceptCluster(AliVCluster* clus, Bool_t acceptMC) const
{
  if (!acceptMC && clus->Chi2() == 100)
    return kFALSE;

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::AcceptTrack(AliVTrack* track, Bool_t acceptMC) const
{
  if (!acceptMC && track->GetLabel() == 100)
    return kFALSE;
  
  if (fAnaType == kFullAcceptance) {
    return kTRUE;
  }
  else if (fAnaType == kEMCAL || fAnaType == kEMCALFiducial) {
    return (Bool_t)(track->Eta() > fMinEta && track->Eta() < fMaxEta && track->Phi() > fMinPhi && track->Phi() < fMaxPhi);
  }
  else {
    AliWarning("Analysis type not recognized! Assuming kFullAcceptance...");
    return kTRUE;
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserExec(Option_t *) 
{
  Init();

  RetrieveEventObjects();

  FillHistograms();
    
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
