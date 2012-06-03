// $Id$
//
// Jet analysis task (S.Aiola).
//
//

#include <TObject.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TParameter.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliVEventHandler.h"
#include "AliLog.h"

#include "AliAnalysisTaskSAJF.h"

ClassImp(AliAnalysisTaskSAJF)

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSAJF"),
  fMinRC2LJ(1.0),
  fEmbJetsName("EmbJets"),
  fRandTracksName("TracksRandomized"),
  fRandCaloName("CaloClustersRandomized"),
  fRhoName("Rho"),
  fSkipEventsWithNoJets(kTRUE),
  fEmbJets(0),
  fRandTracks(0),
  fRandCaloClusters(0),
  fRho(0),
  fHistCentrality(0),
  fHistRejectedEvents(0),
  fHistRhoVSleadJetPt(0),
  fHistRCPhiEta(0),
  fHistRCPtExLJVSDPhiLJ(0),
  fHistRhoVSRCPt(0),
  fHistEmbJetPhiEta(0),
  fHistEmbPartPhiEta(0),
  fHistRhoVSEmbBkg(0)

{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistJetPhiEta[i] = 0;
    fHistJetsPt[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHist2LeadingJetPt[i] = 0;
    fHistJetsNEFvsPt[i] = 0;
    fHistJetsZvsPt[i] = 0;
    fHistTracksPtLJ[i] = 0;
    fHistClusEtLJ[i] = 0;
    fHistTracksPtBkg[i] = 0;
    fHistClusEtBkg[i] = 0;
    fHistRho[i] = 0;
    fHistCorrJetsPt[i] = 0;
    fHistCorrLeadingJetPt[i] = 0;
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistRCPtRand[i] = 0;
    fHistDeltaPtRC[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistDeltaPtRCRand[i] = 0;
    fHistEmbJets[i] = 0;
    fHistEmbPart[i] = 0;
    fHistDeltaPtEmb[i] = 0;
  }

}

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF(const char *name) : 
  AliAnalysisTaskEmcalJet(name),
  fMinRC2LJ(1.0),
  fEmbJetsName("EmbJets"),
  fRandTracksName("TracksRandomized"),
  fRandCaloName("CaloClustersRandomized"),
  fRhoName("Rho"),
  fSkipEventsWithNoJets(kTRUE),
  fEmbJets(0),
  fRandTracks(0),
  fRandCaloClusters(0),
  fRho(0),
  fHistCentrality(0),
  fHistRejectedEvents(0),
  fHistRhoVSleadJetPt(0),
  fHistRCPhiEta(0),
  fHistRCPtExLJVSDPhiLJ(0),
  fHistRhoVSRCPt(0),
  fHistEmbJetPhiEta(0),
  fHistEmbPartPhiEta(0),
  fHistRhoVSEmbBkg(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistJetPhiEta[i] = 0;
    fHistJetsPt[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHist2LeadingJetPt[i] = 0;
    fHistJetsNEFvsPt[i] = 0;
    fHistJetsZvsPt[i] = 0;
    fHistTracksPtLJ[i] = 0;
    fHistClusEtLJ[i] = 0;
    fHistTracksPtBkg[i] = 0;
    fHistClusEtBkg[i] = 0;
    fHistRho[i] = 0;
    fHistCorrJetsPt[i] = 0;
    fHistCorrLeadingJetPt[i] = 0;
    fHistRCPt[i] = 0;
    fHistRCPtExLJ[i] = 0;
    fHistRCPtRand[i] = 0;
    fHistDeltaPtRC[i] = 0;
    fHistDeltaPtRCExLJ[i] = 0;
    fHistDeltaPtRCRand[i] = 0;
    fHistEmbJets[i] = 0;
    fHistEmbPart[i] = 0;
    fHistDeltaPtEmb[i] = 0;
  }

}

//________________________________________________________________________
AliAnalysisTaskSAJF::~AliAnalysisTaskSAJF()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Create histograms

  const Float_t binWidth = (fMaxBinPt - fMinBinPt) / fNbins;
  
  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner(); 
  
  fHistCentrality = new TH1F("fHistCentrality","fHistCentrality", fNbins, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistRejectedEvents = new TH1F("fHistRejectedEvents","fHistRejectedEvents", 6, 0, 6);
  fHistRejectedEvents->GetXaxis()->SetTitle("Reason");
  fHistRejectedEvents->GetYaxis()->SetTitle("counts");
  fHistRejectedEvents->GetXaxis()->SetBinLabel(1, "Rho <= 0");
  fHistRejectedEvents->GetXaxis()->SetBinLabel(2, "Max Jet <= 0");
  fHistRejectedEvents->GetXaxis()->SetBinLabel(3, "Max Jet not found");
  fHistRejectedEvents->GetXaxis()->SetBinLabel(4, "No jets");
  fOutput->Add(fHistRejectedEvents);

  fHistRhoVSleadJetPt = new TH2F("fHistRhoVSleadJetPt","fHistRhoVSleadJetPt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistRhoVSleadJetPt->GetXaxis()->SetTitle("#rho * area [GeV/c]");
  fHistRhoVSleadJetPt->GetYaxis()->SetTitle("Leading jet p_{T} [GeV/c]");
  fOutput->Add(fHistRhoVSleadJetPt);

  fHistRCPhiEta = new TH2F("fHistRCPhiEta","Phi-Eta distribution of rigid cones", 40, -2, 2, 64, 0, 6.4);
  fHistRCPhiEta->GetXaxis()->SetTitle("#eta");
  fHistRCPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistRCPhiEta);

  fHistRCPtExLJVSDPhiLJ = new TH2F("fHistRCPtExLJVSDPhiLJ","fHistRCPtExLJVSDPhiLJ", fNbins, fMinBinPt, fMaxBinPt, 128, -1.6, 4.8);
  fHistRCPtExLJVSDPhiLJ->GetXaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
  fHistRCPtExLJVSDPhiLJ->GetYaxis()->SetTitle("#Delta#phi");
  fOutput->Add(fHistRCPtExLJVSDPhiLJ);

  fHistRhoVSRCPt = new TH2F("fHistRhoVSRCPt","fHistRhoVSRCPt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistRhoVSRCPt->GetXaxis()->SetTitle("#rho * area [GeV/c]");
  fHistRhoVSRCPt->GetYaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
  fOutput->Add(fHistRhoVSRCPt);

  fHistEmbJetPhiEta = new TH2F("fHistEmbJetPhiEta","Phi-Eta distribution of embedded jets", 40, -2, 2, 64, 0, 6.4);
  fHistEmbJetPhiEta->GetXaxis()->SetTitle("#eta");
  fHistEmbJetPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistEmbJetPhiEta);

  fHistEmbPartPhiEta = new TH2F("fHistEmbPartPhiEta","Phi-Eta distribution of embedded particles", 40, -2, 2, 64, 0, 6.4);
  fHistEmbPartPhiEta->GetXaxis()->SetTitle("#eta");
  fHistEmbPartPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistEmbPartPhiEta);

  fHistRhoVSEmbBkg = new TH2F("fHistRhoVSEmbBkg","fHistRhoVSEmbBkg", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistRhoVSEmbBkg->GetXaxis()->SetTitle("rho * area [GeV/c]");
  fHistRhoVSEmbBkg->GetYaxis()->SetTitle("background of embedded track [GeV/c]");
  fOutput->Add(fHistRhoVSEmbBkg);

  TString histname;

  for (Int_t i = 0; i < 4; i++) {
    histname = "fHistJetPhiEta_";
    histname += i;
    fHistJetPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 40, -2, 2, 64, 0, 6.4);
    fHistJetPhiEta[i]->GetXaxis()->SetTitle("#eta");
    fHistJetPhiEta[i]->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistJetPhiEta[i]);

    histname = "fHistJetsPt_";
    histname += i;
    fHistJetsPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
    fHistJetsPt[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistJetsPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistJetsPt[i]);

    histname = "fHistJetsPtArea_";
    histname += i;
    fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 20, 0, fJetRadius * fJetRadius * TMath::Pi() * 1.5);
    fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
    fOutput->Add(fHistJetsPtArea[i]);

    histname = "fHistLeadingJetPt_";
    histname += i;
    fHistLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
    fHistLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
    fHistLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistLeadingJetPt[i]);

    histname = "fHist2LeadingJetPt_";
    histname += i;
    fHist2LeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
    fHist2LeadingJetPt[i]->GetXaxis()->SetTitle("p_{T} [GeV]");
    fHist2LeadingJetPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHist2LeadingJetPt[i]);

    histname = "fHistJetsZvsPt_";
    histname += i;
    fHistJetsZvsPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetsZvsPt[i]->GetXaxis()->SetTitle("Z");
    fHistJetsZvsPt[i]->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    fOutput->Add(fHistJetsZvsPt[i]);

    if (fAnaType == kEMCAL) {
      histname = "fHistJetsNEFvsPt_";
      histname += i;
      fHistJetsNEFvsPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
      fHistJetsNEFvsPt[i]->GetXaxis()->SetTitle("NEF");
      fHistJetsNEFvsPt[i]->GetYaxis()->SetTitle("p_{T} [GeV/c]");
      fOutput->Add(fHistJetsNEFvsPt[i]);

      histname = "fHistClusEtLJ_";
      histname += i;
      fHistClusEtLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt / 5);
      fHistClusEtLJ[i]->GetXaxis()->SetTitle("E_{T} [GeV]");
      fHistClusEtLJ[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistClusEtLJ[i]);
      
      histname = "fHistClusEtBkg_";
      histname += i;
      fHistClusEtBkg[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt / 5);
      fHistClusEtBkg[i]->GetXaxis()->SetTitle("E_{T} [GeV]");
      fHistClusEtBkg[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistClusEtBkg[i]);
    }

    histname = "fHistTracksPtLJ_";
    histname += i;
    fHistTracksPtLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt / 5);
    fHistTracksPtLJ[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistTracksPtLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksPtLJ[i]);

    histname = "fHistTracksPtBkg_";
    histname += i;
    fHistTracksPtBkg[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt / 5);
    fHistTracksPtBkg[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistTracksPtBkg[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksPtBkg[i]);

    histname = "fHistRho_";
    histname += i;
    fHistRho[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRho[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistRho[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRho[i]);

    histname = "fHistCorrJetsPt_";
    histname += i;
    fHistCorrJetsPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistCorrJetsPt[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistCorrJetsPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCorrJetsPt[i]);

    histname = "fHistCorrLeadingJetPt_";
    histname += i;
    fHistCorrLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistCorrLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    fHistCorrLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCorrLeadingJetPt[i]);
    
    histname = "fHistRCPt_";
    histname += i;
    fHistRCPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPt[i]->GetXaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
    fHistRCPt[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPt[i]);

    histname = "fHistRCPtExLJ_";
    histname += i;
    fHistRCPtExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPtExLJ[i]->GetXaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
    fHistRCPtExLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPtExLJ[i]);

    histname = "fHistRCPtRand_";
    histname += i;
    fHistRCPtRand[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt * 2);
    fHistRCPtRand[i]->GetXaxis()->SetTitle("rigid cone p_{T} [GeV/c]");
    fHistRCPtRand[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistRCPtRand[i]);

    histname = "fHistDeltaPtRC_";
    histname += i;
    fHistDeltaPtRC[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt - fMaxBinPt / 2 + binWidth / 2, fMinBinPt + fMaxBinPt / 2 + binWidth / 2);
    fHistDeltaPtRC[i]->GetXaxis()->SetTitle("#deltap_{T} [GeV/c]");
    fHistDeltaPtRC[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRC[i]);

    histname = "fHistDeltaPtRCExLJ_";
    histname += i;
    fHistDeltaPtRCExLJ[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt - fMaxBinPt / 2 + binWidth / 2, fMinBinPt + fMaxBinPt / 2 + binWidth / 2);
    fHistDeltaPtRCExLJ[i]->GetXaxis()->SetTitle("#deltap_{T} [GeV/c]");
    fHistDeltaPtRCExLJ[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRCExLJ[i]);

    histname = "fHistDeltaPtRCRand_";
    histname += i;
    fHistDeltaPtRCRand[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt - fMaxBinPt / 2 + binWidth / 2, fMinBinPt + fMaxBinPt / 2 + binWidth / 2);
    fHistDeltaPtRCRand[i]->GetXaxis()->SetTitle("#deltap_{T} [GeV/c]");
    fHistDeltaPtRCRand[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtRCRand[i]);

    histname = "fHistEmbJets_";
    histname += i;
    fHistEmbJets[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
    fHistEmbJets[i]->GetXaxis()->SetTitle("embedded jet p_{T} [GeV/c]");
    fHistEmbJets[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEmbJets[i]);

    histname = "fHistEmbPart_";
    histname += i;
    fHistEmbPart[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
    fHistEmbPart[i]->GetXaxis()->SetTitle("embedded particle p_{T} [GeV/c]");
    fHistEmbPart[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEmbPart[i]);

    histname = "fHistDeltaPtEmb_";
    histname += i;
    fHistDeltaPtEmb[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt - fMaxBinPt / 2 + binWidth / 2, fMinBinPt + fMaxBinPt / 2 + binWidth / 2);
    fHistDeltaPtEmb[i]->GetXaxis()->SetTitle("#deltap_{T} [GeV/c]");
    fHistDeltaPtEmb[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtEmb[i]);
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::RetrieveEventObjects()
{
  if(!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  fRho = -1;
  
  if (strcmp(fRhoName,"")) {
    TParameter<Double_t> *rho = dynamic_cast<TParameter<Double_t>*>(InputEvent()->FindListObject(fRhoName));
  
    if (rho) {
      fRho = rho->GetVal();
    }
    else {
      AliWarning(Form("Could not retrieve rho %s!", fRhoName.Data()));
    }
  }
  
  if (strcmp(fEmbJetsName,"")) {
    fEmbJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fEmbJetsName));
    if (!fEmbJets) {
      AliWarning(Form("Could not retrieve emb jets %s!", fEmbJetsName.Data())); 
    }
  }

  if (strcmp(fRandCaloName,"") && fAnaType == kEMCAL) {
    fRandCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fRandCaloName));
    if (!fRandCaloClusters) {
      AliWarning(Form("Could not retrieve randomized clusters %s!", fRandCaloName.Data())); 
    }
  }

  if (strcmp(fRandTracksName,"")) {
    fRandTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fRandTracksName));
    if (!fRandTracks) {
      AliWarning(Form("Could not retrieve randomized tracks %s!", fRandTracksName.Data())); 
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::FillHistograms()
{
  if (fRho < 0) {
    fHistRejectedEvents->Fill("Rho <= 0", 1);
    return kFALSE;
  }

  Int_t maxJetIndex  = -1;
  Int_t max2JetIndex = -1;

  // ************
  // General histograms
  // _________________________________

  GetLeadingJets(maxJetIndex, max2JetIndex);
  
  if (fSkipEventsWithNoJets && maxJetIndex < 0) {
    fHistRejectedEvents->Fill("No jets", 1);
    return kFALSE;
  }

  AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fJets->At(maxJetIndex));
  if (fSkipEventsWithNoJets && !jet) {
    fHistRejectedEvents->Fill("Max Jet not found", 1);
    return kFALSE;
  }

  Float_t maxJetCorrPt = 0; 

  if (jet)
    maxJetCorrPt = jet->Pt() - fRho * jet->Area();

  if (fSkipEventsWithNoJets && maxJetCorrPt <= 0) {
    fHistRejectedEvents->Fill("Max Jet <= 0", 1);
    return kFALSE;
  }

  fHistCentrality->Fill(fCent);
  fHistRho[fCentBin]->Fill(fRho);

  if (jet) {
    fHistLeadingJetPt[fCentBin]->Fill(jet->Pt());
    fHistRhoVSleadJetPt->Fill(fRho * jet->Area(), jet->Pt());
    fHistCorrLeadingJetPt[fCentBin]->Fill(maxJetCorrPt);
  }
  
  AliEmcalJet* jet2 = 0;
  if (max2JetIndex >= 0)
    jet2 = dynamic_cast<AliEmcalJet*>(fJets->At(max2JetIndex));

  if (jet2)
    fHist2LeadingJetPt[fCentBin]->Fill(jet2->Pt());

  DoJetLoop();
  
  DoTrackLoop(maxJetIndex);

  if (fAnaType == kEMCAL)
    DoClusterLoop(maxJetIndex);

  // ************
  // Random cones
  // _________________________________
  
  const Float_t rcArea = fJetRadius * fJetRadius * TMath::Pi();

  // Simple random cones
  Float_t RCpt = 0;
  Float_t RCeta = 0;
  Float_t RCphi = 0;
  GetRigidCone(RCpt, RCeta, RCphi, kFALSE, 0);
  if (RCpt > 0) {
    fHistRCPt[fCentBin]->Fill(RCpt / rcArea);
    fHistDeltaPtRC[fCentBin]->Fill(RCpt - rcArea * fRho);
  }
  
  // Random cones far from leading jet
  Float_t RCptExLJ = 0;
  Float_t RCetaExLJ = 0;
  Float_t RCphiExLJ = 0;
  GetRigidCone(RCptExLJ, RCetaExLJ, RCphiExLJ, kFALSE, jet);
  if (RCptExLJ > 0) {
    fHistRCPhiEta->Fill(RCetaExLJ, RCphiExLJ);
    fHistRhoVSRCPt->Fill(fRho, RCptExLJ / rcArea);

    Float_t dphi = RCphiExLJ - jet->Phi();
    if (dphi > 4.8) dphi -= TMath::Pi() * 2;
    if (dphi < -1.6) dphi += TMath::Pi() * 2; 
    fHistRCPtExLJVSDPhiLJ->Fill(RCptExLJ, dphi);
    
    fHistRCPtExLJ[fCentBin]->Fill(RCptExLJ / rcArea);
    fHistDeltaPtRCExLJ[fCentBin]->Fill(RCptExLJ - rcArea * fRho);
  }

  // Random cones with randomized particles
  Float_t RCptRand = 0;
  Float_t RCetaRand = 0;
  Float_t RCphiRand = 0;
  GetRigidCone(RCptRand, RCetaRand, RCphiRand, kTRUE, 0, fRandTracks, fRandCaloClusters);
  if (RCptRand > 0) {
    fHistRCPtRand[fCentBin]->Fill(RCptRand / rcArea);
    fHistDeltaPtRCRand[fCentBin]->Fill(RCptRand - rcArea * fRho);
  }

  // ************
  // Embedding
  // _________________________________

  if (!fEmbJets)
    return kTRUE;

  AliEmcalJet *embJet  = 0;
  TObject     *maxPart = 0;

  DoEmbJetLoop(embJet, maxPart);

  if (embJet) {

    fHistEmbJets[fCentBin]->Fill(embJet->Pt());
    fHistEmbPart[fCentBin]->Fill(embJet->MCPt());
    fHistEmbJetPhiEta->Fill(embJet->Eta(), embJet->Phi());
    fHistDeltaPtEmb[fCentBin]->Fill(embJet->Pt() - embJet->Area() * fRho - embJet->MCPt());
    fHistRhoVSEmbBkg->Fill(embJet->Area() * fRho, embJet->Pt() - embJet->MCPt());

    AliVCluster *cluster = dynamic_cast<AliVCluster*>(maxPart);
    if (cluster) {
      Float_t pos[3];
      cluster->GetPosition(pos);
      TVector3 clusVec(pos);
      
      fHistEmbPartPhiEta->Fill(clusVec.Eta(), clusVec.Phi());
    }
    else {
      AliVTrack *track = dynamic_cast<AliVTrack*>(maxPart);
      if (track) {
	fHistEmbPartPhiEta->Fill(track->Eta(), track->Phi());
      }
      else {
	AliWarning(Form("%s - Embedded particle type not found or not recognized (neither AliVCluster nor AliVTrack)!", GetName()));
	return kTRUE;
      }
    }

  }
  else {
    AliWarning(Form("%s - Embedded jet not found in the event!", GetName()));
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::GetLeadingJets(Int_t &maxJetIndex, Int_t &max2JetIndex)
{
  if (!fJets)
    return;

  Int_t njets = fJets->GetEntriesFast();

  Float_t maxJetPt = 0;
  Float_t max2JetPt = 0;
  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    Float_t corrPt = jet->Pt() - fRho * jet->Area();

    if (maxJetIndex == -1 || maxJetPt < corrPt) {
      max2JetPt = maxJetPt;
      max2JetIndex = maxJetIndex;
      maxJetPt = corrPt;
      maxJetIndex = ij;
    }
    else if (max2JetIndex == -1 || max2JetPt < corrPt) {
      max2JetPt = corrPt;
      max2JetIndex = ij;
    }
  } //jet loop 
}


//________________________________________________________________________
void AliAnalysisTaskSAJF::DoJetLoop()
{
  if (!fJets)
    return;

  Int_t njets = fJets->GetEntriesFast();

  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    Float_t corrPt = jet->Pt() - fRho * jet->Area();

    fHistJetsPt[fCentBin]->Fill(jet->Pt());
    fHistJetsPtArea[fCentBin]->Fill(corrPt, jet->Area());
    fHistCorrJetsPt[fCentBin]->Fill(corrPt);

    fHistJetPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());

    if (fAnaType == kEMCAL)
      fHistJetsNEFvsPt[fCentBin]->Fill(jet->NEF(), jet->Pt());

    if (fTracks) {
      for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
	AliVParticle *track = jet->TrackAt(it, fTracks);
	if (track)
	  fHistJetsZvsPt[fCentBin]->Fill(track->Pt() / jet->Pt(), jet->Pt());
      }
    }

    if (fCaloClusters) {
      for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
	AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);
	
	if (cluster) {
	  TLorentzVector nPart;
	  cluster->GetMomentum(nPart, fVertex);
	fHistJetsZvsPt[fCentBin]->Fill(nPart.Et() / jet->Pt(), jet->Pt());
	}
      }
    }
  } //jet loop 
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::DoEmbJetLoop(AliEmcalJet* &embJet, TObject* &maxPart)
{
  if (!fEmbJets)
    return;

  TLorentzVector maxClusVect;

  Int_t nembjets = fEmbJets->GetEntriesFast();

  for (Int_t ij = 0; ij < nembjets; ij++) {
      
    AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(fEmbJets->At(ij));
      
    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    } 
      
    if (!AcceptJet(jet))
      continue;

    if (!jet->IsMC())
      continue;

    AliVParticle *maxTrack = 0;

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, fTracks);
      
      if (!track) 
	continue;

      if (track->GetLabel() != 100 && !track->InheritsFrom("AliMCParticle"))
	continue;
      
      if (!maxTrack || track->Pt() > maxTrack->Pt())
	maxTrack = track;
    }
    
    AliVCluster *maxClus = 0;

    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);
      
      if (!cluster) 
	continue;

      if (cluster->Chi2() != 100)
	continue;
      
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      
      if (!maxClus || nPart.Et() > maxClusVect.Et()) {
	new (&maxClusVect) TLorentzVector(nPart);
	maxClus = cluster;
      } 
    }

    if ((maxClus && maxTrack && maxClusVect.Et() > maxTrack->Pt()) || (maxClus && !maxTrack)) {
      maxPart = maxClus;
      embJet = jet;
    }
    else if (maxTrack) {
      maxPart = maxTrack;
      embJet = jet;
    }

    return;  // MC jets found, exit
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::DoTrackLoop(Int_t maxJetIndex)
{ 
  if (!fTracks)
    return;

  AliEmcalJet* jet = 0;
  if (maxJetIndex >= 0 && fJets)
    jet = dynamic_cast<AliEmcalJet*>(fJets->At(maxJetIndex));

  Int_t ntracks = fTracks->GetEntriesFast();

  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliVTrack* track = dynamic_cast<AliVTrack*>(fTracks->At(iTracks));         
    if(!track) {
      AliError(Form("Could not retrieve track %d",iTracks)); 
      continue; 
    }
    
    if (!AcceptTrack(track)) continue;

    if (jet && IsJetTrack(jet, iTracks)) {
      fHistTracksPtLJ[fCentBin]->Fill(track->Pt());
    }
    else {
      fHistTracksPtBkg[fCentBin]->Fill(track->Pt());
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::DoClusterLoop(Int_t maxJetIndex)
{
  if (!fCaloClusters)
    return;

  AliEmcalJet* jet = 0;
  if (maxJetIndex >= 0 && fJets)
    jet = dynamic_cast<AliEmcalJet*>(fJets->At(maxJetIndex));

  Int_t nclusters =  fCaloClusters->GetEntriesFast();
  for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
    AliVCluster* cluster = dynamic_cast<AliVCluster*>(fCaloClusters->At(iClusters));
    if (!cluster) {
      AliError(Form("Could not receive cluster %d", iClusters));
      continue;
    }  
    
    if (!(cluster->IsEMCAL())) continue;
    
    if (!AcceptCluster(cluster)) continue;

    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);

    if (jet && IsJetCluster(jet, iClusters)) {
      fHistClusEtLJ[fCentBin]->Fill(nPart.Et());
    }
    else {
      fHistClusEtBkg[fCentBin]->Fill(nPart.Et());
    }
  } //cluster loop 
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::GetRigidCone(Float_t &pt, Float_t &eta, Float_t &phi, Bool_t acceptMC,
				       AliEmcalJet *jet, TClonesArray* tracks, TClonesArray* clusters) const
{
  if (!tracks)
    tracks = fTracks;

  if (!clusters)
    clusters = fCaloClusters;

  if (!tracks && !clusters)
    return;

  eta = 0;
  phi = 0;
  pt = 0;

  Float_t LJeta = 999;
  Float_t LJphi = 999;

  if (jet) {
    LJeta = jet->Eta();
    LJphi = jet->Phi();
  }

  Float_t maxEta = fMaxEta - fJetRadius;
  Float_t minEta = fMinEta + fJetRadius;
  Float_t maxPhi = fMaxPhi - fJetRadius;
  Float_t minPhi = fMinPhi + fJetRadius;

  if (maxPhi > TMath::Pi() * 2) maxPhi = TMath::Pi() * 2;
  if (minPhi < 0) minPhi = 0;
  
  Float_t dLJ = 0;
  Int_t repeats = 0;
  do {
    eta = gRandom->Rndm() * (maxEta - minEta) + minEta;
    phi = gRandom->Rndm() * (maxPhi - minPhi) + minPhi;
    dLJ = TMath::Sqrt((LJeta - eta) * (LJeta - eta) + (LJphi - phi) * (LJphi - phi));
    repeats++;
  } while (dLJ < fMinRC2LJ && repeats < 999);

  if (repeats == 999)
    return;

  if (fAnaType == kEMCAL && clusters) {
    Int_t nclusters =  clusters->GetEntriesFast();
    for (Int_t iClusters = 0; iClusters < nclusters; iClusters++) {
      AliVCluster* cluster = dynamic_cast<AliVCluster*>(clusters->At(iClusters));
      if (!cluster) {
	AliError(Form("Could not receive cluster %d", iClusters));
	continue;
      }  
      
      if (!(cluster->IsEMCAL())) continue;
      
      if (!AcceptCluster(cluster, acceptMC)) continue;
      
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, const_cast<Double_t*>(fVertex));
      
      Float_t pos[3];
      cluster->GetPosition(pos);
      TVector3 clusVec(pos);
      
      Float_t d = TMath::Sqrt((clusVec.Eta() - eta) * (clusVec.Eta() - eta) + (clusVec.Phi() - phi) * (clusVec.Phi() - phi));
      
      if (d <= fJetRadius)
	pt += nPart.Pt();
    }
  }

  if (tracks) {
    Int_t ntracks = tracks->GetEntriesFast();
    for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
      AliVTrack* track = dynamic_cast<AliVTrack*>(tracks->At(iTracks));         
      if(!track) {
	AliError(Form("Could not retrieve track %d",iTracks)); 
	continue; 
      }
      
      if (!AcceptTrack(track, acceptMC)) continue;
      
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
  }
}
//________________________________________________________________________
void AliAnalysisTaskSAJF::Init()
{
  AliAnalysisTaskEmcalJet::Init();

  const Float_t semiDiag = TMath::Sqrt((fMaxPhi - fMinPhi) * (fMaxPhi - fMinPhi) + (fMaxEta - fMinEta) * (fMaxEta - fMinEta)) / 2;
  if (fMinRC2LJ > semiDiag * 0.5) {
    AliWarning(Form("The parameter fMinRC2LJ = %f is too large for the considered acceptance. Will use fMinRC2LJ = %f", fMinRC2LJ, semiDiag * 0.5));
    fMinRC2LJ = semiDiag * 0.5;
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
