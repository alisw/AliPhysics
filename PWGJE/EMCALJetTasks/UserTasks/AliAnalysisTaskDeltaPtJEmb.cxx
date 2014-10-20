// $Id$
//
// Jet embedding deltaPt task.
//
// Author: M. Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskDeltaPtJEmb.h"

ClassImp(AliAnalysisTaskDeltaPtJEmb)

//________________________________________________________________________
AliAnalysisTaskDeltaPtJEmb::AliAnalysisTaskDeltaPtJEmb() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskDeltaPtJEmb", kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  fMinFractionShared(0.5),
  fHistEmbJetsPtArea(0),
  fHistEmbJetsCorrPtArea(0),
  fHistEmbPartPtvsJetPt(0),
  fHistEmbPartPtvsJetCorrPt(0),
  fHistJetPtvsJetCorrPt(0),
  fHistDistLeadPart2JetAxis(0),
  fHistEmbBkgArea(0),
  fHistRhoVSEmbBkg(0),
  fHistDeltaPtEmbArea(0),
  fHistDeltaPtEmbvsEP(0),
  fHistDeltaPtEmbPtProbe(0),
  fHistEmbJetsPhiEta(0),
  fHistLeadPartPhiEta(0)
{
  // Default constructor.

  fHistEmbJetsPtArea = new TH3*[fNcentBins];
  fHistEmbJetsCorrPtArea = new TH3*[fNcentBins];
  fHistEmbPartPtvsJetPt = new TH2*[fNcentBins];
  fHistEmbPartPtvsJetCorrPt = new TH2*[fNcentBins];
  fHistJetPtvsJetCorrPt = new TH2*[fNcentBins];
  fHistDistLeadPart2JetAxis = new TH1*[fNcentBins];
  fHistEmbBkgArea = new TH2*[fNcentBins];
  fHistRhoVSEmbBkg = new TH2*[fNcentBins];
  fHistDeltaPtEmbArea = new TH2*[fNcentBins];
  fHistDeltaPtEmbvsEP = new TH2*[fNcentBins];
  fHistDeltaPtEmbPtProbe = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistEmbJetsPtArea[i] = 0;
    fHistEmbJetsCorrPtArea[i] = 0;
    fHistEmbPartPtvsJetPt[i] = 0;
    fHistEmbPartPtvsJetCorrPt[i] = 0;
    fHistJetPtvsJetCorrPt[i] = 0;
    fHistDistLeadPart2JetAxis[i] = 0;
    fHistEmbBkgArea[i] = 0;
    fHistRhoVSEmbBkg[i] = 0;
    fHistDeltaPtEmbArea[i] = 0;
    fHistDeltaPtEmbvsEP[i] = 0;
    fHistDeltaPtEmbPtProbe[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskDeltaPtJEmb::AliAnalysisTaskDeltaPtJEmb(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  fMinFractionShared(0.5),
  fHistEmbJetsPtArea(0),
  fHistEmbJetsCorrPtArea(0),
  fHistEmbPartPtvsJetPt(0),
  fHistEmbPartPtvsJetCorrPt(0),
  fHistJetPtvsJetCorrPt(0),
  fHistDistLeadPart2JetAxis(0),
  fHistEmbBkgArea(0),
  fHistRhoVSEmbBkg(0),
  fHistDeltaPtEmbArea(0),
  fHistDeltaPtEmbvsEP(0),
  fHistDeltaPtEmbPtProbe(0),
  fHistEmbJetsPhiEta(0),
  fHistLeadPartPhiEta(0)
{
  // Standard constructor.

  fHistEmbJetsPtArea = new TH3*[fNcentBins];
  fHistEmbJetsCorrPtArea = new TH3*[fNcentBins];
  fHistEmbPartPtvsJetPt = new TH2*[fNcentBins];
  fHistEmbPartPtvsJetCorrPt = new TH2*[fNcentBins];
  fHistJetPtvsJetCorrPt = new TH2*[fNcentBins];
  fHistDistLeadPart2JetAxis = new TH1*[fNcentBins];
  fHistEmbBkgArea = new TH2*[fNcentBins];
  fHistRhoVSEmbBkg = new TH2*[fNcentBins];
  fHistDeltaPtEmbArea = new TH2*[fNcentBins];
  fHistDeltaPtEmbvsEP = new TH2*[fNcentBins];
  fHistDeltaPtEmbPtProbe = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistEmbJetsPtArea[i] = 0;
    fHistEmbJetsCorrPtArea[i] = 0;
    fHistEmbPartPtvsJetPt[i] = 0;
    fHistEmbPartPtvsJetCorrPt[i] = 0;
    fHistJetPtvsJetCorrPt[i] = 0;
    fHistDistLeadPart2JetAxis[i] = 0;
    fHistEmbBkgArea[i] = 0;
    fHistRhoVSEmbBkg[i] = 0;
    fHistDeltaPtEmbArea[i] = 0;
    fHistDeltaPtEmbvsEP[i] = 0;
    fHistDeltaPtEmbPtProbe[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPtJEmb::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetsCont = GetJetContainer("Jets");
  fTracksCont = GetParticleContainer("Tracks");
  fCaloClustersCont = GetClusterContainer("CaloClusters");
  
  fHistEmbJetsPhiEta = new TH2F("fHistEmbJetsPhiEta","fHistEmbJetsPhiEta", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
  fHistEmbJetsPhiEta->GetXaxis()->SetTitle("#eta");
  fHistEmbJetsPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistEmbJetsPhiEta);
  
  fHistLeadPartPhiEta = new TH2F("fHistLeadPartPhiEta","fHistLeadPartPhiEta", 100, -1, 1, 201, 0, TMath::Pi() * 2.01);
  fHistLeadPartPhiEta->GetXaxis()->SetTitle("#eta");
  fHistLeadPartPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistLeadPartPhiEta);
  
  TString histname;

  const Int_t nbinsZ = 12;
  Double_t binsZ[nbinsZ+1] = {0,1,2,3,4,5,6,7,8,9,10,20,1000};

  Double_t *binsPt       = GenerateFixedBinArray(fNbins, fMinBinPt, fMaxBinPt);
  Double_t *binsCorrPt   = GenerateFixedBinArray(fNbins*2, -fMaxBinPt, fMaxBinPt);
  Double_t *binsArea     = GenerateFixedBinArray(50, 0, 2);

  for (Int_t i = 0; i < fNcentBins; i++) {
    histname = "fHistEmbJetsPtArea_";
    histname += i;
    fHistEmbJetsPtArea[i] = new TH3F(histname.Data(), histname.Data(), 50, binsArea, fNbins, binsPt, nbinsZ, binsZ);
    fHistEmbJetsPtArea[i]->GetXaxis()->SetTitle("area");
    fHistEmbJetsPtArea[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb,raw} (GeV/#it{c})");
    fOutput->Add(fHistEmbJetsPtArea[i]);

    histname = "fHistEmbJetsCorrPtArea_";
    histname += i;
    fHistEmbJetsCorrPtArea[i] = new TH3F(histname.Data(), histname.Data(), 50, binsArea, fNbins * 2, binsCorrPt, nbinsZ, binsZ);
    fHistEmbJetsCorrPtArea[i]->GetXaxis()->SetTitle("area");
    fHistEmbJetsCorrPtArea[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb,corr} (GeV/#it{c})");
    fOutput->Add(fHistEmbJetsCorrPtArea[i]);

    histname = "fHistEmbPartPtvsJetPt_";
    histname += i;
    fHistEmbPartPtvsJetPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    fHistEmbPartPtvsJetPt[i]->GetXaxis()->SetTitle("#sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
    fHistEmbPartPtvsJetPt[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} (GeV/#it{c})");
    fHistEmbPartPtvsJetPt[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistEmbPartPtvsJetPt[i]);

    histname = "fHistEmbPartPtvsJetCorrPt_";
    histname += i;
    fHistEmbPartPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), 
					    fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
    fHistEmbPartPtvsJetCorrPt[i]->GetXaxis()->SetTitle("#sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
    fHistEmbPartPtvsJetCorrPt[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - A#rho (GeV/#it{c})");
    fHistEmbPartPtvsJetCorrPt[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistEmbPartPtvsJetCorrPt[i]);

    histname = "fHistJetPtvsJetCorrPt_";
    histname += i;
    fHistJetPtvsJetCorrPt[i] = new TH2F(histname.Data(), histname.Data(), 
					fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
    fHistJetPtvsJetCorrPt[i]->GetXaxis()->SetTitle("#it{p}_{T,jet}^{emb} (GeV/#it{c})");
    fHistJetPtvsJetCorrPt[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - A#rho (GeV/#it{c})");
    fHistJetPtvsJetCorrPt[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJetPtvsJetCorrPt[i]);

    histname = "fHistDistLeadPart2JetAxis_";
    histname += i;
    fHistDistLeadPart2JetAxis[i] = new TH1F(histname.Data(), histname.Data(), 50, 0, 0.5);
    fHistDistLeadPart2JetAxis[i]->GetXaxis()->SetTitle("distance");
    fHistDistLeadPart2JetAxis[i]->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistDistLeadPart2JetAxis[i]);

    histname = "fHistEmbBkgArea_";
    histname += i;
    fHistEmbBkgArea[i] = new TH2F(histname.Data(), histname.Data(), 50, 0, 2, fNbins, fMinBinPt, fMaxBinPt);
    fHistEmbBkgArea[i]->GetXaxis()->SetTitle("area");
    fHistEmbBkgArea[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - #sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
    fOutput->Add(fHistEmbBkgArea[i]);

    histname = "fHistRhoVSEmbBkg_";
    histname += i;
    fHistRhoVSEmbBkg[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    fHistRhoVSEmbBkg[i]->GetXaxis()->SetTitle("A#rho (GeV/#it{c})");
    fHistRhoVSEmbBkg[i]->GetYaxis()->SetTitle("#it{p}_{T,jet}^{emb} - #sum#it{p}_{T,const}^{emb} (GeV/#it{c})");
    fOutput->Add(fHistRhoVSEmbBkg[i]);
      
    histname = "fHistDeltaPtEmbArea_";
    histname += i;
    fHistDeltaPtEmbArea[i] = new TH2F(histname.Data(), histname.Data(), 
				      50, 0, 2, fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtEmbArea[i]->GetXaxis()->SetTitle("area");
    fHistDeltaPtEmbArea[i]->GetYaxis()->SetTitle("#delta#it{p}_{T}^{emb} (GeV/#it{c})");
    fHistDeltaPtEmbArea[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtEmbArea[i]);

    histname = "fHistDeltaPtEmbvsEP_";
    histname += i;
    fHistDeltaPtEmbvsEP[i] = new TH2F(histname.Data(), histname.Data(), 
				      402, -TMath::Pi()*1.01, TMath::Pi()*3.01, fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtEmbvsEP[i]->GetXaxis()->SetTitle("#phi_{jet} - #Psi_{EP}");
    fHistDeltaPtEmbvsEP[i]->GetYaxis()->SetTitle("#delta#it{p}_{T}^{emb} (GeV/#it{c})");
    fHistDeltaPtEmbvsEP[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtEmbvsEP[i]);

    histname = "fHistDeltaPtEmbPtProbe_";
    histname += i;
    fHistDeltaPtEmbPtProbe[i] = new TH2F(histname.Data(), histname.Data(), 
					 fNbins, fMinBinPt, fMaxBinPt, fNbins * 2, -fMaxBinPt, fMaxBinPt);
    fHistDeltaPtEmbPtProbe[i]->GetXaxis()->SetTitle("#it{p}_{T,probe} (GeV/#it{c})");
    fHistDeltaPtEmbPtProbe[i]->GetYaxis()->SetTitle("#delta#it{p}_{T}^{emb} (GeV/#it{c})");
    fHistDeltaPtEmbPtProbe[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaPtEmbPtProbe[i]);
  }

  delete[] binsPt;
  delete[] binsCorrPt;
  delete[] binsArea;

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDeltaPtJEmb::FillHistograms()
{
  // Fill histograms.

  // ************
  // Embedding
  // Jets should already be matched, GetClosestJet gives matched partner. Cut on shared fraction to be applied
  // _________________________________

  if (fJetsCont) {
    AliEmcalJet *jet1 = NULL;
    AliEmcalJet *jet2 = NULL;
    fJetsCont->ResetCurrentID();
    while((jet1 = fJetsCont->GetNextAcceptJet())) {

      jet2 = jet1->ClosestJet();
      if(!jet2) continue;

      Double_t fraction = fJetsCont->GetFractionSharedPt(jet1);
      if(fMinFractionShared>0. && fraction<fMinFractionShared) continue;

      TLorentzVector mom;
      fJetsCont->GetLeadingHadronMomentum(mom,jet1);
      
      Double_t distLeading2Jet = TMath::Sqrt((jet1->Eta() - mom.Eta()) * (jet1->Eta() - mom.Eta()) + (jet1->Phi() - mom.Phi()+TMath::Pi()) * (jet1->Phi() - mom.Phi()+TMath::Pi()));
      
      fHistEmbPartPtvsJetPt[fCentBin]->Fill(jet2->Pt(), jet1->Pt());
      fHistEmbPartPtvsJetCorrPt[fCentBin]->Fill(jet2->Pt(), jet1->Pt() - jet1->Area() * fRhoVal);
      fHistLeadPartPhiEta->Fill(mom.Eta(), mom.Phi()+TMath::Pi());
      fHistDistLeadPart2JetAxis[fCentBin]->Fill(distLeading2Jet);
      
      fHistEmbJetsPtArea[fCentBin]->Fill(jet1->Area(), jet1->Pt(), mom.Pt());
      fHistEmbJetsCorrPtArea[fCentBin]->Fill(jet1->Area(), jet1->Pt() - fRhoVal * jet1->Area(), mom.Pt());
      fHistEmbJetsPhiEta->Fill(jet1->Eta(), jet1->Phi());
      fHistJetPtvsJetCorrPt[fCentBin]->Fill(jet1->Pt(), jet1->Pt() - fRhoVal * jet1->Area());
      
      fHistEmbBkgArea[fCentBin]->Fill(jet1->Area(), jet1->Pt() - jet2->Pt());
      fHistRhoVSEmbBkg[fCentBin]->Fill(fRhoVal * jet1->Area(), jet1->Pt() - jet2->Pt());
      fHistDeltaPtEmbArea[fCentBin]->Fill(jet1->Area(), jet1->Pt() - jet1->Area() * fRhoVal - jet2->Pt());
      fHistDeltaPtEmbvsEP[fCentBin]->Fill(jet1->Phi() - fEPV0, jet1->Pt() - jet1->Area() * fRhoVal - jet2->Pt());
      fHistDeltaPtEmbPtProbe[fCentBin]->Fill(jet2->Pt(), jet1->Pt() - jet1->Area() * fRhoVal - jet2->Pt());
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskDeltaPtJEmb::ExecOnce()
{
  // Initialize the analysis.

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;
  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
}
