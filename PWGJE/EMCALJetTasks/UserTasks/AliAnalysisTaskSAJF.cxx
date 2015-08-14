//
// Jet analysis task.
//
// Author: S.Aiola

#include <TClonesArray.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"

#include "AliAnalysisTaskSAJF.h"

ClassImp(AliAnalysisTaskSAJF)

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSAJF", kTRUE),
  fHistoType(1),
  fHistRejectionReason(0),
  fHistTracksJetPt(0),
  fHistClustersJetPt(0),
  fHistTracksPtDist(0),
  fHistClustersPtDist(0),
  fHistTracksZJetPtJetConst(0),
  fHistClustersZJetPtJetConst(0),
  fHistJetObservables(0),
  fHistJetPtEtaPhi(0),
  fHistJetPtArea(0),
  fHistJetPtEP(0),
  fHistJetPtNEF(0),
  fHistJetPtZ(0),
  fHistJetPtLeadingPartPt(0),
  fHistJetCorrPtEtaPhi(0),
  fHistJetCorrPtArea(0),
  fHistJetCorrPtEP(0),
  fHistJetCorrPtNEF(0),
  fHistJetCorrPtZ(0),
  fHistJetCorrPtLeadingPartPt(0),
  fHistJetPtCorrPt(0),
  fHistJetPtMCPt(0),
  fHistJetMCPtCorrPt(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistoType(1),
  fHistRejectionReason(0),
  fHistTracksJetPt(0),
  fHistClustersJetPt(0),
  fHistTracksPtDist(0),
  fHistClustersPtDist(0),
  fHistTracksZJetPtJetConst(0),
  fHistClustersZJetPtJetConst(0),
  fHistJetObservables(0),
  fHistJetPtEtaPhi(0),
  fHistJetPtArea(0),
  fHistJetPtEP(0),
  fHistJetPtNEF(0),
  fHistJetPtZ(0),
  fHistJetPtLeadingPartPt(0),
  fHistJetCorrPtEtaPhi(0),
  fHistJetCorrPtArea(0),
  fHistJetCorrPtEP(0),
  fHistJetCorrPtNEF(0),
  fHistJetCorrPtZ(0),
  fHistJetCorrPtLeadingPartPt(0),
  fHistJetPtCorrPt(0),
  fHistJetPtMCPt(0),
  fHistJetMCPtCorrPt(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::AllocateTHnSparse()
{
  Double_t jetRadius = 0.4;
  AliJetContainer *jets = GetJetContainer(0);
  if (jets) jetRadius = jets->GetJetRadius();
  
  TString title[20]= {""};
  Int_t nbins[20]  = {0};
  Double_t min[20] = {0.};
  Double_t max[20] = {0.};
  Int_t dim = 0;

  if (fForceBeamType != kpp) {
    title[dim] = "Centrality (%)";
    nbins[dim] = 20;
    min[dim] = 0;
    max[dim] = 100;
    dim++;
    
    title[dim] = "#phi_{jet} - #psi_{RP}";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = TMath::Pi();
    dim++;
  }

  title[dim] = "#eta_{jet}";
  nbins[dim] = 50;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "#phi_{jet} (rad)";
  nbins[dim] = 150;
  min[dim] = 0;
  max[dim] = 2*TMath::Pi();
  dim++;

  title[dim] = "p_{T} (GeV/c)";
  nbins[dim] = fNbins;
  min[dim] = fMinBinPt;
  max[dim] = fMaxBinPt;
  dim++;

  if (fIsEmbedded) {
    title[dim] = "p_{T}^{MC} (GeV/c)";
    nbins[dim] = fNbins;
    min[dim] = fMinBinPt;
    max[dim] = fMaxBinPt;
    dim++;
  }

  if (!GetRhoName().IsNull()) {
    title[dim] = "p_{T}^{corr} (GeV/c)";
    nbins[dim] = fNbins*2;
    min[dim] = -fMaxBinPt;
    max[dim] = fMaxBinPt;
    dim++;
  }

  title[dim] = "A_{jet}";
  nbins[dim] = 100;
  min[dim] = 0;
  max[dim] = 2.0*jetRadius*jetRadius*TMath::Pi();
  dim++;

  if (fCaloClusters) {
    title[dim] = "NEF";
    nbins[dim] = 102;
    min[dim] = 0;
    max[dim] = 1.02;
    dim++;
  }

  title[dim] = "Z";
  nbins[dim] = 102;
  min[dim] = 0;
  max[dim] = 1.02;
  dim++;

  if (fForceBeamType != kpp) {
    title[dim] = "No. of constituents";
    nbins[dim] = 250;
    min[dim] = -0.5;
    max[dim] = 249.5;
    dim++;
  }
  else {
    title[dim] = "No. of constituents";
    nbins[dim] = 50;
    min[dim] = -0.5;
    max[dim] = 49.5;
    dim++;
  }

  title[dim] = "p_{T,particle}^{leading} (GeV/c)";
  nbins[dim] = 120;
  min[dim] = 0;
  max[dim] = 120;
  dim++;

  fHistJetObservables = new THnSparseD("fHistJetObservables","fHistJetObservables",dim,nbins,min,max);
  fOutput->Add(fHistJetObservables);
  for (Int_t i = 0; i < dim; i++)
    fHistJetObservables->GetAxis(i)->SetTitle(title[i]);
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::AllocateTHX()
{
  fHistJetPtEtaPhi = new TH3*[fNcentBins];
  fHistJetPtArea = new TH2*[fNcentBins];
  fHistJetPtEP = new TH2*[fNcentBins];
  fHistJetPtNEF = new TH2*[fNcentBins];
  fHistJetPtZ = new TH2*[fNcentBins];
  fHistJetPtLeadingPartPt = new TH2*[fNcentBins];
  fHistJetCorrPtEtaPhi = new TH3*[fNcentBins];
  fHistJetCorrPtArea = new TH2*[fNcentBins];
  fHistJetCorrPtEP = new TH2*[fNcentBins];
  fHistJetCorrPtNEF = new TH2*[fNcentBins];
  fHistJetCorrPtZ = new TH2*[fNcentBins];
  fHistJetCorrPtLeadingPartPt = new TH2*[fNcentBins];
  fHistJetPtCorrPt = new TH2*[fNcentBins];
  fHistJetPtMCPt = new TH2*[fNcentBins];
  fHistJetMCPtCorrPt = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    TString histname;

    histname = "fHistJetPtEtaPhi_";
    histname += i;
    fHistJetPtEtaPhi[i] = new TH3F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 20, -1, 1, 41, 0, 2*TMath::Pi()*41/40);
    fHistJetPtEtaPhi[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistJetPtEtaPhi[i]->GetYaxis()->SetTitle("#eta");
    fHistJetPtEtaPhi[i]->GetZaxis()->SetTitle("#phi_{jet} (rad)");
    fOutput->Add(fHistJetPtEtaPhi[i]);
      
    histname = "fHistJetPtArea_";
    histname += i;
    fHistJetPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 150, 0, 1.5);
    fHistJetPtArea[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistJetPtArea[i]->GetYaxis()->SetTitle("A_{jet}");
    fHistJetPtArea[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJetPtArea[i]);

    histname = "fHistJetPtEP_";
    histname += i;
    fHistJetPtEP[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 100, 0, TMath::Pi());
    fHistJetPtEP[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistJetPtEP[i]->GetYaxis()->SetTitle("#phi_{jet} - #psi_{RP}");
    fHistJetPtEP[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJetPtEP[i]);

    histname = "fHistJetPtNEF_";
    histname += i;
    fHistJetPtNEF[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 102, 0, 1.02);
    fHistJetPtNEF[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistJetPtNEF[i]->GetYaxis()->SetTitle("NEF");
    fHistJetPtNEF[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJetPtNEF[i]);

    histname = "fHistJetPtZ_";
    histname += i;
    fHistJetPtZ[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 102, 0, 1.02);
    fHistJetPtZ[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistJetPtZ[i]->GetYaxis()->SetTitle("z");
    fHistJetPtZ[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJetPtZ[i]);

    histname = "fHistJetPtLeadingPartPt_";
    histname += i;
    fHistJetPtLeadingPartPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 120, 0, 120);
    fHistJetPtLeadingPartPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistJetPtLeadingPartPt[i]->GetYaxis()->SetTitle("p_{T,particle}^{leading} (GeV/c)");
    fHistJetPtLeadingPartPt[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJetPtLeadingPartPt[i]);

    if (!GetRhoName().IsNull()) {
      histname = "fHistJetCorrPtEtaPhi_";
      histname += i;
      fHistJetCorrPtEtaPhi[i] = new TH3F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 20, -1, 1, 41, 0, 2*TMath::Pi()*201/200);
      fHistJetCorrPtEtaPhi[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistJetCorrPtEtaPhi[i]->GetYaxis()->SetTitle("#eta");
      fHistJetCorrPtEtaPhi[i]->GetZaxis()->SetTitle("#phi_{jet} (rad)");
      fOutput->Add(fHistJetCorrPtEtaPhi[i]);
      
      histname = "fHistJetCorrPtArea_";
      histname += i;
      fHistJetCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 150, 0, 1.5);
      fHistJetCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
      fHistJetCorrPtArea[i]->GetYaxis()->SetTitle("A_{jet}");
      fHistJetCorrPtArea[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetCorrPtArea[i]);

      histname = "fHistJetCorrPtEP_";
      histname += i;
      fHistJetCorrPtEP[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 100, 0, TMath::Pi());
      fHistJetCorrPtEP[i]->GetXaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
      fHistJetCorrPtEP[i]->GetYaxis()->SetTitle("#phi_{jet} - #psi_{RP}");
      fHistJetCorrPtEP[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetCorrPtEP[i]);

      histname = "fHistJetCorrPtNEF_";
      histname += i;
      fHistJetCorrPtNEF[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 102, 0, 1.02);
      fHistJetCorrPtNEF[i]->GetXaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
      fHistJetCorrPtNEF[i]->GetYaxis()->SetTitle("NEF");
      fHistJetCorrPtNEF[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetCorrPtNEF[i]);

      histname = "fHistJetCorrPtZ_";
      histname += i;
      fHistJetCorrPtZ[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 102, 0, 1.02);
      fHistJetCorrPtZ[i]->GetXaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
      fHistJetCorrPtZ[i]->GetYaxis()->SetTitle("z");
      fHistJetCorrPtZ[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetCorrPtZ[i]);

      histname = "fHistJetCorrPtLeadingPartPt_";
      histname += i;
      fHistJetCorrPtLeadingPartPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 120, 0, 120);
      fHistJetCorrPtLeadingPartPt[i]->GetXaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
      fHistJetCorrPtLeadingPartPt[i]->GetYaxis()->SetTitle("p_{T,particle}^{leading} (GeV/c)");
      fHistJetCorrPtLeadingPartPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetCorrPtLeadingPartPt[i]);

      histname = "fHistJetPtCorrPt_";
      histname += i;
      fHistJetPtCorrPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
      fHistJetPtCorrPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistJetPtCorrPt[i]->GetYaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
      fHistJetPtCorrPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetPtCorrPt[i]);

      if (fIsEmbedded) {
	histname = "fHistJetMCPtCorrPt_";
	histname += i;
	fHistJetMCPtCorrPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins*2, -fMaxBinPt, fMaxBinPt);
	fHistJetMCPtCorrPt[i]->GetXaxis()->SetTitle("p_{T}^{MC} (GeV/c)");
	fHistJetMCPtCorrPt[i]->GetYaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
	fHistJetMCPtCorrPt[i]->GetZaxis()->SetTitle("counts");
	fOutput->Add(fHistJetMCPtCorrPt[i]);
      }
    }

    if (fIsEmbedded) {
      histname = "fHistJetPtMCPt_";
      histname += i;
      fHistJetPtMCPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
      fHistJetPtMCPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
      fHistJetPtMCPt[i]->GetYaxis()->SetTitle("p_{T}^{MC} (GeV/c)");
      fHistJetPtMCPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetPtMCPt[i]);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Int_t constituentsNbins = 250;
  Double_t constituentsMax = 249.5;
  
  if (fForceBeamType == kpp) {
    constituentsNbins = 50;
    constituentsMax = 49.5;
  }
  
  if (fHistoType == 0) 
    AllocateTHX();
  else
    AllocateTHnSparse();

  fHistTracksJetPt = new TH2*[fNcentBins];
  fHistClustersJetPt = new TH2*[fNcentBins];
  fHistTracksPtDist = new TH2*[fNcentBins];
  fHistClustersPtDist = new TH2*[fNcentBins];
  fHistRejectionReason = new TH2*[fNcentBins];

  fHistTracksZJetPtJetConst = new TH3*[fNcentBins];
  fHistClustersZJetPtJetConst = new TH3*[fNcentBins];
  
  for (Int_t i = 0; i < fNcentBins; i++) {
    TString histname;

    if (fParticleCollArray.GetEntriesFast() > 0) {
      histname = "fHistTracksJetPt_";
      histname += i;
      fHistTracksJetPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2, fNbins, fMinBinPt, fMaxBinPt);
      fHistTracksJetPt[i]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
      fHistTracksJetPt[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
      fHistTracksJetPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistTracksJetPt[i]);
      
      histname = "fHistTracksPtDist_";
      histname += i;
      fHistTracksPtDist[i] = new TH2F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2, 100, 0, 5);
      fHistTracksPtDist[i]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
      fHistTracksPtDist[i]->GetYaxis()->SetTitle("d");
      fHistTracksPtDist[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistTracksPtDist[i]);

      histname = "fHistTracksZJetPtJetConst_";
      histname += i;
      fHistTracksZJetPtJetConst[i] = new TH3F(histname.Data(), histname.Data(), 120, 0.0, 1.2, fNbins, fMinBinPt, fMaxBinPt, constituentsNbins, -0.5, constituentsMax);
      fHistTracksZJetPtJetConst[i]->GetXaxis()->SetTitle("z_{track}");
      fHistTracksZJetPtJetConst[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
      fHistTracksZJetPtJetConst[i]->GetZaxis()->SetTitle("# of constituents");
      fOutput->Add(fHistTracksZJetPtJetConst[i]);
    }

    if (fClusterCollArray.GetEntriesFast() > 0) {
      histname = "fHistClustersJetPt_";
      histname += i;
      fHistClustersJetPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2, fNbins, fMinBinPt, fMaxBinPt);
      fHistClustersJetPt[i]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
      fHistClustersJetPt[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
      fHistClustersJetPt[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistClustersJetPt[i]);

      histname = "fHistClustersPtDist_";
      histname += i;
      fHistClustersPtDist[i] = new TH2F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2, 100, 0, 5);
      fHistClustersPtDist[i]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
      fHistClustersPtDist[i]->GetYaxis()->SetTitle("d");
      fHistClustersPtDist[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistClustersPtDist[i]);

      histname = "fHistClustersZJetPtJetConst_";
      histname += i;
      fHistClustersZJetPtJetConst[i] = new TH3F(histname.Data(), histname.Data(), 120, 0.0, 1.2, fNbins, fMinBinPt, fMaxBinPt, constituentsNbins, -0.5, constituentsMax);
      fHistClustersZJetPtJetConst[i]->GetXaxis()->SetTitle("z_{clus}");
      fHistClustersZJetPtJetConst[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
      fHistClustersZJetPtJetConst[i]->GetZaxis()->SetTitle("# of constituents");
      fOutput->Add(fHistClustersZJetPtJetConst[i]);
    }

    histname = "fHistRejectionReason_";
    histname += i;
    fHistRejectionReason[i] = new TH2F(histname, histname, 32, 0, 32, 100, 0, 250);
    fHistRejectionReason[i]->GetXaxis()->SetTitle("Rejection reason");
    fHistRejectionReason[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
    fHistRejectionReason[i]->GetZaxis()->SetTitle("counts");
    SetRejectionReasonLabels(fHistRejectionReason[i]->GetXaxis());
    fOutput->Add(fHistRejectionReason[i]);
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram

}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::FillHistograms()
{
  // Fill histograms.

  AliJetContainer *jets = GetJetContainer(0);

  if (!jets) return kFALSE;
  
  AliEmcalJet* jet = 0;

  jets->ResetCurrentID();
  while ((jet = jets->GetNextJet())) {
    if (!jet) {
      AliError("Could not receive jet!");
      continue;
    }

    if (!jets->AcceptJet(jet)) {
      fHistRejectionReason[fCentBin]->Fill(jets->GetRejectionReasonBitPosition(), jet->Pt());
      continue;
    }

    Float_t ptLeading = jets->GetLeadingHadronPt(jet);
    Float_t corrPt = jet->Pt() - fRhoVal * jet->Area();

    TLorentzVector leadPart;

    jets->GetLeadingHadronMomentum(leadPart, jet);

    // Fill THnSparse
    Double_t ep = jet->Phi() - fEPV0;
    while (ep < 0) ep += TMath::Pi();
    while (ep >= TMath::Pi()) ep -= TMath::Pi();

    Double_t z = GetParallelFraction(leadPart.Vect(), jet);
    if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1

    FillJetHisto(fCent, ep, jet->Eta(), jet->Phi(), jet->Pt(), jet->MCPt(), corrPt, jet->Area(), 
		 jet->NEF(), z, jet->GetNumberOfConstituents(), ptLeading);

    if (fTracks) {
      for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
	AliVParticle *track = jet->TrackAt(it, fTracks);
	if (track) {
	  fHistTracksJetPt[fCentBin]->Fill(track->Pt(), jet->Pt());
	  Double_t dist = TMath::Sqrt((track->Eta() - jet->Eta()) * (track->Eta() - jet->Eta()) + (track->Phi() - jet->Phi()) * (track->Phi() - jet->Phi()));
	  fHistTracksPtDist[fCentBin]->Fill(track->Pt(), dist);
          fHistTracksZJetPtJetConst[fCentBin]->Fill(GetParallelFraction(track, jet), jet->Pt(), jet->GetNumberOfConstituents());
	}
      }
    }

    if (fCaloClusters) {
      for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
	AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);
	
	if (cluster) {
	  TLorentzVector nPart;
	  cluster->GetMomentum(nPart, fVertex);

	  fHistClustersJetPt[fCentBin]->Fill(nPart.Pt(), jet->Pt());
	  Double_t dist = TMath::Sqrt((nPart.Eta() - jet->Eta()) * (nPart.Eta() - jet->Eta()) + (nPart.Phi() - jet->Phi()) * (nPart.Phi() - jet->Phi()));
	  fHistClustersPtDist[fCentBin]->Fill(nPart.Pt(), dist);
          fHistClustersZJetPtJetConst[fCentBin]->Fill(GetParallelFraction(nPart.Vect(), jet), jet->Pt(), jet->GetNumberOfConstituents());
	}
      }
    }
  } //jet loop 

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::FillJetHisto(Double_t cent, Double_t ep, Double_t eta, Double_t phi, Double_t pt, Double_t MCpt, Double_t corrpt, Double_t area, 
				       Double_t NEF, Double_t z, Int_t n, Double_t leadingpt)
{
  if (fHistoType == 0) {
    fHistJetPtEtaPhi[fCentBin]->Fill(pt,eta,phi);
    fHistJetPtArea[fCentBin]->Fill(pt,area);
    fHistJetPtEP[fCentBin]->Fill(pt,ep);
    fHistJetPtNEF[fCentBin]->Fill(pt,NEF);
    fHistJetPtZ[fCentBin]->Fill(pt,z);
    fHistJetPtLeadingPartPt[fCentBin]->Fill(pt,leadingpt);
    if (fHistJetCorrPtEtaPhi[fCentBin]) {
      fHistJetCorrPtEtaPhi[fCentBin]->Fill(corrpt,eta,phi);
      fHistJetCorrPtArea[fCentBin]->Fill(corrpt,area);
      fHistJetCorrPtEP[fCentBin]->Fill(corrpt,ep);
      fHistJetCorrPtNEF[fCentBin]->Fill(corrpt,NEF);
      fHistJetCorrPtZ[fCentBin]->Fill(corrpt,z);
      fHistJetCorrPtLeadingPartPt[fCentBin]->Fill(corrpt,leadingpt);
      fHistJetPtCorrPt[fCentBin]->Fill(pt,corrpt);
      if (fIsEmbedded)
	fHistJetMCPtCorrPt[fCentBin]->Fill(MCpt,corrpt);
    }
    if (fIsEmbedded)
      fHistJetPtMCPt[fCentBin]->Fill(pt,MCpt);
  }
  else {
 
    Double_t contents[20]={0};

    for (Int_t i = 0; i < fHistJetObservables->GetNdimensions(); i++) {
      TString title(fHistJetObservables->GetAxis(i)->GetTitle());
      if (title=="Centrality (%)")
	contents[i] = cent;
      else if (title=="#phi_{jet} - #psi_{RP}")
	contents[i] = ep;
      else if (title=="#eta_{jet}")
	contents[i] = eta;
      else if (title=="#phi_{jet} (rad)")
	contents[i] = phi;
      else if (title=="p_{T} (GeV/c)")
	contents[i] = pt;
      else if (title=="p_{T}^{MC} (GeV/c)")
	contents[i] = MCpt;
      else if (title=="p_{T}^{corr} (GeV/c)")
	contents[i] = corrpt;
      else if (title=="A_{jet}")
	contents[i] = area;
      else if (title=="NEF")
	contents[i] = NEF;
      else if (title=="Z")
	contents[i] = z;
      else if (title=="No. of constituents")
	contents[i] = n;
      else if (title=="p_{T,particle}^{leading} (GeV/c)")
	contents[i] = leadingpt;
      else 
	AliWarning(Form("Unable to fill dimension %s!",title.Data()));
    }

    fHistJetObservables->Fill(contents);
  }
}
