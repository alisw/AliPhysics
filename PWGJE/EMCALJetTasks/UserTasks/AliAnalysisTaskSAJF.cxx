// $Id$
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

#include "AliAnalysisTaskSAJF.h"

ClassImp(AliAnalysisTaskSAJF)

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSAJF", kTRUE),
  fHistoType(1),
  fHistJetObservables(0)

{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistTracksJetPt[i] = 0;
    fHistClustersJetPt[i] = 0;
    fHistTracksPtDist[i] = 0;
    fHistClustersPtDist[i] = 0;

    fHistJetPtEtaPhi[i] = 0;
    fHistJetPtArea[i] = 0;
    fHistJetPtEP[i] = 0;
    fHistJetPtNEF[i] = 0;
    fHistJetPtZ[i] = 0;
    fHistJetPtLeadingPartPt[i] = 0;
    fHistJetCorrPtEtaPhi[i] = 0;
    fHistJetCorrPtArea[i] = 0;
    fHistJetCorrPtEP[i] = 0;
    fHistJetCorrPtNEF[i] = 0;
    fHistJetCorrPtZ[i] = 0;
    fHistJetCorrPtLeadingPartPt[i] = 0;
    fHistJetPtCorrPt[i] = 0;
    fHistJetPtMCPt[i] = 0;
    fHistJetMCPtCorrPt[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistoType(1),
  fHistJetObservables(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistTracksJetPt[i] = 0;
    fHistClustersJetPt[i] = 0;
    fHistTracksPtDist[i] = 0;
    fHistClustersPtDist[i] = 0;

    fHistJetPtEtaPhi[i] = 0;
    fHistJetPtArea[i] = 0;
    fHistJetPtEP[i] = 0;
    fHistJetPtNEF[i] = 0;
    fHistJetPtZ[i] = 0;
    fHistJetPtLeadingPartPt[i] = 0;
    fHistJetCorrPtEtaPhi[i] = 0;
    fHistJetCorrPtArea[i] = 0;
    fHistJetCorrPtEP[i] = 0;
    fHistJetCorrPtNEF[i] = 0;
    fHistJetCorrPtZ[i] = 0;
    fHistJetCorrPtLeadingPartPt[i] = 0;
    fHistJetPtCorrPt[i] = 0;
    fHistJetPtMCPt[i] = 0;
    fHistJetMCPtCorrPt[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::AllocateTHnSparse()
{
    TString title[20]= {""};
  Int_t nbins[20]  = {0};
  Double_t min[20] = {0.};
  Double_t max[20] = {0.};
  Int_t dim = 0;

  if (fForceBeamType != kpp) {
    title[dim] = "Centrality (%)";
    nbins[dim] = 22;
    min[dim] = -5;
    max[dim] = 105;
    dim++;
    
    title[dim] = "#phi_{jet} - #psi_{RP}";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = TMath::Pi();
    dim++;
  }

  title[dim] = "#eta";
  nbins[dim] = 100;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "#phi_{jet} (rad)";
  nbins[dim] = 201;
  min[dim] = 0;
  max[dim] = 2*TMath::Pi()*nbins[dim]/(nbins[dim]-1);
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
  nbins[dim] = 150;
  min[dim] = 0;
  max[dim] = 1.5;
  dim++;

  title[dim] = "NEF";
  nbins[dim] = 102;
  min[dim] = 0;
  max[dim] = 1.02;
  dim++;

  title[dim] = "Z";
  nbins[dim] = 102;
  min[dim] = 0;
  max[dim] = 1.02;
  dim++;

  title[dim] = "No. of constituents";
  nbins[dim] = 250;
  min[dim] = -0.5;
  max[dim] = 249.5;
  dim++;

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
  for (Int_t i = 0; i < 4; i++) {
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

  if (fHistoType == 0)
    AllocateTHX();
  else
    AllocateTHnSparse();

  for (Int_t i = 0; i < 4; i++) {
    TString histname;

    if (fParticleCollArray.GetEntriesFast()>0) {
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
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
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
    }
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram

}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::FillHistograms()
{
  // Fill histograms.

  if (!fJets) {
    AliError(Form("%s - Jet array not provided, returning...", GetName()));
    return kFALSE;
  }

  for (Int_t ij = 0; ij < fJets->GetEntriesFast(); ij++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }

    if (!AcceptJet(jet))
      continue;

    Float_t ptLeading = GetLeadingHadronPt(jet);
    Float_t corrPt = jet->Pt() - fRhoVal * jet->Area();

    // Fill THnSparse
    Double_t ep = jet->Phi() - fEPV0;
    while (ep < 0) ep += TMath::Pi();
    while (ep >= TMath::Pi()) ep -= TMath::Pi();

    FillJetHisto(fCent, ep, jet->Eta(), jet->Phi(), jet->Pt(), jet->MCPt(), corrPt, jet->Area(), 
		 jet->NEF(), ptLeading/jet->Pt(), jet->GetNumberOfConstituents(), ptLeading);

    if (fTracks) {
      for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
	AliVParticle *track = jet->TrackAt(it, fTracks);
	if (track) {
	  fHistTracksJetPt[fCentBin]->Fill(track->Pt(), jet->Pt());
	  Double_t dist = TMath::Sqrt((track->Eta() - jet->Eta()) * (track->Eta() - jet->Eta()) + (track->Phi() - jet->Phi()) * (track->Phi() - jet->Phi()));
	  fHistTracksPtDist[fCentBin]->Fill(track->Pt(), dist);
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
      else if (title=="#eta")
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
