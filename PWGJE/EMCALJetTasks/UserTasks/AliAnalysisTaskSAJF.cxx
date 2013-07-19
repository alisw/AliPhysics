// $Id$
//
// Jet analysis task.
//
// Author: S.Aiola

#include <TClonesArray.h>
#include <TH2F.h>
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
  fHistJetObservables(0)

{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistTracksJetPt[i] = 0;
    fHistClustersJetPt[i] = 0;
    fHistTracksPtDist[i] = 0;
    fHistClustersPtDist[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistJetObservables(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistTracksJetPt[i] = 0;
    fHistClustersJetPt[i] = 0;
    fHistTracksPtDist[i] = 0;
    fHistClustersPtDist[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  TString title[20]= {""};
  Int_t nbins[20]  = {0};
  Double_t min[20] = {0.};
  Double_t max[20] = {0.};
  Int_t dim = 0;

  if (fForceBeamType != kpp) {
    title[dim] = "Centrality (%)";
    nbins[dim] = 101;
    min[dim] = 0;
    max[dim] = 101;
    dim++;
    
    title[dim] = "#phi_{jet} - #psi_{RP}";
    nbins[dim] = fNbins/2;
    min[dim] = -TMath::Pi()/2;
    max[dim] = TMath::Pi()*5/2;
    dim++;
  }

  title[dim] = "#eta";
  nbins[dim] = fNbins/2;
  min[dim] = -1;
  max[dim] = 1;
  dim++;

  title[dim] = "#phi_{jet} (rad)";
  nbins[dim] = fNbins/2+1;
  min[dim] = 0;
  max[dim] = 2*TMath::Pi()*nbins[dim]/(nbins[dim]-1);
  dim++;

  title[dim] = "p_{T} (GeV/c)";
  nbins[dim] = fNbins;
  min[dim] = 0;
  max[dim] = 250;
  dim++;

  if (fIsEmbedded) {
    title[dim] = "p_{T}^{MC} (GeV/c)";
    nbins[dim] = fNbins;
    min[dim] = 0;
    max[dim] = 250;
    dim++;
  }

  if (!fRhoName.IsNull()) {
    title[dim] = "p_{T}^{corr} (GeV/c)";
    nbins[dim] = fNbins*2;
    min[dim] = -250;
    max[dim] = 250;
    dim++;
  }

  title[dim] = "A_{jet}";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 1;
  dim++;

  title[dim] = "NEF";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 1.25;
  dim++;

  title[dim] = "Z";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 1.25;
  dim++;

  title[dim] = "No. of constituents";
  nbins[dim] = 250;
  min[dim] = -0.5;
  max[dim] = 249.5;
  dim++;

  title[dim] = "p_{T,particle}^{leading} (GeV/c)";
  nbins[dim] = fNbins/2;
  min[dim] = 0;
  max[dim] = 125;
  dim++;

  fHistJetObservables = new THnSparseD("fHistJetObservables","fHistJetObservables",dim,nbins,min,max);
  fOutput->Add(fHistJetObservables);
  for (Int_t i = 0; i < dim; i++)
    fHistJetObservables->GetAxis(i)->SetTitle(title[i]);

  for (Int_t i = 0; i < 4; i++) {
    TString histname;

    if (!fTracksName.IsNull()) {
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

    if (!fCaloName.IsNull()) {
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
    FillJetHisto(fCent, fEPV0, jet->Eta(), jet->Phi(), jet->Pt(), jet->MCPt(), corrPt, jet->Area(), 
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
  Double_t contents[20]={0};

  for (Int_t i = 0; i < fHistJetObservables->GetNdimensions(); i++) {
    TString title(fHistJetObservables->GetAxis(i)->GetTitle());
    if (title=="Centrality (%)")
      contents[i] = cent;
    else if (title=="#phi_{jet} - #psi_{RP}")
      contents[i] = phi - ep;
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
