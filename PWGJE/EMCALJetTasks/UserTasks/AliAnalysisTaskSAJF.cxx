// $Id$
//
// Jet analysis task.
//
// Author: S.Aiola

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"

#include "AliAnalysisTaskSAJF.h"

ClassImp(AliAnalysisTaskSAJF)

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskSAJF", kTRUE),
  fNjetsVsCent(0)

{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistEvents[i] = 0;
    fHistLeadingJetPhiEta[i] = 0;
    fHistLeadingJetPtArea[i] = 0;
    fHistLeadingJetCorrPtArea[i] = 0;
    fHistRhoVSleadJetPt[i] = 0;
    fHistJetPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
    fHistJetsNEFvsPt[i] = 0;
    fHistJetsCEFvsCEFPt[i] = 0;
    fHistJetsZvsPt[i] = 0;
    fHistConstituents[i] = 0;
    fHistTracksJetPt[i] = 0;
    fHistClustersJetPt[i] = 0;
    fHistTracksPtDist[i] = 0;
    fHistClustersPtDist[i] = 0;
    fHistJetNconstVsPt[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAJF::AliAnalysisTaskSAJF(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fNjetsVsCent(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistEvents[i] = 0;
    fHistLeadingJetPhiEta[i] = 0;
    fHistLeadingJetPtArea[i] = 0;
    fHistLeadingJetCorrPtArea[i] = 0;
    fHistRhoVSleadJetPt[i] = 0;
    fHistJetPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
    fHistJetsNEFvsPt[i] = 0;
    fHistJetsCEFvsCEFPt[i] = 0;
    fHistJetsZvsPt[i] = 0;
    fHistConstituents[i] = 0;
    fHistTracksJetPt[i] = 0;
    fHistClustersJetPt[i] = 0;
    fHistTracksPtDist[i] = 0;
    fHistClustersPtDist[i] = 0;
    fHistJetNconstVsPt[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskSAJF::~AliAnalysisTaskSAJF()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fNjetsVsCent = new TH2F("fNjetsVsCent","fNjetsVsCent", 100, 0, 100, 150, -0.5, 149.5);
  fNjetsVsCent->GetXaxis()->SetTitle("Centrality (%)");
  fNjetsVsCent->GetYaxis()->SetTitle("# of jets");
  fOutput->Add(fNjetsVsCent);

  const Int_t nbinsZ = 12;
  Float_t binsZ[nbinsZ+1] = {0,1,2,3,4,5,6,7,8,9,10,20,1000};

  Float_t *binsPt       = GenerateFixedBinArray(fNbins, fMinBinPt, fMaxBinPt);
  Float_t *binsCorrPt   = GenerateFixedBinArray(fNbins*2, -fMaxBinPt, fMaxBinPt);
  Float_t *binsArea     = GenerateFixedBinArray(30, 0, fJetRadius * fJetRadius * TMath::Pi() * 3);
  Float_t *binsEta      = GenerateFixedBinArray(50,-1, 1);
  Float_t *binsPhi      = GenerateFixedBinArray(101, 0, TMath::Pi() * 2.02);
  Float_t *bins120      = GenerateFixedBinArray(120, 0, 1.2);
  Float_t *bins150      = GenerateFixedBinArray(150, -0.5, 149.5);

  TString histname;

  for (Int_t i = 0; i < 4; i++) {
    histname = "fHistEvents_";
    histname += i;
    fHistEvents[i] = new TH1F(histname,histname, 6, 0, 6);
    fHistEvents[i]->GetXaxis()->SetTitle("Event state");
    fHistEvents[i]->GetYaxis()->SetTitle("counts");
    fHistEvents[i]->GetXaxis()->SetBinLabel(1, "No jets");
    fHistEvents[i]->GetXaxis()->SetBinLabel(2, "Max jet not found");
    fHistEvents[i]->GetXaxis()->SetBinLabel(3, "Rho == 0");
    fHistEvents[i]->GetXaxis()->SetBinLabel(4, "Max jet <= 0");
    fHistEvents[i]->GetXaxis()->SetBinLabel(5, "OK");
    fOutput->Add(fHistEvents[i]);

    histname = "fHistLeadingJetPhiEta_";
    histname += i;
    fHistLeadingJetPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 50,-1, 1, 101, 0, TMath::Pi() * 2.02);
    fHistLeadingJetPhiEta[i]->GetXaxis()->SetTitle("#eta");
    fHistLeadingJetPhiEta[i]->GetYaxis()->SetTitle("#phi");
    fOutput->Add(fHistLeadingJetPhiEta[i]);

    histname = "fHistLeadingJetPtArea_";
    histname += i;
    fHistLeadingJetPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 30, 0, fJetRadius * fJetRadius * TMath::Pi() * 3);
    fHistLeadingJetPtArea[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
    fHistLeadingJetPtArea[i]->GetYaxis()->SetTitle("area");
    fOutput->Add(fHistLeadingJetPtArea[i]);

    if (!fRhoName.IsNull()) {
      histname = "fHistLeadingJetCorrPtArea_";
      histname += i;
      fHistLeadingJetCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins * 2, -fMaxBinPt, fMaxBinPt, 30, 0, fJetRadius * fJetRadius * TMath::Pi() * 3);
      fHistLeadingJetCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
      fHistLeadingJetCorrPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistLeadingJetCorrPtArea[i]);
      
      histname = "fHistRhoVSleadJetPt_";
      histname += i;
      fHistRhoVSleadJetPt[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt*2, fNbins, fMinBinPt, fMaxBinPt);
      fHistRhoVSleadJetPt[i]->GetXaxis()->SetTitle("#rho * area (GeV/c)");
      fHistRhoVSleadJetPt[i]->GetYaxis()->SetTitle("Leading jet p_{T} (GeV/c)");
      fOutput->Add(fHistRhoVSleadJetPt[i]);
    }

    histname = "fHistJetPhiEta_";
    histname += i;
    fHistJetPhiEta[i] = new TH3F(histname.Data(), histname.Data(), 
				 50, binsEta, 
				 101, binsPhi, 
				 nbinsZ, binsZ);
    fHistJetPhiEta[i]->GetXaxis()->SetTitle("#eta");
    fHistJetPhiEta[i]->GetYaxis()->SetTitle("#phi");
    fHistJetPhiEta[i]->GetZaxis()->SetTitle("p_{T,lead} (GeV/c)");
    fOutput->Add(fHistJetPhiEta[i]);
    
    histname = "fHistJetsPtArea_";
    histname += i;
    fHistJetsPtArea[i] = new TH3F(histname.Data(), histname.Data(), 
				  fNbins, binsPt, 
				  30, binsArea,
				  nbinsZ, binsZ);
    fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
    fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
    fHistJetsPtArea[i]->GetZaxis()->SetTitle("p_{T,lead} (GeV/c)");
    fOutput->Add(fHistJetsPtArea[i]);

    histname = "fHistJetsZvsPt_";
    histname += i;
    fHistJetsZvsPt[i] = new TH3F(histname.Data(), histname.Data(), 
				 120, bins120, 
				 fNbins, binsPt,
				 nbinsZ, binsZ);
    fHistJetsZvsPt[i]->GetXaxis()->SetTitle("Z");
    fHistJetsZvsPt[i]->GetYaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
    fHistJetsZvsPt[i]->GetZaxis()->SetTitle("p_{T,lead} (GeV/c)");
    fOutput->Add(fHistJetsZvsPt[i]);

    histname = "fHistJetsNEFvsPt_";
    histname += i;
    fHistJetsNEFvsPt[i] = new TH3F(histname.Data(), histname.Data(), 
				   120, bins120, 
				   fNbins, binsPt,
				   nbinsZ, binsZ);
    fHistJetsNEFvsPt[i]->GetXaxis()->SetTitle("NEF");
    fHistJetsNEFvsPt[i]->GetYaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
    fHistJetsNEFvsPt[i]->GetZaxis()->SetTitle("p_{T,lead} (GeV/c)");
    fOutput->Add(fHistJetsNEFvsPt[i]);
    
    histname = "fHistJetsCEFvsCEFPt_";
    histname += i;
    fHistJetsCEFvsCEFPt[i] = new TH3F(histname.Data(), histname.Data(), 
				      120, bins120, 
				      fNbins, binsPt,
				      nbinsZ, binsZ);
    fHistJetsCEFvsCEFPt[i]->GetXaxis()->SetTitle("1-NEF");
    fHistJetsCEFvsCEFPt[i]->GetYaxis()->SetTitle("(1-NEF)*p_{T}^{raw} (GeV/c)");
    fHistJetsCEFvsCEFPt[i]->GetZaxis()->SetTitle("p_{T,lead} (GeV/c)");
    fOutput->Add(fHistJetsCEFvsCEFPt[i]);

    if (!fRhoName.IsNull()) {
      histname = "fHistJetsCorrPtArea_";
      histname += i;
      fHistJetsCorrPtArea[i] = new TH3F(histname.Data(), histname.Data(), 
					fNbins * 2, binsCorrPt, 
					30, binsArea,
					nbinsZ, binsZ);
      fHistJetsCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
      fHistJetsCorrPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsCorrPtArea[i]);
    }

    histname = "fHistConstituents_";
    histname += i;
    fHistConstituents[i] = new TH2F(histname.Data(), histname.Data(), 100, 1, 101, 100, -0.5, 99.5);
    fHistConstituents[i]->GetXaxis()->SetTitle("p_{T,part} (GeV/c)");
    fHistConstituents[i]->GetYaxis()->SetTitle("no. of particles");
    fHistConstituents[i]->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistConstituents[i]);

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

    histname = "fHistJetNconstVsPt_";
    histname += i;
    fHistJetNconstVsPt[i] = new TH3F(histname.Data(), histname.Data(), 150, bins150, fNbins, binsPt, nbinsZ, binsZ);
    fHistJetNconstVsPt[i]->GetXaxis()->SetTitle("# of constituents");
    fHistJetNconstVsPt[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
    fHistJetNconstVsPt[i]->GetZaxis()->SetTitle("p_{T,lead} (GeV/c)");
    fOutput->Add(fHistJetNconstVsPt[i]);
  }

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram

  delete[] binsPt;
  delete[] binsCorrPt;
  delete[] binsArea;
  delete[] binsEta;
  delete[] binsPhi;
  delete[] bins120;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSAJF::FillHistograms()
{
  // Fill histograms.

  if (!fJets) {
    AliError(Form("%s - Jet array not provided, returning...", GetName()));
    return kFALSE;
  }

  if (fJets->GetEntriesFast() < 1) { // no jets in array, skipping
    fHistEvents[fCentBin]->Fill("No jets", 1);
    return kTRUE;
  }

  static Int_t sortedJets[9999] = {-1};
  GetSortedArray(sortedJets, fJets, fRhoVal);
  
  if (sortedJets[0] < 0) { // no accepted jets, skipping
    fHistEvents[fCentBin]->Fill("No jets", 1);
    return kTRUE;
  }

  // OK, event accepted

  if (fRhoVal == 0) 
    fHistEvents[fCentBin]->Fill("Rho == 0", 1);
  else
    fHistEvents[fCentBin]->Fill("OK", 1);

  for (Int_t i = 0; i < fNLeadingJets && i < fJets->GetEntriesFast(); i++) {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(sortedJets[i]));

    if (!jet) {
      AliError(Form("Could not receive jet %d", sortedJets[i]));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    fHistLeadingJetPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());
    fHistLeadingJetPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());

    Float_t corrPt = jet->Pt() - fRhoVal * jet->Area();

    if (fHistLeadingJetCorrPtArea[fCentBin])
      fHistLeadingJetCorrPtArea[fCentBin]->Fill(corrPt, jet->Area());

    if (i==0 && fHistRhoVSleadJetPt[fCentBin]) 
      fHistRhoVSleadJetPt[fCentBin]->Fill(fRhoVal, jet->Pt());
  }

  Int_t njets = DoJetLoop();

  fNjetsVsCent->Fill(fCent, njets);

  return kTRUE;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSAJF::DoJetLoop()
{
  // Do the jet loop.

  if (!fJets)
    return 0;

  const Int_t njets = fJets->GetEntriesFast();
  Int_t nAccJets = 0;

  TH1F constituents("constituents", "constituents", 
		    fHistConstituents[0]->GetNbinsX(), fHistConstituents[0]->GetXaxis()->GetXmin(), fHistConstituents[0]->GetXaxis()->GetXmax()); 

  for (Int_t ij = 0; ij < njets; ij++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(ij));

    if (!jet) {
      AliError(Form("Could not receive jet %d", ij));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    Float_t corrPt = jet->Pt() - fRhoVal * jet->Area();

    Float_t ptLeading = GetLeadingHadronPt(jet);

    fHistJetPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi(), ptLeading);
    fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area(), ptLeading);
    if (fHistJetsCorrPtArea[fCentBin])
      fHistJetsCorrPtArea[fCentBin]->Fill(corrPt, jet->Area(), ptLeading);
    fHistJetNconstVsPt[fCentBin]->Fill(jet->GetNumberOfConstituents(), jet->Pt(), ptLeading);

    fHistJetsNEFvsPt[fCentBin]->Fill(jet->NEF(), jet->Pt(), ptLeading);
    fHistJetsCEFvsCEFPt[fCentBin]->Fill(1-jet->NEF(), (1-jet->NEF())*jet->Pt(), ptLeading);

    if (fTracks) {
      for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
	AliVParticle *track = jet->TrackAt(it, fTracks);
	if (track) {
	  fHistJetsZvsPt[fCentBin]->Fill(track->Pt() / jet->Pt(), jet->Pt(), ptLeading);
	  constituents.Fill(track->Pt());
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
	  fHistJetsZvsPt[fCentBin]->Fill(nPart.Et() / jet->Pt(), jet->Pt(), ptLeading);
	  constituents.Fill(nPart.Pt());
	  fHistClustersJetPt[fCentBin]->Fill(nPart.Pt(), jet->Pt());
	  Double_t dist = TMath::Sqrt((nPart.Eta() - jet->Eta()) * (nPart.Eta() - jet->Eta()) + (nPart.Phi() - jet->Phi()) * (nPart.Phi() - jet->Phi()));
	  fHistClustersPtDist[fCentBin]->Fill(nPart.Pt(), dist);
	}
      }
    }

    for (Int_t i = 1; i <= constituents.GetNbinsX(); i++) {
      fHistConstituents[fCentBin]->Fill(constituents.GetBinCenter(i), constituents.GetBinContent(i));
    }

    constituents.Reset();
    nAccJets++;
  } //jet loop 

  return nAccJets;
}

//________________________________________________________________________
void AliAnalysisTaskSAJF::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
