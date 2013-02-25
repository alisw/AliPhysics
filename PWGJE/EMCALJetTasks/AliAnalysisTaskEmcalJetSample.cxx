// $Id$
//
// Jet sample analysis task.
//
// Author: S.Aiola

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"

#include "AliAnalysisTaskEmcalJetSample.h"

ClassImp(AliAnalysisTaskEmcalJetSample)

//________________________________________________________________________
AliAnalysisTaskEmcalJetSample::AliAnalysisTaskEmcalJetSample() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetSample", kTRUE)

{
  // Default constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistTracksPt[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetSample::AliAnalysisTaskEmcalJetSample(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE)
{
  // Standard constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistTracksPt[i] = 0;
    fHistClustersPt[i] = 0;
    fHistLeadingJetPt[i] = 0;
    fHistJetsPhiEta[i] = 0;
    fHistJetsPtArea[i] = 0;
    fHistJetsPtLeadHad[i] = 0;
    fHistJetsCorrPtArea[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetSample::~AliAnalysisTaskEmcalJetSample()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSample::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  TString histname;

  for (Int_t i = 0; i < 4; i++) {
    if (!fTracksName.IsNull()) {
      histname = "fHistTracksPt_";
      histname += i;
      fHistTracksPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistTracksPt[i]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
      fHistTracksPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistTracksPt[i]);
    }

    if (!fCaloName.IsNull()) {
      histname = "fHistClustersPt_";
      histname += i;
      fHistClustersPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistClustersPt[i]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
      fHistClustersPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistClustersPt[i]);
    }

    if (!fJetsName.IsNull()) {
      histname = "fHistLeadingJetPt_";
      histname += i;
      fHistLeadingJetPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt);
      fHistLeadingJetPt[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistLeadingJetPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistLeadingJetPt[i]);
      
      histname = "fHistJetsPhiEta_";
      histname += i;
      fHistJetsPhiEta[i] = new TH2F(histname.Data(), histname.Data(), 50, -1, 1, 101, 0, TMath::Pi()*2 + TMath::Pi()/200);
      fHistJetsPhiEta[i]->GetXaxis()->SetTitle("#eta");
      fHistJetsPhiEta[i]->GetYaxis()->SetTitle("#phi");
      fOutput->Add(fHistJetsPhiEta[i]);
      
      histname = "fHistJetsPtArea_";
      histname += i;
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 30, 0, 3);
      fHistJetsPtArea[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtArea[i]->GetYaxis()->SetTitle("area");
      fOutput->Add(fHistJetsPtArea[i]);

      histname = "fHistJetsPtLeadHad_";
      histname += i;
      fHistJetsPtLeadHad[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistJetsPtLeadHad[i]->GetXaxis()->SetTitle("p_{T}^{raw} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetYaxis()->SetTitle("p_{T,lead} (GeV/c)");
      fHistJetsPtLeadHad[i]->GetZaxis()->SetTitle("counts");
      fOutput->Add(fHistJetsPtLeadHad[i]);
    
      if (!fRhoName.IsNull()) {
	histname = "fHistJetsCorrPtArea_";
	histname += i;
	fHistJetsCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 30, 0, 3);
	fHistJetsCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
	fHistJetsCorrPtArea[i]->GetYaxis()->SetTitle("area");
	fOutput->Add(fHistJetsCorrPtArea[i]);
      }
    }
  }
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSample::FillHistograms()
{
  // Fill histograms.

  if (fTracks) {
    const Int_t ntracks = fTracks->GetEntriesFast();
    
    for (Int_t it = 0; it < ntracks; it++) {
      AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(it));

      if (!track) {
	AliError(Form("Could not receive track %d", it));
	continue;
      }
     
      if (!AcceptTrack(track))
	continue;

      fHistTracksPt[fCentBin]->Fill(track->Pt()); 
    }
  }
  
  if (fCaloClusters) {
    const Int_t nclusters = fCaloClusters->GetEntriesFast();
    
    for (Int_t ic = 0; ic < nclusters; ic++) {
      AliVCluster *cluster = static_cast<AliVCluster*>(fCaloClusters->At(ic));
      
      if (!cluster) {
	AliError(Form("Could not receive cluster %d", ic));
	continue;
      }

      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      fHistClustersPt[fCentBin]->Fill(nPart.Pt());
    }
  }

  if (fJets) {
    static Int_t sortedJets[9999] = {-1};
    GetSortedArray(sortedJets, fJets);

    if (sortedJets[0]>=0) {
      AliEmcalJet* leadJet = static_cast<AliEmcalJet*>(fJets->At(sortedJets[0]));
      if (leadJet)
	fHistLeadingJetPt[fCentBin]->Fill(leadJet->Pt());
      else
	AliError("Could not retrieve leading jet!");
    }

    const Int_t njets = fJets->GetEntriesFast();
    for (Int_t ij = 0; ij < njets; ij++) {

      AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(ij));
      if (!jet) {
	AliError(Form("Could not receive jet %d", ij));
	continue;
      }  

      if (!AcceptJet(jet))
	continue;

      fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
      fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());

      Float_t ptLeading = GetLeadingHadronPt(jet);
      fHistJetsPtLeadHad[fCentBin]->Fill(jet->Pt(), ptLeading);

      if (fRho) {
	Float_t corrPt = jet->Pt() - fRhoVal * jet->Area();
	fHistJetsCorrPtArea[fCentBin]->Fill(corrPt, jet->Area());
      }
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSample::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSample::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
