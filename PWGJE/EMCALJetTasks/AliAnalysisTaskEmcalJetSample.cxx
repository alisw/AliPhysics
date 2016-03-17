// $Id$
//
// Jet sample analysis task.
//
// Author: S.Aiola

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskEmcalJetSample.h"

ClassImp(AliAnalysisTaskEmcalJetSample)

//________________________________________________________________________
AliAnalysisTaskEmcalJetSample::AliAnalysisTaskEmcalJetSample() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetSample", kTRUE),
  fHistTracksPt(0),
  fHistNTracks(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistClustDx(0),
  fHistClustDz(0),
  fNAccJets(0),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0)

{
  // Default constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistNTracks        = new TH1*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistNTracks[i] = 0;
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
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistTracksPt(0),
  fHistNTracks(0),
  fHistClustersPt(0),
  fHistLeadingJetPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsPtLeadHad(0),
  fHistJetsCorrPtArea(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fHistClustDx(0),
  fHistClustDz(0),
  fNAccJets(0),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0)
{
  // Standard constructor.

  fHistTracksPt       = new TH1*[fNcentBins];
  fHistNTracks        = new TH1*[fNcentBins];
  fHistClustersPt     = new TH1*[fNcentBins];
  fHistLeadingJetPt   = new TH1*[fNcentBins];
  fHistJetsPhiEta     = new TH2*[fNcentBins];
  fHistJetsPtArea     = new TH2*[fNcentBins];
  fHistJetsPtLeadHad  = new TH2*[fNcentBins];
  fHistJetsCorrPtArea = new TH2*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fHistTracksPt[i] = 0;
    fHistNTracks[i] = 0;
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

  fJetsCont           = GetJetContainer(0);
  if(fJetsCont) { //get particles and clusters connected to jets
    fTracksCont       = fJetsCont->GetParticleContainer();
    fCaloClustersCont = fJetsCont->GetClusterContainer();
  } else {        //no jets, just analysis tracks and clusters
    fTracksCont       = GetParticleContainer(0);
    fCaloClustersCont = GetClusterContainer(0);
  }
  if(fTracksCont) fTracksCont->SetClassName("AliVTrack");
  if(fCaloClustersCont) fCaloClustersCont->SetClassName("AliVCluster");

  TString histname;

  for (Int_t i = 0; i < fNcentBins; i++) {
    if (fParticleCollArray.GetEntriesFast()>0) {
      histname = "fHistTracksPt_";
      histname += i;
      fHistTracksPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistTracksPt[i]->GetXaxis()->SetTitle("p_{T,track} (GeV/c)");
      fHistTracksPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistTracksPt[i]);
      
      histname = "fHistNTracks_";
      histname += i;
      fHistNTracks[i] = new TH1F(histname.Data(), histname.Data(), 200, 0., 199.);
      fOutput->Add(fHistNTracks[i]);
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
      histname = "fHistClustersPt_";
      histname += i;
      fHistClustersPt[i] = new TH1F(histname.Data(), histname.Data(), fNbins / 2, fMinBinPt, fMaxBinPt / 2);
      fHistClustersPt[i]->GetXaxis()->SetTitle("p_{T,clus} (GeV/c)");
      fHistClustersPt[i]->GetYaxis()->SetTitle("counts");
      fOutput->Add(fHistClustersPt[i]);
    }

    if (fJetCollArray.GetEntriesFast()>0) {
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
      fHistJetsPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins, fMinBinPt, fMaxBinPt, 50, 0, 1);
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
    
      if (!(GetJetContainer()->GetRhoName().IsNull())) {
	histname = "fHistJetsCorrPtArea_";
	histname += i;
	fHistJetsCorrPtArea[i] = new TH2F(histname.Data(), histname.Data(), fNbins*2, -fMaxBinPt, fMaxBinPt, 50, 0, 1);
	fHistJetsCorrPtArea[i]->GetXaxis()->SetTitle("p_{T}^{corr} [GeV/c]");
	fHistJetsCorrPtArea[i]->GetYaxis()->SetTitle("area");
	fOutput->Add(fHistJetsCorrPtArea[i]);
      }
    }
  }

  histname = "fHistPtDEtaDPhiTrackClus";
  fHistPtDEtaDPhiTrackClus = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{track};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiTrackClus);

  histname = "fHistPtDEtaDPhiClusTrack";
  fHistPtDEtaDPhiClusTrack = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{clus};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiClusTrack);

  fHistClustDx = new TH1F("fHistClustDx","fHistClustDx;Dx",1000,0.,1.);
  fOutput->Add(fHistClustDx);

  fHistClustDz = new TH1F("fHistClustDz","fHistClustDz;Dz",1000,0.,1.);
  fOutput->Add(fHistClustDz);

  fNAccJets = new TH1F("fNAccJets","fNAccJets;N/ev",11,-0.5, 9.5);
  fOutput->Add(fNAccJets);
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetSample::FillHistograms()
{
  // Fill histograms.

  if (fTracksCont) {
     fHistNTracks[fCentBin]->Fill(fTracksCont->GetNAcceptedParticles());
     fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    while(track) {
      fHistTracksPt[fCentBin]->Fill(track->Pt()); 
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
      
    }
  }

  if (fCaloClustersCont) {
    fCaloClustersCont->ResetCurrentID();
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(); 
    while(cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      fHistClustersPt[fCentBin]->Fill(nPart.Pt());
      Double_t dx = cluster->GetTrackDx();
      Double_t dz = cluster->GetTrackDz();
      fHistClustDx->Fill(dx);
      fHistClustDz->Fill(dz);
      cluster = fCaloClustersCont->GetNextAcceptCluster();
    }
  }

  if (fJetsCont) {
     Int_t count = 0;
    fJetsCont->ResetCurrentID();
    AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(); 
    while(jet) {
       count++;
      fHistJetsPtArea[fCentBin]->Fill(jet->Pt(), jet->Area());
      fHistJetsPhiEta[fCentBin]->Fill(jet->Eta(), jet->Phi());

      Float_t ptLeading = fJetsCont->GetLeadingHadronPt(jet);
      fHistJetsPtLeadHad[fCentBin]->Fill(jet->Pt(), ptLeading);

      if (fHistJetsCorrPtArea[fCentBin]) {
	Float_t corrPt = jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
	fHistJetsCorrPtArea[fCentBin]->Fill(corrPt, jet->Area());
      }
      jet = fJetsCont->GetNextAcceptJet(); 
    }
    fNAccJets->Fill(count);
    jet = fJetsCont->GetLeadingJet();
    if(jet) fHistLeadingJetPt[fCentBin]->Fill(jet->Pt());
  }

  CheckClusTrackMatching();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSample::CheckClusTrackMatching()
{
  
  if(!fTracksCont || !fCaloClustersCont)
    return;

  Double_t deta = 999;
  Double_t dphi = 999;

  //Get closest cluster to track
  fTracksCont->ResetCurrentID();
  AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()); 
  while(track) {
    //Get matched cluster
    Int_t emc1 = track->GetEMCALcluster();
    if(fCaloClustersCont && emc1>=0) {
      AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
      if(clusMatch) {
	AliPicoTrack::GetEtaPhiDiff(track, clusMatch, dphi, deta);
	fHistPtDEtaDPhiTrackClus->Fill(track->Pt(),deta,dphi);
      }
    }
    track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
  }
  
  //Get closest track to cluster
  fCaloClustersCont->ResetCurrentID();
  AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(); 
  while(cluster) {
    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);
    fHistClustersPt[fCentBin]->Fill(nPart.Pt());
    
    //Get matched track
    AliVTrack *mt = NULL;      
    AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
    if(acl) {
      if(acl->GetNTracksMatched()>1)
	mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
    }
    else {
      AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
      Int_t im = ecl->GetTrackMatchedIndex();
      if(fTracksCont && im>=0) {
	mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
      }
    }
    if(mt) {
      AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
      fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);
      
      /* //debugging
	 if(mt->IsEMCAL()) {
	 Int_t emc1 = mt->GetEMCALcluster();
	 Printf("current id: %d  emc1: %d",fCaloClustersCont->GetCurrentID(),emc1);
	 AliVCluster *clm = fCaloClustersCont->GetCluster(emc1);
	 AliPicoTrack::GetEtaPhiDiff(mt, clm, dphi, deta);
	 Printf("deta: %f dphi: %f",deta,dphi);
	 }
      */
    }
    cluster = fCaloClustersCont->GetNextAcceptCluster();
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetSample::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

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
