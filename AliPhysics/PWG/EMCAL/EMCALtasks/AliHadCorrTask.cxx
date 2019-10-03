//
// Hadronic correction task.
//
// Author: R.Reed, C.Loizides, S.Aiola

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliVEventHandler.h"
#include "AliEMCALGeometry.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliHadCorrTask.h"

ClassImp(AliHadCorrTask)

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask() : 
  AliAnalysisTaskEmcal("AliHadCorrTask", kFALSE),
  fOutCaloName(),
  fPhiMatch(0.05),
  fEtaMatch(0.025),
  fDoTrackClus(kTRUE),
  fHadCorr(0),
  fEexclCell(0),
  fDoExact(kFALSE),
  fEsdMode(kTRUE),
  fOutClusters(0),
  fHistMatchEtaPhiAll(0),
  fHistMatchEtaPhiAllTr(0),
  fHistMatchEtaPhiAllCl(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNClusMatchCent(0)
{
  // Default constructor.

  for(Int_t i=0; i<10; i++) {
    fHistEsubPch[i]    = 0;
    fHistEsubPchRat[i] = 0;

    if (i<5) {
      fHistEsubPchRatAll[i] = 0;

      fHistMatchEvsP[i]    = 0;
      fHistMatchdRvsEP[i]  = 0;
      fHistNMatchEnergy[i] = 0;

      fHistEmbTrackMatchesOversub[i] = 0;
      fHistNonEmbTrackMatchesOversub[i] = 0;
      fHistOversubMCClusters[i] = 0;
      fHistOversubNonMCClusters[i] = 0;
      fHistOversub[i] = 0;

      for(Int_t j=0; j<4; j++)
        fHistNCellsEnergy[i][j] = 0;
    }
    
    for(Int_t j=0; j<9; j++) {
      for(Int_t k=0; k<2; k++) 
        fHistMatchEtaPhi[i][j][k] = 0;
    }
  } 
}

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fOutCaloName("CaloClustersCorr"),
  fPhiMatch(0.05),
  fEtaMatch(0.025),
  fDoTrackClus(kTRUE),
  fHadCorr(0),
  fEexclCell(0),
  fDoExact(kFALSE),
  fEsdMode(kTRUE),
  fOutClusters(0),
  fHistMatchEtaPhiAll(0),
  fHistMatchEtaPhiAllTr(0),
  fHistMatchEtaPhiAllCl(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNClusMatchCent(0)
{
  // Standard constructor.

  for(Int_t i=0; i<10; i++) {
    fHistEsubPch[i]    = 0;
    fHistEsubPchRat[i] = 0;

    if (i<5) {
      fHistEsubPchRatAll[i] = 0;

      fHistMatchEvsP[i]    = 0;
      fHistMatchdRvsEP[i]  = 0;
      fHistNMatchEnergy[i] = 0;

      fHistEmbTrackMatchesOversub[i] = 0;
      fHistNonEmbTrackMatchesOversub[i] = 0;
      fHistOversubMCClusters[i] = 0;
      fHistOversubNonMCClusters[i] = 0;
      fHistOversub[i] = 0;

      for(Int_t j=0; j<4; j++)
        fHistNCellsEnergy[i][j] = 0;
    }
    
    for(Int_t j=0; j<9; j++) {
      for(Int_t k=0; k<2; k++) 
        fHistMatchEtaPhi[i][j][k] = 0;
    }
  } 
  
  SetMakeGeneralHistograms(histo);

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}

//________________________________________________________________________
AliHadCorrTask::~AliHadCorrTask()
{
  // Destructor
}

//________________________________________________________________________
UInt_t AliHadCorrTask::GetMomBin(Double_t p) const
{
  // Get momenum bin.

  UInt_t pbin=0;
  if (p<0.5) 
    pbin=0;
  else if (p>=0.5 && p<1.0) 
    pbin=1;
  else if (p>=1.0 && p<1.5) 
    pbin=2;
  else if (p>=1.5 && p<2.) 
    pbin=3;
  else if (p>=2. && p<3.) 
    pbin=4;
  else if (p>=3. && p<4.) 
    pbin=5;
  else if (p>=4. && p<5.) 
    pbin=6;
  else if (p>=5. && p<8.) 
    pbin=7;
  else 
    pbin=8;

  return pbin;
}

//________________________________________________________________________
Double_t AliHadCorrTask::GetEtaSigma(Int_t pbin) const
{
  // Get sigma in eta.

  Double_t EtaSigma[9]={0.0097,0.0075,0.0059,0.0055,0.0053,0.005,0.005,0.0045,0.0042};
  return 2.0*EtaSigma[pbin];
}

//________________________________________________________________________
Double_t AliHadCorrTask::GetPhiMean(Int_t pbin, Int_t centbin) const
{
  // Get phi mean.

  if (centbin==0) { 
    Double_t PhiMean[9]={0.036,
			 0.021,
			 0.0121,
			 0.0084,
			 0.0060,
			 0.0041,
			 0.0031,
			 0.0022,
			 0.001};
    return PhiMean[pbin];
  } else if (centbin==1) { 
    Double_t PhiMean[9]={0.036,
			 0.021,
			 0.0121,
			 0.0084,
			 0.0060,
			 0.0041,
			 0.0031,
			 0.0022,
			 0.001};
    return PhiMean[pbin];
  } else if (centbin==2) { 
    Double_t PhiMean[9]={0.036,
			 0.021,
			 0.0121,
			 0.0084,
			 0.0060,
			 0.0041,
			 0.0031,
			 0.0022,
			 0.001};
    return PhiMean[pbin];
  } else if (centbin==3) { 
    Double_t PhiMean[9]={0.036,
			 0.021,
			 0.0121,
			 0.0084,
			 0.0060,
			 0.0041,
			 0.0031,
			 0.0022,
			 0.001};  
    return PhiMean[pbin];
  } else if (centbin==4) { 
    Double_t PhiMean[9]={0.036,
			 0.021,
			 0.0127,
			 0.0089,
			 0.0068,
			 0.0049,
			 0.0038,
			 0.0028,
			 0.0018};
    return PhiMean[pbin]*(-1.);
  } else if (centbin==5) { 
    Double_t PhiMean[9]={0.036,
			 0.021,
			 0.0127,
			 0.0089,
			 0.0068,
			 0.0048,
			 0.0038,
			 0.0028,
			 0.0018};
    return PhiMean[pbin]*(-1.);
  } else if (centbin==6) { 
    Double_t PhiMean[9]={0.036,
			 0.021,
			 0.0127,
			 0.0089,
			 0.0068,
			 0.0045,
			 0.0035,
			 0.0028,
			 0.0018};
    return PhiMean[pbin]*(-1.);
  } else if (centbin==7) { 
    Double_t PhiMean[9]={0.036,
			 0.021,
			 0.0127,
			 0.0089,
			 0.0068,
			 0.0043,
			 0.0035,
			 0.0028,
			 0.0018};
    return PhiMean[pbin]*(-1.);
  }

  return 0;
}

//________________________________________________________________________
Double_t AliHadCorrTask::GetPhiSigma(Int_t pbin, Int_t centbin) const
{
  // Get phi sigma.

  if (centbin==0) { 
    Double_t PhiSigma[9]={0.0221,
			  0.0128,
			  0.0074,
			  0.0064,
			  0.0059,
			  0.0055,
			  0.0052,
			  0.0049,
			  0.0045};
    return 2.*PhiSigma[pbin];
  } else if (centbin==1) { 
    Double_t PhiSigma[9]={0.0217,
			  0.0120,
			  0.0076,
			  0.0066,
			  0.0062,
			  0.0058,
			  0.0054,
			  0.0054,
0.0045};
    return 2.*PhiSigma[pbin];
  } else if (centbin==2) { 
    Double_t PhiSigma[9]={0.0211,
			  0.0124,
			  0.0080,
			  0.0070,
			  0.0067,
			  0.0061,
			  0.0059,
			  0.0054,
			  0.0047};
    return 2.*PhiSigma[pbin];
  } else if (centbin==3) { 
    Double_t PhiSigma[9]={0.0215,
			  0.0124,
			  0.0082,
			  0.0073,
			  0.0069,
			  0.0064,
			  0.0060,
			  0.0055,
			  0.0047};  
    return 2.*PhiSigma[pbin];
  } else if (centbin==4) { 
    Double_t PhiSigma[9]={0.0199,
			  0.0108,
			  0.0072,
			  0.0071,
			  0.0060,
			  0.0055,
			  0.0052,
			  0.0049,
			  0.0045};
    return 2.*PhiSigma[pbin];
  } else if (centbin==5) { 
    Double_t PhiSigma[9]={0.0200,
			  0.0110,
			  0.0074,
			  0.0071,
			  0.0064,
			  0.0059,
			  0.0055,
			  0.0052,
			  0.0045};
    return 2.*PhiSigma[pbin];
  } else if (centbin==6) { 
    Double_t PhiSigma[9]={0.0202,
			  0.0113,
			  0.0077,
			  0.0071,
			  0.0069,
			  0.0064,
			  0.0060,
			  0.0055,
			  0.0050};
    return 2.*PhiSigma[pbin];
  } else if (centbin==7) { 
    Double_t PhiSigma[9]={0.0205,
			  0.0113,
			  0.0080,
			  0.0074,
			  0.0078,
			  0.0067,
			  0.0062,
			  0.0055,
			  0.0050};
    return 2.*PhiSigma[pbin];
  }

  return 0;
}

//________________________________________________________________________
void AliHadCorrTask::UserCreateOutputObjects()
{
  // Create my user objects.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  if (!fCreateHisto) return;

  TString name;
  TString temp;

  const Int_t nCentChBins = fNcentBins * 2;

  fHistMatchEtaPhiAll = new TH2F("fHistMatchEtaPhiAll", "fHistMatchEtaPhiAll;#Delta#eta;#Delta#phi", fNbins, -0.1, 0.1, fNbins, -0.1, 0.1);
  fOutput->Add(fHistMatchEtaPhiAll);

  fHistMatchEtaPhiAllTr = new TH2F("fHistMatchEtaPhiAllTr", "fHistMatchEtaPhiAllTr;#Delta#eta;#Delta#phi", fNbins, -0.1, 0.1, fNbins, -0.1, 0.1);
  fOutput->Add(fHistMatchEtaPhiAllTr);

  fHistMatchEtaPhiAllCl = new TH2F("fHistMatchEtaPhiAllCl", "fHistMatchEtaPhiAllCl;#Delta#eta;#Delta#phi", fNbins, -0.1, 0.1, fNbins, -0.1, 0.1);
  fOutput->Add(fHistMatchEtaPhiAllCl);

  for(Int_t icent=0; icent<nCentChBins; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {
        name = Form("fHistMatchEtaPhi_%i_%i_%i",icent,ipt,ieta);
        fHistMatchEtaPhi[icent][ipt][ieta] = new TH2F(name, name, fNbins, -0.1, 0.1, fNbins, -0.1, 0.1);
        fHistMatchEtaPhi[icent][ipt][ieta]->SetXTitle("#Delta#eta");
        fHistMatchEtaPhi[icent][ipt][ieta]->SetYTitle("#Delta#phi");
        fOutput->Add(fHistMatchEtaPhi[icent][ipt][ieta]);
      }
    }

    name = Form("fHistEsubPch_%i",icent);
    temp = Form("%s (Nmatches==1)",name.Data());
    fHistEsubPch[icent]=new TH1F(name, temp, fNbins, fMinBinPt, fMaxBinPt);
    fHistEsubPch[icent]->SetXTitle("#sum p (GeV) weighted with E_{sub}");
    fOutput->Add(fHistEsubPch[icent]);
    
    name = Form("fHistEsubPchRat_%i",icent);
    temp = Form("%s (Nmatches==1)",name.Data());
    fHistEsubPchRat[icent]=new TH2F(name, temp, fNbins, fMinBinPt, fMaxBinPt, fNbins*2, 0., 10.);
    fHistEsubPchRat[icent]->SetXTitle("#Sigma p (GeV)");
    fHistEsubPchRat[icent]->SetYTitle("E_{sub} / #sum p");
    fOutput->Add(fHistEsubPchRat[icent]);

    if (icent<fNcentBins) {
      for(Int_t itrk=0; itrk<4; ++itrk) {
        name = Form("fHistNCellsEnergy_%i_%i",icent,itrk);
        temp = Form("%s (Nmatches==%d);N_{cells};E_{clus} (GeV)",name.Data(),itrk);
        fHistNCellsEnergy[icent][itrk]  = new TH2F(name, temp, fNbins, fMinBinPt, fMaxBinPt, 101, -0.5, 100.5);
        fOutput->Add(fHistNCellsEnergy[icent][itrk]);
      }

      name = Form("fHistEsubPchRatAll_%i",icent);
      temp = Form("%s (all Nmatches)",name.Data());
      fHistEsubPchRatAll[icent]=new TH2F(name, temp, fNbins, fMinBinPt, fMaxBinPt, fNbins*2, 0., 10.);
      fHistEsubPchRatAll[icent]->SetXTitle("#Sigma p (GeV)");
      fHistEsubPchRatAll[icent]->SetYTitle("E_{sub} / #sum p");
      fOutput->Add(fHistEsubPchRatAll[icent]);

      name = Form("fHistMatchEvsP_%i",icent);
      temp = Form("%s (all Nmatches)",name.Data());
      fHistMatchEvsP[icent] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt, fNbins*2, 0., 10.);
      fHistMatchEvsP[icent]->SetXTitle("E_{clus} (GeV)");
      fHistMatchEvsP[icent]->SetYTitle("E_{clus} / #sum p");
      fOutput->Add(fHistMatchEvsP[icent]);

      name = Form("fHistMatchdRvsEP_%i",icent);
      temp = Form("%s (all Nmatches)",name.Data());
      fHistMatchdRvsEP[icent] = new TH2F(name, temp, fNbins, 0., 0.2, fNbins*2, 0., 10.);
      fHistMatchdRvsEP[icent]->SetXTitle("#Delta R between track and cluster");
      fHistMatchdRvsEP[icent]->SetYTitle("E_{clus} / p");
      fOutput->Add(fHistMatchdRvsEP[icent]);

      name = Form("fHistNMatchEnergy_%i",icent);
      fHistNMatchEnergy[icent] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt, 101, -0.5, 100.5);
      fHistNMatchEnergy[icent]->SetXTitle("E_{clus} (GeV)");
      fHistNMatchEnergy[icent]->SetYTitle("N_{matches}");
      fOutput->Add(fHistNMatchEnergy[icent]);
    }
  }

  fHistNclusvsCent      = new TH1F("Nclusvscent",      "NclusVsCent; Cent (%)",                     100, 0, 100);
  fHistNclusMatchvsCent = new TH1F("NclusMatchvscent", "NclusMatchVsCent (all Nmatches); Cent (%)", 100, 0, 100);
  fHistEbefore          = new TH1F("Ebefore",          "Ebefore; Cent (%); E_{clus} (GeV)",         100, 0, 100);
  fHistEafter           = new TH1F("Eafter",           "Eafter;  Cent (%); E_{clus} (GeV)",         100, 0, 100);
  fHistEoPCent          = new TH2F("EoPCent",          "EoPCent; Cent (%); E_{clus} / #sum p",      100, 0, 100, fNbins*2, 0, 10);
  fHistNMatchCent       = new TH2F("NMatchesCent",     "NMatchesCent; Cent (%); Nmatches",          100, 0, 100, 11, -0.5,  10.5);
  fHistNClusMatchCent   = new TH2F("NClusMatchesCent", "NClusMatchesCent; Cent (%); Nmatches",      100, 0, 100, 11, -0.5,  10.5);

  fOutput->Add(fHistNclusMatchvsCent);
  fOutput->Add(fHistNclusvsCent);
  fOutput->Add(fHistEbefore);
  fOutput->Add(fHistEafter);
  fOutput->Add(fHistEoPCent);
  fOutput->Add(fHistNMatchCent);
  fOutput->Add(fHistNClusMatchCent);

  if (fIsEmbedded) {
    for(Int_t icent=0; icent<fNcentBins; ++icent) {
      name = Form("fHistEmbTrackMatchesOversub_%d",icent);
      fHistEmbTrackMatchesOversub[icent] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt, fNbins, 0, 1.2);
      fHistEmbTrackMatchesOversub[icent]->GetXaxis()->SetTitle("E_{clus}^{raw} (GeV)");
      fHistEmbTrackMatchesOversub[icent]->GetYaxis()->SetTitle("E_{oversub} / E_{clus}^{raw}");
      fOutput->Add(fHistEmbTrackMatchesOversub[icent]);

      name = Form("fHistNonEmbTrackMatchesOversub_%d",icent);
      fHistNonEmbTrackMatchesOversub[icent] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt, fNbins, 0, 1.2);
      fHistNonEmbTrackMatchesOversub[icent]->GetXaxis()->SetTitle("E_{clus}^{raw} (GeV)");
      fHistNonEmbTrackMatchesOversub[icent]->GetYaxis()->SetTitle("E_{oversub} / E_{clus}^{raw}");
      fOutput->Add(fHistNonEmbTrackMatchesOversub[icent]);

      name = Form("fHistOversubMCClusters_%d",icent);
      fHistOversubMCClusters[icent] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt, fNbins, 0, 1.2);
      fHistOversubMCClusters[icent]->GetXaxis()->SetTitle("E_{clus}^{raw} (GeV)");
      fHistOversubMCClusters[icent]->GetYaxis()->SetTitle("E_{oversub} / E_{clus}^{raw}");
      fOutput->Add(fHistOversubMCClusters[icent]);

      name = Form("fHistOversubNonMCClusters_%d",icent);
      fHistOversubNonMCClusters[icent] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt, fNbins, 0, 1.2);
      fHistOversubNonMCClusters[icent]->GetXaxis()->SetTitle("E_{clus}^{raw} (GeV)");
      fHistOversubNonMCClusters[icent]->GetYaxis()->SetTitle("E_{oversub} / E_{clus}^{raw}");
      fOutput->Add(fHistOversubNonMCClusters[icent]);

      name = Form("fHistOversub_%d",icent);
      fHistOversub[icent] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt, fNbins, 0, 1.2);
      fHistOversub[icent]->GetXaxis()->SetTitle("E_{clus}^{raw} (GeV)");
      fHistOversub[icent]->GetYaxis()->SetTitle("E_{oversub} / E_{clus}^{raw}");
      fOutput->Add(fHistOversub[icent]);
    }
  }

  PostData(1, fOutput);
}

//________________________________________________________________________
void AliHadCorrTask::ExecOnce() 
{
  // Initialize the analysis.

  AliParticleContainer *tracks = GetParticleContainer(0);
  if (!tracks) {
    AliError(Form("%s: This task needs a particle container!", GetName()));
    return;
  }
  
  TClass trackClass(tracks->GetClassName());
  if (!trackClass.InheritsFrom("AliVTrack")) {
    tracks->SetClassName("AliVTrack"); // enforce only AliVTrack and derived classes
  }

  AliClusterContainer *clusters = GetClusterContainer(0);
  if (!clusters) {
    AliError(Form("%s: This task needs a cluster container!", GetName()));
    return;
  }

  TClass clusterClass(clusters->GetClassName());
  if (!clusterClass.InheritsFrom("AliVCluster")) {
    clusters->SetClassName("AliVCluster"); // enforce only AliVCluster and derived classes
  }

  // Do base class initializations and if it fails -> bail out
  AliAnalysisTaskEmcal::ExecOnce();
  if (!fLocalInitialized) return;

  if (dynamic_cast<AliAODEvent*>(InputEvent())) fEsdMode = kFALSE;

  if (!fOutCaloName.IsNull()) {
    if (fEsdMode) {
      fOutClusters = new TClonesArray("AliESDCaloCluster");
    }
    else {
      fOutClusters = new TClonesArray("AliAODCaloCluster");
    }

    fOutClusters->SetName(fOutCaloName);
    AddObjectToEvent(fOutClusters);
  }
}

//________________________________________________________________________
void AliHadCorrTask::DoMatchedTracksLoop(Int_t icluster,
                                         Double_t &totalTrkP, Int_t &Nmatches, Double_t &trkPMCfrac, Int_t &NMCmatches) 
{
  // Do the loop over matched tracks for the cluster.

  AliParticleContainer* tracks = GetParticleContainer(0);
  AliClusterContainer* clusters = GetClusterContainer(0);

  AliVCluster* cluster = clusters->GetCluster(icluster);

  if (!cluster) return;
  
  // loop over matched tracks
  Int_t Ntrks = cluster->GetNTracksMatched();
  for (Int_t i = 0; i < Ntrks; ++i) {
    AliVTrack* track = 0;

    if (fEsdMode) {
      Int_t itrack = cluster->GetTrackMatchedIndex(i);
      if (itrack >= 0) track = static_cast<AliVTrack*>(tracks->GetAcceptParticle(itrack));
    }
    else {
      track = static_cast<AliVTrack*>(cluster->GetTrackMatched(i));
      UInt_t rejectionReason = 0;
      if (!tracks->AcceptParticle(track, rejectionReason)) track = 0;
    }

    if (!track) continue;

    Double_t etadiff = 999;
    Double_t phidiff = 999;
    GetEtaPhiDiff(track, cluster, phidiff, etadiff);
    if (fCreateHisto) fHistMatchEtaPhiAllCl->Fill(etadiff, phidiff);

    // check if track also points to cluster
    if (fDoTrackClus && (track->GetEMCALcluster() != icluster)) continue;

    Double_t mom = track->P();
    UInt_t mombin = GetMomBin(mom); 
    Int_t centbinch = fCentBin;
    
    if (track->Charge() < 0) centbinch += fNcentBins;

    if (fCreateHisto) {
      Int_t etabin = 0;
      if (track->Eta() > 0) etabin=1;
      fHistMatchEtaPhi[centbinch][mombin][etabin]->Fill(etadiff, phidiff);
      fHistMatchEtaPhiAll->Fill(etadiff, phidiff);
    }
    
    Double_t etaCut   = 0.0;
    Double_t phiCutlo = 0.0;
    Double_t phiCuthi = 0.0;

    if (fPhiMatch > 0) {
      phiCutlo = -fPhiMatch;
      phiCuthi = +fPhiMatch;
    }
    else {
      phiCutlo = GetPhiMean(mombin, centbinch) - GetPhiSigma(mombin, fCentBin);
      phiCuthi = GetPhiMean(mombin, centbinch) + GetPhiSigma(mombin, fCentBin);
    }

    if (fEtaMatch > 0) {
      etaCut = fEtaMatch;
    }
    else {
      etaCut = GetEtaSigma(mombin);
    }

    if ((phidiff < phiCuthi && phidiff > phiCutlo) && TMath::Abs(etadiff) < etaCut) {
      if (track->GetLabel() > fMinMCLabel) {
	++NMCmatches;
	trkPMCfrac += mom;
      }
      ++Nmatches;
      totalTrkP += mom;

      if (fCreateHisto) {
        if (fHadCorr > 1) {
          Double_t dphi = 0;
          Double_t deta = 0;
          GetEtaPhiDiff(track, cluster, dphi, deta);
          Double_t dR = TMath::Sqrt(dphi*dphi + deta*deta);
          Double_t energyclus = cluster->GetNonLinCorrEnergy();
          fHistMatchdRvsEP[fCentBin]->Fill(dR, energyclus / mom);
        }
      }
    }
  }

  if (totalTrkP > 0) trkPMCfrac /= totalTrkP;
}

//________________________________________________________________________
Bool_t AliHadCorrTask::Run() 
{
  // Run the hadronic correction

  AliClusterContainer *clusters = GetClusterContainer(0);

  // delete output
  if (fOutClusters) fOutClusters->Delete();

  Int_t clusCount = 0;

   // loop over all clusters
  clusters->ResetCurrentID();
  AliVCluster *cluster = 0;
  while ((cluster = clusters->GetNextAcceptCluster())) {

    Double_t energyclus = 0;
    if (fCreateHisto) {
      fHistEbefore->Fill(fCent, cluster->GetNonLinCorrEnergy());
      fHistNclusvsCent->Fill(fCent);
    }
  
    // apply correction / subtraction
    // to subtract only the closest track set fHadCor to a %
    // to subtract all tracks within the cut set fHadCor to %+1
    if (fHadCorr > 1) {
      energyclus = ApplyHadCorrAllTracks(clusters->GetCurrentID(), fHadCorr - 1);
    }
    else if (fHadCorr > 0) {
      energyclus = ApplyHadCorrOneTrack(clusters->GetCurrentID(), fHadCorr);
    }
    else {
      energyclus = cluster->GetNonLinCorrEnergy();
    }

    if (energyclus < 0) energyclus = 0;

    cluster->SetHadCorrEnergy(energyclus);

    if (fCreateHisto) fHistEafter->Fill(fCent, energyclus);

    if (fOutClusters && energyclus > 0) { // create corrected cluster
      AliVCluster *oc;
      if (fEsdMode) {
	AliESDCaloCluster *ec = dynamic_cast<AliESDCaloCluster*>(cluster);
	if (!ec) continue;
	oc = new ((*fOutClusters)[clusCount]) AliESDCaloCluster(*ec);
      }
      else { 
	AliAODCaloCluster *ac = dynamic_cast<AliAODCaloCluster*>(cluster);
	if (!ac) continue;
	oc = new ((*fOutClusters)[clusCount]) AliAODCaloCluster(*ac);
      }
      ++clusCount;
      oc->SetE(energyclus);
      oc->SetNonLinCorrEnergy(cluster->GetNonLinCorrEnergy()); //get the non linearity corrected energy of the clusters. Works only if the clustermaker was ran before
      oc->SetHadCorrEnergy(energyclus); //same as the default energy field of this specific copy of the cluster container
    }
  }
  
  return kTRUE;
}

//________________________________________________________________________
Double_t AliHadCorrTask::ApplyHadCorrOneTrack(Int_t icluster, Double_t hadCorr) 
{
  // Apply the hadronic correction with one track only.

  AliParticleContainer *tracks = GetParticleContainer(0);
  AliClusterContainer *clusters = GetClusterContainer(0);

  AliVCluster* cluster = clusters->GetCluster(icluster);
  
  Double_t energyclus = cluster->GetNonLinCorrEnergy();
  
  AliVTrack* track = 0;

  if (cluster->GetNTracksMatched() > 0) {
    if (fEsdMode) {
      Int_t itrack = cluster->GetTrackMatchedIndex(0);
      if (itrack >= 0) track = static_cast<AliVTrack*>(tracks->GetAcceptParticle(itrack));
    }
    else {
      track = static_cast<AliVTrack*>(cluster->GetTrackMatched(0));
      UInt_t rejectionReason = 0;
      if (!tracks->AcceptParticle(track, rejectionReason)) track = 0;
    }
  }
  
  if (!track || track->P() < 1e-6) return energyclus;

  Double_t mom = track->P();

  Double_t dEtaMin = 1e9;
  Double_t dPhiMin = 1e9;
  GetEtaPhiDiff(track, cluster, dPhiMin, dEtaMin);
  
  if (fCreateHisto) fHistMatchEtaPhiAllCl->Fill(dEtaMin, dPhiMin);

  // check if track also points to cluster
  Int_t cid = track->GetEMCALcluster();
  if (fDoTrackClus && (cid != icluster)) return energyclus;

  UInt_t mombin = GetMomBin(mom);
  Int_t centbinch = fCentBin;
  if (track->Charge() < 0) centbinch += fNcentBins;

  // plot some histograms if switched on
  if (fCreateHisto) {
    Int_t etabin = 0;
    if(track->Eta() > 0) etabin = 1;
	    
    fHistMatchEtaPhi[centbinch][mombin][etabin]->Fill(dEtaMin, dPhiMin);
    fHistMatchEtaPhiAll->Fill(dEtaMin, dPhiMin);
    
    if (mom > 0) {
      Double_t etadiff = 0;
      Double_t phidiff = 0;
      GetEtaPhiDiff(track, cluster, phidiff, etadiff);
      Double_t dRmin = TMath::Sqrt(etadiff*etadiff + phidiff*phidiff);
      fHistMatchEvsP[fCentBin]->Fill(energyclus, energyclus / mom);
      fHistEoPCent->Fill(fCent, energyclus / mom);
      fHistMatchdRvsEP[fCentBin]->Fill(dRmin, energyclus / mom);
    }
  }
	   
  // define eta/phi cuts
  Double_t etaCut   = 0.0;
  Double_t phiCutlo = 0.0;
  Double_t phiCuthi = 0.0;
  if (fPhiMatch > 0) {
    phiCutlo = -fPhiMatch;
    phiCuthi = +fPhiMatch;
  }
  else {
    phiCutlo = GetPhiMean(mombin, centbinch) - GetPhiSigma(mombin, fCentBin);
    phiCuthi = GetPhiMean(mombin, centbinch) + GetPhiSigma(mombin, fCentBin);
  }
  if (fEtaMatch > 0) {
    etaCut = fEtaMatch;
  }
  else {
    etaCut = GetEtaSigma(mombin);
  }
  
  // apply the correction if the track is in the eta/phi window
  if ((dPhiMin < phiCuthi && dPhiMin > phiCutlo) && TMath::Abs(dEtaMin) < etaCut) {
    energyclus -= hadCorr * mom;
  }

  return energyclus;
}

//________________________________________________________________________
Double_t AliHadCorrTask::ApplyHadCorrAllTracks(Int_t icluster, Double_t hadCorr) 
{
  // Apply the hadronic correction with all tracks.

  AliClusterContainer *clusters = GetClusterContainer(0);
  AliParticleContainer *tracks = GetParticleContainer(0);
  
  AliVCluster* cluster = clusters->GetCluster(icluster);
  
  Double_t energyclus = cluster->GetNonLinCorrEnergy();
  Double_t cNcells = cluster->GetNCells();
  
  Double_t totalTrkP = 0.0;  // count total track momentum
  Int_t Nmatches = 0;        // count total number of matches

  Double_t trkPMCfrac = 0.0; // count total track momentum
  Int_t NMCmatches = 0;      // count total number of matches
  
  // do the loop over the matched tracks and get the number of matches and the total momentum
  DoMatchedTracksLoop(icluster, totalTrkP, Nmatches, trkPMCfrac, NMCmatches);

  Double_t Esub = hadCorr * totalTrkP;
	
  if (Esub > energyclus) Esub = energyclus;
	
  // applying Peter's proposed algorithm
  // never subtract the full energy of the cluster 
  Double_t clusEexcl = fEexclCell * cNcells;
  if (energyclus < clusEexcl) clusEexcl = energyclus;
  if ((energyclus - Esub) < clusEexcl) Esub = (energyclus - clusEexcl);

  // embedding
  Double_t EsubMC       = 0;
  Double_t EsubBkg      = 0;
  Double_t EclusMC      = 0;
  Double_t EclusBkg     = 0;
  Double_t EclusCorr    = 0;
  Double_t EclusMCcorr  = 0;
  Double_t EclusBkgcorr = 0;
  Double_t overSub = 0;
  if (fIsEmbedded) {
    EsubMC   = hadCorr * totalTrkP * trkPMCfrac;
    EsubBkg  = hadCorr * totalTrkP - EsubMC;
    EclusMC  = energyclus * cluster->GetMCEnergyFraction();
    EclusBkg = energyclus - EclusMC;
 
    if (energyclus > Esub)
      EclusCorr = energyclus - Esub;

    if (EclusMC > EsubMC)
      EclusMCcorr = EclusMC - EsubMC;
  
    if (EclusBkg > EsubBkg)
      EclusBkgcorr = EclusBkg - EsubBkg;
  
    overSub = EclusMCcorr + EclusBkgcorr - EclusCorr;
  }

  // plot some histograms if switched on
  if (fCreateHisto) {
    fHistNMatchCent->Fill(fCent, Nmatches);
    fHistNMatchEnergy[fCentBin]->Fill(energyclus, Nmatches);
    
    if (Nmatches > 0) fHistNclusMatchvsCent->Fill(fCent);

    if (Nmatches < 3) { 
      fHistNCellsEnergy[fCentBin][Nmatches]->Fill(energyclus, cNcells);
    }
    else {
      fHistNCellsEnergy[fCentBin][3]->Fill(energyclus, cNcells);
    }
      
    if (totalTrkP > 0) {
      Double_t EoP = energyclus / totalTrkP;
      fHistEoPCent->Fill(fCent, EoP);
      fHistMatchEvsP[fCentBin]->Fill(energyclus, EoP);

      fHistEsubPchRatAll[fCentBin]->Fill(totalTrkP, Esub / totalTrkP);

      if (Nmatches == 1) {
        AliVTrack* track = 0;
        if (fEsdMode) {
          Int_t itrack = cluster->GetTrackMatchedIndex(0);
          if (itrack >= 0) track = static_cast<AliVTrack*>(tracks->GetAcceptParticle(itrack));
        }
        else {
          track = static_cast<AliVTrack*>(cluster->GetTrackMatched(0));
          UInt_t rejectionReason = 0;
          if (!tracks->AcceptParticle(track, rejectionReason)) track = 0;
        }
        if (track) {
          Int_t centbinchm = fCentBin;
          if (track->Charge() < 0) centbinchm += fNcentBins;
          fHistEsubPchRat[centbinchm]->Fill(totalTrkP, Esub / totalTrkP);
          fHistEsubPch[centbinchm]->Fill(totalTrkP, Esub);
        }
      }

      if (fIsEmbedded) {
	fHistOversub[fCentBin]->Fill(energyclus, overSub / energyclus);

	if (cluster->GetMCEnergyFraction() > 0.95)
	  fHistOversubMCClusters[fCentBin]->Fill(energyclus, overSub / energyclus);
	else if (cluster->GetMCEnergyFraction() < 0.05)
	  fHistOversubNonMCClusters[fCentBin]->Fill(energyclus, overSub / energyclus);

	if (trkPMCfrac < 0.05)
	  fHistNonEmbTrackMatchesOversub[fCentBin]->Fill(energyclus, overSub / energyclus);
	else if (trkPMCfrac > 0.95)
	  fHistEmbTrackMatchesOversub[fCentBin]->Fill(energyclus, overSub / energyclus);
      }
    }
  }

  if (fIsEmbedded && fDoExact) {
    Esub -= overSub;
    if (EclusBkgcorr + EclusMCcorr > 0) {
      Double_t newfrac = EclusMCcorr / (EclusBkgcorr + EclusMCcorr);
      cluster->SetMCEnergyFraction(newfrac);
    }
  }

  // apply the correction
  energyclus -= Esub;

  return energyclus;
}
