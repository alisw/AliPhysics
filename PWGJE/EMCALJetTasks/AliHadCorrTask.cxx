// $Id$
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
#include "AliPicoTrack.h"
#include "AliVEventHandler.h"
#include "AliEmcalParticle.h"
#include "AliEMCALGeometry.h"

#include "AliHadCorrTask.h"

ClassImp(AliHadCorrTask)

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask() : 
  AliAnalysisTaskEmcal("AliHadCorrTask", kFALSE),
  fOutCaloName(),
  fPhiMatch(0.05),
  fEtaMatch(0.025),
  fDoTrackClus(0),
  fHadCorr(0),
  fEexclCell(0),
  fEsdMode(kTRUE),
  fOutClusters(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNClusMatchCent(0)
{
  // Default constructor.

  for(Int_t i=0; i<8; i++) {
    fHistEsubPch[i]    = 0;
    fHistEsubPchRat[i] = 0;
    fHistEsubPchRatAll[i] = 0;

    if (i<4) {
      fHistMatchEvsP[i]    = 0;
      fHistMatchdRvsEP[i]  = 0;
      fHistNMatchEnergy[i] = 0;

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
  fDoTrackClus(1),
  fHadCorr(0),
  fEexclCell(0),
  fEsdMode(kTRUE),
  fOutClusters(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNClusMatchCent(0)
{
  // Standard constructor.

  for(Int_t i=0; i<8; i++) {
    fHistEsubPch[i]    = 0;
    fHistEsubPchRat[i] = 0;
    fHistEsubPchRatAll[i] = 0;

    if (i<4) {
      fHistMatchEvsP[i]    = 0;
      fHistMatchdRvsEP[i]  = 0;
      fHistNMatchEnergy[i] = 0;

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
Int_t AliHadCorrTask::GetMomBin(Double_t p) const
{
  // Get momenum bin.

  Int_t pbin=-1;
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
  else if (p>=8.) 
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

  if (!fCreateHisto)
    return;

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  TString name;

  for(Int_t icent=0; icent<8; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {
	name = Form("fHistMatchEtaPhi_%i_%i_%i",icent,ipt,ieta);
	fHistMatchEtaPhi[icent][ipt][ieta] = new TH2F(name, name, 200, -0.1, 0.1, 400, -0.2, 0.2);
	fOutput->Add(fHistMatchEtaPhi[icent][ipt][ieta]);
      }
    }

    name = Form("fHistEsubPch_%i",icent);
    fHistEsubPch[icent]=new TH1F(name, name, 400, 0., 100.);
    fOutput->Add(fHistEsubPch[icent]);
    
    name = Form("fHistEsubPchRat_%i",icent);
    fHistEsubPchRat[icent]=new TH2F(name, name, 400, 0., 200., 1000, 0., 10.);
    fOutput->Add(fHistEsubPchRat[icent]);

    name = Form("fHistEsubPchRatAll_%i",icent);
    fHistEsubPchRatAll[icent]=new TH2F(name, name, 400, 0., 200., 1000, 0., 10.);
    fOutput->Add(fHistEsubPchRatAll[icent]);
    
    if (icent<4) {
      for(Int_t itrk=0; itrk<4; ++itrk) {
	name = Form("fHistNCellsEnergy_%i_%i",icent,itrk);
	fHistNCellsEnergy[icent][itrk]  = new TH2F(name, name, 1000, 0, 100, 101, -0.5, 100.5);
	fOutput->Add(fHistNCellsEnergy[icent][itrk]);
      }

      name = Form("fHistMatchEvsP_%i",icent);
      fHistMatchEvsP[icent] = new TH2F(name, name, 400, 0., 200., 1000, 0., 10.);
      fOutput->Add(fHistMatchEvsP[icent]);

      name = Form("fHistMatchdRvsEP_%i",icent);
      fHistMatchdRvsEP[icent] = new TH2F(name, name, 1000, 0., 1., 1000, 0., 10.);
      fOutput->Add(fHistMatchdRvsEP[icent]);

      name = Form("fHistNMatchEnergy_%i",icent);
      fHistNMatchEnergy[icent]  = new TH2F(name, name, 1000, 0, 100, 101, -0.5, 100.5);
      fOutput->Add(fHistNMatchEnergy[icent]);
    }
  }

  fHistNclusvsCent      = new TH1F("Nclusvscent",      "NclusVsCent",      100, 0, 100);
  fHistNclusMatchvsCent = new TH1F("NclusMatchvscent", "NclusMatchVsCent", 100, 0, 100);
  fHistEbefore          = new TH1F("Ebefore",          "Ebefore",          100, 0, 100);
  fHistEafter           = new TH1F("Eafter",           "Eafter",           100, 0, 100);
  fHistEoPCent          = new TH2F("EoPCent",          "EoPCent",          100, 0, 100, 1000, 0,   10);
  fHistNMatchCent       = new TH2F("NMatchesCent",     "NMatchesCent",     100, 0, 100, 11, -0.5, 10.5);
  fHistNClusMatchCent   = new TH2F("NClusMatchesCent", "NClusMatchesCent", 100, 0, 100, 11, -0.5, 10.5);

  fOutput->Add(fHistNclusMatchvsCent);
  fOutput->Add(fHistNclusvsCent);
  fOutput->Add(fHistEbefore);
  fOutput->Add(fHistEafter);
  fOutput->Add(fHistEoPCent);
  fOutput->Add(fHistNMatchCent);
  fOutput->Add(fHistNClusMatchCent);

  PostData(1, fOutput);
}

//________________________________________________________________________
void AliHadCorrTask::ExecOnce() 
{
  // Do base class initializations and if it fails -> bail out
  if (!fInitialized)
    AliAnalysisTaskEmcal::ExecOnce();
  if (!fInitialized)
    return;

  // Init the analysis.

  AliAnalysisTaskEmcal::ExecOnce();

  if (!fInitialized)
    return;

  if (dynamic_cast<AliAODEvent*>(InputEvent()))
    fEsdMode = kFALSE;

  if (fEsdMode) { // optimization in case autobranch loading is off
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    if (fCaloName == "CaloClusters")
      am->LoadBranch("CaloClusters");
    if (fTracksName == "Tracks")
      am->LoadBranch("Tracks");
    am->LoadBranch("Centrality.");      
  }

  if (fEsdMode) 
    fOutClusters = new TClonesArray("AliESDCaloCluster");
  else 
    fOutClusters = new TClonesArray("AliAODCaloCluster");

  fOutClusters->SetName(fOutCaloName);

  // post output in event if not yet present
  if (!(InputEvent()->FindListObject(fOutCaloName))) {
    InputEvent()->AddObject(fOutClusters);
  }
  else {
    fInitialized = kFALSE;
    AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutCaloName.Data()));
    fInitialized = kFALSE;
    return;
  }
}

//________________________________________________________________________
void AliHadCorrTask::DoTrackLoop() 
{
  const Int_t Ntracks = fTracks->GetEntries();
  for (Int_t iTrk = 0; iTrk < Ntracks; ++iTrk) {
    Int_t NmatchClus = 0;

    AliEmcalParticle *emctrack = static_cast<AliEmcalParticle*>(fTracks->At(iTrk));
    if (!emctrack)
      continue;
    if (!emctrack->IsEMCAL())
      continue;

    AliVTrack *track = emctrack->GetTrack();
    if (!track)
      continue;
    if (!AcceptTrack(track))
      continue;
    
    if (track->GetTrackEtaOnEMCal() < fGeom->GetArm1EtaMin() + fEtaMatch ||
	track->GetTrackEtaOnEMCal() > fGeom->GetArm1EtaMax() - fEtaMatch || 
	track->GetTrackPhiOnEMCal() < fGeom->GetArm1PhiMin() * TMath::DegToRad() + fPhiMatch || 
	track->GetTrackPhiOnEMCal() > fGeom->GetArm1PhiMax() * TMath::DegToRad() - fPhiMatch)
      continue;
    
    Int_t Nclus = emctrack->GetNumberOfMatchedObj();

    for (Int_t iClus = 0; iClus < Nclus; ++iClus) {
      AliEmcalParticle *emccluster = static_cast<AliEmcalParticle*>(fCaloClusters->At(iClus));
      if (!emccluster)
	continue;

      AliVCluster *cluster = emccluster->GetCluster();
      if (!cluster)
	continue;
      if (!AcceptCluster(cluster))
	continue;

      Double_t etadiff = 999;
      Double_t phidiff = 999;
      AliPicoTrack::GetEtaPhiDiff(track, cluster, phidiff, etadiff);

      if (TMath::Abs(phidiff) < fPhiMatch && TMath::Abs(etadiff) < fEtaMatch) 
	NmatchClus++;
    }
    fHistNClusMatchCent->Fill(fCent, NmatchClus);
  }
}

//________________________________________________________________________
void AliHadCorrTask::DoMatchedTracksLoop(AliEmcalParticle *emccluster, Double_t &totalTrkP, Int_t &Nmatches) 
{
  // Do the loop over matched tracks for cluster emccluster.

  AliVCluster *cluster = emccluster->GetCluster();
  Int_t iClus = emccluster->IdInCollection();
  Double_t energyclus = cluster->E();

  // loop over matched tracks
  const Int_t Ntrks = emccluster->GetNumberOfMatchedObj();
  for (Int_t i = 0; i < Ntrks; ++i) {
    Int_t    iTrack = emccluster->GetMatchedObjId(i);
    Double_t dR     = emccluster->GetMatchedObjDistance(i);
    
    AliEmcalParticle *emctrack = static_cast<AliEmcalParticle*>(fTracks->At(iTrack));
    if (!emctrack)
      continue;
    AliVTrack *track = emctrack->GetTrack();
    if (!track)
      continue;
    if (!AcceptTrack(track))
      continue;

    // check if track also points to cluster
    if (fDoTrackClus && (track->GetEMCALcluster()) != iClus)
      continue;

    Double_t etadiff = 999;
    Double_t phidiff = 999;
    AliPicoTrack::GetEtaPhiDiff(track, cluster, phidiff, etadiff);
    
    Double_t mom       = track->P();
    Int_t    mombin    = GetMomBin(mom); 
    Int_t    centbinch = fCentBin;
    if (track->Charge()<0) 
      centbinch += 4;

    if (fCreateHisto) {
      Int_t etabin = 0;
      if(track->Eta() > 0) etabin=1;
      fHistMatchEtaPhi[centbinch][mombin][etabin]->Fill(etadiff, phidiff);
    }
    
    Double_t etaCut   = 0.0;
    Double_t phiCutlo = 0.0;
    Double_t phiCuthi = 0.0;

    if (fPhiMatch > 0) {
      phiCutlo = -fPhiMatch;
      phiCuthi = +fPhiMatch;
    } else {
      phiCutlo = GetPhiMean(mombin, centbinch) - GetPhiSigma(mombin, fCentBin);
      phiCuthi = GetPhiMean(mombin, centbinch) + GetPhiSigma(mombin, fCentBin);
    }

    if (fEtaMatch > 0) {
      etaCut = fEtaMatch;
    } else {
      etaCut = GetEtaSigma(mombin);
    }

    if ((phidiff < phiCuthi && phidiff > phiCutlo) && TMath::Abs(etadiff) < etaCut) {
      ++Nmatches;
      totalTrkP += track->P();

      if (fCreateHisto) {
        if (fHadCorr > 1 && mombin > -1) {
          fHistMatchdRvsEP[fCentBin]->Fill(dR, energyclus / mom);
        }
      }
    }
  }
}

//________________________________________________________________________
Bool_t AliHadCorrTask::Run() 
{
  // Run the hadronic correction
  
  if (fCreateHisto)
    DoTrackLoop();

  // delete output
  fOutClusters->Delete();

   // loop over all clusters
  const Int_t Nclus = fCaloClusters->GetEntries();
  for (Int_t iClus = 0, clusCount=0; iClus < Nclus; ++iClus) {
    AliEmcalParticle *emccluster = static_cast<AliEmcalParticle*>(fCaloClusters->At(iClus));
    if (!emccluster)
      continue;

    AliVCluster *cluster = emccluster->GetCluster();
    if (!cluster)
      continue;
    if (!AcceptCluster(cluster))
      continue;

    Double_t energyclus = 0;
    if (fCreateHisto)
      fHistEbefore->Fill(fCent, cluster->E());
  
    // apply correction / subtraction
    // to subtract only the closest track set fHadCor to a %
    // to subtract all tracks within the cut set fHadCor to %+1
    if (fHadCorr > 1)
      energyclus = ApplyHadCorrAllTracks(emccluster, fHadCorr - 1);	
    else if (fHadCorr > 0)
      energyclus = ApplyHadCorrOneTrack(emccluster, fHadCorr);	
    else 
      energyclus = cluster->E();

    if (energyclus < 0) 
      energyclus = 0;

    if (energyclus > 0) { // create corrected cluster
      AliVCluster *oc;
      if (fEsdMode) {
        AliESDCaloCluster *ec = dynamic_cast<AliESDCaloCluster*>(cluster);
        if (!ec)
          continue;
	oc = new ((*fOutClusters)[clusCount]) AliESDCaloCluster(*ec);
      } else { 
        AliAODCaloCluster *ac = dynamic_cast<AliAODCaloCluster*>(cluster);
        if (!ac)
          continue;
	oc = new ((*fOutClusters)[clusCount]) AliAODCaloCluster(*ac);
      }
      oc->SetE(energyclus);

      ++clusCount;

      if (fCreateHisto)
	fHistEafter->Fill(fCent, energyclus);
    }
  }
  
  return kTRUE;
}

//________________________________________________________________________
Double_t AliHadCorrTask::ApplyHadCorrOneTrack(AliEmcalParticle *emccluster, Double_t hadCorr) 
{
  // Apply the hadronic correction with one track only.

  AliVCluster *cluster = emccluster->GetCluster();
  Double_t energyclus  = cluster->E();
  Int_t    iMin        = emccluster->GetMatchedObjId();
  if (iMin < 0)
    return energyclus;

  AliEmcalParticle *emctrack = static_cast<AliEmcalParticle*>(fTracks->At(iMin));
  if (!emctrack)
    return energyclus;

  // check if track also points to cluster
  Int_t cid = emctrack->GetMatchedObjId();
  if (fDoTrackClus && (cid!=emccluster->IdInCollection())) 
    return energyclus;

  AliVTrack *track = emctrack->GetTrack();
  if (!track)
    return energyclus;

  Double_t mom = track->P();
  if (mom == 0)
    return energyclus;

  Double_t dRmin      = emccluster->GetMatchedObjDistance();
  Double_t dEtaMin    = 1e9;
  Double_t dPhiMin    = 1e9;
  AliPicoTrack::GetEtaPhiDiff(track, cluster, dPhiMin, dEtaMin);

  Int_t mombin = GetMomBin(mom);
  Int_t centbinch = fCentBin;
  if (track->Charge()<0) 
    centbinch += 4;

  // plot some histograms if switched on
  if (fCreateHisto) {
    Int_t etabin = 0;
    if(track->Eta() > 0) 
      etabin = 1;
	    
    fHistMatchEtaPhi[centbinch][mombin][etabin]->Fill(dEtaMin, dPhiMin);
    
    if (mom > 0) {
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
  if(fEtaMatch > 0) {
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
Double_t AliHadCorrTask::ApplyHadCorrAllTracks(AliEmcalParticle *emccluster, Double_t hadCorr) 
{
  // Apply the hadronic correction with all tracks.

  AliVCluster *cluster = emccluster->GetCluster();
  
  Double_t energyclus = cluster->E();
  Double_t cNcells = cluster->GetNCells();
  
  Double_t totalTrkP  = 0.0; // count total track momentum
  Int_t    Nmatches   = 0;   // count total number of matches
  
  // do the loop over the matched tracks and get the number of matches and the total momentum
  DoMatchedTracksLoop(emccluster, totalTrkP, Nmatches);

  Double_t Esub = hadCorr * totalTrkP;
	
  if (Esub > energyclus) 
    Esub = energyclus;
	
  // applying Peter's proposed algorithm
  // never subtract the full energy of the cluster 
  Double_t clusEexcl = fEexclCell * cNcells;
  if (energyclus < clusEexcl) 
    clusEexcl = energyclus;
  if ((energyclus - Esub) < clusEexcl) 
    Esub = (energyclus - clusEexcl);

  Double_t EoP = 0;
  if (totalTrkP>0)
    EoP = energyclus / totalTrkP;

  // plot some histograms if switched on
  if (fCreateHisto) {
    fHistNclusvsCent->Fill(fCent);
    fHistNMatchCent->Fill(fCent, Nmatches);
    fHistNMatchEnergy[fCentBin]->Fill(energyclus, Nmatches);
    
    if(Nmatches > 0) 
      fHistNclusMatchvsCent->Fill(fCent);

    if(Nmatches < 3)  
      fHistNCellsEnergy[fCentBin][Nmatches]->Fill(energyclus, cNcells);
    else
      fHistNCellsEnergy[fCentBin][3]->Fill(energyclus, cNcells);

    if (EoP > 0) {
      fHistEoPCent->Fill(fCent, EoP);
      fHistMatchEvsP[fCentBin]->Fill(energyclus, EoP);
    }

    if (Nmatches == 1 && totalTrkP > 0) {
      Int_t iMin = emccluster->GetMatchedObjId();
      AliEmcalParticle *emctrack = static_cast<AliEmcalParticle*>(fTracks->At(iMin));
      AliVTrack *track = emctrack->GetTrack();
      Int_t centbinchm = fCentBin;
      if (track->Charge()<0) 
	centbinchm += 4;
      
      fHistEsubPchRat[centbinchm]->Fill(totalTrkP, Esub / totalTrkP);
      fHistEsubPch[centbinchm]->Fill(totalTrkP, Esub);
    }
    if (totalTrkP > 0) {
      fHistEsubPchRatAll[fCentBin]->Fill(totalTrkP, Esub / totalTrkP);
    }
  }

  // apply the correction
  energyclus -= Esub;

  return energyclus;
}

