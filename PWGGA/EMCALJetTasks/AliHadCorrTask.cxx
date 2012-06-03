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
  fOutClusters(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNMatchCent_trk(0),
  fHistCentrality(0),
  fHistNoMatchEtaPhi(0)
{
  // Default constructor.

  for(Int_t i=0; i<8; i++) {
      fHistEsubPch[i]    = 0;
      fHistEsubPchRat[i]    = 0;
    for(Int_t j=0; j<4; j++) {
      fHistNCellsEnergy[i][j] = 0;
    }
    if (i<4) {
      fHistMatchEvsP[i]    = 0;
      fHistMatchdRvsEP[i]  = 0;
      fHistNMatchEnergy[i] = 0;
    }
    for(Int_t j=0; j<9; j++) {
      for(Int_t k=0; k<2; k++) {
	fHistMatchEtaPhi[i][j][k] = 0;
      }
    }
  } 
}

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask(const char *name) : 
  AliAnalysisTaskEmcal(name, kFALSE),
  fOutCaloName("CaloClustersCorr"),
  fPhiMatch(0.05),
  fEtaMatch(0.025),
  fDoTrackClus(0),
  fHadCorr(0),
  fEexclCell(0),
  fOutClusters(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNMatchCent_trk(0),
  fHistCentrality(0),
  fHistNoMatchEtaPhi(0)

{
  // Standard constructor.

  for(Int_t i=0; i<8; i++) {
      fHistEsubPch[i] = 0;
      fHistEsubPchRat[i] = 0;
    for(Int_t j=0; j<3; j++) {
    }
    
    if (i<4) {
      fHistMatchEvsP[i]   = 0;
      fHistMatchdRvsEP[i] = 0;
    }
    
    for(Int_t j=0; j<9; j++) {
      for(Int_t k=0; k<2; k++) {
	fHistMatchEtaPhi[i][j][k] = 0;
      }
    } 
  }
  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask(const char *name, Bool_t histo) : 
  AliAnalysisTaskEmcal(name, histo),
  fOutCaloName("CaloClustersCorr"),
  fPhiMatch(0.05),
  fEtaMatch(0.025),
  fDoTrackClus(0),
  fHadCorr(0),
  fEexclCell(0),
  fOutClusters(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNMatchCent_trk(0),
  fHistCentrality(0),
  fHistNoMatchEtaPhi(0)

{
  // Standard constructor.

   for(Int_t i=0; i<8; i++) {
    if (i<4) {
      fHistMatchEvsP[i]   = 0;
      fHistMatchdRvsEP[i] = 0;
      for(Int_t j=0; j<3; j++) {
	fHistNCellsEnergy[i][j] = 0;
      }
    }
    fHistEsubPch[i] = 0;
    fHistEsubPchRat[i] = 0;
    for(Int_t j=0; j<9; j++) {
      for(Int_t k=0; k<2; k++) {
	fHistMatchEtaPhi[i][j][k] = 0;
      }
    }
  } 

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
  if (p>=0 && p<0.5) 
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
  Double_t EtaSigma[9]={0.0097,0.0075,0.0059,0.0055,0.0053,0.005,0.005,0.0045,0.0042};
  return 2.0*EtaSigma[pbin];
}

//________________________________________________________________________
Double_t AliHadCorrTask::GetPhiMean(Int_t pbin, Int_t centbin) const
{
  if (centbin==0){ 
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
  }else if(centbin==1){ 
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
  }else if(centbin==2){ 
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
  }else if(centbin==3){ 
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
  }else if(centbin==4){ 
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
  }else if(centbin==5){ 
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
  }else if(centbin==6){ 
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
  }else if(centbin==7){ 
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
  if (centbin==0){ 
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
  }else if(centbin==1){ 
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
  }else if(centbin==2){ 
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
  }else if(centbin==3){ 
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
  }else if(centbin==4){ 
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
  }else if(centbin==5){ 
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
  }else if(centbin==6){ 
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
  }else if(centbin==7){ 
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

  AliVEventHandler* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!handler) {
    AliError("Input handler not available!");
    return;
  }
 
  if (handler->InheritsFrom("AliESDInputHandler")) {
    fOutClusters = new TClonesArray("AliESDCaloCluster");
  } else if (handler->InheritsFrom("AliAODInputHandler")) {
    fOutClusters = new TClonesArray("AliAODCaloCluster");
  } else {
    AliError("Input handler not recognized!");
    return;
  }
  fOutClusters->SetName(fOutCaloName);

  if (!fCreateHisto)
    return;

  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

  TString name;

  for(Int_t icent=0; icent<8; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      for(Int_t ieta=0; ieta<2; ++ieta) {

	name = Form("fHistMatchEtaPhi_%i_%i_%i",icent,ipt,ieta);
	fHistMatchEtaPhi[icent][ipt][ieta] = new TH2F(name, name, 200, -0.1, 0.1, 400, -0.2, 0.2);
	fOutput->Add(fHistMatchEtaPhi[icent][ipt][ieta]);
      }
    }


    for(Int_t itrk=0; itrk<4; ++itrk) {
      
      name = Form("fHistNCellsEnergy_%i_%i",icent,itrk);
      fHistNCellsEnergy[icent][itrk]  = new TH2F(name, name, 1000, 0, 100, 101, -0.5, 100.5);
      fOutput->Add(fHistNCellsEnergy[icent][itrk]);
    }    
  
    name = Form("fHistEsubPch_%i",icent);
    fHistEsubPch[icent]=new TH1F(name, name, 400, 0., 100.);
    fOutput->Add(fHistEsubPch[icent]);
    
    name = Form("fHistEsubPchRat_%i",icent);
    fHistEsubPchRat[icent]=new TH2F(name, name, 400, 0., 200., 1000, 0., 10.);
    fOutput->Add(fHistEsubPchRat[icent]);
  
    
    if(icent<4){
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

  fHistCentrality       = new TH1F("fHistCentrality",  "Centrality",       100, 0, 100);

  fHistNclusvsCent      = new TH1F("Nclusvscent",      "NclusVsCent",      100, 0, 100);
  fHistNclusMatchvsCent = new TH1F("NclusMatchvscent", "NclusMatchVsCent", 100, 0, 100);
  fHistEbefore          = new TH1F("Ebefore",          "Ebefore",          100, 0, 100);
  fHistEafter           = new TH1F("Eafter",           "Eafter",           100, 0, 100);
  fHistEoPCent          = new TH2F("EoPCent",          "EoPCent",          100, 0, 100, 1000, 0,   10);
  fHistNMatchCent       = new TH2F("NMatchesCent",     "NMatchesCent",     100, 0, 100, 11, -0.5, 10.5);
  fHistNMatchCent_trk   = new TH2F("NMatchesCent_trk", "NMatchesCent_trk", 100, 0, 100, 11, -0.5, 10.5);

  fHistNoMatchEtaPhi   = new TH2F("NoMatchEtaPhi","NoMatchEtaPhi",200,-1.0,1.0,90,1.0,4.0);

  fOutput->Add(fHistNclusMatchvsCent);
  fOutput->Add(fHistNclusvsCent);
  fOutput->Add(fHistEbefore);
  fOutput->Add(fHistEafter);
  fOutput->Add(fHistEoPCent);
  fOutput->Add(fHistNMatchCent);
  fOutput->Add(fHistNMatchCent_trk);
  fOutput->Add(fHistCentrality);
  fOutput->Add(fHistNoMatchEtaPhi);

  PostData(1, fOutput);
}

//________________________________________________________________________
void AliHadCorrTask::DoTrackClusLoop() 
{
  const Int_t Ntrks = fTracks->GetEntries();

  // loop over all tracks
  for (Int_t t = 0; t < Ntrks; ++t) {
    AliEmcalParticle *emctrack = dynamic_cast<AliEmcalParticle*>(fTracks->At(t));
    if (!emctrack)
      continue;

    AliVTrack *track = emctrack->GetTrack();
    if (!track)
      continue;
    if (!AcceptTrack(track))
      continue;

    Int_t     Nclus    = emctrack->GetNumberOfMatchedObj();
    Int_t     Nmatches = 0;

    // loop over matched clusters
    for (Int_t i = 0; i < Nclus; ++i) {
      Int_t c = emctrack->GetMatchedObjId(i);
      AliEmcalParticle *emccluster = dynamic_cast<AliEmcalParticle*>(fCaloClusters->At(c));
      if (!emccluster)
	continue;

      AliVCluster *cluster = emccluster->GetCluster();
      if (!cluster)
	continue;

      Double_t etadiff = 999;
      Double_t phidiff = 999;
      AliPicoTrack::GetEtaPhiDiff(track, cluster, phidiff, etadiff);

      if (TMath::Abs(phidiff) < 0.050 && TMath::Abs(etadiff) < 0.025) // pp cuts!!!
	Nmatches++;
    }

    fHistNMatchCent_trk->Fill(fCent, Nmatches);

    if(Nmatches == 0 && track->Pt() > 2.0)
      fHistNoMatchEtaPhi->Fill(track->GetTrackEtaOnEMCal(), track->GetTrackPhiOnEMCal());
  }
}

//________________________________________________________________________
void AliHadCorrTask::DoMatchedTracksLoop(AliEmcalParticle *emccluster, Double_t &totalTrkP, Int_t &Nmatches) 
{
  // Do the loop over matched tracks for cluster emccluster

  AliVCluster *cluster = emccluster->GetCluster();
  Int_t iClus = emccluster->IdInCollection();
  Double_t energyclus = cluster->E();

  // loop over matched tracks
  Int_t Ntrks = emccluster->GetNumberOfMatchedObj();
  for (Int_t i = 0; i < Ntrks; ++i) {
    Int_t    iTrack = emccluster->GetMatchedObjId(i);
    Double_t dR     = emccluster->GetMatchedObjDistance(i);
    
    AliEmcalParticle *emctrack = dynamic_cast<AliEmcalParticle*>(fTracks->At(iTrack));
    if (!emctrack)
      continue;

    AliVTrack *track = emctrack->GetTrack();
    if (!track)
      continue;
    if (!AcceptTrack(track))
      continue;

    Double_t etadiff = 999;
    Double_t phidiff = 999;
    AliPicoTrack::GetEtaPhiDiff(track, cluster, phidiff, etadiff);
    
    Double_t mom       = track->P();
    Int_t    mombin    = GetMomBin(mom); 
    Int_t    centbinch = fCentBin;
    if (track->Charge() == -1 || track->Charge() == 255) 
      centbinch += 4;

    if (fCreateHisto) {
      Int_t etabin = 0;
      if(track->Eta() > 0) etabin=1;
	
      fHistMatchEtaPhi[centbinch][mombin][etabin]->Fill(etadiff, phidiff);
    }
    
    Double_t EtaCut   = 0.0;
    Double_t PhiCutlo = 0.0;
    Double_t PhiCuthi = 0.0;

    if (fPhiMatch > 0){
      PhiCutlo = -fPhiMatch;
      PhiCuthi = fPhiMatch;
    }
    else {
      PhiCutlo = GetPhiMean(mombin, centbinch) - GetPhiSigma(mombin, fCentBin);
      PhiCuthi = GetPhiMean(mombin, centbinch) + GetPhiSigma(mombin, fCentBin);
    }

    if (fEtaMatch > 0) {
      EtaCut = fEtaMatch;
    }
    else {
      EtaCut = GetEtaSigma(mombin);
    }

    if ((phidiff < PhiCuthi && phidiff > PhiCutlo) && TMath::Abs(etadiff) < EtaCut) {
      if ((fDoTrackClus && (track->GetEMCALcluster()) == iClus) || !fDoTrackClus) {
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
}

//________________________________________________________________________
Bool_t AliHadCorrTask::Run() 
{
  // Run the hadronic correction

  // post output in event if not yet present
  if (!(InputEvent()->FindListObject(fOutCaloName)))
    InputEvent()->AddObject(fOutClusters);
  
  // delete output
  fOutClusters->Delete();

  // esd or aod mode
  Bool_t esdMode = kTRUE;
  if (dynamic_cast<AliAODEvent*>(InputEvent()))
      esdMode = kFALSE;

  // optimization in case autobranch loading is off
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (fCaloName == "CaloClusters")
    am->LoadBranch("CaloClusters");
  if (fTracksName == "Tracks")
    am->LoadBranch("Tracks");
  am->LoadBranch("Centrality");      

  if (fCreateHisto)
    fHistCentrality->Fill(fCent);

  const Int_t Nclus = fCaloClusters->GetEntries();
 
  if (fDoTrackClus && fCreateHisto)
    DoTrackClusLoop();

  // loop over all clusters
  for (Int_t iClus = 0, clusCount=0; iClus < Nclus; ++iClus) {
    AliEmcalParticle *emccluster = dynamic_cast<AliEmcalParticle*>(fCaloClusters->At(iClus));
    if (!emccluster)
      continue;

    AliVCluster *cluster = emccluster->GetCluster();
    if (!cluster)
      continue;
    if (!AcceptCluster(cluster))
      continue;

    Double_t energyclus = 0;
    fHistEbefore->Fill(fCent, cluster->E());
  
    // apply correction / subtraction
    if (fHadCorr > 0) {
      // to subtract only the closest track set fHadCor to a %
      // to subtract all tracks within the cut set fHadCor to %+1
      if (fHadCorr > 1)
	energyclus = ApplyHadCorrAllTracks(emccluster, fHadCorr - 1);	
      else 
	energyclus = ApplyHadCorrOneTrack(emccluster, fHadCorr);	
    }

    if (energyclus < 0) 
      energyclus = 0;

    if (energyclus > 0) { // create corrected cluster
      AliVCluster *oc;
      if (esdMode) 
	oc = new ((*fOutClusters)[clusCount]) AliESDCaloCluster(*(dynamic_cast<AliESDCaloCluster*>(cluster)));
      else 
	oc = new ((*fOutClusters)[clusCount]) AliAODCaloCluster(*(dynamic_cast<AliAODCaloCluster*>(cluster)));
      
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
  AliVCluster *cluster = emccluster->GetCluster();
  
  Double_t energyclus = cluster->E();
  Int_t    iMin       = emccluster->GetMatchedObjId();

  if (iMin < 0)
    return energyclus;

  AliEmcalParticle *emctrack = dynamic_cast<AliEmcalParticle*>(fTracks->At(iMin));
  if (!emctrack)
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
  if (track->Charge() == -1 || track->Charge() == 255) 
    centbinch += 4;

  // Plot some histograms if switched on
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
	   
  // Define eta/phi cuts
  Double_t EtaCut   = 0.0;
  Double_t PhiCutlo = 0.0;
  Double_t PhiCuthi = 0.0;
  if (fPhiMatch > 0) {
    PhiCutlo = -fPhiMatch;
    PhiCuthi = fPhiMatch;
  }
  else {
    PhiCutlo = GetPhiMean(mombin, centbinch) - GetPhiSigma(mombin, fCentBin);
    PhiCuthi = GetPhiMean(mombin, centbinch) + GetPhiSigma(mombin, fCentBin);
  }
  if(fEtaMatch > 0) {
    EtaCut = fEtaMatch;
  }
  else {
    EtaCut = GetEtaSigma(mombin);
  }
  
  // Apply the correction if the track is in the eta/phi window
  if ((dPhiMin < PhiCuthi && dPhiMin > PhiCutlo) && TMath::Abs(dEtaMin) < EtaCut) {

    if ((fDoTrackClus && (track->GetEMCALcluster()) == emccluster->IdInCollection()) || !fDoTrackClus){
      energyclus -= hadCorr * mom;
    }
  }

  return energyclus;
}

//________________________________________________________________________
Double_t AliHadCorrTask::ApplyHadCorrAllTracks(AliEmcalParticle *emccluster, Double_t hadCorr) 
{
  AliVCluster *cluster = emccluster->GetCluster();
  
  Double_t energyclus = cluster->E();
  Double_t cNcells = cluster->GetNCells();
  
  Double_t totalTrkP  = 0.0; // count total track momentum
  Int_t    Nmatches   = 0;   // count total number of matches
  
  // Do the loop over the matched tracks and get the number of matches and the total momentum
  DoMatchedTracksLoop(emccluster, totalTrkP, Nmatches);

  Double_t Esub = hadCorr * totalTrkP;

  Double_t EoP = -1;
  if (totalTrkP > 0)
    EoP = energyclus / totalTrkP;

  // Plot some histograms if switched on
  if(fCreateHisto) {
    fHistNclusvsCent->Fill(fCent);
    fHistNMatchCent->Fill(fCent, Nmatches);
    fHistNMatchEnergy[fCentBin]->Fill(energyclus, Nmatches);
    
    if(Nmatches > 0) 
      fHistNclusMatchvsCent->Fill(fCent);
      
    if(Nmatches == 0)
      fHistNCellsEnergy[fCentBin][0]->Fill(energyclus, cNcells);
    else if(Nmatches == 1)
      fHistNCellsEnergy[fCentBin][1]->Fill(energyclus, cNcells);	
    else if(Nmatches == 2)
      fHistNCellsEnergy[fCentBin][2]->Fill(energyclus, cNcells);
    else
      fHistNCellsEnergy[fCentBin][3]->Fill(energyclus, cNcells);

    if (EoP > 0) {
      fHistEoPCent->Fill(fCent, EoP);
      fHistMatchEvsP[fCentBin]->Fill(energyclus, EoP);
    }

    if (Nmatches == 1) {
      Int_t iMin = emccluster->GetMatchedObjId();
      AliEmcalParticle *emctrack = dynamic_cast<AliEmcalParticle*>(fTracks->At(iMin));
      AliVTrack *track = emctrack->GetTrack();
      Int_t centbinchm = fCentBin;
      if (track->Charge() == -1 || track->Charge() == 255) 
	centbinchm += 4;
      
      if (totalTrkP > 0)
	fHistEsubPchRat[centbinchm]->Fill(totalTrkP, Esub / totalTrkP);

      fHistEsubPch[centbinchm]->Fill(totalTrkP, Esub);
    } 
  }

  if (totalTrkP <= 0)
    return energyclus;
 
  Double_t clusEexcl = fEexclCell * cNcells;

  if (energyclus < clusEexcl) 
    clusEexcl = energyclus;
	
  if (Esub > energyclus) 
    Esub = energyclus;
	
  //applying Peter's proposed algorithm
  //Never subtract the full energy of the cluster 
  if ((energyclus - Esub) < clusEexcl) 
    Esub = (energyclus - clusEexcl);

  //apply the correction
  energyclus -= Esub;

  return energyclus;
}

void AliHadCorrTask::Terminate(Option_t *)
{
  // Nothing to be done.
}
