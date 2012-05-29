// $Id$
//
// Hadronic correction task.
//
// Author: R.Reed, C.Loizides

#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAODCaloCluster.h"
#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDCaloCluster.h"
#include "AliESDtrack.h"
#include "AliPicoTrack.h"
#include "AliVEventHandler.h"

#include "AliHadCorrTask.h"

ClassImp(AliHadCorrTask)

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask() : 
  AliAnalysisTaskSE("AliHadCorrTask"),
  fTracksName(),
  fCaloName(),
  fOutCaloName(),
  fPhiMatch(0.05),
  fEtaMatch(0.025),
  fDoTrackClus(0),
  fHadCorr(0),
  fMinPt(0.15),
  fCreateHisto(kFALSE),
  fOutClusters(0),
  fOutputList(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNMatchCent_trk(0),
  fHistCentrality(0)
{
  // Default constructor.

  for(Int_t i=0; i<8; i++) {
    if (i<4) {
      fHistMatchEvsP[i]   = 0;
      fHistMatchdRvsEP[i] = 0;
      for(Int_t j=0; j<3; j++) {
	fHistEsubPch[i][j]    = 0;
	fHistEsubPchRat[i][j]    = 0;
      }
    }
    for(Int_t j=0; j<9; j++) {
      fHistMatchEtaPhi[i][j] = 0;
    }

  } 
}

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask(const char *name) : 
  AliAnalysisTaskSE(name),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fOutCaloName("CaloClustersCorr"),
  fPhiMatch(0.05),
  fEtaMatch(0.025),
  fDoTrackClus(0),
  fHadCorr(0),
  fMinPt(0.15),
  fCreateHisto(kFALSE),
  fOutClusters(0),
  fOutputList(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNMatchCent_trk(0),
  fHistCentrality(0)

{
  // Standard constructor.

   for(Int_t i=0; i<8; i++) {
    if (i<4) {
      fHistMatchEvsP[i]   = 0;
      fHistMatchdRvsEP[i] = 0;
      for(Int_t j=0; j<3; j++) {
	fHistEsubPch[i][j] = 0;
	fHistEsubPchRat[i][j] = 0;
      }
    }
    for(Int_t j=0; j<9; j++) {
      fHistMatchEtaPhi[i][j] = 0;
    }
  } 

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";
}

//________________________________________________________________________
AliHadCorrTask::AliHadCorrTask(const char *name, Bool_t histo) : 
  AliAnalysisTaskSE(name),
  fTracksName("Tracks"),
  fCaloName("CaloClusters"),
  fOutCaloName("CaloClustersCorr"),
  fPhiMatch(0.05),
  fEtaMatch(0.025),
  fDoTrackClus(0),
  fHadCorr(0),
  fMinPt(0.15),
  fCreateHisto(histo),
  fOutClusters(0),
  fOutputList(0),
  fHistNclusvsCent(0),
  fHistNclusMatchvsCent(0),
  fHistEbefore(0),
  fHistEafter(0),
  fHistEoPCent(0),
  fHistNMatchCent(0),
  fHistNMatchCent_trk(0),
  fHistCentrality(0)

{
  // Standard constructor.

   for(Int_t i=0; i<8; i++) {
    if (i<4) {
      fHistMatchEvsP[i]   = 0;
      fHistMatchdRvsEP[i] = 0;
      for(Int_t j=0; j<3; j++) {
	fHistEsubPch[i][j] = 0;
	fHistEsubPchRat[i][j] = 0;
      }
    }
    for(Int_t j=0; j<9; j++) {
      fHistMatchEtaPhi[i][j] = 0;
    }
  } 

  fBranchNames="ESD:AliESDRun.,AliESDHeader.,PrimaryVertex.";

  if (fCreateHisto)
    DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliHadCorrTask::~AliHadCorrTask()
{
  // Destructor
}

//________________________________________________________________________
Int_t AliHadCorrTask::GetCentBin(Double_t cent) const
{
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10) 
    centbin=0;
  else if (cent>=10 && cent<30) 
    centbin=1;
  else if (cent>=30 && cent<50) 
    centbin=2;
  else if (cent>=50 && cent<=100) 
    centbin=3;

  return centbin;
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

  Double_t EtaSigma[9]={0.0097,0.0075,0.0059,0.0055,0.0053,0.005,0.005,0.045,0.042};
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
  fOutputList = new TList();
  fOutputList->SetOwner();

  if (!fCreateHisto) {
    PostData(1, fOutputList);
    return;
  }

  for(Int_t icent=0; icent<8; ++icent) {
    for(Int_t ipt=0; ipt<9; ++ipt) {
      TString name(Form("fHistMatchEtaPhi_%i_%i",icent,ipt));
      fHistMatchEtaPhi[icent][ipt] = new TH2F(name, name, 400, -0.2, 0.2, 1600, -0.8, 0.8);
      fOutputList->Add(fHistMatchEtaPhi[icent][ipt]);
    }

    if(icent<4){
      TString name(Form("fHistMatchEvsP_%i",icent));
      fHistMatchEvsP[icent] = new TH2F(name, name, 400, 0., 200., 1000, 0., 10.);
      fOutputList->Add(fHistMatchEvsP[icent]);

      name = Form("fHistMatchdRvsEP_%i",icent);
      fHistMatchdRvsEP[icent] = new TH2F(name, name, 1000, 0., 1., 1000, 0., 10.);
      fOutputList->Add(fHistMatchdRvsEP[icent]);

      for(Int_t itrk=0; itrk<3; ++itrk){
	name = Form("fHistEsubPch_%i_%i",icent,itrk);
	fHistEsubPch[icent][itrk]=new TH1F(name, name, 400, 0., 100.);
	fHistEsubPch[icent][itrk]->Sumw2();
	fOutputList->Add(fHistEsubPch[icent][itrk]);

	name = Form("fHistEsubPchRat_%i_%i",icent,itrk);
	fHistEsubPchRat[icent][itrk]=new TH2F(name, name, 400, 0., 200., 1000, 0., 10.);
	fOutputList->Add(fHistEsubPchRat[icent][itrk]);
      }
    }
  }

  fHistCentrality       = new TH1F("fHistCentrality",  "Centrality",       100, 0, 100);

  fHistNclusvsCent      = new TH1F("Nclusvscent",      "NclusVsCent",      100, 0, 100);
  fHistNclusMatchvsCent = new TH1F("NclusMatchvscent", "NclusMatchVsCent", 100, 0, 100);
  fHistEbefore          = new TH1F("Ebefore",          "Ebefore",          100, 0, 100);
  fHistEafter           = new TH1F("Eafter",           "Eafter",           100, 0, 100);
  fHistEoPCent          = new TH2F("EoPCent",          "EoPCent",          100, 0, 100, 1000, 0,   10);
  fHistNMatchCent       = new TH2F("NMatchesCent",     "NMatchesCent",     100, 0, 100, 101, -0.5, 100.5);
  fHistNMatchCent_trk   = new TH2F("NMatchesCent_trk", "NMatchesCent_trk", 100, 0, 100, 101, -0.5, 100.5);
  fOutputList->Add(fHistNclusMatchvsCent);
  fOutputList->Add(fHistNclusvsCent);
  fOutputList->Add(fHistEbefore);
  fOutputList->Add(fHistEafter);
  fOutputList->Add(fHistEoPCent);
  fOutputList->Add(fHistNMatchCent);
  fOutputList->Add(fHistNMatchCent_trk);
  fOutputList->Add(fHistCentrality);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliHadCorrTask::UserExec(Option_t *) 
{
  // Execute per event.

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
  
  TList *l = InputEvent()->GetList();
  
  // get centrality 
  Double_t cent = -1; 
 
  AliCentrality *centrality = InputEvent()->GetCentrality() ;

  if (centrality)
    cent = centrality->GetCentralityPercentile("V0M");
  else
    cent=99; // probably pp data
  
  if (cent<0) {
    AliWarning(Form("Centrality negative: %f, assuming 99", cent));
    cent = 99;
  }
  
  Int_t centbin = GetCentBin(cent);

  if (fCreateHisto)
    fHistCentrality->Fill(cent);

  // get input collections
  TClonesArray *tracks = 0;
  TClonesArray *clus   = 0;
 
  tracks = dynamic_cast<TClonesArray*>(l->FindObject(fTracksName));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracksName.Data() ));
    return;
  }
  clus = dynamic_cast<TClonesArray*>(l->FindObject(fCaloName));
  if (!clus) {
    AliError(Form("Pointer to clus %s == 0", fCaloName.Data() ));
    return;
  }

  // get primary vertex
  Double_t vertex[3] = {0, 0, 0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  // loop over clusters and tracks
  const Int_t Nclus = clus->GetEntries();
  const Int_t Ntrks = tracks->GetEntries();
 
  if (fDoTrackClus) {
    for(Int_t t = 0; t<Ntrks; ++t) {
      AliVTrack *track = dynamic_cast<AliVTrack*>(tracks->At(t));
      if (!track)
        continue;
      Int_t Nmatches = 0;
      Double_t dEtaMin  = 1e9;
      Double_t dPhiMin  = 1e9;
      Int_t    imin     = -1;
      for(Int_t i=0; i < Nclus; ++i) {
        AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(i));
        if (!c)
          continue;
        Double_t etadiff=999;
        Double_t phidiff=999;
        AliPicoTrack::GetEtaPhiDiff(track,c,phidiff,etadiff);
        Double_t dR = TMath::Sqrt(etadiff*etadiff+phidiff*phidiff);
        Double_t dRmin = TMath::Sqrt(dEtaMin*dEtaMin+dPhiMin*dPhiMin);
        if(dR > 25) 
          continue;
	if(dR<dRmin){
          dEtaMin = etadiff;
          dPhiMin = phidiff;
          imin = i;
        }
	if (TMath::Abs(phidiff)<0.05 && TMath::Abs(etadiff)<0.025) { // pp cuts!!!
	  Nmatches++;
	}
      }

      if (fCreateHisto)
	fHistNMatchCent_trk->Fill(cent,Nmatches);

      track->SetEMCALcluster(imin);
    }
         
  }

  for (Int_t iClus = 0, clusCount=0; iClus < Nclus; ++iClus) {
    AliVCluster *c = dynamic_cast<AliVCluster*>(clus->At(iClus));
    if (!c)
      continue;
    if (!c->IsEMCAL())
      continue;

    // make primary particle out of cluster
    TLorentzVector nPart;
    c->GetMomentum(nPart, vertex);

    Double_t etclus = nPart.Pt();
    if (etclus<fMinPt) 
      continue;

    Double_t energyclus = nPart.P();

    // reset cluster/track matching
    c->SetEmcCpvDistance(-1);
    c->SetTrackDistance(999,999);

    // run cluster-track matching
    Int_t    imin       = -1;
    Double_t dEtaMin    = 1e9;
    Double_t dPhiMin    = 1e9;
    Double_t dRmin      = 1e9;
    Double_t totalTrkP  = 0.0; // count total track momentum
    Int_t    Nmatches   = 0;   // count total number of matches
    for (Int_t t = 0; t<Ntrks; ++t) {

      AliVTrack *track = dynamic_cast<AliVTrack*>(tracks->At(t));

      if (!track)
        continue;
      if (!track->IsEMCAL())
        continue;
      if (track->Pt() < fMinPt)
        continue;

      Double_t etadiff = 999;
      Double_t phidiff = 999;
      AliPicoTrack::GetEtaPhiDiff(track, c, phidiff, etadiff);
      Double_t dR = TMath::Sqrt(etadiff * etadiff + phidiff * phidiff);
      if(dR<dRmin){
        dEtaMin = etadiff;
        dPhiMin = phidiff;
        dRmin = dR;
        imin = t;
      }
    
      Double_t mom = track->P();
      Int_t mombin = GetMomBin(mom);
      Int_t centbinch = centbin;
      if (track->Charge() == -1 || track->Charge() == 255) 
	centbinch += 4;

      if (fCreateHisto) {
	if (fHadCorr > 1 && mombin > -1) {
	  fHistMatchEtaPhi[centbinch][mombin]->Fill(etadiff, phidiff);
	  fHistMatchdRvsEP[centbin]->Fill(dR, energyclus / mom);
	}
      }
    
      Double_t EtaCut = 0.0;
      Double_t PhiCutlo = 0.0;
      Double_t PhiCuthi = 0.0;
      if(fPhiMatch > 0){
	PhiCutlo = -fPhiMatch;
	PhiCuthi = fPhiMatch;
      }
      else {
	PhiCutlo = GetPhiMean(mombin,centbinch)-GetPhiSigma(mombin,centbin);
	PhiCuthi = GetPhiMean(mombin,centbinch)+GetPhiSigma(mombin,centbin);
      }

      if(fEtaMatch > 0){
	EtaCut = fEtaMatch;
      }
      else {
	EtaCut = GetEtaSigma(mombin);
      }

      if ((phidiff < PhiCuthi && phidiff > PhiCutlo) && TMath::Abs(etadiff) < EtaCut) {
	if((fDoTrackClus && (track->GetEMCALcluster()) == iClus) || !fDoTrackClus){
	  ++Nmatches;
	  totalTrkP += track->P();
	}
      }
    }

    // store closest track
    c->SetEmcCpvDistance(imin);
    c->SetTrackDistance(dPhiMin, dEtaMin);

    if(fCreateHisto) {
      fHistNclusvsCent->Fill(cent);
      fHistEbefore->Fill(cent, energyclus);
      fHistNMatchCent->Fill(cent, Nmatches);

      if(Nmatches > 0) 
	fHistNclusMatchvsCent->Fill(cent);
    }
  
    // apply correction / subtraction
    if (fHadCorr > 0) {
      // to subtract only the closest track set fHadCor to a %
      // to subtract all tracks within the cut set fHadCor to %+1
      if (fHadCorr > 1) {
        if (totalTrkP > 0) {
          Double_t EoP  = energyclus / totalTrkP;
	  Double_t Esub = (fHadCorr-1)*totalTrkP;
	  if (Esub > energyclus) 
            Esub = energyclus;

	  if(fCreateHisto) {
	    fHistEoPCent->Fill(cent, EoP);
	    fHistMatchEvsP[centbin]->Fill(energyclus, EoP);
	    
	    if (Nmatches == 1) fHistEsubPchRat[centbin][0]->Fill(totalTrkP,Esub/totalTrkP);
	    else if(Nmatches == 2) fHistEsubPchRat[centbin][1]->Fill(totalTrkP,Esub/totalTrkP);
	    else fHistEsubPchRat[centbin][2]->Fill(totalTrkP,Esub/totalTrkP);
	    
	    if(Nmatches==1) fHistEsubPch[centbin][0]->Fill(totalTrkP,Esub);
	    else if(Nmatches==2) fHistEsubPch[centbin][1]->Fill(totalTrkP,Esub);
	    else fHistEsubPch[centbin][2]->Fill(totalTrkP,Esub);
	  }
        }
        energyclus -= (fHadCorr - 1) * totalTrkP;
      } 
      else if (imin >= 0) {
        AliVTrack *t = dynamic_cast<AliVTrack*>(tracks->At(imin));
        if (t) {
          Double_t mom = t->P();
          Int_t mombin = GetMomBin(mom);
          Int_t centbinch = centbin;
          if (t->Charge() == -1 || t->Charge() == 255) 
            centbinch += 4;

	  if (fCreateHisto) {
	    fHistMatchEtaPhi[centbinch][mombin]->Fill(dEtaMin, dPhiMin);
	    
	    if (mom > 0) {
	      fHistMatchEvsP[centbin]->Fill(energyclus, energyclus / mom);
	      fHistEoPCent->Fill(cent, energyclus / mom);
	      fHistMatchdRvsEP[centbin]->Fill(dRmin, energyclus / mom);
	    }
	  }

	  Double_t EtaCut = 0.0;
	  Double_t PhiCutlo = 0.0;
	  Double_t PhiCuthi = 0.0;
	  if(fPhiMatch > 0) {
	    PhiCutlo = -fPhiMatch;
	    PhiCuthi = fPhiMatch;
	  }
	  else {
	    PhiCutlo = GetPhiMean(mombin, centbinch) - GetPhiSigma(mombin, centbin);
	    PhiCuthi = GetPhiMean(mombin, centbinch) + GetPhiSigma(mombin, centbin);
	  }
	  
	  if(fEtaMatch > 0) {
	    EtaCut = fEtaMatch;
	  }
	  else {
	    EtaCut = GetEtaSigma(mombin);
	  }
	  
	  if ((dPhiMin<PhiCuthi && dPhiMin>PhiCutlo) && TMath::Abs(dEtaMin) < EtaCut) {
	    
	    if((fDoTrackClus && (t->GetEMCALcluster()) == iClus) || !fDoTrackClus){
	      energyclus -= fHadCorr * t->P();
	    }
          }
        }
      }
    }

    if (energyclus < 0) 
      energyclus = 0;

    if (fCreateHisto)
      fHistEafter->Fill(cent, energyclus);

    if (energyclus > 0) { // create corrected cluster
      AliVCluster *oc;
      if (esdMode) {
        oc = new ((*fOutClusters)[clusCount]) AliESDCaloCluster(*(dynamic_cast<AliESDCaloCluster*>(c)));
      } else {
        oc = new ((*fOutClusters)[clusCount]) AliAODCaloCluster(*(dynamic_cast<AliAODCaloCluster*>(c)));
      }
      oc->SetE(energyclus);
      ++clusCount;
    }
  }

  if (fCreateHisto)
    PostData(1, fOutputList);
}

//________________________________________________________________________
void AliHadCorrTask::Terminate(Option_t *) 
{
  // Nothing to be done.
}

