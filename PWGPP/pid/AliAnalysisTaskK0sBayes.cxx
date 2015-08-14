#include "AliAnalysisTaskK0sBayes.h"

// ROOT includes
#include <TMath.h>

// AliRoot includes
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliCentrality.h"
#include "AliVHeader.h"
#include "AliAODVZERO.h"
#include "TFile.h"
#include "AliOADBContainer.h"
#include "TH2F.h"
#include "TF1.h"
#include "AliGenHijingEventHeader.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "TChain.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliEventplane.h"
#include "AliAnalysisManager.h"
#include "TRandom.h"


Float_t AliAnalysisTaskK0sBayes::fPtKsMin[AliAnalysisTaskK0sBayes::nPtBin] = {0.,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4.,5.,6.};// ptmin bin
Float_t AliAnalysisTaskK0sBayes::fPtKsMax[AliAnalysisTaskK0sBayes::nPtBin] = {0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4.,5.,6.,10.};// ptmax bin

// STL includes
//#include <iostream>
//using namespace std;

ClassImp(AliAnalysisTaskK0sBayes)
//_____________________________________________________________________________
AliAnalysisTaskK0sBayes::AliAnalysisTaskK0sBayes():
  AliAnalysisTaskSE(),
  fVtxCut(10.0),  // cut on |vertex| < fVtxCut
  fEtaCut(0.8),   // cut on |eta| < fEtaCut
  fMinPt(0.15),   // cut on pt > fMinPt
  fIsMC(kFALSE),
  fQAsw(kFALSE),
  fNcluster(70),
  fFilterBit(1),
  fList(new TList()),
  fList2(new TList()),
  fList3(new TList()),
  fCentrality(-1),
  fPsi(0),
  fPtKs(0.),
  fPhiKs(0.),
  fEtaKs(0.),
  fMassV0(0.),
  fPtKp(0.),
  fPhiKp(0.),
  fEtaKp(0.),
  fPtKn(0.),
  fPhiKn(0.),
  fEtaKn(0.),
  fPidKp(0),
  fPidKn(0),
  fTOFTPCsignal(0),
  fCombsignal(0),
  fCutsDaughter(NULL),
  fPIDCombined(NULL),
  fContPid(NULL),
  fContPid2(NULL),
  fContUser(NULL),
  fContUser2(NULL),
  fNK0s(0),
  fNpiPos(0),
  fNpiNeg(0),
  fHmismTOF(0),
  fHchannelTOFdistr(0),
  fTypeCol(2),
  fPIDuserCut(NULL),
  fToEP(kFALSE),
  fSpeciesRef(2)
{
  // Default constructor (should not be used)
  fList->SetName("contKsBayes1");
  fList2->SetName("contKsBayes2");
  fList3->SetName("contKsBayes3");

  fList->SetOwner(kTRUE); 
  fList2->SetOwner(kTRUE); 
  fList3->SetOwner(kTRUE); 

  TFile *fmism = new TFile("$ALICE_ROOT/TOF/data/TOFmismatchDistr.root");
  fHmismTOF = (TH1F *) fmism->Get("TOFmismDistr");

  TFile *fchDist = new TFile("$ALICE_ROOT/TOF/data/TOFchannelDist.root");
  fHchannelTOFdistr = (TH1D *) fchDist->Get("hTOFchanDist"); 

  for(Int_t i=0;i < nCentrBin;i++){
    fElTOF[i] = NULL; 
    fElTPC[i] = NULL; 
    fPiTOF[i] = NULL; 
    fPiTPC[i] = NULL; 
    fKaTOF[i] = NULL; 
    fKaTPC[i] = NULL; 
    fPrTOF[i] = NULL; 
    fPrTPC[i] = NULL; 
  }
  for(Int_t i=0;i < 4;i++){
    hMatching[i] = NULL;
    hTracking[i] = NULL;
  }
  for(Int_t i=0;i < 1000;i++){
    fPhiK0s[i] = 0.0;
    fPtK0s[i] = 0.0;
    fIPiPos[i] = 0;
    fIPiNeg[i] = 0;
    fIpiP[i] = 0;
    fIpiN[i] = 0;
    fMassKs[i] = 0.0;
  }
}

//______________________________________________________________________________
AliAnalysisTaskK0sBayes::AliAnalysisTaskK0sBayes(const char *name):
  AliAnalysisTaskSE(name),
  fVtxCut(10.0),  // cut on |vertex| < fVtxCut
  fEtaCut(0.8),   // cut on |eta| < fEtaCut
  fMinPt(0.15),   // cut on pt > fMinPt
  fIsMC(kFALSE),
  fQAsw(kFALSE),
  fNcluster(70),
  fFilterBit(1),
  fList(new TList()),
  fList2(new TList()),
  fList3(new TList()),
  fCentrality(-1),
  fPsi(0),
  fPtKs(0.),
  fPhiKs(0.),
  fEtaKs(0.),
  fMassV0(0.),
  fPtKp(0.),
  fPhiKp(0.),
  fEtaKp(0.),
  fPtKn(0.),
  fPhiKn(0.),
  fEtaKn(0.),
  fPidKp(0),
  fPidKn(0),
  fTOFTPCsignal(0),
  fCombsignal(0),
  fCutsDaughter(NULL),
  fPIDCombined(NULL),
  fContPid(NULL),
  fContPid2(NULL),
  fContUser(NULL),
  fContUser2(NULL),
  fNK0s(0),
  fNpiPos(0),
  fNpiNeg(0),
  fHmismTOF(0),
  fHchannelTOFdistr(0),
  fTypeCol(2),
  fPIDuserCut(NULL),
  fToEP(kFALSE),
  fSpeciesRef(2)
{

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());

  // Output slot #1 writes into a TTree
  fList->SetName("contKsBayes1");
  fList2->SetName("contKsBayes2");
  fList3->SetName("contKsBayes3");

  fList->SetOwner(kTRUE); 
  fList2->SetOwner(kTRUE); 
  fList3->SetOwner(kTRUE); 

  TFile *fmism = new TFile("$ALICE_ROOT/TOF/data/TOFmismatchDistr.root");
  fHmismTOF = (TH1F *) fmism->Get("TOFmismDistr");

  TFile *fchDist = new TFile("$ALICE_ROOT/TOF/data/TOFchannelDist.root");
  fHchannelTOFdistr = (TH1D *) fchDist->Get("hTOFchanDist"); 

  for(Int_t i=0;i < nCentrBin;i++){
    fElTOF[i] = NULL; 
    fElTPC[i] = NULL; 
    fPiTOF[i] = NULL; 
    fPiTPC[i] = NULL; 
    fKaTOF[i] = NULL; 
    fKaTPC[i] = NULL; 
    fPrTOF[i] = NULL; 
    fPrTPC[i] = NULL; 
  }
  for(Int_t i=0;i < 4;i++){
    hMatching[i] = NULL;
    hTracking[i] = NULL;
  }
  for(Int_t i=0;i < 1000;i++){
    fPhiK0s[i] = 0.0;
    fPtK0s[i] = 0.0;
    fIPiPos[i] = 0;
    fIPiNeg[i] = 0;
    fIpiP[i] = 0;
    fIpiN[i] = 0;
    fMassKs[i] = 0.0;
  }
}
//_____________________________________________________________________________
AliAnalysisTaskK0sBayes::~AliAnalysisTaskK0sBayes()
{

}

//______________________________________________________________________________
void AliAnalysisTaskK0sBayes::UserCreateOutputObjects()
{

  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);

  Float_t invMmin = 0.497-0.005*8;
  Float_t invMmax = 0.497+0.005*8;

  const Int_t nBinPid = 14;

  Int_t binPid[nBinPid] = {1/*ptKs*/,1+7*(!fToEP)/*EtaPiP*/,20/*pt+*/,1/*pt-*/,5/*P+*/,1/*P-*/,2/*TOFmatch+*/,1/*TOFmatch-*/,2/*istrue*/,4/*Nsigma+*/,1/*Nsigma-*/,1+4*(fToEP)/*DeltaPhi+*/,1/*DeltaPhi-*/,1+4*(fToEP)/*Psi*/};

  Int_t binPid2[nBinPid] = {1/*ptKs*/,1+7*(!fToEP)/*EtaPiN*/,1/*pt+*/,20/*pt-*/,1/*P+*/,5/*P-*/,1/*TOFmatch+*/,2/*TOFmatch-*/,2/*istrue*/,1/*Nsigma+*/,4/*Nsigma-*/,1/*DeltaPhi+*/,1+4*(fToEP)/*DeltaPhi-*/,1+4*(fToEP)/*Psi*/};

  fContPid = new AliPIDperfContainer("contPID",nBinPid,binPid);
  fContPid->SetTitleX("M_{K^{0}_{s}}");
  fContPid->SetTitleY("centrality (%)");
  fContPid->SetVarName(0,"p_{T}^{K^{0}_{s}}");
  fContPid->SetVarRange(0,0,10);
  fContPid->SetVarName(1,"#eta^{K^{0}_{s}}");
  fContPid->SetVarRange(1,-0.8,0.8);
  fContPid->SetVarName(2,"p_{T}^{Kp}");
  fContPid->SetVarRange(2,0.3,4.3);
  fContPid->SetVarName(3,"p_{T}^{Kn}");
  fContPid->SetVarRange(3,0.3,4.3);
  fContPid->SetVarName(4,"BayesProb^{Kp}");
  fContPid->SetVarRange(4,0,1.);
  fContPid->SetVarName(5,"BayesProb^{Kn}");
  fContPid->SetVarRange(5,0,1.);
  fContPid->SetVarName(6,"isTOFmatch^{Kp}");
  fContPid->SetVarRange(6,-0.5,1.5);
  fContPid->SetVarName(7,"isTOFmatch^{Kn}");
  fContPid->SetVarRange(7,-0.5,1.5);
  fContPid->SetVarName(8,"isKsTrue^{Kn}");
  fContPid->SetVarRange(8,-0.5,1.5);
  fContPid->SetVarName(9,"N#sigma^{Kp}");
  fContPid->SetVarRange(9,1.25,6.25);
  fContPid->SetVarName(10,"N#sigma^{Kn}");
  fContPid->SetVarRange(10,1.25,6.25);
  fContPid->SetVarName(11,"#Delta#phi^{Kp}");
  fContPid->SetVarRange(11,-TMath::Pi(),TMath::Pi());
  fContPid->SetVarName(12,"#Delta#phi^{Kn}");
  fContPid->SetVarRange(12,-TMath::Pi(),TMath::Pi());
  fContPid->SetVarName(13,"#Psi");
  fContPid->SetVarRange(13,-TMath::Pi()/2,TMath::Pi()/2);

  fContPid2 = new AliPIDperfContainer("contPID2",nBinPid,binPid2);
  fContPid2->SetTitleX("M_{K^{0}_{s}}");
  fContPid2->SetTitleY("centrality (%)");
  fContPid2->SetVarName(0,"p_{T}^{K^{0}_{s}}");
  fContPid2->SetVarRange(0,0,10);
  fContPid2->SetVarName(1,"#eta^{K^{0}_{s}}");
  fContPid2->SetVarRange(1,-0.8,0.8);
  fContPid2->SetVarName(2,"p_{T}^{Kp}");
  fContPid2->SetVarRange(2,0.3,4.3);
  fContPid2->SetVarName(3,"p_{T}^{Kn}");
  fContPid2->SetVarRange(3,0.3,4.3);
  fContPid2->SetVarName(4,"BayesProb^{Kp}");
  fContPid2->SetVarRange(4,0,1.);
  fContPid2->SetVarName(5,"BayesProb^{Kn}");
  fContPid2->SetVarRange(5,0,1.);
  fContPid2->SetVarName(6,"isTOFmatch^{Kp}");
  fContPid2->SetVarRange(6,-0.5,1.5);
  fContPid2->SetVarName(7,"isTOFmatch^{Kn}");
  fContPid2->SetVarRange(7,-0.5,1.5);
  fContPid2->SetVarName(8,"isKsTrue^{Kn}");
  fContPid2->SetVarRange(8,-0.5,1.5);
  fContPid2->SetVarName(9,"N#sigma^{Kp}");
  fContPid2->SetVarRange(9,1.25,6.25);
  fContPid2->SetVarName(10,"N#sigma^{Kn}");
  fContPid2->SetVarRange(10,1.25,6.25);
  fContPid2->SetVarName(11,"#Delta#phi^{Kp}");
  fContPid2->SetVarRange(11,-TMath::Pi(),TMath::Pi());
  fContPid2->SetVarName(12,"#Delta#phi^{Kn}");
  fContPid2->SetVarRange(12,-TMath::Pi(),TMath::Pi());
  fContPid2->SetVarName(13,"#Psi");
  fContPid2->SetVarRange(13,-TMath::Pi()/2,TMath::Pi()/2);

  const Int_t nDETsignal = 100; // mass
  Double_t binDETsignal[nDETsignal+1];
  for(Int_t i=0;i<nDETsignal+1;i++){
    binDETsignal[i] = invMmin + i*(invMmax - invMmin) / nDETsignal;
  }
  const Int_t nDETsignal2 = 10; // centrality
  Double_t binDETsignal2[nDETsignal2+1];
  for(Int_t i=0;i<nDETsignal2+1;i++){
    binDETsignal2[i] = i*100./ nDETsignal2;
  }
  fContPid->AddSpecies("K0s",nDETsignal,binDETsignal,nDETsignal2,binDETsignal2);
  fContPid2->AddSpecies("K0s2",nDETsignal,binDETsignal,nDETsignal2,binDETsignal2);

  fList->Add(fContPid);
  fList->Add(fContPid2);

  const Int_t nBinUser = 6;
  Int_t binUser[nBinUser] = {8/*Eta*/,20/*pt*/,2/*istrue*/,4/*whatSelection*/,1/*DeltaPhi*/,1/*Psi*/};
  fContUser = new AliPIDperfContainer("contUserPID",nBinUser,binUser);
  fContUser->SetTitleX("M_{K^{0}_{s}}");
  fContUser->SetTitleY("centrality (%)");
  fContUser->SetVarName(0,"#eta^{K^{0}_{s}}");
  fContUser->SetVarRange(0,-0.8,0.8);
  fContUser->SetVarName(1,"p_{T}");
  fContUser->SetVarRange(1,0.3,4.3);
  fContUser->SetVarName(2,"isKsTrue^{Kn}");
  fContUser->SetVarRange(2,-0.5,1.5);
  fContUser->SetVarName(3,"whatSelected"); // 0=no, 1=pi, 2=K, 3=p
  fContUser->SetVarRange(3,-0.5,3.5);
  fContUser->SetVarName(4,"#Delta#phi");
  fContUser->SetVarRange(4,-TMath::Pi(),TMath::Pi());
  fContUser->SetVarName(5,"#Psi");
  fContUser->SetVarRange(5,-TMath::Pi()/2,TMath::Pi()/2);

  fContUser2 = new AliPIDperfContainer("contUserPID2",nBinUser,binUser);
  fContUser2->SetTitleX("M_{K^{0}_{s}}");
  fContUser2->SetTitleY("centrality (%)");
  fContUser2->SetVarName(0,"#eta^{K^{0}_{s}}");
  fContUser2->SetVarRange(0,-0.8,0.8);
  fContUser2->SetVarName(1,"p_{T}");
  fContUser2->SetVarRange(1,0.3,4.3);
  fContUser2->SetVarName(2,"isKsTrue^{Kn}");
  fContUser2->SetVarRange(2,-0.5,1.5);
  fContUser2->SetVarName(3,"whatSelected");
  fContUser2->SetVarRange(3,-0.5,3.5);
  fContUser2->SetVarName(4,"#Delta#phi");
  fContUser2->SetVarRange(4,-TMath::Pi(),TMath::Pi());
  fContUser2->SetVarName(5,"#Psi");
  fContUser2->SetVarRange(5,-TMath::Pi()/2,TMath::Pi()/2);

  fContUser->AddSpecies("K0s",nDETsignal,binDETsignal,nDETsignal2,binDETsignal2);
  fContUser2->AddSpecies("K0s2",nDETsignal,binDETsignal,nDETsignal2,binDETsignal2);

  fList->Add(fContUser);
  fList->Add(fContUser2);

  hMatching[0] = new TH2F("hMatchAll","TOF matched (all);p_{T} (GeV/#it{c});centrality (%)",50,0,10,nDETsignal2,0,100);
  hMatching[1] = new TH2F("hMatchPi","TOF matched (#pi);p_{T} (GeV/#it{c});centrality (%)",50,0,10,nDETsignal2,0,100);
  hMatching[2] = new TH2F("hMatchKa","TOF matched (K);p_{T} (GeV/#it{c});centrality (%)",50,0,10,nDETsignal2,0,100);
  hMatching[3] = new TH2F("hMatchPr","TOF matched (p);p_{T} (GeV/#it{c});centrality (%)",50,0,10,nDETsignal2,0,100);

  hTracking[0] = new TH2F("hTrackingAll","TPC tracks (all);p_{T} (GeV/#it{c});centrality (%)",50,0,10,nDETsignal2,0,100);
  hTracking[1] = new TH2F("hTrackingPi","TPC tracks (#pi);p_{T} (GeV/#it{c});centrality (%)",50,0,10,nDETsignal2,0,100);
  hTracking[2] = new TH2F("hTrackingKa","TPC tracks (K);p_{T} (GeV/#it{c});centrality (%)",50,0,10,nDETsignal2,0,100);
  hTracking[3] = new TH2F("hTrackingPr","TPC tracks (p);p_{T} (GeV/#it{c});centrality (%)",50,0,10,nDETsignal2,0,100);

  fList2->Add(hMatching[0]);
  fList2->Add(hMatching[1]);
  fList2->Add(hMatching[2]);
  fList2->Add(hMatching[3]);
  fList2->Add(hTracking[0]);
  fList2->Add(hTracking[1]);
  fList2->Add(hTracking[2]);
  fList2->Add(hTracking[3]);

  fTOFTPCsignal = new TH2F("hTOFTPCsignal","TOF-TPC signal for pions (0.8 < p_{T} < 0.9 GeV/#it{c}, cent. = 0-20);N#sigma_{TOF};N#sigma_{TPC}",100,-5,5,100,-5,5);
  fList2->Add(fTOFTPCsignal);
  fCombsignal = new TH1F("hCombsignal","Combined signal for pions (0.8 < p_{T} < 0.9 GeV/#it{c}, cent. = 0-20);N#sigma",100,0,5);
  fList2->Add(fCombsignal);

  // QA plots
  char name[100],title[400];
  for(Int_t i=0;i < nCentrBin;i++){
    snprintf(name,100,"hPiTPCQA_%i",i);
    snprintf(title,400,"TPC signal for #pi cent=%i-%i%c;p_{T} (GeV/#it{c});N#sigma",i*10,(i+1)*10,'%');
    fPiTPC[i] = new TH2F(name,title,50,0,4,200,-10,10);
    fList3->Add(fPiTPC[i]);

    snprintf(name,100,"hKaTPCQA_%i",i);
    snprintf(title,400,"TPC signal for K cent=%i-%i%c;p_{T} (GeV/#it{c});N#sigma",i*10,(i+1)*10,'%');
    fKaTPC[i] = new TH2F(name,title,50,0,4,200,-10,10);
    fList3->Add(fKaTPC[i]);

    snprintf(name,100,"hPrTPCQA_%i",i);
    snprintf(title,400,"TPC signal for p cent=%i-%i%c;p_{T} (GeV/#it{c});N#sigma",i*10,(i+1)*10,'%');
    fPrTPC[i] = new TH2F(name,title,50,0,4,200,-10,10);
    fList3->Add(fPrTPC[i]);

    snprintf(name,100,"hElTPCQA_%i",i);
    snprintf(title,400,"TPC signal for e cent=%i-%i%c;p_{T} (GeV/#it{c});N#sigma",i*10,(i+1)*10,'%');
    fElTPC[i] = new TH2F(name,title,50,0,4,200,-10,10);
    fList3->Add(fElTPC[i]);

    snprintf(name,100,"hPiTOFQA_%i",i);
    snprintf(title,400,"TOF signal for #pi cent=%i-%i%c;p_{T} (GeV/#it{c});N#sigma",i*10,(i+1)*10,'%');
    fPiTOF[i] = new TH2F(name,title,50,0,4,200,-10,10);
    fList3->Add(fPiTOF[i]);

    snprintf(name,100,"hKaTOFQA_%i",i);
    snprintf(title,400,"TOF signal for K cent=%i-%i%c;p_{T} (GeV/#it{c});N#sigma",i*10,(i+1)*10,'%');
    fKaTOF[i] = new TH2F(name,title,50,0,4,200,-10,10);
    fList3->Add(fKaTOF[i]);

    snprintf(name,100,"hPrTOFQA_%i",i);
    snprintf(title,400,"TOF signal for p cent=%i-%i%c;p_{T} (GeV/#it{c});N#sigma",i*10,(i+1)*10,'%');
    fPrTOF[i] = new TH2F(name,title,50,0,4,200,-10,10);
    fList3->Add(fPrTOF[i]);

    snprintf(name,100,"hElTOFQA_%i",i);
    snprintf(title,400,"TOF signal for e cent=%i-%i%c;p_{T} (GeV/#it{c});N#sigma",i*10,(i+1)*10,'%');
    fElTOF[i] = new TH2F(name,title,50,0,4,200,-10,10);
    fList3->Add(fElTOF[i]);

  }

  // Post output data.
  PostData(1, fList);
  PostData(2, fList2);
  PostData(3, fList3);
}

//______________________________________________________________________________
void AliAnalysisTaskK0sBayes::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
    
  if(fIsMC){
    AliVEvent *vEvent = InputEvent();
    for (Int_t itrack=0; itrack<vEvent->GetNumberOfTracks(); ++itrack) {
      AliAODTrack *aodTrack = const_cast<AliAODTrack*>(static_cast<AliAODTrack*>(vEvent->GetTrack(itrack)));
      if (!aodTrack) continue;
      AliAODPid *aodPid = const_cast<AliAODPid*>(aodTrack->GetDetPid());
      if (!aodPid) continue;
      aodPid->SetTPCsignalN(136.5/126.5 * aodPid->GetTPCsignalN());
    } 
  }

    fOutputAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fOutputAOD){
	Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
	this->Dump();
	return;
    }
    
    Float_t zvtx = GetVertex(fOutputAOD);



    //Get the MC object
    if(fIsMC){
      AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(fOutputAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
      if (!mcHeader) {
	AliError("Could not find MC Header in AOD");
	return;
      }
    }

    if (TMath::Abs(zvtx) < fVtxCut) {

      SelectK0s();

      //Centrality
      Float_t v0Centr  = -10.;
      Float_t trkCentr  = -10.;
      AliCentrality *centrality = fOutputAOD->GetCentrality();
      if (centrality){
	trkCentr  = centrality->GetCentralityPercentile("V0M");
	v0Centr = centrality->GetCentralityPercentile("TRK"); 
      }

      if(!fTypeCol){
	v0Centr=100./(fOutputAOD->GetNumberOfTracks()/12.+1);
	trkCentr=v0Centr;
      }

      if((TMath::Abs(v0Centr - trkCentr) < 5.0 || (fTypeCol!=2)) && v0Centr>0){ // consistency cut on centrality selection
        fCentrality = v0Centr;
	Analyze(fOutputAOD); // Do analysis!!!!

      }
    }
    
}

//________________________________________________________________________
void AliAnalysisTaskK0sBayes::Analyze(AliAODEvent* aodEvent)
{

  Int_t ntrack = aodEvent->GetNumberOfTracks();

  fPtKp=0.,fPhiKp=0.,fEtaKp=0.;
  fPtKn=0.,fPhiKn=0.,fEtaKn=0.;
  fPidKp=0,fPidKn=0;
  fMassV0=-1;
  
  TLorentzVector ksV;
  
  Double_t px,py,pz,E;

  Float_t invMmin = 0.497-0.005*8;
  Float_t invMmax = 0.497+0.005*8;
  
  Int_t icentr = 8;
  if(fCentrality < 0) icentr = 8;
  else if(fCentrality < 10) icentr = 0;
  else if(fCentrality < 20) icentr = 1;
  else if(fCentrality < 30) icentr = 2;
  else if(fCentrality < 40) icentr = 3;
  else if(fCentrality < 50) icentr = 4;
  else if(fCentrality < 60) icentr = 5;
  else if(fCentrality < 70) icentr = 6;
  else if(fCentrality < 80) icentr = 7;

  Float_t addMismatchForMC = 0.005;
  if(fCentrality < 50) addMismatchForMC += 0.005;
  if(fCentrality < 20) addMismatchForMC += 0.02;

  if(fTypeCol == 0 || fTypeCol == 1) addMismatchForMC = 0.005;

  fPsi = 0;
  /* Compute TPC EP */
  Double_t Qx2 = 0, Qy2 = 0;
  Double_t Qx3 = 0, Qy3 = 0;
  for(Int_t iT = 0; iT < ntrack; iT++) {
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iT));
    if(!aodTrack) AliFatal("Not a standard AOD");
    
    if (!aodTrack){
      continue;
    }
    
    Bool_t trkFlag = aodTrack->TestFilterBit(fFilterBit);
    
    if ((TMath::Abs(aodTrack->Eta()) > 0.8) || (aodTrack->Pt() < 0.2) || (aodTrack->GetTPCNcls() < fNcluster)  || !trkFlag) 
      continue;
    
    Double_t b[2] = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
    if (!aodTrack->PropagateToDCA(aodEvent->GetPrimaryVertex(), aodEvent->GetMagneticField(), 100., b, bCov))
      continue;
    
    if ((TMath::Abs(b[0]) > 3.0) || (TMath::Abs(b[1]) > 2.4))
      continue;
    
    Qx2 += TMath::Cos(2*aodTrack->Phi()); 
    Qy2 += TMath::Sin(2*aodTrack->Phi());
    Qx3 += TMath::Cos(3*aodTrack->Phi()); 
    Qy3 += TMath::Sin(3*aodTrack->Phi());
    
  }

  fPsi = TMath::ATan2(Qy2, Qx2)/2.;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  AliPIDResponse *PIDResponse=inputHandler->GetPIDResponse();

  PIDResponse->SetTOFResponse(aodEvent,AliPIDResponse::kTOF_T0);

  PIDResponse->GetTOFResponse().SetTOFtailAllPara(-23,1.1);

  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);

  Double_t probP[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t probN[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t nSigmaTPC,nSigmaTOF=6,nSigmaTPC2,nSigmaTOF2=6,nSigmaComb,nSigmaComb2;
  Double_t nSigmaTPCRef,nSigmaTOFRef=6,nSigmaTPC2Ref,nSigmaTOF2Ref=6,nSigmaCombRef,nSigmaComb2Ref;

  
  AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
  TClonesArray *mcArray = NULL;
  if (mcHeader)
    mcArray = (TClonesArray*)aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());

  Int_t nmc = 0;
  if(mcArray)
    nmc = mcArray->GetEntries();

  for(Int_t i=0;i < ntrack;i++){ // loop on tracks
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
    if(!track) AliFatal("Not a standard AOD");
    
    AliAODMCParticle *mcp = NULL;
    Int_t pdg = 0;
    
    if (!track){
      continue;
    }
    
    Int_t tofMatch = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    
    Int_t label = -1;
    if(mcArray){
      label = track->GetLabel();
      if(label != -1 && label < nmc){
	label = TMath::Abs(label);
	mcp = (AliAODMCParticle*)mcArray->At(label);
	pdg = TMath::Abs(mcp->GetPdgCode());
      }
      else
	label = -1;
    }
    else{
      /*UInt_t detUsed =*/ fPIDCombined->ComputeProbabilities(track, PIDResponse, probP);
    }
    
    if(track->TestFilterBit(fFilterBit) && TMath::Abs(track->Eta()) < 0.8 && track->Charge() > 0){
      hTracking[0]->Fill(track->Pt(),fCentrality);
      if(pdg == 211)
	hTracking[1]->Fill(track->Pt(),fCentrality);
      else if(pdg == 321)
	hTracking[2]->Fill(track->Pt(),fCentrality);
      else if(pdg == 2212)
	hTracking[3]->Fill(track->Pt(),fCentrality);
      else if(! mcp){ // Fill matching histo with the prob
	hTracking[1]->Fill(track->Pt(),fCentrality,probP[2]);
	hTracking[2]->Fill(track->Pt(),fCentrality,probP[3]);
	hTracking[3]->Fill(track->Pt(),fCentrality,probP[4]);
      }
    }
    
    if(!tofMatch) continue;
    
    if(track->TestFilterBit(fFilterBit) && TMath::Abs(track->Eta()) < 0.8 && track->Charge() > 0){
      hMatching[0]->Fill(track->Pt(),fCentrality);
      if(pdg == 211)
	hMatching[1]->Fill(track->Pt(),fCentrality);
      else if(pdg == 321)
	hMatching[2]->Fill(track->Pt(),fCentrality);
      else if(pdg == 2212)
	hMatching[3]->Fill(track->Pt(),fCentrality);
      else if(! mcp){ // Fill matching histo with the prob
	hMatching[1]->Fill(track->Pt(),fCentrality,probP[2]);
	hMatching[2]->Fill(track->Pt(),fCentrality,probP[3]);
	hMatching[3]->Fill(track->Pt(),fCentrality,probP[4]);
      }
    }
  }
  
//   Int_t pdg1 = -1;
//   Int_t pdg2 = -1;


  // start analysis K0s
  for(Int_t i=0;i < ntrack;i++){ // loop on positive tracks
    AliAODTrack* KpTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
    if(!KpTrack) AliFatal("Not a standard AOD");
        
    if (!KpTrack){
      continue;
    }
    
    if(!(KpTrack->Charge() > 0 && KpTrack->Pt() > 0.3  && TMath::Abs(KpTrack->Eta()) < 0.8)) continue;

    nSigmaComb=5;
    nSigmaCombRef=5;
    nSigmaTOF = 5;
    nSigmaTOFRef = 5;
    fPtKp=KpTrack->Pt(),fPhiKp=KpTrack->Phi(),fEtaKp=KpTrack->Eta();
    fPidKp=0;

    UInt_t detUsedP = fPIDCombined->ComputeProbabilities(KpTrack, PIDResponse, probP);

    Double_t oldpP[10];
    fPIDCombined->GetPriors(KpTrack, oldpP, PIDResponse, detUsedP);

    nSigmaTPC = PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kPion);
    nSigmaTPCRef = PIDResponse->NumberOfSigmasTPC(KpTrack,(AliPID::EParticleType) fSpeciesRef);

    if(! (TMath::Abs(nSigmaTPC) < 5)) continue;

    Int_t tofMatch1 = (KpTrack->GetStatus() & AliVTrack::kTOFout) && (KpTrack->GetStatus() & AliVTrack::kTIME);

    /*
    if(mcArray){
      Int_t labelK = TMath::Abs(KpTrack->GetLabel());
      AliAODMCParticle *mcp1 = (AliAODMCParticle*)mcArray->At(labelK);
      pdg1 = TMath::Abs(mcp1->GetPdgCode());
    }
    */

    fPidKp = Int_t(probP[2]*100);

    if(tofMatch1){
      if(!IsChannelValid(TMath::Abs(KpTrack->Eta()))){
	// remove this tof hit
	fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
	detUsedP = fPIDCombined->ComputeProbabilities(KpTrack, PIDResponse, probP);
	fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
	fPidKp = Int_t(probP[4]*100);
	tofMatch1=0;
      }
      else{
	if(probP[2] > probP[3] && probP[2] > probP[4] && probP[2] > probP[0]) fPidKp += 128; // max prob
	
	nSigmaTOF = PIDResponse->NumberOfSigmasTOF(KpTrack,AliPID::kProton);
	if(TMath::Abs(PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kProton))<1) fPrTOF[icentr]->Fill(fPtKp,nSigmaTOF);
	if(TMath::Abs(nSigmaTOF)<1) fPrTPC[icentr]->Fill(fPtKp,PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kProton));
	nSigmaTOF = PIDResponse->NumberOfSigmasTOF(KpTrack,AliPID::kElectron);
	if(TMath::Abs(PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kElectron))<1) fElTOF[icentr]->Fill(fPtKp,nSigmaTOF);
	if(TMath::Abs(nSigmaTOF)<1) fElTPC[icentr]->Fill(fPtKp,PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kElectron));
	nSigmaTOF = PIDResponse->NumberOfSigmasTOF(KpTrack,AliPID::kKaon);
	if(TMath::Abs(PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kKaon))<1) fKaTOF[icentr]->Fill(fPtKp,nSigmaTOF);
	if(TMath::Abs(nSigmaTOF)<1) fKaTPC[icentr]->Fill(fPtKp,PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kKaon));
	nSigmaTOF = PIDResponse->NumberOfSigmasTOF(KpTrack,AliPID::kPion);
	if(TMath::Abs(PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kPion))<1) fPiTOF[icentr]->Fill(fPtKp,nSigmaTOF);
	if(TMath::Abs(nSigmaTOF)<1) fPiTPC[icentr]->Fill(fPtKp,PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kPion));
	nSigmaTOFRef = PIDResponse->NumberOfSigmasTOF(KpTrack,(AliPID::EParticleType) fSpeciesRef);
		
	if(fIsMC){
	  Float_t mismAdd = addMismatchForMC;
	  if(KpTrack->Pt() < 1) mismAdd = addMismatchForMC/KpTrack->Pt();
	  
	  if(gRandom->Rndm() < mismAdd){
	    Float_t etaAbs = TMath::Abs(KpTrack->Eta());
	    Int_t channel = Int_t(4334.09 - 4758.36 * etaAbs -1989.71 * etaAbs*etaAbs + 1957.62*etaAbs*etaAbs*etaAbs);
	    channel = channel % 8736;
	    Float_t distIP = fHchannelTOFdistr->GetBinContent(channel);
	    
	    // generate random time
	    Float_t timeRandom = fHmismTOF->GetRandom() + distIP*3.35655419905265973e+01;
	    Double_t times[AliPID::kSPECIESC];
	    KpTrack->GetIntegratedTimes(times,AliPID::kSPECIESC);
	    nSigmaTOF = TMath::Abs(timeRandom - times[2])/PIDResponse->GetTOFResponse().GetExpectedSigma(KpTrack->P(), times[2], AliPID::kPion);
	    nSigmaTOFRef = TMath::Abs(timeRandom - times[fSpeciesRef])/PIDResponse->GetTOFResponse().GetExpectedSigma(KpTrack->P(), times[fSpeciesRef], fSpeciesRef);
	  }
	}

	if(fCentrality < 20 && KpTrack->Pt() < 0.9 && KpTrack->Pt() > 0.8)fTOFTPCsignal->Fill(nSigmaTOF,nSigmaTPC);
        nSigmaTOF = TMath::Abs(nSigmaTOF);

	if(nSigmaTOF < 2) fPidKp += 256;
	else if(nSigmaTOF < 3) fPidKp += 512;
      }
    }
    
    if(tofMatch1){
      nSigmaComb = TMath::Sqrt(0.5*(nSigmaTOF*nSigmaTOF + nSigmaTPC*nSigmaTPC));
      nSigmaCombRef = TMath::Sqrt(0.5*(nSigmaTOFRef*nSigmaTOFRef + nSigmaTPCRef*nSigmaTPCRef));
      if(nSigmaTOF < 5 && fCentrality < 20 && KpTrack->Pt() < 0.9 && KpTrack->Pt() > 0.8){
	fCombsignal->Fill(nSigmaComb);
      }
    } else {
      nSigmaComb = TMath::Abs(nSigmaTPC);
      nSigmaCombRef = TMath::Abs(nSigmaTPCRef);
    }

    nSigmaTOFRef = TMath::Abs(nSigmaTOFRef);

    // use sigmaTOF instead of sigmaComb
    //if(tofMatch1){
    //nSigmaComb = nSigmaTOF;
      //      nSigmaCombRef = nSigmaTOFRef;
    //}

    if(nSigmaComb < 2) nSigmaComb = 2;
    else if(nSigmaComb < 3) nSigmaComb = 3;
    else if(nSigmaComb < 5) nSigmaComb = 4.99;
    else nSigmaComb = 6;

    if(nSigmaCombRef < 2) nSigmaCombRef = 2;
    else if(nSigmaCombRef < 3) nSigmaCombRef = 3;
    else if(nSigmaCombRef < 5) nSigmaCombRef = 4.99;
    else nSigmaCombRef = 6;

    Int_t iks=-1;
    for(Int_t k=0;k < fNK0s;k++){ // find the K0s which contains the positive track
      if(i == fIpiP[k]) iks = k;
    }

    if(fPtKp > 4.299) fPtKp = 4.299;

    if(iks > -1 && fIpiN[iks] > -1){
      //for(Int_t j=0;j < ntrack;j++){ // loop on negative tracks
      Int_t j = fIpiN[iks];
      AliAODTrack* KnTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(j));
      if(!KnTrack) AliFatal("Not a standard AOD");
      
      if (!KnTrack){
	continue;
      }

      if(!(KnTrack->Charge() < 0 && KnTrack->Pt() > 0.3 && TMath::Abs(KnTrack->Eta()) < 0.8)) continue;

      fPtKn=KnTrack->Pt(),fPhiKn=KnTrack->Phi(),fEtaKn=KnTrack->Eta();
      fPidKn=0;

      UInt_t detUsedN = fPIDCombined->ComputeProbabilities(KnTrack, PIDResponse, probN);
      Double_t oldpN[10];
      fPIDCombined->GetPriors(KnTrack, oldpN, PIDResponse, detUsedN);

      nSigmaTPC2 = PIDResponse->NumberOfSigmasTPC(KnTrack,AliPID::kPion);
      nSigmaTPC2Ref = PIDResponse->NumberOfSigmasTPC(KnTrack,(AliPID::EParticleType) fSpeciesRef);
      
      if(! (TMath::Abs(nSigmaTPC2) < 5)) continue;
      
      nSigmaComb2=5;
      nSigmaTOF2=5;
      nSigmaComb2Ref=5;
      nSigmaTOF2Ref=5;

      Int_t tofMatch2 = (KnTrack->GetStatus() & AliVTrack::kTOFout) && (KnTrack->GetStatus() & AliVTrack::kTIME);
      /*
      if(mcArray){
	Int_t labelK = TMath::Abs(KnTrack->GetLabel());
	AliAODMCParticle *mcp2 = (AliAODMCParticle*)mcArray->At(labelK);
 	pdg2 = TMath::Abs(mcp2->GetPdgCode());
      }
      */

      fPidKn = Int_t(probN[2]*100);

      if(tofMatch2){
	if(!IsChannelValid(TMath::Abs(KnTrack->Eta()))){
	  // remove this tof hit
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
	  detUsedP = fPIDCombined->ComputeProbabilities(KnTrack, PIDResponse, probN);
	  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
	  fPidKn = Int_t(probN[4]*100);
	  tofMatch2=0;
	}
	else{
	  if(probN[2] > probN[3] && probN[2] > probN[4] && probN[2] > probN[0]) fPidKn += 128; // max prob
	  
	  nSigmaTOF2 = PIDResponse->NumberOfSigmasTOF(KnTrack,AliPID::kPion);	  	  
	  nSigmaTOF2 = TMath::Abs(nSigmaTOF2);
	  nSigmaTOF2Ref = PIDResponse->NumberOfSigmasTOF(KnTrack,(AliPID::EParticleType) fSpeciesRef);	  	  
	  nSigmaTOF2Ref = TMath::Abs(nSigmaTOF2Ref);
	  
	  if(fIsMC){
	    Float_t mismAdd = addMismatchForMC;
	    if(KnTrack->Pt() < 1) mismAdd = addMismatchForMC/KnTrack->Pt();
	    
	    if(gRandom->Rndm() < mismAdd){
	      Float_t etaAbs = TMath::Abs(KnTrack->Eta());
	      Int_t channel = Int_t(4334.09 - 4758.36 * etaAbs -1989.71 * etaAbs*etaAbs + 1957.62*etaAbs*etaAbs*etaAbs);
	      channel = channel % 8736;
	      Float_t distIP = fHchannelTOFdistr->GetBinContent(channel);
	      
	      // generate random time
	      Float_t timeRandom = fHmismTOF->GetRandom() + distIP*3.35655419905265973e+00;
	      Double_t times[AliPID::kSPECIESC];
	      KnTrack->GetIntegratedTimes(times,AliPID::kSPECIESC);
	      nSigmaTOF2 = TMath::Abs(timeRandom - times[2])/PIDResponse->GetTOFResponse().GetExpectedSigma(KnTrack->P(), times[2], AliPID::kPion);
	      nSigmaTOF2Ref = TMath::Abs(timeRandom - times[fSpeciesRef])/PIDResponse->GetTOFResponse().GetExpectedSigma(KnTrack->P(), times[fSpeciesRef], fSpeciesRef);
	    }
	  }

	  if(fCentrality < 20 && KnTrack->Pt() < 1.2 && KnTrack->Pt() > 1) fTOFTPCsignal->Fill(nSigmaTOF2,nSigmaTPC2);

	  if(nSigmaTOF2 < 2) fPidKn += 256;
	  else if(nSigmaTOF2 < 3) fPidKn += 512;
	}
      }

      px = KpTrack->Px() + KnTrack->Px();
      py = KpTrack->Py() + KnTrack->Py();
      pz = KpTrack->Pz() + KnTrack->Pz();
      E = TMath::Sqrt(KpTrack->P()*KpTrack->P() + 1.39e-01*1.39e-01);
      E += TMath::Sqrt(KnTrack->P()*KnTrack->P()+ 1.39e-01*1.39e-01);

      ksV.SetPxPyPzE(px,py,pz,E);
      fMassV0 = fMassKs[iks];
      
      if(fMassV0 <  invMmin || fMassV0 > invMmax) continue;


      fPtKs = ksV.Pt();
      fEtaKs = ksV.Eta();
      fPhiKs = ksV.Phi();

      if(tofMatch2){
	nSigmaComb2 = TMath::Sqrt(0.5*(nSigmaTOF2*nSigmaTOF2+ nSigmaTPC2*nSigmaTPC2));
	nSigmaComb2Ref = TMath::Sqrt(0.5*(nSigmaTOF2Ref*nSigmaTOF2Ref+ nSigmaTPC2Ref*nSigmaTPC2Ref));
	if(nSigmaTOF2 < 5 && fCentrality < 20 && KnTrack->Pt() < 1.2 && KnTrack->Pt() > 1){
	  fCombsignal->Fill(nSigmaComb2);
	}
      } else {
	nSigmaComb2 = TMath::Abs(nSigmaTPC2);
	nSigmaComb2Ref = TMath::Abs(nSigmaTPC2Ref);
      }

      // use sigmaTOF instead of sigmaComb
      //if(tofMatch2){
	//nSigmaComb2 = nSigmaTOF2;
	//nSigmaComb2Ref = nSigmaTOF2Ref;
      //}

      if(nSigmaComb2 < 2) nSigmaComb2 = 2;
      else if(nSigmaComb2 < 3) nSigmaComb2 = 3;
      else if(nSigmaComb2 < 5) nSigmaComb2 = 4.99;
      else nSigmaComb2 = 6;  
      if(nSigmaComb2Ref < 2) nSigmaComb2Ref = 2;
      else if(nSigmaComb2Ref < 3) nSigmaComb2Ref = 3;
      else if(nSigmaComb2Ref < 5) nSigmaComb2Ref = 4.99;
      else nSigmaComb2Ref = 6;  

      Bool_t isTrue = kFALSE;

      if(mcArray){
      	Int_t labelP = TMath::Abs(KpTrack->GetLabel());
      	Int_t labelN = TMath::Abs(KnTrack->GetLabel());

	if(labelP > -1 && labelN > -1){
	  AliAODMCParticle *partP = (AliAODMCParticle*)mcArray->At(labelP);
	  AliAODMCParticle *partN = (AliAODMCParticle*)mcArray->At(labelN);

	  Int_t mP = partP->GetMother();
	  Int_t mN = partN->GetMother();
	  if(mP == mN && mP > -1){
	    AliAODMCParticle *partM = (AliAODMCParticle*)mcArray->At(mP);
	    Int_t pdgM = partM->GetPdgCode();
	    if(pdgM == 310) isTrue = kTRUE;
	  }
	}

      }

      Double_t deltaphi1 = KpTrack->Phi() - fPsi;
      Double_t deltaphi2 = KnTrack->Phi() - fPsi;

      if(gRandom->Rndm() < 0.5){
	deltaphi1 += TMath::Pi();
	deltaphi2 += TMath::Pi();
      }

      while(deltaphi1 > TMath::Pi()) deltaphi1 -= TMath::Pi()*2;
      while(deltaphi1 < -TMath::Pi()) deltaphi1 += TMath::Pi()*2;
      while(deltaphi2 > TMath::Pi()) deltaphi2 -= TMath::Pi()*2;
      while(deltaphi2 < -TMath::Pi()) deltaphi2 += TMath::Pi()*2;

      if(fPtKn > 4.299) fPtKn = 4.299;

      Float_t xTOfill[] = {static_cast<Float_t>(fPtKs),static_cast<Float_t>(KpTrack->Eta()),static_cast<Float_t>(fPtKp),static_cast<Float_t>(fPtKn),static_cast<Float_t>(probP[2]),static_cast<Float_t>(probN[2]),static_cast<Float_t>(tofMatch1),static_cast<Float_t>(tofMatch2),static_cast<Float_t>(isTrue),static_cast<Float_t>(nSigmaComb),static_cast<Float_t>(nSigmaComb2),static_cast<Float_t>(deltaphi1),static_cast<Float_t>(deltaphi2),static_cast<Float_t>(fPsi)};
      Float_t xTOfill2[] = {static_cast<Float_t>(fPtKs),static_cast<Float_t>(KpTrack->Eta()),static_cast<Float_t>(fPtKp),static_cast<Float_t>(fPtKn),static_cast<Float_t>(probP[2]),static_cast<Float_t>(probN[2]),static_cast<Float_t>(tofMatch1),static_cast<Float_t>(tofMatch2),static_cast<Float_t>(isTrue),static_cast<Float_t>(nSigmaComb),static_cast<Float_t>(nSigmaComb2),static_cast<Float_t>(deltaphi1),static_cast<Float_t>(deltaphi2),static_cast<Float_t>(fPsi)};
      
      Int_t ipt = 0;
      while(ipt < nPtBin && fPtKsMin[ipt] < fPtKs){
	ipt++;
      }
      ipt--;
      if(ipt < 0) ipt = 0; // just to be sure

      if(TMath::Abs(fEtaKs) < 0.8 && fPtKp > 0.3 && fPtKn > 0.3){
	if(fSpeciesRef != 2){
	  xTOfill[4] = probP[fSpeciesRef];
	  xTOfill2[5] = probN[fSpeciesRef];

	  xTOfill[9] = nSigmaCombRef;
	  xTOfill2[10] = nSigmaComb2Ref;

	}

	fContPid->Fill(0,fMassV0,fCentrality,xTOfill);
        xTOfill[1] = KnTrack->Eta();
	fContPid2->Fill(0,fMassV0,fCentrality,xTOfill2);

	if(fPIDuserCut){
	  Float_t xUser[] = {static_cast<Float_t>(KpTrack->Eta()),static_cast<Float_t>(fPtKp),static_cast<Float_t>(isTrue),0,static_cast<Float_t>(deltaphi1),static_cast<Float_t>(fPsi)};
	  Float_t xUser2[] = {static_cast<Float_t>(KnTrack->Eta()),static_cast<Float_t>(fPtKn),static_cast<Float_t>(isTrue),0,static_cast<Float_t>(deltaphi2),static_cast<Float_t>(fPsi)};

	  if(fPIDuserCut->IsSelected(KpTrack,AliPID::kPion)){ // to be filled for positive
	    xUser[3] = 1;
	  } else if(fPIDuserCut->IsSelected(KpTrack,AliPID::kKaon)){
	    xUser[3] = 2;
	  } else if(fPIDuserCut->IsSelected(KpTrack,AliPID::kProton)){
	    xUser[3] = 3;
	  }
	  if(fPIDuserCut->IsSelected(KnTrack,AliPID::kPion)){ // to be filled for negative
	    xUser2[3] = 1;
	  } else if(fPIDuserCut->IsSelected(KnTrack,AliPID::kKaon)){
	    xUser2[3] = 2;
	  } else if(fPIDuserCut->IsSelected(KnTrack,AliPID::kProton)){
	    xUser2[3] = 3;
	  }
	  fContUser->Fill(0,fMassV0,fCentrality,xUser);
	  fContUser2->Fill(0,fMassV0,fCentrality,xUser2);
	}

      }


    }
  } // end analysi K0s
}

//_____________________________________________________________________________
Float_t AliAnalysisTaskK0sBayes::GetVertex(AliAODEvent* aod) const
{

  Float_t zvtx = -999;

  const AliAODVertex* vtxAOD = aod->GetPrimaryVertex();
  if (!vtxAOD)
    return zvtx;
  if(vtxAOD->GetNContributors()>0)
    zvtx = vtxAOD->GetZ();
  
  return zvtx;
}
//_____________________________________________________________________________
void AliAnalysisTaskK0sBayes::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
//=======================================================================
void AliAnalysisTaskK0sBayes::SelectK0s(){
  fNK0s=0;
  fNpiPos=0;
  fNpiNeg=0;

  Int_t nV0s = fOutputAOD->GetNumberOfV0s();
  AliAODv0 *myV0;
  Double_t dMASS=0.0;
  for (Int_t i=0; i!=nV0s; ++i) {
    myV0 = (AliAODv0*) fOutputAOD->GetV0(i);
    if(!myV0) continue;
    if(myV0->Pt()<0.1 || TMath::Abs(myV0->Eta()) > 0.8) continue; // skipping low momentum
    Int_t pass = PassesAODCuts(myV0,fOutputAOD,0); // check for K0s
    if(pass) {
      dMASS = myV0->MassK0Short();
      Float_t massLambda = myV0->MassLambda();
      Float_t massAntiLambda = myV0->MassAntiLambda();

      if(TMath::Abs(dMASS-0.497)/0.005 < 8 && TMath::Abs(massLambda-1.115)/0.005 > 8 && TMath::Abs(massAntiLambda-1.115)/0.005 > 8){
	AliAODTrack *iT=(AliAODTrack*) myV0->GetDaughter(0); // positive
	AliAODTrack *jT=(AliAODTrack*) myV0->GetDaughter(1); // negative
	if(iT->Charge()<0){
	  iT=(AliAODTrack*) myV0->GetDaughter(1); // positive
	  jT=(AliAODTrack*) myV0->GetDaughter(0); // negative
	}
	fPhiK0s[fNK0s] = myV0->Phi();
	fPtK0s[fNK0s] = myV0->Pt();
	fIpiP[fNK0s] = FindDaugheterIndex(iT);
	fIpiN[fNK0s] = FindDaugheterIndex(jT);
	fMassKs[fNK0s] = dMASS;
	if(fIpiP[fNK0s] > -1 && fIpiN[fNK0s] > -1)
	  fNK0s++;
      }
    }
  }

  /* My V0 code
  // fill pion stacks
  Int_t nAODTracks = fOutputAOD->GetNumberOfTracks();
  for(Int_t iT = 0; iT < nAODTracks; iT++) { // loop on the tracks
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(fOutputAOD->GetTrack(iT));
    if(!aodTrack) AliFatal("Not a standard AOD");
    
    if (!aodTrack){
      continue;
    }
    
    Bool_t trkFlag = aodTrack->TestFilterBit(fFilterBit);

    if ((TMath::Abs(aodTrack->Eta()) > fEtaCut) || (aodTrack->Pt() < fMinPt) || (aodTrack->GetTPCNcls() < fNcluster) || !trkFlag){
      continue;
    }

    Double_t b[2] = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
    if (!aodTrack->PropagateToDCA(fOutputAOD->GetPrimaryVertex(), fOutputAOD->GetMagneticField(), 100., b, bCov))
      continue;
    
    if(TMath::Abs(b[0]) < 0.3/aodTrack->Pt()) continue;


    Int_t charge = aodTrack->Charge();
    if(charge > 0){
      fIPiPos[fNpiPos] = iT;
      fNpiPos++;
    }
    else{
      fIPiNeg[fNpiNeg] = iT;
      fNpiNeg++;
    }     
  }
  
  for(Int_t i=0;i < fNpiPos;i++){
    AliAODTrack *pip = fOutputAOD->GetTrack(fIPiPos[i]);
    AliESDtrack pipE(pip);

    for(Int_t j=0;j < fNpiNeg;j++){
      AliAODTrack *pin = fOutputAOD->GetTrack(fIPiNeg[j]);
      AliESDtrack pinE(pin);

      Double_t xn, xp, mindist=pinE.GetDCA(&pipE,fOutputAOD->GetMagneticField(),xn,xp);

      Double_t pPos[3];
      Double_t pNeg[3];
      pipE.GetPxPyPzAt(xp,fOutputAOD->GetMagneticField(),pPos);
      pinE.GetPxPyPzAt(xn,fOutputAOD->GetMagneticField(),pNeg);

      Float_t length = (xp+xn)*0.5;

      Float_t pxs = pPos[0] + pNeg[0];
      Float_t pys = pPos[1] + pNeg[1];
      Float_t pzs = pPos[2] + pNeg[2];
      Float_t es = TMath::Sqrt(pPos[0]*pPos[0] + pPos[1]*pPos[1] + pPos[2]*pPos[2] + 0.13957*0.13957) + TMath::Sqrt(pNeg[0]*pNeg[0] + pNeg[1]*pNeg[1] + pNeg[2]*pNeg[2] + 0.13957*0.13957);

      Float_t pt = TMath::Sqrt(pxs*pxs + pys*pys);
      Float_t phi = TMath::ATan2(pys,pxs);
      Float_t mass = TMath::Sqrt(es*es - pt*pt - pzs*pzs);
      
      //      if(length > 1) printf("length = %f - distance = %f - mass= %f\n",length,mindist,mass);

      if(mindist < 0.4&& length > 0.7 && length < 25){
        Float_t esL = TMath::Sqrt(pPos[0]*pPos[0] + pPos[1]*pPos[1] + pPos[2]*pPos[2] + 0.938*0.938) + TMath::Sqrt(pNeg[0]*pNeg[0] + pNeg[1]*pNeg[1] + pNeg[2]*pNeg[2] + 0.13957*0.13957);
        Float_t esAL = TMath::Sqrt(pPos[0]*pPos[0] + pPos[1]*pPos[1] + pPos[2]*pPos[2] + 0.13957*0.13957) + TMath::Sqrt(pNeg[0]*pNeg[0] + pNeg[1]*pNeg[1] + pNeg[2]*pNeg[2] + 0.938*0.938);

        Float_t massaL = TMath::Sqrt(esL*esL - pt*pt - pzs*pzs);
        Float_t massaAL = TMath::Sqrt(esAL*esAL - pt*pt - pzs*pzs);

        if(TMath::Abs(mass-0.497)/0.005 < 8 && massaL > 1.15 && massaAL > 1.15){
          fPhiK0s[fNK0s] = phi;
          fPtK0s[fNK0s] = pt;
	  fIpiP[fNK0s] =fIPiPos[i] ;
	  fIpiN[fNK0s] = fIPiNeg[j];
          fMassKs[fNK0s] = mass;
          fNK0s++;
        }
      }
    }
  }
  */
}

//=======================================================================
Int_t AliAnalysisTaskK0sBayes::PassesAODCuts(AliAODv0 *myV0, AliAODEvent *tAOD,Int_t specie)
{ 
  if (myV0->GetOnFlyStatus() ) return 0;
  
  //the following is needed in order to evualuate track-quality
  AliAODTrack *iT, *jT;
  AliAODVertex *vV0s = myV0->GetSecondaryVtx();
  Double_t pos[3],cov[6];
  vV0s->GetXYZ(pos);
  vV0s->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  
  // TESTING CHARGE
  int iPos, iNeg;
  iT=(AliAODTrack*) myV0->GetDaughter(0);
  if(iT->Charge()>0) {
    iPos = 0; iNeg = 1;
  } else {
    iPos = 1; iNeg = 0;
  }
  // END OF TEST

  iT=(AliAODTrack*) myV0->GetDaughter(iPos); // positive

  jT=(AliAODTrack*) myV0->GetDaughter(iNeg); // negative

  Bool_t trkFlag = iT->TestFilterBit(fFilterBit);
  if(!trkFlag) return 0;
  Bool_t trkFlag2 = jT->TestFilterBit(fFilterBit);
  if(!trkFlag2) return 0;

  Double_t pvertex[3];
  pvertex[0]=tAOD->GetPrimaryVertex()->GetX();
  pvertex[1]=tAOD->GetPrimaryVertex()->GetY();
  pvertex[2]=tAOD->GetPrimaryVertex()->GetZ();

  Double_t dDL=myV0->DecayLengthV0( pvertex );
  if(dDL  < 0.5 || dDL > 25) return 0;

  Double_t dDCA=myV0->DcaV0Daughters();
  if(dDCA >0.5) return 0;

  Double_t dCTP=myV0->CosPointingAngle( pvertex );
  if(dCTP < -1) return 0;

//   AliESDtrack ieT( iT );
//   Double_t dD0P=ieT.GetD(pvertex[0],pvertex[1],tAOD->GetMagneticField());
//   if(TMath::Abs(dD0P) < 0]) return 0;

//   AliESDtrack jeT( jT );
//   Double_t dD0M=jeT.GetD(pvertex[0],pvertex[1],tAOD->GetMagneticField());
//   if(TMath::Abs(dD0M) < 0) return 0;

//   Double_t dD0D0=dD0P*dD0M;
//   if(dD0D0>0) return 0;

//   Double_t dETA=myV0->Eta();
//   if(dETA <-0.8) return 0;
//   if(dETA >0.8) return 0;

//   Double_t dQT=myV0->PtArmV0();
//   if(specie==0) if(dQT<???) return 0;

  Double_t dALPHA=myV0->AlphaV0(); // AlphaV0 -> AODRecoDecat::Alpha -> return 1.-2./(1.+QlProng(0)/QlProng(1));
  if(myV0->ChargeProng(iPos)<0) dALPHA = -dALPHA; // protects for a change in convention

  if(specie==1 && dALPHA<0) return 2; // antilambda
  return 1; // K0s or lambda
}
//-------------------------------------------------------------------------
Int_t AliAnalysisTaskK0sBayes::FindDaugheterIndex(AliAODTrack *trk){
  Int_t ntrack = fOutputAOD->GetNumberOfTracks();

  for(Int_t i=0;i < ntrack;i++){ // loop on tracks
    AliAODTrack* track = dynamic_cast<AliAODTrack*>(fOutputAOD->GetTrack(i));
    if(!track) AliFatal("Not a standard AOD");
    if(track == trk) return i;
  }
  
  printf("Daughter for %p not found\n",trk);
  return -1;
}
//-------------------------------------------------------------------------
Int_t AliAnalysisTaskK0sBayes::IsChannelValid(Float_t etaAbs){
  if(!fIsMC) return 1; // used only on MC

  if(fTypeCol == 2){ // LHC10h or LHC11h because of TOF matching window at 3 cm
    Int_t channel = Int_t(4334.09 - 4758.36 * etaAbs -1989.71 * etaAbs*etaAbs + 1957.62*etaAbs*etaAbs*etaAbs); 
  
    if(!(channel%20)) return 0; // 5% additional loss in MC
  }

  return 1;
}
