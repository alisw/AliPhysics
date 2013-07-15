#include "AliAnalysisTaskPhiBayes.h"

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


Float_t AliAnalysisTaskPhiBayes::fPtPhiMin[AliAnalysisTaskPhiBayes::nPtBin] = {0.,0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4.,5.,6.};// ptmin bin
Float_t AliAnalysisTaskPhiBayes::fPtPhiMax[AliAnalysisTaskPhiBayes::nPtBin] = {0.25,0.5,0.75,1,1.5,2,2.5,3,3.5,4.,5.,6.,10.};// ptmax bin

// STL includes
//#include <iostream>
//using namespace std;

ClassImp(AliAnalysisTaskPhiBayes)
//_____________________________________________________________________________
AliAnalysisTaskPhiBayes::AliAnalysisTaskPhiBayes():
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
  fPtPhi(0.),
  fPhiPhi(0.),
  fEtaPhi(0.),
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
  fHmismTOF(0),
  fHchannelTOFdistr(0),
  fTypeCol(2)
{
  // Default constructor (should not be used)
  fList->SetName("contPhiBayes1");
  fList2->SetName("contPhiBayes2");
  fList3->SetName("contPhiBayes3");

  fList->SetOwner(kTRUE); 
  fList2->SetOwner(kTRUE); 
  fList3->SetOwner(kTRUE); 

  TFile *fmism = new TFile("$ALICE_ROOT/TOF/data/TOFmismatchDistr.root");
  fHmismTOF = (TH1F *) fmism->Get("TOFmismDistr");

  TFile *fchDist = new TFile("$ALICE_ROOT/TOF/data/TOFchannelDist.root");
  fHchannelTOFdistr = (TH1D *) fchDist->Get("hTOFchanDist"); 
}

//______________________________________________________________________________
AliAnalysisTaskPhiBayes::AliAnalysisTaskPhiBayes(const char *name):
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
  fPtPhi(0.),
  fPhiPhi(0.),
  fEtaPhi(0.),
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
  fHmismTOF(0),
  fHchannelTOFdistr(0),
  fTypeCol(2)
{

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());

  // Output slot #1 writes into a TTree
  fList->SetName("contPhiBayes1");
  fList2->SetName("contPhiBayes2");
  fList3->SetName("contPhiBayes3");

  fList->SetOwner(kTRUE); 
  fList2->SetOwner(kTRUE); 
  fList3->SetOwner(kTRUE); 

  TFile *fmism = new TFile("$ALICE_ROOT/TOF/data/TOFmismatchDistr.root");
  fHmismTOF = (TH1F *) fmism->Get("TOFmismDistr");

  TFile *fchDist = new TFile("$ALICE_ROOT/TOF/data/TOFchannelDist.root");
  fHchannelTOFdistr = (TH1D *) fchDist->Get("hTOFchanDist"); 
}
//_____________________________________________________________________________
AliAnalysisTaskPhiBayes::~AliAnalysisTaskPhiBayes()
{

}

//______________________________________________________________________________
void AliAnalysisTaskPhiBayes::UserCreateOutputObjects()
{

  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF);

  Float_t invMmin = 0.985;
  Float_t invMmax = 1.045;

  const Int_t nBinPid = 14;

  Int_t binPid[nBinPid] = {1/*ptPhi*/,8/*EtaKP*/,20/*pt+*/,1/*pt-*/,5/*P+*/,1/*P-*/,2/*TOFmatch+*/,1/*TOFmatch-*/,2/*istrue*/,4/*Nsigma+*/,1/*Nsigma-*/,1/*DeltaPhi+*/,1/*DeltaPhi-*/,1/*Psi*/};

  Int_t binPid2[nBinPid] = {1/*ptPhi*/,8/*EtaKN*/,1/*pt+*/,20/*pt-*/,1/*P+*/,5/*P-*/,1/*TOFmatch+*/,2/*TOFmatch-*/,2/*istrue*/,1/*Nsigma+*/,4/*Nsigma-*/,1/*DeltaPhi+*/,1/*DeltaPhi-*/,1/*Psi*/};

  fContPid = new AliPIDperfContainer("contPID",nBinPid,binPid);
  fContPid->SetTitleX("M_{#phi}");
  fContPid->SetTitleY("centrality (%)");
  fContPid->SetVarName(0,"p_{T}^#phi}");
  fContPid->SetVarRange(0,0,10);
  fContPid->SetVarName(1,"#eta^{#phi}");
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
  fContPid->SetVarName(8,"isPhiTrue^{Kn}");
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
  fContPid2->SetTitleX("M_{#phi}");
  fContPid2->SetTitleY("centrality (%)");
  fContPid2->SetVarName(0,"p_{T}^{#phi}");
  fContPid2->SetVarRange(0,0,10);
  fContPid2->SetVarName(1,"#eta^{#phi}");
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
  fContPid2->SetVarName(8,"isPhiTrue^{Kn}");
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
  fContPid->AddSpecies("phi",nDETsignal,binDETsignal,nDETsignal2,binDETsignal2);
  fContPid2->AddSpecies("phi2",nDETsignal,binDETsignal,nDETsignal2,binDETsignal2);

  fList->Add(fContPid);
  fList->Add(fContPid2);

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
void AliAnalysisTaskPhiBayes::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
    
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

      //Centrality
      Float_t v0Centr  = -10.;
      Float_t trkCentr  = -10.;
      AliCentrality *centrality = fOutputAOD->GetCentrality();
      if (centrality){
	trkCentr  = centrality->GetCentralityPercentile("V0M");
	v0Centr = centrality->GetCentralityPercentile("TRK"); 
      }

      if((TMath::Abs(v0Centr - trkCentr) < 5.0 || (fTypeCol!=2)) && v0Centr>0){ // consistency cut on centrality selection
        fCentrality = v0Centr;
	Analyze(fOutputAOD); // Do analysis!!!!

      }
    }
    
}

//________________________________________________________________________
void AliAnalysisTaskPhiBayes::Analyze(AliAODEvent* aodEvent)
{
  Int_t ntrack = aodEvent->GetNumberOfTracks();

  fPtKp=0.,fPhiKp=0.,fEtaKp=0.;
  fPtKn=0.,fPhiKn=0.,fEtaKn=0.;
  fPidKp=0,fPidKn=0;
  fMassV0=-1;
  
  TLorentzVector phiV;
  
  Double_t px,py,pz,E;

  Float_t invMmin = 0.985;
  Float_t invMmax = 1.045;
  
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
  if(fCentrality < 50) addMismatchForMC += 0.01;
  if(fCentrality < 20) addMismatchForMC += 0.02;

  if(fTypeCol == 0) addMismatchForMC = 0.005;
  else if(fTypeCol == 1) addMismatchForMC *= 0.5;

  fPsi = 0;
  /* Compute TPC EP */
  Double_t Qx2 = 0, Qy2 = 0;
  Double_t Qx3 = 0, Qy3 = 0;
  for(Int_t iT = 0; iT < ntrack; iT++) {
    AliAODTrack* aodTrack = aodEvent->GetTrack(iT);
    
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

  PIDResponse->GetTOFResponse().SetTrackParameter(0,0.);
  PIDResponse->GetTOFResponse().SetTrackParameter(1,0.);
  PIDResponse->GetTOFResponse().SetTrackParameter(2,0.018);
  PIDResponse->GetTOFResponse().SetTrackParameter(3,50.0);

  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);

  Double_t probP[10],probN[10];
  Double_t nSigmaTPC,nSigmaTOF=6,nSigmaTPC2,nSigmaTOF2=6,nSigmaComb,nSigmaComb2;

  
  AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
  TClonesArray *mcArray = NULL;
  if (mcHeader)
    mcArray = (TClonesArray*)aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());

  Int_t nmc = 0;
  if(mcArray)
    nmc = mcArray->GetEntries();

  for(Int_t i=0;i < ntrack;i++){ // loop on tracks
    AliAODTrack* track = aodEvent->GetTrack(i);
    
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


  // start analysis phi
  for(Int_t i=0;i < ntrack;i++){ // loop on positive tracks
    AliAODTrack* KpTrack = aodEvent->GetTrack(i);
        
    if (!KpTrack){
      continue;
    }
    
    if(!(KpTrack->Charge() > 0 && KpTrack->Pt() > 0.3  && TMath::Abs(KpTrack->Eta()) < 0.8)) continue;

    nSigmaComb=5;
    fPtKp=KpTrack->Pt(),fPhiKp=KpTrack->Phi(),fEtaKp=KpTrack->Eta();
    fPidKp=0;

    UInt_t detUsedP = fPIDCombined->ComputeProbabilities(KpTrack, PIDResponse, probP);

    Double_t oldpP[10];
    fPIDCombined->GetPriors(KpTrack, oldpP, PIDResponse, detUsedP);

    nSigmaTPC = PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kProton);
    fPrTPC[icentr]->Fill(fPtKp,nSigmaTPC);
    nSigmaTPC = PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kElectron);
    fElTPC[icentr]->Fill(fPtKp,nSigmaTPC);
    nSigmaTPC = PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kPion);
    fPiTPC[icentr]->Fill(fPtKp,nSigmaTPC);
    nSigmaTPC = PIDResponse->NumberOfSigmasTPC(KpTrack,AliPID::kKaon);
    fKaTPC[icentr]->Fill(fPtKp,nSigmaTPC);

    if(! (TMath::Abs(nSigmaTPC) < 5)) continue;

    Int_t tofMatch1 = (KpTrack->GetStatus() & AliVTrack::kTOFout) && (KpTrack->GetStatus() & AliVTrack::kTIME);
    /*
    if(mcArray){
      Int_t labelK = TMath::Abs(KpTrack->GetLabel());
      AliAODMCParticle *mcp1 = (AliAODMCParticle*)mcArray->At(labelK);
      pdg1 = TMath::Abs(mcp1->GetPdgCode());
    }
    */

    fPidKp = Int_t(probP[3]*100);

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
	if(probP[3] > probP[2] && probP[3] > probP[4] && probP[3] > probP[0]) fPidKp += 128; // max prob
	
	nSigmaTOF = PIDResponse->NumberOfSigmasTOF(KpTrack,AliPID::kProton);
	fPrTOF[icentr]->Fill(fPtKp,nSigmaTOF);
	nSigmaTOF = PIDResponse->NumberOfSigmasTOF(KpTrack,AliPID::kElectron);
	fElTOF[icentr]->Fill(fPtKp,nSigmaTOF);
	nSigmaTOF = PIDResponse->NumberOfSigmasTOF(KpTrack,AliPID::kPion);
	fPiTOF[icentr]->Fill(fPtKp,nSigmaTOF);
	nSigmaTOF = PIDResponse->NumberOfSigmasTOF(KpTrack,AliPID::kKaon);
	fKaTOF[icentr]->Fill(fPtKp,nSigmaTOF);
		
	nSigmaTOF = TMath::Abs(nSigmaTOF);
	
	if(fIsMC){
	  Float_t mismAdd = addMismatchForMC;
	  if(KpTrack->Pt() < 1) mismAdd = addMismatchForMC/KpTrack->Pt();
	  
	  if(gRandom->Rndm() < mismAdd){
	    Float_t etaAbs = TMath::Abs(KpTrack->Eta());
	    Int_t channel = Int_t(4334.09 - 4758.36 * etaAbs -1989.71 * etaAbs*etaAbs + 1957.62*etaAbs*etaAbs*etaAbs);
	    channel = channel % 8736;
	    Float_t distIP = fHchannelTOFdistr->GetBinContent(channel);
	    
	    // generate random time
	    Float_t timeRandom = fHmismTOF->GetRandom() + distIP*3.35655419905265973e+00;
	    Double_t times[10];
	    KpTrack->GetIntegratedTimes(times);
	    nSigmaTOF = TMath::Abs(timeRandom - times[3])/PIDResponse->GetTOFResponse().GetExpectedSigma(KpTrack->P(), times[3], AliPID::kKaon);
	  }
	}

	if(fCentrality < 20 && KpTrack->Pt() < 0.9 && KpTrack->Pt() > 0.8)fTOFTPCsignal->Fill(nSigmaTOF,nSigmaTPC);

	if(nSigmaTOF < 2) fPidKp += 256;
	else if(nSigmaTOF < 3) fPidKp += 512;
      }
    }
    
    if(tofMatch1){
      nSigmaComb = TMath::Sqrt(nSigmaTOF*nSigmaTOF + nSigmaTPC*nSigmaTPC);
      if(nSigmaTOF < 5 && fCentrality < 20 && KpTrack->Pt() < 0.9 && KpTrack->Pt() > 0.8){
	fCombsignal->Fill(nSigmaComb);
      }
    } else {
      nSigmaComb = TMath::Abs(nSigmaTPC);
    }

    // use sigmaTOF instead of sigmaComb
    //nSigmaComb = nSigmaTOF;

    if(nSigmaComb < 2) nSigmaComb = 2;
    else if(nSigmaComb < 3) nSigmaComb = 3;
    else if(nSigmaComb < 5) nSigmaComb = 4.99;
    else nSigmaComb = 6;

    if(fPtKp > 4.299) fPtKp = 4.299;

    for(Int_t j=0;j < ntrack;j++){ // loop on negative tracks
      AliAODTrack* KnTrack = aodEvent->GetTrack(j);
      
      if (!KnTrack){
	continue;
      }

      if(!(KnTrack->Charge() < 0 && KnTrack->Pt() > 0.3 && TMath::Abs(KnTrack->Eta()) < 0.8)) continue;

      px = KpTrack->Px() + KnTrack->Px();
      py = KpTrack->Py() + KnTrack->Py();
      pz = KpTrack->Pz() + KnTrack->Pz();
      E = TMath::Sqrt(KpTrack->P()*KpTrack->P() + 4.93676999999999977e-01*4.93676999999999977e-01);
      E += TMath::Sqrt(KnTrack->P()*KnTrack->P()+ 4.93676999999999977e-01*4.93676999999999977e-01);

      phiV.SetPxPyPzE(px,py,pz,E);
      fMassV0 = phiV.M();
      
      if(fMassV0 <  invMmin || fMassV0 > invMmax) continue;

      fPtKn=KnTrack->Pt(),fPhiKn=KnTrack->Phi(),fEtaKn=KnTrack->Eta();
      fPidKn=0;

      UInt_t detUsedN = fPIDCombined->ComputeProbabilities(KnTrack, PIDResponse, probN);
      Double_t oldpN[10];
      fPIDCombined->GetPriors(KnTrack, oldpN, PIDResponse, detUsedN);

      nSigmaTPC2 = PIDResponse->NumberOfSigmasTPC(KnTrack,AliPID::kKaon);
      
      if(! (TMath::Abs(nSigmaTPC2) < 5)) continue;
      
      nSigmaComb2=5;

      Int_t tofMatch2 = (KnTrack->GetStatus() & AliVTrack::kTOFout) && (KnTrack->GetStatus() & AliVTrack::kTIME);
      /*
      if(mcArray){
	Int_t labelK = TMath::Abs(KnTrack->GetLabel());
	AliAODMCParticle *mcp2 = (AliAODMCParticle*)mcArray->At(labelK);
 	pdg2 = TMath::Abs(mcp2->GetPdgCode());
      }
      */

      fPidKn = Int_t(probN[3]*100);

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
	  if(probN[3] > probN[2] && probN[3] > probN[4] && probN[3] > probN[0]) fPidKn += 128; // max prob
	  
	  nSigmaTOF2 = PIDResponse->NumberOfSigmasTOF(KnTrack,AliPID::kKaon);
	  	  
	  nSigmaTOF2 = TMath::Abs(nSigmaTOF2);
	  
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
	      Double_t times[10];
	      KnTrack->GetIntegratedTimes(times);
	      nSigmaTOF = TMath::Abs(timeRandom - times[3])/PIDResponse->GetTOFResponse().GetExpectedSigma(KnTrack->P(), times[3], AliPID::kKaon);
	    }
	  }

	  if(fCentrality < 20 && KnTrack->Pt() < 1.2 && KnTrack->Pt() > 1) fTOFTPCsignal->Fill(nSigmaTOF2,nSigmaTPC2);

	  if(nSigmaTOF2 < 2) fPidKn += 256;
	  else if(nSigmaTOF2 < 3) fPidKn += 512;
	}
      }

      fPtPhi = phiV.Pt();
      fEtaPhi = phiV.Eta();
      fPhiPhi = phiV.Phi();

      if(tofMatch2){
	nSigmaComb2 = TMath::Sqrt(nSigmaTOF2*nSigmaTOF2+ nSigmaTPC2*nSigmaTPC2);
	if(nSigmaTOF2 < 5 && fCentrality < 20 && KnTrack->Pt() < 1.2 && KnTrack->Pt() > 1){
	  fCombsignal->Fill(nSigmaComb2);
	}
      } else {
	nSigmaComb2 = TMath::Abs(nSigmaTPC2);
      }

      // use sigmaTOF instead of sigmaComb
      //nSigmaComb2 = nSigmaTOF2;

      if(nSigmaComb2 < 2) nSigmaComb2 = 2;
      else if(nSigmaComb2 < 3) nSigmaComb2 = 3;
      else if(nSigmaComb2 < 5) nSigmaComb2 = 4.99;
      else nSigmaComb2 = 6;  

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
	    if(pdgM == 333) isTrue = kTRUE;
	  }
	}

      }

      Double_t deltaphi1 = KpTrack->Phi() - fPsi;
      Double_t deltaphi2 = KnTrack->Phi() - fPsi;

      while(deltaphi1 > TMath::Pi()) deltaphi1 -= TMath::Pi()*2;
      while(deltaphi1 < -TMath::Pi()) deltaphi1 += TMath::Pi()*2;
      while(deltaphi2 > TMath::Pi()) deltaphi2 -= TMath::Pi()*2;
      while(deltaphi2 < -TMath::Pi()) deltaphi2 += TMath::Pi()*2;

      if(fPtKn > 4.299) fPtKn = 4.299;

      Float_t xTOfill[] = {fPtPhi,KpTrack->Eta(),fPtKp,fPtKn,(fPidKp%128)*0.01,(fPidKn%128)*0.01,tofMatch1,tofMatch2,isTrue,nSigmaComb,nSigmaComb2,deltaphi1,deltaphi2,fPsi};
      

      Int_t ipt = 0;
      while(fPtPhiMin[ipt] < fPtPhi && ipt < nPtBin){
	ipt++;
      }
      ipt--;
      if(ipt < 0) ipt = 0; // just to be sure

      if(TMath::Abs(fEtaPhi) < 0.8 && fPtKp > 0.3 && fPtKn > 0.3){
	if((fPidKn%128) > 80) fContPid->Fill(0,fMassV0,fCentrality,xTOfill); // use tagging on negative track
        xTOfill[1] = KnTrack->Eta();
	if((fPidKp%128) > 80) fContPid2->Fill(0,fMassV0,fCentrality,xTOfill);// use tagging on positive track

      }


    }
  } // end analysi phi

}

//_____________________________________________________________________________
Float_t AliAnalysisTaskPhiBayes::GetVertex(AliAODEvent* aod) const
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
void AliAnalysisTaskPhiBayes::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
//-------------------------------------------------------------------------
Int_t AliAnalysisTaskPhiBayes::IsChannelValid(Float_t etaAbs){
  if(!fIsMC) return 1; // used only on MC

  if(fTypeCol == 2){ // LHC10h or LHC11h because of TOF matching window at 3 cm
    Int_t channel = Int_t(4334.09 - 4758.36 * etaAbs -1989.71 * etaAbs*etaAbs + 1957.62*etaAbs*etaAbs*etaAbs); 
  
    if(!(channel%20)) return 0; // 5% additional loss in MC
  }

  return 1;
}
