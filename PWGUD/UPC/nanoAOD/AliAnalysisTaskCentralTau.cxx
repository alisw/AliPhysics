/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TClonesArray.h"
#include "TColor.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1I.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliDataFile.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDZDC.h"
#include "AliInputEventHandler.h"
#include "AliOADBContainer.h"
#include "AliPIDResponse.h"
#include "AliTimeRangeCut.h"
#include "AliTOFTriggerMask.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"

// my headers
#include "AliAnalysisTaskCentralTau.h"

ClassImp(AliAnalysisTaskCentralTau);

using std::cout;
using std::endl;

//analysis skeleton of UPC nano AODs,

//_____________________________________________________________________________
AliAnalysisTaskCentralTau::AliAnalysisTaskCentralTau()
  : AliAnalysisTaskSE(),
    fPIDResponse(0), fTrackCutsBit0(0), fTrackCutsBit1(0), fTrackCutsBit4(0), cutEta(0.9), fOutputList(0),tTwoTracks(0), hTriggerCounter(0), hParticleTypeCounter(0), fPt(0), fY(0), fM(0), fPhi(0),
    fZNAenergy(0), fZNCenergy(0), fChannel(0), fSign(0), fRunNumber(0), fADAdecision(0), fADCdecision(0), fV0Adecision(0), fV0Cdecision(0)
{

//Dummy constructor

}//AliAnalysisTaskCentralTau


//_____________________________________________________________________________
AliAnalysisTaskCentralTau::AliAnalysisTaskCentralTau(const char *name)
  : AliAnalysisTaskSE(name),
    fPIDResponse(0), fTrackCutsBit0(0), fTrackCutsBit1(0), fTrackCutsBit4(0), cutEta(0.9), fOutputList(0), tTwoTracks(0), hTriggerCounter(0), hParticleTypeCounter(0), fPt(0), fY(0), fM(0), fPhi(0),
    fZNAenergy(0), fZNCenergy(0), fChannel(0), fSign(0), fRunNumber(0), fADAdecision(0), fADCdecision(0), fV0Adecision(0), fV0Cdecision(0)
{
  for(Int_t i = 0; i<(sizeof(fTriggers)/sizeof(fTriggers[0])); i++)          fTriggers[i] = kFALSE;
  for(Int_t i = 0; i<(sizeof(fTriggerClass)/sizeof(fTriggerClass[0]));  i++) fTriggerClass[i] = kFALSE;
  DefineOutput(1, TList::Class());

}//AliAnalysisTaskCentralTau

//_____________________________________________________________________________
AliAnalysisTaskCentralTau::~AliAnalysisTaskCentralTau()
{
  // Destructor
  
  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
     delete fOutputList;
     fOutputList = 0x0;
  }

}//~AliAnalysisTaskCentralTau


//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::UserCreateOutputObjects()
{

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
  fTrackCutsBit0 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fTrackCutsBit0->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  fTrackCutsBit1 = AliESDtrackCuts::GetStandardITSSATrackCuts2010(kFALSE,kTRUE);
  fTrackCutsBit4 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);

  fOutputList = new TList();
  fOutputList ->SetOwner(); // @suppress("Ambiguous problem")
  
  tTwoTracks = new TTree("tTwoTracks", "tTwoTracks");
  tTwoTracks ->Branch("fPt", &fPt, "fPt/F");
  tTwoTracks ->Branch("fPtDaughter", &fPtDaughter[0], "fPtDaughter[2]/F");
  tTwoTracks ->Branch("fVectElectron", &fVectDaughter[0]);
  tTwoTracks ->Branch("fVectOtherTrack", &fVectDaughter[1]);
  tTwoTracks ->Branch("fSignDaughter", &fSignDaughter[0], "fSignDaughter[2]/I");
  tTwoTracks ->Branch("fY", &fY, "fY/F");
  tTwoTracks ->Branch("fM", &fM, "fM/F");
  tTwoTracks ->Branch("fPhi", &fPhi, "fPhi/F");
  tTwoTracks ->Branch("fChannel", &fChannel, "fChannel/I");
  tTwoTracks ->Branch("fSign", &fSign, "fSign/I");
  tTwoTracks ->Branch("fZNAenergy", &fZNAenergy,"fZNAenergy/F");
  tTwoTracks ->Branch("fZNCenergy", &fZNCenergy,"fZNCenergy/F");
  tTwoTracks ->Branch("fZNAtime", &fZNAtime[0],"fZNAtime[4]/F");
  tTwoTracks ->Branch("fZNCtime", &fZNCtime[0],"fZNCtime[4]/F");
  tTwoTracks ->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
  tTwoTracks ->Branch("fTriggers", &fTriggers, Form("fTriggers[%i]/O",(sizeof(fTriggers)/sizeof(fTriggers[0]))));
  tTwoTracks ->Branch("fADAdecision", &fADAdecision, "fADAdecision/I");
  tTwoTracks ->Branch("fADCdecision", &fADCdecision, "fADCdecision/I");
  tTwoTracks ->Branch("fV0Adecision", &fV0Adecision, "fV0Adecision/I");
  tTwoTracks ->Branch("fV0Cdecision", &fV0Cdecision, "fV0Cdecision/I");
  tTwoTracks ->Branch("fPIDpt", &fPIDpt[0], "fPIDpt[2]/F");
  tTwoTracks ->Branch("fTPCsignal", &fTPCsignal[0], "fTPCsignal[2]/F");
  tTwoTracks ->Branch("fTOFsignal", &fTOFsignal[0], "fTOFsignal[2]/F");
  tTwoTracks ->Branch("fTPCmostProbableTrackID", &fTPCmostProbableTrackType[0], "fTPCmostProbableTrackID/I");
  tTwoTracks ->Branch("fTOFmostProbableTrackID", &fTOFmostProbableTrackType[0], "fTOFmostProbableTrackID/I");
  fOutputList->Add(tTwoTracks);

  hTriggerCounter = new TH2I("hTriggerCounter","Number of analyzed UPC triggers per run",3,1,4,3000,295000,298000);
  fOutputList->Add(hTriggerCounter);
  hParticleTypeCounter = new TH1I("hParticleTypeCounter","Electron, Muon, Pion",4,-0.5,3.5);
  fOutputList->Add(hParticleTypeCounter);
     
  PostData(1, fOutputList);

}//UserCreateOutputObjects


//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::UserExec(Option_t *)
{
  //
  // META
  //
  AliVEvent *fEvent = InputEvent();
  if(!fEvent) return;

  fRunNumber = fEvent->GetRunNumber();

  AliTimeRangeCut fTimeRangeCut;
  fTimeRangeCut.InitFromEvent(InputEvent());
  if(fTimeRangeCut.CutEvent(InputEvent()))return;

  //
  // TRIGGER
  //
  TString trigger = fEvent->GetFiredTriggerClasses();
  
  if(fRunNumber>=295881 && fRunNumber<296594)	if(!trigger.Contains("CCUP29-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP30-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP31-B-SPD2-CENTNOTRD"))return;
  if(fRunNumber>=296594) 					          	if(!trigger.Contains("CCUP29-U-SPD2-CENTNOTRD") && !trigger.Contains("CCUP30-B-SPD2-CENTNOTRD") && !trigger.Contains("CCUP31-B-SPD2-CENTNOTRD"))return;
  if(fRunNumber< 295881)  						        if(!trigger.Contains("CCUP29-B-NOPF-CENTNOTRD") && !trigger.Contains("CCUP30-B-NOPF-CENTNOTRD") && !trigger.Contains("CCUP31-B-NOPF-CENTNOTRD"))return;


  UInt_t fL0inputs = fEvent->GetHeader()->GetL0TriggerInputs();
  fTriggers[0] = trigger.Contains("CCUP29-B-SPD2-CENTNOTRD");
  fTriggers[1] = trigger.Contains("CCUP29-B-NOPF-CENTNOTRD");
  fTriggers[2] = trigger.Contains("CCUP29-U-SPD2-CENTNOTRD");
  fTriggers[3] = trigger.Contains("CCUP30-B-NOPF-CENTNOTRD");
  fTriggers[4] = trigger.Contains("CCUP30-B-SPD2-CENTNOTRD");
  fTriggers[5] = trigger.Contains("CCUP31-B-NOPF-CENTNOTRD");
  fTriggers[6] = trigger.Contains("CCUP31-B-SPD2-CENTNOTRD");
  fTriggers[7] =  fL0inputs & (1 << 11);//OM2
  fTriggers[8] =  fL0inputs & (1 << 12);//OMU
  
  if(trigger.Contains("CCUP29-B") || trigger.Contains("CCUP29-U"))  hTriggerCounter->Fill(1,fRunNumber);
  if(trigger.Contains("CCUP30-B"))                                  hTriggerCounter->Fill(2,fRunNumber);
  if(trigger.Contains("CCUP31-B"))                                  hTriggerCounter->Fill(3,fRunNumber);

  //
  // VERTEX CUT
  //
  const AliVVertex *fVertex = fEvent->GetPrimaryVertex();
  if(fVertex->GetNContributors()<2) return;
  if(TMath::Abs(fVertex->GetZ())>15)return;
  
  //
  // OFFLINE VETOS INFORMATION
  //
  AliVVZERO *fV0data = fEvent->GetVZEROData();
  AliVAD *fADdata = fEvent->GetADData();
  
  AliESDZDC *fZDCdata = (AliESDZDC*)fEvent->GetZDCData();
  fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
  fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
  Int_t detChZNA  = fZDCdata->GetZNATDCChannel();
  Int_t detChZNC  = fZDCdata->GetZNCTDCChannel();
  if (fEvent->GetRunNumber()>=245726 && fEvent->GetRunNumber()<=245793) detChZNA = 10;
  for (Int_t i=0;i<4;i++){
    fZNAtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNA,i);
    fZNCtime[i] = fZDCdata->GetZDCTDCCorrected(detChZNC,i);
  }

  fV0Adecision = fV0data->GetV0ADecision();
  fV0Cdecision = fV0data->GetV0CDecision();
  
  fADAdecision = fADdata->GetADADecision();
  fADCdecision = fADdata->GetADCDecision();

  //
  // TRACKS INFO
  //
  UInt_t nGoodTracksTPC=0;
  UInt_t nGoodTracksSPD=0;
  UInt_t nGoodTracksTOF=0;
  Int_t TrackIndexTPC[5] = {-1,-1,-1,-1,-1};

  //
  // LOOP OVER TRACKS
  //
  if (fEvent->GetNumberOfTracks() < 1) return;
  for(Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
    Bool_t goodTPCTrack = kTRUE;
    Bool_t goodITSTrack = kTRUE;
    AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(iTrack));
    if( !trk ) continue;
    if(!fTrackCutsBit1->AcceptTrack(trk)) goodITSTrack = kFALSE;
    if(!fTrackCutsBit4->AcceptTrack(trk)) goodTPCTrack = kFALSE;
    else{ if(trk->HasPointOnITSLayer(0) && trk->HasPointOnITSLayer(1))  nGoodTracksSPD++;}
    if(goodTPCTrack){
    	TrackIndexTPC[nGoodTracksTPC] = iTrack;
    	nGoodTracksTPC++;
    }
    if (nGoodTracksTPC == 5) {
      Printf("***** WARNING: Event has 5 TPC good tracks. Loop over tracks terminated.");
      break;
    }
  }//Track loop

  //
  // SETUP STG
  //
  Int_t crossedFO[4];
  TBits fFOCrossedChips(1200);
  const AliVMultiplicity *mult = fEvent->GetMultiplicity();
  TBits fFOFiredChips = mult->GetFastOrFiredChips();

  //
  // TRACKS SETTING
  //
  TDatabasePDG *pdgdat = TDatabasePDG::Instance();
  Float_t muonMass = pdgdat->GetParticle( 13 )->Mass();
  Float_t electronMass = pdgdat->GetParticle( 11 )->Mass();
  Float_t pionMass = pdgdat->GetParticle( 211 )->Mass();

  Short_t qTrack[2];
  TLorentzVector vMuon[2], vElectron[2], vPion[2];
  TLorentzVector vDiTauCandidate;

  Float_t nSigmaMuon[2], nSigmaElectron[2], nSigmaPion[2], dEdx[2];

  //
  // SELECT TWO TRACKS
  //
  if(nGoodTracksTPC == 2 && nGoodTracksSPD == 2){
  	fFOCrossedChips.ResetAllBits(kFALSE);

    //
    // LOOP OVER TWO TRACKS
    //
  	Bool_t isOneTrackElectron = kFALSE;
  	Int_t electronTrackID = -1;
  	for(Int_t iTrack=0; iTrack<2; iTrack++) {

  	  //
  	  // Fast-OR chips crossed STG
  	  //
      AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(TrackIndexTPC[iTrack]));
      if( !trk ) continue;
      crossedFO[0] = trk->GetITSModuleIndex(0);
      crossedFO[1] = trk->GetITSModuleIndex(1);
      crossedFO[2] = trk->GetITSModuleIndex(6);
      crossedFO[3] = trk->GetITSModuleIndex(7);
      SetCrossed(crossedFO, fFOCrossedChips);

  	  //
  	  // PID AND KINEMATIC INFO
  	  //

      fPtDaughter[iTrack] = trk->Pt();
      fSignDaughter[iTrack] = trk->Charge();

      Float_t PIDTPCElectron = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
      Float_t PIDTPCMuon = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
      Float_t PIDTPCPion = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
      Float_t PIDTPCKaon = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
      Float_t PIDTPCProton = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);
      Float_t PIDTOFElectron = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
      Float_t PIDTOFMuon = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
      Float_t PIDTOFPion = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
      Float_t PIDTOFKaon = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
      Float_t PIDTOFProton = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);

      qTrack[iTrack] = trk->Charge();

      vElectron[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), electronMass);
      vMuon[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), muonMass);
      vPion[iTrack].SetPtEtaPhiM(trk->Pt(), trk->Eta(), trk->Phi(), pionMass);

      nSigmaElectron[iTrack] = PIDTPCElectron;
      nSigmaMuon[iTrack] = PIDTPCMuon;
      nSigmaPion[iTrack] = PIDTPCPion;

      dEdx[iTrack] = trk->GetTPCsignal();

      Int_t trackPIDid = TestPIDTPChypothesis(nSigmaElectron[iTrack], nSigmaMuon[iTrack], nSigmaPion[iTrack]);
      hParticleTypeCounter->Fill(trackPIDid);
      if (trackPIDid == 1) {
        isOneTrackElectron = kTRUE;
        electronTrackID = iTrack;
      }
    }

  	if (!isOneTrackElectron) return;

    //
    // SAVE STG DECISION
    //
    TBits fFOCrossFiredChips = fFOCrossedChips & fFOFiredChips;
    fTriggers[9] = IsSTGFired(fFOCrossFiredChips,fRunNumber >= 295753 ? 9 : 3);

    //
    // CHARGE? (OPPOSITE, LIKE)
    //
    if(qTrack[0]*qTrack[1]<0)fSign = -1;
    if(qTrack[0]*qTrack[1]>0)fSign = 1;

    //
    // SAVE KINEMATICS
    //
    Int_t otherTrackID = -1;
    if (electronTrackID == 0) otherTrackID = 1;
    else otherTrackID = 0;
    if(TestPIDTPChypothesis(nSigmaElectron[otherTrackID], nSigmaMuon[otherTrackID], nSigmaPion[otherTrackID]) == 1) {
      vDiTauCandidate = vElectron[electronTrackID]+vElectron[otherTrackID];
      fVectDaughter[0] = vElectron[electronTrackID];
      fVectDaughter[1] = vElectron[otherTrackID];
      fChannel = 1;
      FillTree(tTwoTracks,vDiTauCandidate);
    }
    if(TestPIDTPChypothesis(nSigmaElectron[otherTrackID], nSigmaMuon[otherTrackID], nSigmaPion[otherTrackID]) == 2) {
      vDiTauCandidate = vElectron[electronTrackID]+vMuon[otherTrackID];
      fVectDaughter[0] = vElectron[electronTrackID];
      fVectDaughter[1] = vMuon[otherTrackID];
      fChannel = 2;
      FillTree(tTwoTracks,vDiTauCandidate);
    }
    if(TestPIDTPChypothesis(nSigmaElectron[otherTrackID], nSigmaMuon[otherTrackID], nSigmaPion[otherTrackID]) == 3) {
      vDiTauCandidate = vElectron[electronTrackID]+vPion[otherTrackID];
      fVectDaughter[0] = vElectron[electronTrackID];
      fVectDaughter[1] = vPion[otherTrackID];
      fChannel = 3;
      FillTree(tTwoTracks,vDiTauCandidate);
    }
  }//Two track loop
  
  Bool_t isCUP29, isCUP30, isCUP31;
  isCUP29 = fTriggers[0] || fTriggers[1] || fTriggers[2];
  isCUP30 = fTriggers[3] || fTriggers[4];
  isCUP31 = fTriggers[5] || fTriggers[6];
  fTriggerClass[0] = isCUP29; fTriggerClass[1] = isCUP30; fTriggerClass[2] = isCUP31;

  PostData(1, fOutputList);

}//UserExec

//_____________________________________________________________________________
Int_t AliAnalysisTaskCentralTau::TestPIDTPChypothesis(Float_t e, Float_t m, Float_t p){
  // Choose, which particle it is accroding to PID
  if      (e < m && e < p) return 1; // It is an electron
  else if (m < e && m < p) return 2; // It is a muon
  else if (p < e && p < m) return 3; // It is a pion
  else return 0;
}

//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::TPCsignalInfo(AliESDtrack *trk, Int_t trkID){


//  Float_tfTPCsignal[2],fTOFsignal[2];
//  Int_t fTPCmostProbableTrackType[2], fTOFmostProbableTrackType[2];

  fPIDpt[trkID] = trk->Pt();
  fTPCsignal[trkID] = trk->GetTPCsignal();
  fTOFsignal[trkID] = trk->GetTOFsignal();

  Float_t PIDTPC[5];
  PIDTPC[0] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kElectron);
  PIDTPC[1] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kMuon);
  PIDTPC[2] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kPion);
  PIDTPC[3] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kKaon);
  PIDTPC[4] = fPIDResponse->NumberOfSigmasTPC(trk,AliPID::kProton);

  fTPCmostProbableTrackType[trkID] = *std::min_element(PIDTPC, PIDTPC + 5);

  Float_t PIDTOF[5];
  PIDTOF[0] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kElectron);
  PIDTOF[1] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kMuon);
  PIDTOF[2] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kPion);
  PIDTOF[3] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kKaon);
  PIDTOF[4] = fPIDResponse->NumberOfSigmasTOF(trk,AliPID::kProton);

  fTOFmostProbableTrackType[trkID] = *std::min_element(PIDTOF, PIDTOF + 5);

}

//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::SetCrossed(Int_t spd[4], TBits &crossed){

  Int_t chipId2;
  for(Int_t iLayer = 0; iLayer<4 ;iLayer++)
    if(spd[iLayer]>0) { crossed.SetBitNumber(GetChipId(spd[iLayer],chipId2)); crossed.SetBitNumber(chipId2); }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskCentralTau::GetChipId(Int_t index, Int_t &chipId2, Bool_t debug){
  Int_t status   = (index%1000000)/100000;
  Int_t iModule  = index/1000000;           // 0 - 239
  Int_t iPhi     = iModule/4;               // 0-19 - inner, 20-59 outer
  Int_t iModuleZ = iModule%4;               // 0-3
  Int_t iSign    = (index%100000)/10000;    // 1-4
  Int_t signZ    = iPhi<20 ? (iSign%2==1 ? 1 : -1) : (iSign%2==0 ? 1 : -1); // 1 or -1
  Int_t iX       = (index%10000)/100;       // ??
  Int_t iZ       = index%100;               // 0-36 [mm]
  Int_t signZiZ  = (36-signZ*iZ);
  Int_t chipId   = iModule*5+signZiZ*5/72;
  if (chipId<0) return 1200;
  if (chipId>=1200) return 1201;
  if (signZiZ<0) return 1202;
  if (signZiZ>72) return 1203;
  if (signZiZ==72 && chipId%20==0 && chipId>=400) return 1204;
  chipId2=chipId;
  
  if (signZiZ==0  && chipId%20!=0)  chipId2=chipId-1;
  if (signZiZ==72 && chipId%20!=19) chipId2=chipId+1;
  if (signZiZ==13)  chipId2=chipId+1;
  if (signZiZ==14)  chipId2=chipId+1;
  if (signZiZ==15)  chipId2=chipId-1;
  if (signZiZ==16)  chipId2=chipId-1;
  if (signZiZ==27)  chipId2=chipId+1;
  if (signZiZ==28)  chipId2=chipId+1;
  if (signZiZ==29)  chipId2=chipId-1;
  if (signZiZ==30)  chipId2=chipId-1;
  if (signZiZ==42)  chipId2=chipId+1;
  if (signZiZ==43)  chipId2=chipId+1;
  if (signZiZ==44)  chipId2=chipId-1;
  if (signZiZ==45)  chipId2=chipId-1;
  if (signZiZ==56)  chipId2=chipId+1;
  if (signZiZ==57)  chipId2=chipId+1;
  if (signZiZ==58)  chipId2=chipId-1;
  if (signZiZ==59)  chipId2=chipId-1;
  if (debug) printf("%4i %4i %3i %3i %3i\n",chipId,chipId2,iX,signZiZ,iSign);
  return chipId;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCentralTau::IsSTGFired(TBits bits, Int_t dphiMin, Int_t dphiMax, Bool_t tolerance){
  Int_t n1 = bits.CountBits(400);
  Int_t n0 = bits.CountBits()-n1;
  //cout<<n0<<" "<<n1<<endl;
  if (n0<1 || n1<1) return 0;
  Bool_t stg = 0;
  Bool_t l0[20]={0};
  Bool_t l1[40]={0};
  Bool_t phi[20]={0};
  for (Int_t i=0;   i< 400; ++i) if (bits.TestBitNumber(i)) l0[      i/20] = 1;
  for (Int_t i=400; i<1200; ++i) if (bits.TestBitNumber(i)) l1[(i-400)/20] = 1;
  for (Int_t i=0; i<20; ++i) {
    if (tolerance) phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40] | l1[(2*i+2)%40] | l1[(2*i+39)%40]);
    else           phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40]);
  }
  for (Int_t dphi=dphiMin;dphi<=dphiMax;dphi++)
    for (Int_t i=0; i<20; ++i) stg |= phi[i] & phi[(i+dphi)%20];
  return stg;
}




//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::Terminate(Option_t *)
{
  cout<<"Analysis complete."<<endl;
}//Terminate

//_____________________________________________________________________________
void AliAnalysisTaskCentralTau::FillTree(TTree *t, TLorentzVector v) {

  fPt      = v.Pt();
  if(v.E() != v.Pz())fY = v.Rapidity();
  else fY = -999;
  fM       = v.M();
  fPhi     = v.Phi();

  t->Fill();

}
