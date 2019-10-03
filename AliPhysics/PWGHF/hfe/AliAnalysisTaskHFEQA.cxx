/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//
// QA task
// 
// Authors:
//   Raphaelle Bailhache <R.Bailhache@gsi.de>
//
#include "TROOT.h"
#include "TChain.h"
#include "TMath.h"
#include <TString.h>
#include <TBits.h>
#include <TH1F.h>

#include <TDirectory.h>
#include <TTreeStream.h>

#include "AliVEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliESDUtils.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODTrack.h"
#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"

#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliAnalysisTaskHFEQA.h"
#include "AliAODMCHeader.h"


ClassImp(AliAnalysisTaskHFEQA)

//____________________________________________________________________
AliAnalysisTaskHFEQA::AliAnalysisTaskHFEQA() :
  AliAnalysisTaskSE(),
  fListHist(0x0), 
  fAODAnalysis(kFALSE),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fHFECuts(0),
  fPIDTPConly(0),
  fPIDTRDonly(0),
  fPIDTOFTPC(0),
  fPIDTPCTRD(0),
  fPIDTPCEMCal(0),
  fPIDqaTRDonly(0),
  fPIDqaTOFTPC(0),
  fPIDqaTPCTRD(0),
  fPIDqaTPCEMCal(0),
  fCentralityEstimator("V0M"),
  fCollisionSystem(3),
  fNbEvent(0),
  fTPConly(0),
  fTOFTPC(0),
  fTPCTRD(0),
  fTPCEMCal(0),
  fTPConlydo(kFALSE),
  fTRDonlydo(kFALSE),
  fTOFTPCdo(kFALSE),
  fTPCTRDdo(kFALSE),
  fTPCEMCaldo(kFALSE)
{
  // Constructor
   
}
//______________________________________________________________________________
AliAnalysisTaskHFEQA:: AliAnalysisTaskHFEQA(const char *name) :
  AliAnalysisTaskSE(name),
  fListHist(0x0),
  fAODAnalysis(kFALSE),
  fAODMCHeader(NULL),
  fAODArrayMCInfo(NULL),
  fHFECuts(0),
  fPIDTPConly(0),
  fPIDTRDonly(0),
  fPIDTOFTPC(0),
  fPIDTPCTRD(0),
  fPIDTPCEMCal(0),
  fPIDqaTRDonly(0),
  fPIDqaTOFTPC(0),
  fPIDqaTPCTRD(0),
  fPIDqaTPCEMCal(0),
  fCentralityEstimator("V0M"),
  fCollisionSystem(3),
  fNbEvent(0),
  fTPConly(0),
  fTOFTPC(0),
  fTPCTRD(0),
  fTPCEMCal(0),
  fTPConlydo(kFALSE),
  fTRDonlydo(kFALSE),
  fTOFTPCdo(kFALSE),
  fTPCTRDdo(kFALSE),
  fTPCEMCaldo(kFALSE)
{
  //
  // named ctor
  //
 
  fPIDTPConly = new AliHFEpid("hfePidTPConly");
  fPIDTRDonly = new AliHFEpid("hfePidTRDonly");
  fPIDTOFTPC = new AliHFEpid("hfePidTOFTPC");
  fPIDTPCTRD = new AliHFEpid("hfePidTPCTRD");
  fPIDTPCEMCal = new AliHFEpid("hfePidTPCEMCal");

  fPIDqaTRDonly = new AliHFEpidQAmanager;
  fPIDqaTOFTPC = new AliHFEpidQAmanager;
  fPIDqaTPCTRD = new AliHFEpidQAmanager;
  fPIDqaTPCEMCal = new AliHFEpidQAmanager;

  SetPbPbAnalysis();

  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
    
}
//____________________________________________________________
AliAnalysisTaskHFEQA::AliAnalysisTaskHFEQA(const AliAnalysisTaskHFEQA &ref):
  AliAnalysisTaskSE(ref),
  fListHist(NULL),
  fAODAnalysis(ref.fAODAnalysis), 
  fAODMCHeader(ref.fAODMCHeader),
  fAODArrayMCInfo(ref.fAODArrayMCInfo),
  fHFECuts(NULL),
  fPIDTPConly(0),
  fPIDTRDonly(0),
  fPIDTOFTPC(0),
  fPIDTPCTRD(0),
  fPIDTPCEMCal(0),
  fPIDqaTRDonly(0),
  fPIDqaTOFTPC(0),
  fPIDqaTPCTRD(0),
  fPIDqaTPCEMCal(0),
  fCentralityEstimator(ref.fCentralityEstimator),
  fCollisionSystem(ref.fCollisionSystem),
  fNbEvent(ref.fNbEvent),
  fTPConly(ref.fTPConly),
  fTOFTPC(ref.fTOFTPC),
  fTPCTRD(ref.fTPCTRD),
  fTPCEMCal(ref.fTPCEMCal),
  fTPConlydo(ref.fTPConlydo),
  fTRDonlydo(ref.fTRDonlydo),
  fTOFTPCdo(ref.fTOFTPCdo),
  fTPCTRDdo(ref.fTPCTRDdo),
  fTPCEMCaldo(ref.fTPCEMCaldo)
{
  //
  // Copy Constructor
  //

  ref.Copy(*this);
}

//____________________________________________________________
AliAnalysisTaskHFEQA &AliAnalysisTaskHFEQA::operator=(const AliAnalysisTaskHFEQA &ref){
  //
  // Assignment operator
  //
  if(this == &ref) 
    ref.Copy(*this);
  return *this;
}

//____________________________________________________________
void AliAnalysisTaskHFEQA::Copy(TObject &o) const {
  // 
  // Copy into object o
  //
  AliAnalysisTaskHFEQA &target = dynamic_cast<AliAnalysisTaskHFEQA &>(o);
  target.fListHist = fListHist;
  target.fAODAnalysis = fAODAnalysis;
  target.fAODMCHeader = fAODMCHeader;
  target.fAODArrayMCInfo = fAODArrayMCInfo;
  target.fHFECuts = fHFECuts;
  target.fPIDTPConly = fPIDTPConly;
  target.fPIDTRDonly = fPIDTRDonly;
  target.fPIDTOFTPC = fPIDTOFTPC;
  target.fPIDTPCTRD = fPIDTPCTRD;
  target.fPIDTPCEMCal = fPIDTPCEMCal;
  target.fPIDqaTRDonly = fPIDqaTRDonly;
  target.fPIDqaTOFTPC = fPIDqaTOFTPC;
  target.fPIDqaTPCTRD = fPIDqaTPCTRD;
  target.fPIDqaTPCEMCal = fPIDqaTPCEMCal;
  target.fCentralityEstimator = fCentralityEstimator;
  target.fCollisionSystem = fCollisionSystem;
  target.fNbEvent = fNbEvent;
  target.fTPConly = fTPConly;
  target.fTOFTPC = fTOFTPC;
  target.fTPCTRD = fTPCTRD;
  target.fTPCEMCal = fTPCEMCal;
  target.fTPConlydo = fTPConlydo;
  target.fTRDonlydo = fTRDonlydo;  
  target.fTOFTPCdo = fTOFTPCdo;
  target.fTPCTRDdo = fTPCTRDdo;
  target.fTPCEMCaldo = fTPCEMCaldo;
  	 

}
//____________________________________________________________
AliAnalysisTaskHFEQA::~AliAnalysisTaskHFEQA(){
  //
  // Destructor
  //
  

  if(fListHist) delete fListHist;
  if(fHFECuts) delete fHFECuts;
  if(fPIDTPConly) delete fPIDTPConly;
  if(fPIDTRDonly) delete fPIDTRDonly;
  if(fPIDTOFTPC) delete fPIDTOFTPC;
  if(fPIDTPCTRD) delete fPIDTPCTRD;
  if(fPIDTPCEMCal) delete fPIDTPCEMCal;

   
}
//________________________________________________________________________
void AliAnalysisTaskHFEQA::UserCreateOutputObjects()
{

  //********************
  // Create histograms
  //********************
  AliDebug(2,"AliAnalysisTaskHFEQA: User create output objects");


  // AOD or ESD
  AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
    SetAODAnalysis(kTRUE);
    AliDebug(2,"Put AOD analysis on");
  } else {
    SetAODAnalysis(kFALSE);
  }

  AliDebug(2,"AliAnalysisTaskHFEQA: AOD ESD");

  // HFE cuts

  if(!fHFECuts){
    fHFECuts = new AliHFEcuts;
    fHFECuts->CreateStandardCuts();
  }
  fHFECuts->Initialize();
  if(fAODAnalysis) {
    fHFECuts->SetAOD();
  }  

  AliDebug(2,"AliAnalysisTaskHFEQA: HFE cuts");


  // PIDTPConly HFE
  if(!fPIDTPConly) {
    fPIDTPConly =new AliHFEpid("hfePidTPConly");
  }
  if(!fPIDTPConly->GetNumberOfPIDdetectors()) fPIDTPConly->AddDetector("TPC", 0);
  fPIDTPConly->InitializePID();
  fPIDTPConly->SortDetectors();

  // PIDTRDonly HFE
  if(!fPIDTRDonly) {
    fPIDTRDonly =new AliHFEpid("hfePidTRDonly");
  }
  if(!fPIDTRDonly->GetNumberOfPIDdetectors()) fPIDTRDonly->AddDetector("TRD", 0);
  fPIDTRDonly->InitializePID();
  fPIDqaTRDonly->Initialize(fPIDTRDonly);
  fPIDTRDonly->SortDetectors();

  // PIDTOFTPC HFE
  if(!fPIDTOFTPC) {
    fPIDTOFTPC =new AliHFEpid("hfePidTOFTPC");
  }
  if(!fPIDTOFTPC->GetNumberOfPIDdetectors()) {
    fPIDTOFTPC->AddDetector("TOF", 0);
    fPIDTOFTPC->AddDetector("TPC", 1);
  }
  fPIDTOFTPC->InitializePID();
  fPIDqaTOFTPC->Initialize(fPIDTOFTPC);
  fPIDTOFTPC->SortDetectors();


  // PIDTPCTRD HFE
  if(!fPIDTPCTRD) {
    fPIDTPCTRD =new AliHFEpid("hfePidTPCTRD");
  }
  if(!fPIDTPCTRD->GetNumberOfPIDdetectors()) {
    fPIDTPCTRD->AddDetector("TPC", 0);
    fPIDTPCTRD->AddDetector("TRD", 1);
  }
  fPIDTPCTRD->InitializePID();
  fPIDqaTPCTRD->Initialize(fPIDTPCTRD);
  fPIDTPCTRD->SortDetectors();

  // PIDTPCEMCal HFE
  if(!fPIDTPCEMCal) {
    fPIDTPCEMCal =new AliHFEpid("hfePidTPCEMCal");
  }
  if(!fPIDTPCEMCal->GetNumberOfPIDdetectors()) {
    fPIDTPCEMCal->AddDetector("TPC", 0);
    fPIDTPCEMCal->AddDetector("EMCal", 1);
  }
  fPIDTPCEMCal->InitializePID();
  fPIDqaTPCEMCal->Initialize(fPIDTPCEMCal);
  fPIDTPCEMCal->SortDetectors();
 
  // Histograms
  fNbEvent = new TH1F("NbEvent", "",11,0,11);
  fNbEvent->Sumw2();
  Double_t ptbinning[36] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};
  fTPConly = new TH1F("TPCOnly", "",35,&ptbinning[0]);
  fTPConly->Sumw2();
  fTOFTPC = new TH1F("TOFTPC", "",35,&ptbinning[0]);
  fTOFTPC->Sumw2();
  fTPCTRD = new TH1F("TPCTRD", "",35,&ptbinning[0]);
  fTPCTRD->Sumw2();
  fTPCEMCal = new TH1F("TPCEMCal", "",35,&ptbinning[0]);
  fTPCEMCal->Sumw2();


  // List
    
  fListHist = new TList();
  fListHist->SetOwner();

  fListHist->Add(fPIDqaTRDonly->MakeList("HFEpidQATRDonly"));
  fListHist->Add(fPIDqaTOFTPC->MakeList("HFEpidQATOFTPC"));
  fListHist->Add(fPIDqaTPCTRD->MakeList("HFEpidQATPCTRD"));
  fListHist->Add(fPIDqaTPCEMCal->MakeList("HFEpidQATPCEMCal"));

  fListHist->Add(fNbEvent);
  fListHist->Add(fTPConly);
  fListHist->Add(fTOFTPC);
  fListHist->Add(fTPCTRD);
  fListHist->Add(fTPCEMCal);

  AliDebug(2,"AliAnalysisTaskHFEQA: list");


  fListHist->Print();

  PostData(1, fListHist);

  AliDebug(2,"AliAnalysisTaskHFEQA: post");


}
   
//________________________________________________________________________
void AliAnalysisTaskHFEQA::UserExec(Option_t */*option*/)
{
  //
  // Loop over event
  //
   
  Double_t binct = 11.5;

  AliMCEvent *mcEvent = MCEvent();
 

  AliDebug(2,"MC info");
  // MC info
  Bool_t mcthere = kTRUE;
  if(fAODAnalysis) {
    AliAODEvent *aodE = dynamic_cast<AliAODEvent *>(fInputEvent);
    if(!aodE){
      AliError("No AOD Event");
      return;
    }
    fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fInputEvent->FindListObject(AliAODMCHeader::StdBranchName()));
    if(!fAODMCHeader){ 
      mcthere = kFALSE;
    }
    fAODArrayMCInfo = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fAODArrayMCInfo){ 
      mcthere = kFALSE;
    }
    else {
      fHFECuts->SetMCEvent(aodE);
    }
  }
  else {
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH){
      mcthere=kFALSE;
    }
    if(mcEvent) fHFECuts->SetMCEvent(mcEvent);
  }


  ////////////////////////////////////
  // Number of contributors
  ///////////////////////////////////
  AliDebug(2,"Number of contributors");
  Int_t ncontribVtx = 0;
  if(fAODAnalysis) {
    AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
    if(!fAOD){
      AliError("AOD Event required for AOD Analysis");
      return;
    }  
    AliAODVertex *priVtx = fAOD->GetPrimaryVertex();
    if(priVtx){
      ncontribVtx = priVtx->GetNContributors();
    }
  }
  else {
    AliESDEvent *fESD = dynamic_cast<AliESDEvent *>(fInputEvent);
    if(!fESD){
      AliError("ESD Event required for ESD Analysis");
      return;
    }
    const AliESDVertex *priVtx = fESD->GetPrimaryVertexTracks();
    if(priVtx){
      ncontribVtx = priVtx->GetNContributors();
    }
  }
  AliDebug(2,Form("Number of contributors %d",ncontribVtx));


  /////////////////////////////////
  // centrality
  ////////////////////////////////

  //printf("Centrality \n");
  AliCentrality *centrality = fInputEvent->GetCentrality();
  AliDebug(2,"Got the centrality");
  Float_t cntr = 0.;
  if(centrality && (! Ispp())) { 
    cntr = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
    if((0.0< cntr) && (cntr<5.0)) binct = 0.5;
    if((5.0< cntr) && (cntr<10.0)) binct = 1.5;
    if((10.0< cntr) && (cntr<20.0)) binct = 2.5;
    if((20.0< cntr) && (cntr<30.0)) binct = 3.5;
    if((30.0< cntr) && (cntr<40.0)) binct = 4.5;
    if((40.0< cntr) && (cntr<50.0)) binct = 5.5;
    if((50.0< cntr) && (cntr<60.0)) binct = 6.5;
    if((60.0< cntr) && (cntr<70.0)) binct = 7.5;
    if((70.0< cntr) && (cntr<80.0)) binct = 8.5;
    if((80.0< cntr) && (cntr<90.0)) binct = 9.5;
    if((90.0< cntr) && (cntr<100.0)) binct = 10.5;
    if(binct > 11.0) return;
  }
  else binct = 0.5;
  AliDebug(2,Form("Centrality %f with %s",binct,fCentralityEstimator.Data()));
 
  //////////////////////
  // run number
  //////////////////////

  Int_t runnumber = fInputEvent->GetRunNumber();
  AliDebug(2,Form("Run number %d",runnumber));
   
  if(!fPIDTPConly->IsInitialized()){
    fPIDTPConly->InitializePID(runnumber);
  }
  if(!fPIDTRDonly->IsInitialized()){
    fPIDTRDonly->InitializePID(runnumber);
  }
  if(!fPIDTOFTPC->IsInitialized()){
    fPIDTOFTPC->InitializePID(runnumber);
  }
  if(!fPIDTPCTRD->IsInitialized()){
    fPIDTPCTRD->InitializePID(runnumber);
  }
  if(!fPIDTPCEMCal->IsInitialized()){
    fPIDTPCEMCal->InitializePID(runnumber);
  }

  //
  fHFECuts->SetRecEvent(fInputEvent);
  


  //////////
  // PID
  //////////
  AliDebug(2,"PID response");
  AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
  if(!pidResponse){
    AliDebug(2,"No PID response set");
    return;
  }
  fPIDTPConly->SetPIDResponse(pidResponse);
  fPIDTRDonly->SetPIDResponse(pidResponse);
  fPIDTOFTPC->SetPIDResponse(pidResponse);
  fPIDTPCTRD->SetPIDResponse(pidResponse);
  fPIDTPCEMCal->SetPIDResponse(pidResponse);
   
  //////////////////
  // Event cut
  //////////////////
  AliDebug(2,"Event cut");
  if(!fHFECuts->CheckEventCuts("fEvRecCuts", fInputEvent)) {
    AliDebug(2,"Does not pass the event cut");
    PostData(1, fListHist);
    return;
  }
  fNbEvent->Fill(binct);
  
  //////////////////////////
  // Loop over track
  //////////////////////////
  Int_t nbtracks = fInputEvent->GetNumberOfTracks();
  AliDebug(2,Form("Number of tracks %d",nbtracks));
  for(Int_t k = 0; k < nbtracks; k++){
      
    AliVTrack *track = (AliVTrack *) fInputEvent->GetTrack(k);
    if(!track) continue;
    Double_t pt = track->Pt();     

    AliDebug(2,"test 0\n");
    
    // RecKine: ITSTPC cuts  
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    AliDebug(2,"test 1\n");

    // RecPrim
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    AliDebug(2,"test 2\n");

    // HFEcuts: ITS layers cuts
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    AliDebug(2,"test 3\n");

    // HFE cuts: TOF and mismatch flag
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTOF + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    AliDebug(2,"test 4\n");

    // HFE cuts: TPC PID cleanup
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTPC + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    AliDebug(2,"test 5\n");

    // HFEcuts: Nb of tracklets TRD0
    if(!fHFECuts->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutStepsMCTrack, (TObject *)track)) continue;
    
    AliDebug(2,"Survived");
    

    ////////////////////////
    // Apply PID
    ////////////////////////
    AliHFEpidObject hfetrack;
    if(!fAODAnalysis) hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
    else hfetrack.SetAnalysisType(AliHFEpidObject::kAODanalysis);
    hfetrack.SetRecTrack(track);
    hfetrack.SetCentrality((Int_t)binct);
    hfetrack.SetMulitplicity(ncontribVtx); // for correction
    if(IsPbPb()) hfetrack.SetPbPb();
    else{
      if(IspPb()) hfetrack.SetpPb();
      else {
	hfetrack.SetPP();
	//printf("pp\n");
      }
    }
    AliDebug(2,Form("centrality %f and %d",binct,hfetrack.GetCentrality()));
   
    //printf("test 7\n");

    // Complete PID TPC alone
    if(fTPConlydo) {
      if(fPIDTPConly->IsSelected(&hfetrack,0x0,"recTrackCont",0x0)) {
	fTPConly->Fill(pt);
      }
    }
    AliDebug(2,"TPC only PID\n");
	
    // Complete PID TRD alone
    if(fTRDonlydo) {
      if(fPIDTRDonly->IsSelected(&hfetrack,0x0,"recTrackCont",fPIDqaTRDonly)) {
	AliDebug(2,"Passed TRD only PID\n");
      }
    }
    AliDebug(2,"TRD only PID\n");
	    
    
    // Complete PID TPC TOF 
    if(fTOFTPCdo) {
      if(fPIDTOFTPC->IsSelected(&hfetrack,0x0,"recTrackCont",fPIDqaTOFTPC)) {
	fTOFTPC->Fill(pt);
      }
    }
    AliDebug(2,"TOF TPC PID\n");
    
    // Complete PID TPC TRD 
    if(fTPCTRDdo) {
      if(fPIDTPCTRD->IsSelected(&hfetrack,0x0,"recTrackCont",fPIDqaTPCTRD)) {
	fTPCTRD->Fill(pt);
      }
    }
    AliDebug(2,"TPC TRD PID\n");


    if(fTPCEMCaldo) {
      if(!fAODAnalysis) {
	// Complete PID TPC TRD 
	if(fPIDTPCEMCal->IsSelected(&hfetrack,0x0,"recTrackCont",fPIDqaTPCEMCal)) {
	  fTPCEMCal->Fill(pt);
	}
      }
    }
    AliDebug(2,"TPC EMCal PID\n");
    
    
  }
  
  PostData(1, fListHist);
  
}
