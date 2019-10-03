/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ...                                                            *
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


//=========================================================================//
//             AliEbyE Analysis for Net Charge Fluctuation for RunII       //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//            Version 1:  PbPb Added with AliVTrack (30/11/2014)           //
//            Version 2:  pp Added with AliVTrack   (1/12/2014) - fixme    //            
//            Version 3:  pA Added with AliVTrack   (2/12/2014) - fixme    //    
//            Version 4:  Array Bug Fix | reduced THn-Bins 11/01/2015      //    
//            Version 5:  Introduting 3D mapping       (14/2/2015)         //    
//=========================================================================//

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "THnSparse.h"
#include "AliESD.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliAODHeader.h"
#include "AliVHeader.h"
#include "AliAODpidUtil.h"
#include "AliHelperPID.h"
#include "TClonesArray.h"

#include "TMath.h"
#include "TAxis.h"
#include "TSystem.h" 
#include "TFile.h" 
#include "TPRegexp.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliTracker.h"

#include "AliEbyENetChargeFluctuationTask.h"

ClassImp(AliEbyENetChargeFluctuationTask)

//-----------------------------------------------------------------------
AliEbyENetChargeFluctuationTask::AliEbyENetChargeFluctuationTask(const char *name) : 
  AliAnalysisTaskSE( name ), 

  fInputEventHandler(NULL),
  fESD(NULL),
  fAOD(NULL),
  fMCEvent(NULL),
  fStack(NULL), 
  fAODHandler(NULL),
  fESDHandler(NULL),
  fMCStack(NULL),
  fArrayMC(NULL),
  fESDtrackCuts(NULL),

  fQaList(NULL),
  fPhyList(NULL),
  fDcaList(NULL),
  fEffList(NULL),
 
  fSystemType(2), 
  fCentralityEstimator("V0M"), 

  fVxMax(3.), 
  fVyMax(3.), 
  fVzMax(10.), 
  fPhiMin(-1000.),
  fPhiMax(1000.),
  fPtMin(0.2), 
  fPtMax(3.), 
  fEtaMin(-1.), 
  fEtaMax(1.),  
  fRapMin(-0.5),
  fRapMax(0.5), 
  fDcaXy(10.),
  fDcaZ(10.),  
  fCentralityBin(-1.),    
  fCentralityPercentile(-1.),
// fNp(NULL),
// fMCNp(NULL),
// fRedFactp(NULL),           
   
  fMinTrackLengthMC(80), 
  fSelectBit(AliVEvent::kMB),
  fAODtrackCutBit(768),   
  fNSubSamples(10),     
  fSubSampleIdx(0),
  fOrder(8),
  fNTriggers(5),
  fHEventStatMax(8),
  fNCentralityBins(11),
  fCentralityBinMax(11),
  fNTracks(0),
  fNbwcBin(100),

  fIsMC(kFALSE),
  fIsRatio(kFALSE),
  fIsAOD(kFALSE),
  fIsSub(kFALSE),
  fIsBS(kFALSE),
  fIsPer(kFALSE),
  fIsEff(kFALSE),
  fDebug(kFALSE),
  fIsQa(kFALSE),
  fNeedQa(kFALSE),
  fIsPhy(kFALSE),
  fIsDca(kFALSE),
  fIsNu(kFALSE),
  fIsTen(kFALSE),
  fIs3D(kFALSE),

  fRan(0),             
  fRanIdx(0),                        
  fHelperPID(0x0) { 

  Printf("Task is Initialized");
  

  //-- -  - - -  - - ---- -- - --- -- --- -- - --
  for (Int_t ii = 0; ii <= fOrder; ++ii)
    for (Int_t jj = 0; jj < 2; ++jj)
      fRedFactp[ii][jj] = 1.;
  
  //-- -  - - -  - - ---- -- - --- -- --- -- - -- 
  for (Int_t kk = 0; kk < 4; ++kk)
    for (Int_t jj = 0; jj < 2; ++jj)
      fNp[kk][jj] = 0;
  
  //-- -  - - -  - - ---- -- - --- -- --- -- - --
  for (Int_t kk = 0; kk < 4; ++kk)
    for (Int_t jj = 0; jj < 2; ++jj)
      fMCNp[kk][jj] = 0;  
   

  DefineOutput(1, TList::Class()); 
  DefineOutput(2, TList::Class()); 
  DefineOutput(3, TList::Class()); 
  DefineOutput(4, TList::Class()); 
}

AliEbyENetChargeFluctuationTask::~AliEbyENetChargeFluctuationTask() {
  //!   Cleaning up
  if (fQaList)    delete fQaList;
  if (fPhyList)   delete fPhyList;
  if (fEffList)   delete fEffList;
  if (fDcaList)   delete fDcaList;
  return;
}


//________________ Static Variables _____________________
const Float_t fGBwRap     = 0.1;
const Float_t fGBwPt      = 0.1; 
const Float_t fGRngEta[]  = {-0.8, 0.8};
const Int_t   fGNBinsEta  = Int_t((fGRngEta[1] - fGRngEta[0])/fGBwRap);

const Float_t fGRngRap[]  = {-0.5, 0.5};
const Int_t   fGNBinsRap  = Int_t((fGRngRap[1] - fGRngRap[0])/fGBwRap);
const Float_t fGRngPhi[]  = {0.0, 6.3};
const Int_t   fGNBinsPhi  = 63;

const Float_t fGRngPt[]   = {0.2, 3.3};
const Int_t   fGNBinsPt   = Int_t((fGRngPt[1] - fGRngPt[0])/fGBwPt); 

const Int_t   fGNBinsSign =  2;
const Float_t fGRngSign[] = {-0.5, 1.5};


//---------------------------------------------------------------------------------
void AliEbyENetChargeFluctuationTask::UserCreateOutputObjects() {
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fQaList = new TList();
  fQaList->SetOwner(kTRUE);

  fPhyList = new TList();
  fPhyList->SetOwner(kTRUE);

  fEffList = new TList();
  fEffList->SetOwner(kTRUE);

  fDcaList = new TList();
  fDcaList->SetOwner(kTRUE);
  
  Printf(" >>>================================================================");
  if (!fIsAOD) {
    if(!fESDtrackCuts)
      fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    else 
      Printf(" >>>>  User Track Cuts <<<< ");
    fESDtrackCuts->Print();
    
    Printf(" >>> DCAxy in TC [%8.4f:%8.4f]", 
	   fESDtrackCuts->GetMinDCAToVertexXY(), fESDtrackCuts->GetMaxDCAToVertexXY());
    Printf(" >>> DCAz in TC  [%8.4f:%8.4f]", 
	   fESDtrackCuts->GetMinDCAToVertexZ(), fESDtrackCuts->GetMaxDCAToVertexZ());
	       
    Float_t r1,r2;
    fESDtrackCuts->GetPtRange(r1,r2);
    Printf(" >>> Pt in TC  [%10.4f:%10.4f]",r1,r2);

    fESDtrackCuts->GetRapRange(r1,r2);
    Printf(" >>> Rap in TC [%10.4f:%10.4f]",r1,r2);

    fESDtrackCuts->GetEtaRange(r1,r2);
    Printf(" >>> Eta in TC [%10.4f:%10.4f]",r1,r2);
  }     

  fRan = new TRandom3();
  fRan->SetSeed();
  
  fRanIdx = new TRandom3();
  fRanIdx->SetSeed();
 
  Printf(" >>> MC%d RA%d AO%d SS%d BS%d PE%d EF%d DB%d QA%d NQ%d PY%d DC%d NU%d TE%d", 
	 fIsMC,  fIsRatio, fIsAOD,  fIsSub, fIsBS, fIsPer, fIsEff, fDebug,   
	 fIsQa,  fNeedQa,  fIsPhy,  fIsDca, fIsNu, fIsTen);
  
  Printf(" >>> Centrality: %s System Type %d ",fCentralityEstimator.Data(),fSystemType); 
  Printf(" >>> Vx %.1f Vy %.1f Vz %.1f",fVxMax, fVyMax, fVzMax);
  Printf(" >>> Phi Range [%.2f:%.2f]",fPhiMin,fPhiMax);   
  Printf(" >>> Pt  Range [%.2f:%.2f]",fPtMin,fPtMax);   
  Printf(" >>> Eta Range [%.2f:%.2f]",fEtaMin,fEtaMax);   
  Printf(" >>> Rap Range [%.2f:%.2f]",fRapMin,fRapMax);   
  Printf(" >>> DCA Range [%.2f:%.2f]",fDcaXy,fDcaZ);   
  Printf(" >>> CBin:Cper [%.2f:%.2f]",fCentralityBin,fCentralityPercentile);   
  Printf(" >>> TLengthMC:%10.5f TrackN:%d TrackBit:%d Trigg:%d",
	 fMinTrackLengthMC, fNTracks, fAODtrackCutBit,fNTriggers);
  Printf(" >>> N.Samples:%d SsIdx:%d Order:%d",fNSubSamples,fSubSampleIdx,fOrder);
  Printf(" >>> EventStatBin:%d CentralityBin:%d CentBinMax:%d BWCbin:%d",
	 fHEventStatMax,fNCentralityBins,fCentralityBinMax, fNbwcBin);

  if (fNeedQa) CreateBasicQA();
  if (fIsQa && fNeedQa)  CreateQA();
  if (fIsPhy) InitPhy();
  if (fIsMC && fIsEff) CreateCE();  
  if (fIsMC && fIsDca) CreateDEM();  
  if (fIsDca) CreateDED();  

  //  if (fIsQa) if (fESDtrackCuts) fQaList->Add(fESDtrackCuts);
  TH1::AddDirectory(oldStatus);
  PostData(1, fPhyList); 
  PostData(2, fQaList);
  PostData(3, fDcaList);
  PostData(4, fEffList);
}

//----------------------------------------------------------------------------------
void AliEbyENetChargeFluctuationTask::UserExec( Option_t * ){
  if (SetupEvent() < 0) {

    //   Printf("%d %d", SetupEvent(),fSystemType);

    PostData(1, fPhyList); 
    PostData(2, fQaList);
    PostData(3, fDcaList);
    PostData(4, fEffList);
    return;
  }
    
  //-- -  - - -  - - ---- -- - --- -- --- -- - --
  for (Int_t ii = 0; ii <= fOrder; ++ii)
    for (Int_t jj = 0; jj < 2; ++jj)
      fRedFactp[ii][jj] = 1.;
  
  //-- -  - - -  - - ---- -- - --- -- --- -- - -- 
  for (Int_t kk = 0; kk < 4; ++kk)
    for (Int_t jj = 0; jj < 2; ++jj)
      fNp[kk][jj] = 0;
  
  //-- -  - - -  - - ---- -- - --- -- --- -- - --
  for (Int_t kk = 0; kk < 4; ++kk)
    for (Int_t jj = 0; jj < 2; ++jj)
      fMCNp[kk][jj] = 0;  
  
   // Printf("Number of Track %d",fNTracks);

  fSubSampleIdx = fRanIdx->Integer(fNSubSamples);
  fNTracks  = (fESD) ? fESD->GetNumberOfTracks() : fAOD->GetNumberOfTracks();  

  if (fIsMC && fIsAOD) {
    fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fArrayMC)
      AliFatal("No array of MC particles found !!!"); 
  }

  //-- -  - - -  - - ---- -- - --- -- --- -- - --
  if(fSystemType == 0)      { 
    Execpp(); 

  } 
  //-- -  - - -  - - ---- -- - --- -- --- -- - --
  else if(fSystemType == 1) { 
    ExecpA(); 
  } 
  //-- -  - - -  - - ---- -- - --- -- --- -- - --
  else if(fSystemType == 2) { 
    if (fIsPhy || fIsQa) ExecAA(); 
    if (fIsMC && fIsEff) CalEC();
    if (fIsDca) CalED();
  } 
  else 
    return;
  

  //-- -  - - -  - - ---- -- - --- -- --- -- - --
  PostData(1, fPhyList); 
  PostData(2, fQaList);
  PostData(3, fDcaList);
  PostData(4, fEffList);
  
}


//________________________________________________________________________
Int_t AliEbyENetChargeFluctuationTask::SetupEvent() {
  ResetCurrentEvent();
  if (!fIsAOD && (SetupESD() < 0)) {
    AliError("ESD Event failed");
    return -1;
  }
  if (fIsAOD && (SetupAOD() < 0)) {
    AliError("AOD Event failed");
    return -1;
  }
  if ((fIsMC && !fIsAOD) && (SetupMC() < 0)) {
      AliError("MC Event failed");
      return -1;
  }
  
  if(SetupEventCR(fESDHandler, fAODHandler, fMCEvent) < 0) {
    AliError("MC Event failed");
    return -1;
  }

  return RejectedEvent() ? -2 : 0;
}

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::ResetCurrentEvent() { 
  fESD = NULL; 
  fAOD = NULL;
  if (fIsMC && !fIsAOD) 
    fMCEvent = NULL; 
  else if(fIsMC && fIsAOD) 
    fArrayMC = NULL; 
  return;
}


//________________________________________________________________________
Int_t AliEbyENetChargeFluctuationTask::SetupMC() {
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
  if (!mcH) {
    AliError("MC event handler not available");
    return -1;
  }

  fMCEvent = mcH->MCEvent();
  if (!fMCEvent) {
    AliError("MC event not available");
    return -1;
  }

  AliHeader* header = fMCEvent->Header();
  if (!header) {
    AliError("MC header not available");
    return -1;
  }

  fMCStack = fMCEvent->Stack(); 
  if (!fMCStack) {
    AliError("MC stack not available");
    return -1;
  }
 
  if (!header->GenEventHeader()) {
    AliError("Could not retrieve genHeader from header");
    return -1;
  }

  if (!fMCEvent->GetPrimaryVertex()){
    AliError("Could not get MC vertex");
    return -1;
  }

  return 0;
}


//________________________________________________________________________
Int_t AliEbyENetChargeFluctuationTask::SetupAOD() {
  fAODHandler= dynamic_cast<AliAODInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fAODHandler) {
    AliError("Could not get AOD input handler");
    return -1;
  } 

  fAOD = fAODHandler->GetEvent();
  if (!fAOD) {
    AliError("Could not get AOD event");
    return -1;
  }

  if (!fAODHandler->GetPIDResponse()) {
    AliError("Could not get PID response");
    return -1;
  } 

  if (!fAOD->GetPrimaryVertex()) {
    AliError("Could not get primary vertex");
    return -1;
  }

  if (!((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP()) {
    AliError("Could not get centrality");
    return -1;
  }

  return 0;
}


//________________________________________________________________________
Int_t AliEbyENetChargeFluctuationTask::SetupESD() {

  fESDHandler= dynamic_cast<AliESDInputHandler*> 
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!fESDHandler) {
    AliError("Could not get ESD input handler");
    return -1;
  } 

  fESD = (AliESDEvent*)fESDHandler->GetEvent();
  if (!fESD) {
    AliError("Could not get ESD event");
    return -1;
  }

  if (!fESDHandler->GetPIDResponse()) {
    AliError("Could not get PID response");
    return -1;
  } 

  if (!fESD->GetPrimaryVertexTracks()) {
    AliError("Could not get vertex from tracks");
    return -1;
  }

  if (!fESD->GetCentrality()) {
    AliError("Could not get centrality");
    return -1;
  }

  return 0;
}

//----------------------------------------------------------------------------------
void AliEbyENetChargeFluctuationTask::Execpp(){
  // Printf("Place");
}

//----------------------------------------------------------------------------------
void AliEbyENetChargeFluctuationTask::ExecpA(){
  
}


void AliEbyENetChargeFluctuationTask::FillQAThnRec(AliVTrack *track, Int_t gPid, Double_t rap) {
  Double_t charge = track->Charge() < 0 ? 0. : 1.;
  Double_t  rapp = (gPid == 0) ? track->Eta() : rap;
  Double_t rec[5] = {fCentralityBin+1,charge,rapp,track->Phi(),track->Pt()};
  
  if (gPid == 0) (static_cast<THnSparseF*>(fQaList->FindObject("fHnNchTrackUnCorr")))->Fill(rec);
  else if (gPid == 1) (static_cast<THnSparseF*>(fQaList->FindObject("fHnNpiTrackUnCorr")))->Fill(rec);
  else if (gPid == 2) (static_cast<THnSparseF*>(fQaList->FindObject("fHnNkaTrackUnCorr")))->Fill(rec);
  else if (gPid == 3) (static_cast<THnSparseF*>(fQaList->FindObject("fHnNprTrackUnCorr")))->Fill(rec);
}

void AliEbyENetChargeFluctuationTask::FillQAThnMc(AliVParticle *particle, Int_t gPid, Double_t rap) {
  Double_t charge = (particle->PdgCode() < 0) ? 0. : 1.;
  Double_t  rapp = (gPid == 0) ? particle->Eta() : rap;
  Double_t rec[5] = {fCentralityBin+1,charge,rapp,particle->Phi(),particle->Pt()};
  if (gPid == 0) (static_cast<THnSparseF*>(fQaList->FindObject("fHnNchTrackMc")))->Fill(rec);
  else if (gPid == 1) (static_cast<THnSparseF*>(fQaList->FindObject("fHnNpiTrackMc")))->Fill(rec);
  else if (gPid == 2) (static_cast<THnSparseF*>(fQaList->FindObject("fHnNkaTrackMc")))->Fill(rec);
  else if (gPid == 3) (static_cast<THnSparseF*>(fQaList->FindObject("fHnNprTrackMc")))->Fill(rec);
}

//----------------------------------------------------------------------------------
void AliEbyENetChargeFluctuationTask::ExecAA(){
  
  // Printf("Number of Track %d",fNTracks);

  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliVTrack *track = (fESD) ? 
      static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : 
      static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 
    if(!AcceptTrackL(track)) continue;
    if(!AcceptTrackLDCA(track)) continue;
    Int_t icharge = track->Charge() < 0 ? 0 : 1;
    fNp[0][icharge] += 1.; 
    if (!fIsNu && fIsQa) FillQAThnRec(track,0,0);
    
    Int_t a = fHelperPID->GetParticleSpecies(track,kTRUE);
    if(a < 0 || a > 2) continue;
    Int_t b = a + 1;

    Double_t rap;
    if (!fIsNu && !TrackRapidity(track,rap,b)) continue;
    fNp[b][icharge] += 1.; 
    if (!fIsNu && fIsQa) FillQAThnRec(track, b,rap);

  }
  
  //---- - -- - - - - -   -  -- - - - ---- - - - ---
  if(fIsPhy) {
    if (fIs3D) {FillSourceHistos(0); }
    else {
    FillBasicHistos(0,0);
    FillBasicHistos(0,1);
    if (fIsRatio) {
      FillRatioHistos(0,0);
      if(fIsPer)  FillRatioHistos(0,1);
    }
    if (fIsSub) {
      FillGroupHistos("PhyBinSS",fSubSampleIdx,kFALSE,0);
      if(fIsPer) FillGroupHistos("PhyPerSS",fSubSampleIdx,kFALSE,1);
    }
    if (fIsBS)  {
      for (Int_t i = 0; i < fNSubSamples; ++i)  {
	FillGroupHistos("PhyBinBS",fRan->Integer(fNSubSamples),kFALSE,0);
	if(fIsPer)  FillGroupHistos("PhyPerBS",fRan->Integer(fNSubSamples),kFALSE,1);
      }
    }
    }
  }
  //---- - -- - - - - -   -  -- - - - ---- - - - ---
  if (fIsMC) {
    if (fIsAOD) {

      //   Printf("Number of Track %d",fArrayMC->GetEntries());

      for (Int_t idxMC = 0; idxMC < fArrayMC->GetEntries(); idxMC++) {
	AliAODMCParticle *particle = static_cast<AliAODMCParticle*>(fArrayMC->At(idxMC));
	if (!particle) 
	  continue;
	if (!AcceptTrackLMC((AliVParticle*)particle, idxMC)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;
	fMCNp[0][icharge]    += 1.;    	   
	if(!fIsNu && fIsQa)FillQAThnMc((AliVParticle*)particle,0,0);
	Int_t iPid = 0;  
	if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
	else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
	else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
	else    iPid = 0;

	if(iPid < 1 || iPid > 3) continue;

	Double_t rap;
	if (!fIsNu && !ParticleRapidity((AliVParticle*)particle, rap, iPid)) continue;
	fMCNp[iPid][icharge] += 1.;    
	if(!fIsNu && fIsQa)FillQAThnMc((AliVParticle*)particle,iPid,rap);   
      }

    } else if (fESD) {
      //---- - -- - - - - -   -  -- - - - ---- - - - --- 
         for (Int_t idxMC = 0; idxMC < fStack->GetNprimary(); ++idxMC) {
	   AliVParticle* particle = fMCEvent->GetTrack(idxMC);
	   if (!particle) 
	     continue;
	   if (!AcceptTrackLMC(particle, idxMC)) continue;
	   Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;
	   fMCNp[0][icharge]  += 1.;    	   
	   if(!fIsNu && fIsQa)FillQAThnMc(particle, 0, 0);

	   Int_t iPid = 0;  
	   if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
	   else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
	   else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
	   else  iPid = 0;
	   if(iPid < 1 || iPid > 3) continue;
	   Double_t rap;
	   if (!fIsNu && !ParticleRapidity(particle, rap, iPid)) continue;
	   fMCNp[iPid][icharge] += 1.;    
	   if(!fIsNu && fIsQa)FillQAThnMc(particle, iPid, rap);    
	 }
    }
    
    //---- - -- - - - - -   -  -- - - - ---- - - - --- 
    if (fIsPhy) {
      if (fIs3D) { FillSourceHistos(1); }
      else {
      FillBasicHistos(kTRUE,0);
      FillBasicHistos(kTRUE,1);
      if (fIsRatio) {
	FillRatioHistos(kTRUE,0);
	if(fIsPer)  FillRatioHistos(kTRUE,1);
      }
      if (fIsSub) {
	FillGroupHistos("MCBinSS",fSubSampleIdx,kFALSE,0);
	if(fIsPer) FillGroupHistos("MCPerSS",fSubSampleIdx,kFALSE,1);
      }
      if (fIsBS)  {
	for (Int_t i = 0; i < fNSubSamples; ++i)  {
	  FillGroupHistos("MCBinBS",fRan->Integer(fNSubSamples),kFALSE,0);
	  if(fIsPer)  FillGroupHistos("MCPerBS",fRan->Integer(fNSubSamples),kFALSE,1);
	}
      }
      }
    }
    //---- - -- - - - - -   -  -- - - - ---- - - - --- 
  }
 }

//___________________________________________________________
void AliEbyENetChargeFluctuationTask::Terminate( Option_t * ){
  Info("AliEbyENetChargeFluctuationTask"," Task Successfully finished");
}

//___________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::AcceptTrackLDCA(AliVTrack *track) const {
  Float_t dca[2], cov[3]; // 
  if (fESD)
    (dynamic_cast<AliESDtrack*>(track))->GetImpactParameters(dca, cov);
  else  {
    AliAODTrack* clone = dynamic_cast<AliAODTrack*>(track->Clone("trk_clone"));
    if(clone->TestBit(AliAODTrack::kIsDCA)){
      dca[0] = clone->DCA();
      dca[1] = clone->ZAtDCA();
    } else {
      Double_t dcaa[2], cova[3];
      Double_t gBzKg = fAOD->GetMagneticField();
      Bool_t propagate = clone->PropagateToDCA(fAOD->GetPrimaryVertex(),gBzKg,1000000.,dcaa,cova);
      if (!propagate) { dca[0] = -999; dca[1] = -999;}
      dca[0] = Float_t(dcaa[0]);
      dca[1] = Float_t(dcaa[1]);
    }
    delete clone;
  }

  if ( TMath::Abs(dca[0]) > fDcaXy ) return kFALSE;
  if ( TMath::Abs(dca[1]) > fDcaZ )  return kFALSE;

  return kTRUE;

}

//___________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::AcceptTrackL(AliVTrack *track) const {
 if (!track) 
   return kFALSE; 
 if (track->Charge() == 0) 
   return kFALSE; 
  
  
 if (fIsAOD) {  // AOD
 AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
 if (!trackAOD) {
   AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
   return kFALSE; 
 }
 
 if (!trackAOD->TestFilterBit(fAODtrackCutBit)) return kFALSE;}
 else { if(!fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))  return kFALSE; }

 if(track->Pt() < fPtMin || track->Pt() > fPtMax )  return kFALSE; 
 if (TMath::Abs(track->Eta()) > fEtaMax) return kFALSE; 
 
  
 if (track->Phi() > fPhiMax) return kFALSE;   
 return kTRUE;
}


//___________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::AcceptTrackLMC(AliVParticle *particle, Int_t idxMC) const {
  if(!particle) return kFALSE;

  if (particle->Charge() == 0.0) 
    return kFALSE;
  
  if (fIsAOD) {
    if(!(static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary()) return kFALSE;
  } else {
    if(!fStack->IsPhysicalPrimary(idxMC)) return kFALSE;
  }
 
  if (particle->Pt() < fPtMin || particle->Pt() > fPtMax) return kFALSE;
  if (TMath::Abs(particle->Eta()) > fEtaMax) return kFALSE;
  if (particle->Phi() > fPhiMax) return kFALSE; 

  return kTRUE;
}



//________________________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::ParticleRapidity(AliVParticle *particle, Double_t &rap, Int_t gCurPid) {
  
 if (gCurPid == 0) {
   rap = particle->Eta(); 
   if (TMath::Abs(rap) > fEtaMax) return kFALSE;
   return kTRUE;
 }
  
 Double_t              mass = AliPID::ParticleMass(AliPID::kPion);
 if(gCurPid == 1)      mass = AliPID::ParticleMass(AliPID::kPion);
 else if(gCurPid == 2) mass = AliPID::ParticleMass(AliPID::kKaon);
 else if(gCurPid == 3) mass = AliPID::ParticleMass(AliPID::kProton);
 
  Double_t p  = particle->P();
  Double_t pz = particle->Pz();
  Double_t eP = TMath::Sqrt(p*p + mass*mass);
  rap          = 0.5 * TMath::Log((eP + pz) / (eP - pz));  

  if (TMath::Abs(rap) > fRapMax) return kFALSE;
  return kTRUE;
}

 
//________________________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::TrackRapidity(AliVTrack *track, Double_t &rap, Int_t gCurPid) {
  if (gCurPid == 0) { 
    rap = track->Eta(); 
    if (TMath::Abs(rap) > fEtaMax) return kFALSE;
    return kTRUE; 
  }
  
  Double_t              mass = AliPID::ParticleMass(AliPID::kPion);
  if(gCurPid == 1)      mass = AliPID::ParticleMass(AliPID::kPion);
  else if(gCurPid == 2) mass = AliPID::ParticleMass(AliPID::kKaon);
  else if(gCurPid == 3) mass = AliPID::ParticleMass(AliPID::kProton);

  Double_t pvec[3]; track->GetPxPyPz(pvec);
  Double_t p  = track->P();
  Double_t eP = TMath::Sqrt(p*p + mass*mass);
  rap = 0.5 * TMath::Log((eP + pvec[2]) / (eP - pvec[2]));
	   
  if (TMath::Abs(rap) > fRapMax) return kFALSE;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::IsFindableInTPC(Int_t label) {
  AliMCParticle *mcParticle = static_cast<AliMCParticle*>(fMCEvent->GetTrack(label));
  if(!mcParticle) return kFALSE;
  Int_t counter; 
  Float_t tpcTrackLength = mcParticle->GetTPCTrackLength(AliTracker::GetBz(), 0.05, counter, 3.0); 
  return (tpcTrackLength > fMinTrackLengthMC);    
}

const Char_t* fGEvtNames[] = {"All", "IsTriggered", "HasVertex", "Vx<Vx_{Max}", "Vy<Vy_{Max}", "Vz<Vz_{Max}", "Centrality [0,100]%", "Centrality [<0,>100]%"};
const Char_t* fGCMxNames[] = {"5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};
const Char_t* fGCMxNames10[] = {"10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};
const Char_t* fGTrgNames[] = {"kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" }; 
const Char_t* fGCntNames[] = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%","50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};
const Char_t* fGCntNames10[] = {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%","50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::CreateQA() {
 
  Int_t    bhuc[5]  = {fCentralityBinMax, fGNBinsSign, fGNBinsRap,fGNBinsPhi,fGNBinsPt };      
  Double_t mnhuc[5] = {0.5,fGRngSign[0],fGRngRap[0],fGRngPhi[0],fGRngPt[0]};  
  Double_t mxhuc[5] = {fCentralityBinMax+0.5,fGRngSign[1],fGRngRap[1],fGRngPhi[1],fGRngPt[1]};  
  const Char_t *ctname = "cent:sign:y:phi:pt";
			     
  fQaList->Add(new THnSparseF("fHnNpiTrackUnCorr",ctname,  5, bhuc, mnhuc, mxhuc));
  fQaList->Add(new THnSparseF("fHnNkaTrackUnCorr", ctname, 5, bhuc, mnhuc, mxhuc));
  fQaList->Add(new THnSparseF("fHnNprTrackUnCorr", ctname, 5, bhuc, mnhuc, mxhuc));

  if (fIsMC) {
    fQaList->Add(new THnSparseF("fHnNpiTrackMc", ctname, 5, bhuc, mnhuc, mxhuc));
    fQaList->Add(new THnSparseF("fHnNkaTrackMc", ctname, 5, bhuc, mnhuc, mxhuc));
    fQaList->Add(new THnSparseF("fHnNprTrackMc", ctname, 5, bhuc, mnhuc, mxhuc));

  }

  ctname = "cent:sign:eta:phi:pt";
  bhuc[2]  = fGNBinsEta;
  mnhuc[2] = fGRngEta[0];
  mxhuc[2] = fGRngEta[1]; 
  
  fQaList->Add(new THnSparseF("fHnNchTrackUnCorr", ctname, 5, bhuc, mnhuc, mxhuc));
  if (fIsMC) {
    fQaList->Add(new THnSparseF("fHnNchTrackMc", ctname, 5, bhuc, mnhuc, mxhuc));
  }
}

const Char_t* fgkPidLatex[4][2]= {{"N_{-}","N_{+}"}, {"N_{#pi^{-}}","N_{#pi^{+}}"},{"N_{K^{-}}","N_{K^{+}}"}, {"N_{#bar{p}}","N_{p}"}};
const Char_t* fgkPidTitles[4][2]= {{"Negative","Positive"},{"Anti-Pions","Pions"},{"Anti-Kaons","Kaons"}, {"Anti-Protons","Protons"}};

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::CreateBasicQA() {
  // -- Basic QA

 
  fQaList->Add(new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5));
  fQaList->Add(new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5));
  
  for ( Int_t ii=0; ii < fHEventStatMax; ii++ ) {
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat0"))))->GetXaxis()->SetBinLabel(ii+1, fGEvtNames[ii]);
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat1"))))->GetXaxis()->SetBinLabel(ii+1, fGEvtNames[ii]);
  }
  if (fIsTen) {
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat0"))))->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", fGCMxNames10[fCentralityBinMax-1]));
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat1"))))->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", fGCMxNames10[fCentralityBinMax-1]));
  }  else {
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat0"))))->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", fGCMxNames[fCentralityBinMax-1]));
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat1"))))->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", fGCMxNames[fCentralityBinMax-1]));
  }

  fQaList->Add(new TH1F("hTriggerStat","Trigger statistics;Trigger;Events", fNTriggers,-0.5,fNTriggers-0.5));

  for ( Int_t ii=0; ii < fNTriggers; ii++ ) {
    (static_cast<TH1F*>(fQaList->FindObject(Form("hTriggerStat"))))->GetXaxis()->SetBinLabel(ii+1, fGTrgNames[ii]);
  }

  
  // -- Initialize trigger statistics histograms
  
  fQaList->Add(new TH1F("hCentralityStat",Form("Centrality statistics (%d);Centrality Bins;Events", fNCentralityBins),
					   fNCentralityBins,0.5,fNCentralityBins+0.5));
  
  for ( Int_t ii=0; ii < fNCentralityBins; ii++ ) {
    if (fIsTen) {
      (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityStat"))))->GetXaxis()->SetBinLabel(ii+1, fGCntNames10[ii]);
    } else {
      (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityStat"))))->GetXaxis()->SetBinLabel(ii+1, fGCntNames[ii]);
    }
  }
  
  fQaList->Add(new TH1F("hCentralityPercentileAccepted",Form("Centrality Percentile statistics (%d);Centrality Bins;Events",fNbwcBin), 
			fNbwcBin,0.5,fNbwcBin + 0.5));
  
  fQaList->Add(new TH1F("hCentralityPercentileAll",Form("Centrality Percentile statistics (%d);Centrality Bins;Events",fNbwcBin), 
			fNbwcBin,0.5,fNbwcBin + 0.5));


  fQaList->Add(new TH2F("fHistQAvx",  Form("Histo Vx Selected;Centrality (%d);Vx",fNbwcBin), fNbwcBin,0.5,fNbwcBin + 0.5, 500, -5., 5.));
  fQaList->Add(new TH2F("fHistQAvy",  Form("Histo Vy Selected;Centrality (%d);Vy",fNbwcBin), fNbwcBin,0.5,fNbwcBin + 0.5, 500, -5., 5.));
  fQaList->Add(new TH2F("fHistQAvz",  Form("Histo Vz Selected;Centrality (%d);Vz",fNbwcBin), fNbwcBin,0.5,fNbwcBin + 0.5, 500, -25., 25.)); 
  fQaList->Add(new TH2F("fHistQAvxA", Form("Histo Vx;Centrality (%d);Vx",fNbwcBin), fNbwcBin,0.5,fNbwcBin + 0.5, 500, -5., 5.));
  fQaList->Add(new TH2F("fHistQAvyA", Form("Histo Vy;Centrality (%d);Vy",fNbwcBin), fNbwcBin,0.5,fNbwcBin + 0.5, 500, -5., 5.));
  fQaList->Add(new TH2F("fHistQAvzA", Form("Histo Vz;Centrality (%d);Vz",fNbwcBin), fNbwcBin,0.5,fNbwcBin + 0.5, 500, -25., 25.));
  
  if (fHelperPID) {
    fQaList->Add(new TList);
    TList *list =  static_cast<TList*>(fQaList->Last());
    list->SetName("HelperPID");
    list->SetOwner(kTRUE);
    TList *ll = (TList*)fHelperPID->GetOutputList();
    for (Int_t ikey = 0; ikey < ll->GetEntries(); ikey++) {
      list->Add(ll->At(ikey));
    }
  }
  
}


//----------------------------------------------------------------------------------
void AliEbyENetChargeFluctuationTask::SetAnal(Int_t i, Int_t j){

  if      (i == 0 ) { fIsQa  = 1;                                                                 }
  else if (i == 1 ) { fIsDca = 1;                                                                 }
  else if (i == 2 ) { fIsEff = 1;                                                                 } 
  else if (i == 3 ) { fIsDca = 1; fIsEff  = 1;                                                    }  
  else if (i == 4 ) { fIsPhy = 1; fIsRatio= 1;                                                    } 
  else if (i == 5 ) { fIsPhy = 1; fIsRatio= 1; fIsPer = 1;                                        }   
  else if (i == 6 ) { fIsPhy = 1; fIsSub  = 1;                                                    }
  else if (i == 7 ) { fIsPhy = 1; fIsSub  = 1; fIsPer = 1;                                        }
  else if (i == 8 ) { fIsPhy = 1; fIsBS   = 1;                                                    }   
  else if (i == 9 ) { fIsPhy = 1; fIsBS   = 1; fIsPer = 1;                                        }    
  else if (i == 10) { fIsPhy = 1; fIsBS   = 1; fIsSub = 1;                                        }    
  else if (i == 11) { fIsPhy = 1; fIsBS   = 1; fIsSub = 1; fIsPer = 1;                            }    
  else if (i == 12) { fIsPhy = 1; fIsBS   = 1; fIsSub = 1; fIsQa  = 1;                            }    
  else if (i == 13) { fIsPhy = 1; fIsBS   = 1; fIsSub = 1; fIsQa  = 1; fIsPer   = 1;              }    
  else if (i == 14) { fIsPhy = 1; fIsBS   = 1; fIsSub = 1; fIsQa  = 1; fIsPer   = 1; fIsRatio = 1;}    
  else if (i == 15) { fIsNu  = 1; fIsPhy  = 1; fIsBS  = 1; fIsSub = 1;                            }     
  else if (i == 16) { fIsNu  = 1; fIsPhy  = 1; fIsBS  = 1; fIsSub = 1; fIsPer   = 1;              }     
  else if (i == 17) { fIsNu  = 1; fIsPhy  = 1; fIsBS  = 1; fIsSub = 1; fIsRatio = 1;              }     
  else if (i == 18) { fIsNu  = 1; fIsPhy  = 1; fIsBS  = 1; fIsSub = 1; fIsRatio = 1; fIsPer = 1;  }     
  else if (i == 19) { fIsPhy = 1; fIsBS   = 1; fIsSub = 1; fIsPer = 1; fIsRatio = 1;              }    
  else if (i == 20) { fIsPhy = 1; fIs3D   = 1;                                                    }    
  else {fIsPhy= 0;    fIsEff = 0; fIsDca  = 0; fIsQa  = 0; fIsSub = 0; fIsBS    = 0; fIsPer = 0; fIsRatio = 0;}
   
  if      (j ==  0) {fIsTen = 0; }
  else if (j ==  1) {fIsTen = 0; fNbwcBin = 100;}
  else if (j ==  2) {fIsTen = 0; fNbwcBin = 50; }
  else if (j ==  3) {fIsTen = 0; fNbwcBin = 40; }
  else if (j ==  4) {fIsTen = 0; fNbwcBin = 25; }
  else if (j ==  5) {fIsTen = 0; fNbwcBin = 20; }
  else if (j ==  6) {fIsTen = 1; fNbwcBin = 100;}
  else if (j ==  7) {fIsTen = 1; fNbwcBin = 50; }
  else if (j ==  8) {fIsTen = 1; fNbwcBin = 40; }
  else if (j ==  9) {fIsTen = 1; fNbwcBin = 25; }
  else if (j == 10) {fIsTen = 1; fNbwcBin = 20; }
  else if (j == 11) {fIsTen = 0; fNbwcBin = 100;fNeedQa  = 1;}
  else if (j == 12) {fIsTen = 0; fNbwcBin = 50; fNeedQa  = 1;}
  else if (j == 13) {fIsTen = 0; fNbwcBin = 40; fNeedQa  = 1;}
  else if (j == 14) {fIsTen = 0; fNbwcBin = 25; fNeedQa  = 1;}
  else if (j == 15) {fIsTen = 0; fNbwcBin = 20; fNeedQa  = 1;}
  else if (j == 16) {fIsTen = 1; fNbwcBin = 100;fNeedQa  = 1;}
  else if (j == 17) {fIsTen = 1; fNbwcBin = 50; fNeedQa  = 1;}
  else if (j == 18) {fIsTen = 1; fNbwcBin = 40; fNeedQa  = 1;}
  else if (j == 19) {fIsTen = 1; fNbwcBin = 25; fNeedQa  = 1;}
  else if (j == 20) {fIsTen = 1; fNbwcBin = 20; fNeedQa  = 1;}
  
  if (fIsTen) {
    fNCentralityBins   = 10;
    fCentralityBinMax  = 10;
  } 


  Printf(" >>> I%d PH%d EF%d DC%d QA%d RA%d SS%d BS%d PE%d NU%d NQ%d TE%d",i, 
	 fIsPhy, fIsEff, fIsDca, fIsQa, fIsRatio, fIsSub, fIsBS, fIsPer, fIsNu, fNeedQa, fIsTen);
  if (fIsEff || fIsDca) {fPtMin = 0.2; fPtMax = 3.3;}  
 
}

const Char_t*      fgkPidName[4] = {"Nch","Npi","Nka","Npr"};
const Char_t*   fgkPidShLatex[4] = {"N","#pi","K","p"};
const Char_t*  fgkNetHistName[4] = {"","Plus","Minus","Net"};
const Char_t* fgkNetHistLatex[4] = {"+ + +","+","-","+ - -"};

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::InitPhy() {

  TString name = Form("#it{p}_{T} [%.1f,%.1f] : #eta [%.1f,%.1f]",fPtMin, fPtMax, fEtaMin, fEtaMax);

  if (fIs3D) {CreateSourceHistos(name.Data(),0); if (fIsMC) CreateSourceHistos(name.Data(),1); }
  else{
    CreateBasicHistos(name.Data(),0,0);
    CreateBasicHistos(name.Data(),0,1);
    
    if (fIsRatio) {
      CreateRatioHistos(name.Data(),0,0);
      if (fIsPer) 
	CreateRatioHistos(name.Data(),0,1);
    }
    
    if (fIsSub) {
      CreateGroupHistos("PhyBinSS",name.Data(),fNSubSamples,0);  
      if (fIsPer) 
	CreateGroupHistos("PhyPerSS",name.Data(),fNSubSamples,1);  
    }
    
    if (fIsBS){
      CreateGroupHistos("PhyBinBS",name.Data(),fNSubSamples,0);  
      if (fIsPer) CreateGroupHistos("PhyPerBS",name.Data(),fNSubSamples,1);  
    }
    
    if (fIsMC) {
      CreateBasicHistos(name.Data(),1,0);
      CreateBasicHistos(name.Data(),1,1);
      
      if (fIsRatio){
	CreateRatioHistos(name.Data(),1,0);
	if (fIsPer)  
	  CreateRatioHistos(name.Data(),1,1);
      }
      
      if (fIsSub){
	CreateGroupHistos("MCBinSS", name.Data(),fNSubSamples,0);  
	if (fIsPer)  
	  CreateGroupHistos("MCPerSS", name.Data(),fNSubSamples,1);  
      }
      
      if (fIsBS){
	CreateGroupHistos("MCBinBS",name.Data(),fNSubSamples,0);  
	if (fIsPer)  
	  CreateGroupHistos("MCPerBS", name.Data(),fNSubSamples,1);  
      }
    }
  }
}

//________________________________________________________________________
void  AliEbyENetChargeFluctuationTask::CreateBasicHistos(const Char_t *title, Bool_t isMC, Bool_t isPer)  {

  TString nmc  = (isMC) ? "MC" : "Phy";
  TString name = (isPer) ? Form("%sPer",nmc.Data()) : Form("%sBin",nmc.Data());

 
  TString strTit(title);
  Double_t etaRange[2] = {fEtaMin,fEtaMax};
  fPhyList->Add(new TList);
  TList *list =  static_cast<TList*>(fPhyList->Last());
  list->SetName(Form("f%s",name.Data()));
  list->SetOwner(kTRUE);
  

  Int_t nBinsCent  =  (isPer) ? fNbwcBin : fCentralityBinMax;
  Double_t centBinRange[2];  
  centBinRange[0]  =  0.5;
  centBinRange[1]  =  (isPer) ?  fNbwcBin + 0.5 : fCentralityBinMax+0.5;

  fQaList->Add(new TH1F(Form("h%s%sEventsNch",nmc.Data(), name.Data()),
			"EventStat;Centrality Bins;Events",
			nBinsCent,centBinRange[0],centBinRange[1]));

  fQaList->Add(new TH1F(Form("h%s%sEventsNpi",nmc.Data(), name.Data()),
			"EventStat;Centrality Bins;Events",
			nBinsCent,centBinRange[0],centBinRange[1]));
  
  fQaList->Add(new TH1F(Form("h%s%sEventsNka",nmc.Data(), name.Data()),
			"EventStat;Centrality Bins;Events",
			nBinsCent,centBinRange[0],centBinRange[1]));
  
  fQaList->Add(new TH1F(Form("h%s%sEventsNpr",nmc.Data(), name.Data()),
			"EventStat;Centrality Bins;Events",
			nBinsCent,centBinRange[0],centBinRange[1]));
  
  
 

  if (!fIsNu) {
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
      TString sNetTitle(Form("%s - %s", fgkPidLatex[iPid][1], fgkPidLatex[iPid][0]));
      strTit = (iPid != 0 ) ? Form(" |y|<%.1f", fRapMax) : Form(" |#eta| < %.1f", etaRange[1]);
      
      list->Add(new TProfile(Form("fProfTot%sPlus%s", fgkPidName[iPid],name.Data()), 
			     Form("(%s) : %s;Centrality(%d);(%s)",fgkPidName[iPid], 
				  strTit.Data(), nBinsCent,sNetTitle.Data()),
			     nBinsCent, centBinRange[0], centBinRange[1]));
      
      list->Add(new TProfile(Form("fProfTot%sMinus%s", fgkPidName[iPid],name.Data()), 
			     Form("(%s) : %s;Centrality(%d);(%s)",fgkPidName[iPid], 
				  strTit.Data(), nBinsCent,sNetTitle.Data()),
			     nBinsCent, centBinRange[0], centBinRange[1]));
      
      
      
      for (Int_t idx = 1; idx <= fOrder; ++idx) {
	list->Add(new TProfile(Form("fProf%s%sNet%dM", fgkPidName[iPid],name.Data(), idx), 
			       Form("(%s)^{%d} : %s;Centrality(%d);(%s)^{%d}", 
				    sNetTitle.Data(), idx, strTit.Data(), nBinsCent, sNetTitle.Data(), idx),
			       nBinsCent, centBinRange[0], centBinRange[1]));
      }
      
      for (Int_t ii = 0; ii <= fOrder; ++ii) {
	for (Int_t kk = 0; kk <= fOrder; ++kk) {
	  list->Add(new TProfile(Form("fProf%s%sNetF%02d%02d", fgkPidName[iPid], name.Data(), ii, kk),
				 Form("f_{%02d%02d} : %s;Centrality(%d);f_{%02d%02d}", 
				      ii, kk, strTit.Data(), nBinsCent,ii, kk),
				 nBinsCent, centBinRange[0], centBinRange[1]));
	}
      }
    }  
  }
  
  for (Int_t iPhy = 0; iPhy < 58; ++iPhy) { 
    list->Add(new TProfile(Form("fProf%sNu%02d",name.Data(),iPhy),
			   Form("Physics Variable for index %d | %s ; Centrality(%d);",
				iPhy,name.Data(),nBinsCent),nBinsCent, 
			   centBinRange[0], centBinRange[1]));
  }
  
  return;
}


//________________________________________________________________________
void  AliEbyENetChargeFluctuationTask::CreateSourceHistos(const Char_t *title, Bool_t isMC)  {

  TString nmc  = (isMC)  ? "MC" : "Phy";
  TString name = "Source";
    
  fPhyList->Add(new TList);
  TList *list =  static_cast<TList*>(fPhyList->Last());
  list->SetName(Form("f%s%s",name.Data(),nmc.Data()));
  list->SetOwner(kTRUE);
  
  list->Add(new TH3I(Form("fHist%sNch",name.Data()),Form("%s = Ch;Cent;0;1",title),
		     100,0.5,100.5,1000,0.5,1000.5,1000,0.5,1000.5));
  list->Add(new TH3I(Form("fHist%sNpi",name.Data()),Form("%s = #pi;Cent;0;1",title),
		     100,0.5,100.5,800,0.5,800.5,800,0.5,800.5));
  list->Add(new TH3I(Form("fHist%sNka",name.Data()),Form("%s = K;Cent;0;1",title),
		     100,0.5,100.5,500,0.5,500.5,500,0.5,500.5));
  list->Add(new TH3I(Form("fHist%sNpr",name.Data()),Form("%s = P;Cent;0;1",title),
		     100,0.5,100.5,200,0.5,200.5,200,0.5,200.5));

}

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::FillSourceHistos(Bool_t isMC)  {
  Double_t np[4][2];
  if (isMC) {
    for (Int_t i = 0; i < 4; i++) {
      np[i][0] = fMCNp[i][0];
      np[i][1] = fMCNp[i][1];
    }
  } else {
    for (Int_t i = 0; i < 4; i++) {
      np[i][0] = fNp[i][0];
      np[i][1] = fNp[i][1];
    }
  }

  Bool_t isZero = kTRUE;
  Bool_t isZeroPid[4] = {1,1,1,1};
  if ((np[0][0] == 0) || (np[0][1] == 0))
    {isZero = kFALSE; isZeroPid[0] = kFALSE;}
  else if ((np[1][0] == 0) || (np[1][1] == 0))
    {isZero = kFALSE; isZeroPid[1] = kFALSE;} 
  else if ((np[2][0] == 0) || (np[2][1] == 0))
    {isZero = kFALSE; isZeroPid[2] = kFALSE;}
  else if ((np[3][0] == 0) || (np[3][1] == 0))
    {isZero = kFALSE; isZeroPid[3] = kFALSE;}
  else isZero = kTRUE;

  TString nmc  = (isMC) ? "MC" : "Phy";
  TString name = "Source";
  
  TList *list = static_cast<TList*>(fPhyList->FindObject(Form("f%s%s",name.Data(),nmc.Data())));
  
  if (isZeroPid[0]) (static_cast<TH3I*>(list->FindObject(Form("fHist%sNch",name.Data()))))->Fill(fCentralityPercentile,np[0][0],np[0][1]);
  if (isZeroPid[1]) (static_cast<TH3I*>(list->FindObject(Form("fHist%sNpi",name.Data()))))->Fill(fCentralityPercentile,np[1][0],np[1][1]);
  if (isZeroPid[2]) (static_cast<TH3I*>(list->FindObject(Form("fHist%sNka",name.Data()))))->Fill(fCentralityPercentile,np[2][0],np[2][1]);
  if (isZeroPid[3]) (static_cast<TH3I*>(list->FindObject(Form("fHist%sNpr",name.Data()))))->Fill(fCentralityPercentile,np[3][0],np[3][1]);
  
}


//________________________________________________________________________
void  AliEbyENetChargeFluctuationTask::CreateRatioHistos(const Char_t *title, Bool_t isMC, Bool_t isPer)  {

  TString nmc  = (isMC)  ? "MC" : "Phy";
  TString name = (isPer) ? Form("%sPer",nmc.Data()) : Form("%sBin",nmc.Data());
 
  TString strTit(title);
  
  fPhyList->Add(new TList);
  TList *list =  static_cast<TList*>(fPhyList->Last());
  list->SetName(Form("fRatio%s",name.Data()));
  list->SetOwner(kTRUE);
  
  Int_t    nRbin  = 10000;
  Double_t mRat[] = {0,2.0};
    
  Int_t nBinsCent  =  (isPer) ? fNbwcBin : fCentralityBinMax;
  Double_t centBinRange[2];  
  centBinRange[0]  =  0.5;
  centBinRange[1]  =  (isPer) ?  fNbwcBin + 0.5 : fCentralityBinMax+0.5;

  TString xyz = Form("|y| < %.1f | #eta [%3.1f-%3.1f]",fRapMax,fEtaMin,fEtaMax); 

  if (fIsNu) {
    for (Int_t i = 0; i < 22; i++) {
      list->Add(new TH2F(Form("fHist%sRatio%02d",name.Data(),i), 
			 Form("(%s %s);Centrality(%d);Ratio Idx %d", xyz.Data(), strTit.Data(),nBinsCent,i),
			 nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
    }
  }
  
  Int_t bin[4] = {2800,2200,1200,600}; 
  Int_t bd[] = {1,2,2,2};
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    Int_t bb = bin[iPid];

    for (Int_t iNet = 0; iNet < 4; ++iNet) {
      Int_t bn     = (iNet == 3) ?  500   : bb/bd[iNet]; 
      Float_t blow = (iNet == 3) ? -250.5 : -0.5;
      Float_t bup  = (iNet == 3) ?  249.5 : bn-0.5;
      
      list->Add(new TH2F(Form("fHist%sDist%s%s",name.Data(), fgkPidName[iPid], 
			      fgkNetHistName[iNet]), 
			 Form("(%s %s) : %s Distribution;Centrality(%d);%s_{(%s)}", 
			      xyz.Data(), strTit.Data(),
			      fgkPidShLatex[iPid], nBinsCent,fgkPidShLatex[iPid],
			      fgkNetHistLatex[iNet]),
			 nBinsCent, centBinRange[0], centBinRange[1], bn, blow,bup));    
    }
  }

return;
}

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::FillBasicHistos(Bool_t isMC, Bool_t isPer)  {
  Double_t np[4][2];
  if (isMC) {
    for (Int_t i = 0; i < 4; i++) {
      np[i][0] = fMCNp[i][0];
      np[i][1] = fMCNp[i][1];
    }
  } else {
    for (Int_t i = 0; i < 4; i++) {
      np[i][0] = fNp[i][0];
      np[i][1] = fNp[i][1];
    }
  }

  Bool_t isZero = kTRUE;
  Bool_t isZeroPid[4] = {1,1,1,1};
  if ((np[0][0] == 0) || (np[0][1] == 0))
    {isZero = kFALSE; isZeroPid[0] = kFALSE;}
  else if ((np[1][0] == 0) || (np[1][1] == 0))
    {isZero = kFALSE; isZeroPid[1] = kFALSE;} 
  else if ((np[2][0] == 0) || (np[2][1] == 0))
    {isZero = kFALSE; isZeroPid[2] = kFALSE;}
  else if ((np[3][0] == 0) || (np[3][1] == 0))
    {isZero = kFALSE; isZeroPid[3] = kFALSE;}
  else isZero = kTRUE;

  TString nmc  = (isMC) ? "MC" : "Phy";
  TString name = (isPer) ? Form("%sPer",nmc.Data()) : Form("%sBin",nmc.Data());
  
  Float_t centralityBin = (isPer) ? (fCentralityPercentile + 1) : (fCentralityBin + 1);
  TList *list = static_cast<TList*>(fPhyList->FindObject(Form("f%s",name.Data())));
  
  if (isZeroPid[0]) (static_cast<TH1F*>(fQaList->FindObject(Form("h%s%sEventsNch",nmc.Data(), name.Data()))))->Fill(centralityBin);
  if (isZeroPid[1]) (static_cast<TH1F*>(fQaList->FindObject(Form("h%s%sEventsNpi",nmc.Data(), name.Data()))))->Fill(centralityBin);
  if (isZeroPid[2]) (static_cast<TH1F*>(fQaList->FindObject(Form("h%s%sEventsNka",nmc.Data(), name.Data()))))->Fill(centralityBin);
  if (isZeroPid[3]) (static_cast<TH1F*>(fQaList->FindObject(Form("h%s%sEventsNpr",nmc.Data(), name.Data()))))->Fill(centralityBin);
  
  if (!fIsNu) {
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
      if (isZeroPid[iPid]) {
	Int_t deltaNp = np[iPid][1]-np[iPid][0];  
	Double_t delta = 1.;
	
	(static_cast<TProfile*>(list->FindObject(Form("fProfTot%sPlus%s", fgkPidName[iPid], name.Data()))))->Fill(centralityBin, np[iPid][1]);
	(static_cast<TProfile*>(list->FindObject(Form("fProfTot%sMinus%s", fgkPidName[iPid], name.Data()))))->Fill(centralityBin, np[iPid][0]);
	
	for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
	  delta *= deltaNp;
	  
	  (static_cast<TProfile*>(list->FindObject(Form("fProf%s%sNet%dM", fgkPidName[iPid], name.Data(), idxOrder))))->Fill(centralityBin, delta);
	}
	
	for (Int_t idxOrder = 0; idxOrder <= fOrder; ++ idxOrder) {
	  fRedFactp[idxOrder][0]  = 1.;
	  fRedFactp[idxOrder][1]  = 1.;
	}
	
	for (Int_t idxOrder = 1; idxOrder <= fOrder; ++ idxOrder) {
	  fRedFactp[idxOrder][0]  = fRedFactp[idxOrder-1][0]  * Double_t(np[iPid][0]-(idxOrder-1));
	  fRedFactp[idxOrder][1]  = fRedFactp[idxOrder-1][1]  * Double_t(np[iPid][1]-(idxOrder-1));
	}
	
	for (Int_t ii = 0; ii <= fOrder; ++ii) {  
	  for (Int_t kk = 0; kk <= fOrder; ++kk) { 
	    Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   
	    (static_cast<TProfile*>(list->FindObject(Form("fProf%s%sNetF%02d%02d", fgkPidName[iPid], name.Data(), ii, kk))))->Fill(centralityBin, fik);
	  }
	}
      } // zero
    } // pid
  } // nu
    
  //Printf("%6d %20s %6.2f %6d %6d %6d %6d  %6d %6d %6d %6d", idx, name, centralityBin,
  //	 np[0][1],  np[0][0], 
  //	 np[1][1],  np[1][0], 
  //	 np[2][1],  np[2][0], 
  //	 np[3][1],  np[3][0]);
  //

  if (isZero) {
    Double_t a[10][4]; Double_t b[18];
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
      a[0][iPid] = np[iPid][0];                   // 0  n-               
      a[1][iPid] = np[iPid][1];                   // 1  n+               
      
      a[2][iPid] = np[iPid][0] * np[iPid][0];     // 2  n-^2            
      a[3][iPid] = np[iPid][1] * np[iPid][1];     // 3  n+^2            
      
      a[4][iPid] = np[iPid][0] * (np[iPid][0]-1); // 4  n- . (n- - 1)
      a[5][iPid] = np[iPid][1] * (np[iPid][1]-1); // 5  n+ . (n+ - 1)
      a[6][iPid] = np[iPid][1] * np[iPid][0];     // 6  n+ . n-
      
      a[7][iPid] = np[iPid][1]+np[iPid][0];       // 7  n+ + n-          
      a[8][iPid] = TMath::Power((np[iPid][1]+np[iPid][0]),2); // 8  n+ + n-           
      a[9][iPid] = (np[iPid][1]+np[iPid][0]) * ((np[iPid][1]+np[iPid][0])-1); // 9  (N(N-1))
    }
    
    b[0]  = a[7][0]*a[7][2];       // 40 N   K
    b[1]  = a[7][1]*a[7][2];       // 41 Pi  K
    b[2]  = a[1][1]*a[1][2];       // 42 pi+ k+
    b[3]  = a[1][1]*a[0][2];       // 43 pi+ k-
    b[4]  = a[0][1]*a[1][2];       // 28 pi- k+  
    b[5]  = a[0][1]*a[0][2];       // 29 pi- k-
    
    b[6]  = a[7][0]*a[7][3];       // 30 N   P
    b[7]  = a[7][2]*a[7][3];       // 31 K   P
    b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
    b[9]  = a[1][2]*a[0][3];       // 33 k+  p-
    b[10] = a[0][2]*a[1][3];       // 34 k-  p+
    b[11] = a[0][2]*a[0][3];       // 35 k-  p-
    
    b[12] = a[7][0]*a[7][1];       // 36 N  Pi
    b[13] = a[7][3]*a[7][1];       // 37 P  Pi
    b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
    b[15] = a[1][3]*a[0][1];       // 39 p+ pi-
    b[16] = a[0][3]*a[1][1];       // 40 p- pi+
    b[17] = a[0][3]*a[0][1];       // 41 p- pi-
      
    Int_t k = 0;
    for (Int_t j = 0; j < 4; j++) {
      for (Int_t i = 0; i < 10; i++) {
	(static_cast<TProfile*>(list->FindObject(Form("fProf%sNu%02d", name.Data(),k))))->Fill(centralityBin,a[i][j]); 
	k++;
      }
    }
    
    for (Int_t j = 0; j < 18; j++) {
      (static_cast<TProfile*>(list->FindObject(Form("fProf%sNu%02d", name.Data(),j + 40))))->Fill(centralityBin,b[j]); 
    }
  }
  return;
}


//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::FillRatioHistos(Bool_t isMC,Bool_t isPer)  {
   
  Double_t np[4][2];
  if (isMC) {
    for (Int_t i = 0; i < 4; i++) {
      np[i][0] = fMCNp[i][0];
      np[i][1] = fMCNp[i][1];
    }
  } else {
    for (Int_t i = 0; i < 4; i++) {
      np[i][0] = fNp[i][0];
      np[i][1] = fNp[i][1];
    }
  }

  Bool_t isZero = kTRUE;  
  if ((np[0][0] == 0) || (np[0][1] == 0))
    isZero = kFALSE; 
  else if ((np[1][0] == 0) || (np[1][1] == 0))
    isZero = kFALSE;
  else if ((np[2][0] == 0) || (np[2][1] == 0))
    isZero = kFALSE;
  else if ((np[3][0] == 0) || (np[3][1] == 0))
    isZero = kFALSE;
  else isZero = kTRUE;
  


  TString nmc  = (isMC) ? "MC" : "Phy";
  TString name = (isPer) ? Form("%sPer",nmc.Data()) : Form("%sBin",nmc.Data());
  Float_t centralityBin = (isPer) ? (fCentralityPercentile + 1) : (fCentralityBin + 1);
  TList *list = static_cast<TList*>(fPhyList->FindObject(Form("fRatio%s",name.Data())));
   

  if (!isZero) return;

  if (fIsNu) {
    Double_t a[22]; Double_t b[4];
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
      a[iPid] = np[iPid][0]/np[iPid][1];
      b[iPid] = np[iPid][0] + np[iPid][1];
    }
    
    a[4] = b[1]/b[0];
    a[5] = b[2]/b[0];
    a[6] = b[2]/b[1];
    
    a[7] = b[3]/b[0];
    a[8] = b[3]/b[1];
    a[9] = b[3]/b[2];
    
    a[10] = np[2][0]/np[1][0]; // k-pi-
    a[11] = np[2][0]/np[1][1]; // k-pi+
    a[12] = np[2][1]/np[1][0]; // k+pi-
    a[13] = np[2][1]/np[1][1]; // k+pi+
    
    a[14] = np[3][0]/np[1][0]; // p-pi-
    a[15] = np[3][0]/np[1][1]; // p-pi+
    a[16] = np[3][1]/np[1][0]; // p+pi-
    a[17] = np[3][1]/np[1][1]; // p+pi+
    
    a[18] = np[3][0]/np[2][0]; // p-k-
    a[19] = np[3][0]/np[2][1]; // p-k+
    a[20] = np[3][1]/np[2][0]; // p+k-
    a[21] = np[3][1]/np[2][1]; // p+k+
    
    for (Int_t i = 0; i < 22; i++) {
      (static_cast<TH2F*>(list->FindObject(Form("fHist%sRatio%02d",name.Data(),i))))->Fill(centralityBin, a[i]);
    }
  }
  
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    (static_cast<TH2F*>(list->FindObject(Form("fHist%sDist%s%s",name.Data(), fgkPidName[iPid], fgkNetHistName[0]))))->Fill(centralityBin, np[iPid][1]+np[iPid][0]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHist%sDist%s%s",name.Data(), fgkPidName[iPid], fgkNetHistName[1]))))->Fill(centralityBin, np[iPid][1]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHist%sDist%s%s",name.Data(), fgkPidName[iPid], fgkNetHistName[2]))))->Fill(centralityBin,                  np[iPid][0]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHist%sDist%s%s",name.Data(), fgkPidName[iPid], fgkNetHistName[3]))))->Fill(centralityBin, np[iPid][1]-np[iPid][0]);   }

  return;
}

//________________________________________________________________________
void  AliEbyENetChargeFluctuationTask::CreateGroupHistos(const Char_t *name, const Char_t *title, Int_t nSample, Bool_t isPer)  {


  TString strTit(title);
  Float_t etaRange[2] = {(Float_t)fEtaMin, (Float_t)fEtaMax};
  
  //TList *list[4];
  fPhyList->Add(new TList);
  TList *list =  static_cast<TList*>(fPhyList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);
  
  TString tname    = Form("%s", name);
  Int_t nBinsCent  = (isPer) ? fNbwcBin : fCentralityBinMax;
  Double_t centBinRange[2];  
  centBinRange[0]  = 0.5;
  centBinRange[1]  = (isPer) ?  fNbwcBin + 0.5 : fCentralityBinMax+0.5;
 
  for (Int_t iSub = 0; iSub <= nSample; ++iSub) {
    
    list->Add(new TList);
    TList *listSub = static_cast<TList*>(list->Last());
    listSub->SetName(Form("%s%02d",name, iSub));
    listSub->SetOwner(kTRUE);
    if (!fIsNu) {
      for (Int_t iPid = 0; iPid < 4; ++iPid) {
	
	TString sNetTitle(Form("%s - %s", fgkPidLatex[iPid][1], fgkPidLatex[iPid][0]));
	strTit = (iPid != 0 ) ? Form("|y| < %.1f", fRapMax) : Form(" |#eta|<%.1f", etaRange[1]);
	
	listSub->Add(new TProfile(Form("fProfS%02dTot%sPlus%s", iSub, fgkPidName[iPid],tname.Data()), 
				  Form("(%s) : %s;Centrality(%d);(%s)",fgkPidName[iPid], 
				       strTit.Data(), nBinsCent,sNetTitle.Data()),
				  nBinsCent, centBinRange[0], centBinRange[1]));
	
	listSub->Add(new TProfile(Form("fProfS%02dTot%sMinus%s",iSub, fgkPidName[iPid],tname.Data()), 
				  Form("(%s) : %s;Centrality(%d);(%s)",fgkPidName[iPid], 
				       strTit.Data(),nBinsCent, sNetTitle.Data()),
				  nBinsCent, centBinRange[0], centBinRange[1]));
	
	
	
	for (Int_t idx = 1; idx <= fOrder; ++idx) {
	  listSub->Add(new TProfile(Form("fProfS%02d%s%sNet%dM",iSub, fgkPidName[iPid],
					 tname.Data(), idx), 
				    Form("(%s)^{%d} : %s;Centrality(%d);(%s)^{%d}", 
					 sNetTitle.Data(), idx, strTit.Data(), nBinsCent,
					 sNetTitle.Data(), idx),
				    nBinsCent, centBinRange[0], centBinRange[1]));
	}
	
	for (Int_t ii = 0; ii <= fOrder; ++ii) {
	  for (Int_t kk = 0; kk <= fOrder; ++kk) {
	    listSub->Add(new TProfile(Form("fProfS%02d%s%sNetF%02d%02d",iSub, fgkPidName[iPid], 
					   tname.Data(), ii, kk),
				      Form("f_{%02d%02d} : %s;Centrality(%d);f_{%02d%02d}", 
					   ii, kk, strTit.Data(),nBinsCent, ii, kk),
				      nBinsCent, centBinRange[0], centBinRange[1]));
	  }
	}
	
      }  // iPid
      //-------------------------------------------
    } 
    
    for (Int_t iPhy = 0; iPhy < 58; ++iPhy) { 
      listSub->Add(new TProfile(Form("fProfS%02d%sNu%02d",iSub,tname.Data(),iPhy),
				Form("Physics Variable for index %d | %s | Sub S%02d; Centrality(%d);",
				     iPhy,tname.Data(), iSub,nBinsCent),
				nBinsCent, centBinRange[0], centBinRange[1]));
    }
    
  }// isub
  
  return;
}

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::FillGroupHistos(const Char_t *name,Int_t iSub, Bool_t isMC,Bool_t isPer)  {

  Float_t centralityBin = (isPer) ? (fCentralityPercentile + 1) : (fCentralityBin + 1);
  TString tname = Form("%s", name);
  
  TList *list    = static_cast<TList*>(fPhyList->FindObject(Form("f%s",name)));
  TList *listSub = static_cast<TList*>(list->FindObject(Form("%s%02d",name, iSub)));
     
  Double_t np[4][2];
  if (isMC) {
    for (Int_t i = 0; i < 4; i++) {
      np[i][0] = fMCNp[i][0];
      np[i][1] = fMCNp[i][1];
    }
  } else {
    for (Int_t i = 0; i < 4; i++) {
      np[i][0] = fNp[i][0];
      np[i][1] = fNp[i][1];
    }
  }

  Bool_t isZero = kTRUE;
  Bool_t isZeroPid[4] = {1,1,1,1};
  if ((np[0][0] == 0) || (np[0][1] == 0))
    {isZero = kFALSE; isZeroPid[0] = kFALSE;}
  else if ((np[1][0] == 0) || (np[1][1] == 0))
    {isZero = kFALSE; isZeroPid[1] = kFALSE;} 
  else if ((np[2][0] == 0) || (np[2][1] == 0))
    {isZero = kFALSE; isZeroPid[2] = kFALSE;}
  else if ((np[3][0] == 0) || (np[3][1] == 0))
    {isZero = kFALSE; isZeroPid[3] = kFALSE;}
  else isZero = kTRUE;

  if (!fIsNu) {
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
      if (isZeroPid[iPid]) { 
	Int_t deltaNp = np[iPid][1]-np[iPid][0];  
	Double_t delta = 1.;
	
	(static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02dTot%sPlus%s",iSub, fgkPidName[iPid], tname.Data()))))->Fill(centralityBin, np[iPid][1]);
	(static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02dTot%sMinus%s",iSub, fgkPidName[iPid], tname.Data()))))->Fill(centralityBin, np[iPid][0]);
	
	for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
	  delta *= deltaNp;
	  (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%s%sNet%dM",iSub, fgkPidName[iPid], tname.Data(), idxOrder))))->Fill(centralityBin, delta);
	}
	
	for (Int_t idxOrder = 0; idxOrder <= fOrder; ++idxOrder) {
	  fRedFactp[idxOrder][0]  = 1.;
	  fRedFactp[idxOrder][1]  = 1.;
	}
	
	for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
	  fRedFactp[idxOrder][0]  = fRedFactp[idxOrder-1][0]  * Double_t(np[iPid][0]-(idxOrder-1));
	  fRedFactp[idxOrder][1]  = fRedFactp[idxOrder-1][1]  * Double_t(np[iPid][1]-(idxOrder-1));
	}
	
	for (Int_t ii = 0; ii <= fOrder; ++ii) {  
	  for (Int_t kk = 0; kk <= fOrder; ++kk) { 
	    Double_t fik = fRedFactp[ii][1] * fRedFactp[kk][0];   
	    (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%s%sNetF%02d%02d",iSub, fgkPidName[iPid], tname.Data(), ii, kk))))->Fill(centralityBin, fik);
	  }
	}
      } // zero
    } //pid
  } //nu
  // Printf("%6.2f   %6.2f %6.2f   %6.2f %6.2f  %6.2f %6.2f   %6.2f %6.2f", centralityBin,
  // 	 np[0][1],  np[0][0], 
  //	 np[1][1],  np[1][0], 
  //	 np[2][1],  np[2][0], 
  //	 np[3][1],  np[3][0]);
  //
  
  if (isZero) {
    
    Double_t a[10][4]; Double_t b[18];
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
      a[0][iPid] = np[iPid][0];                   // 0  n-               
      a[1][iPid] = np[iPid][1];                   // 1  n+               
      
      a[2][iPid] = np[iPid][0] * np[iPid][0];     // 2  n-^2            
      a[3][iPid] = np[iPid][1] * np[iPid][1];     // 3  n+^2            
      
      a[4][iPid] = np[iPid][0] * (np[iPid][0]-1); // 4  n- . (n- - 1)
      a[5][iPid] = np[iPid][1] * (np[iPid][1]-1); // 5  n+ . (n+ - 1)
      a[6][iPid] = np[iPid][1] * np[iPid][0];     // 6  n+ . n-
      
      a[7][iPid] = np[iPid][1]+np[iPid][0];       // 7  n+ + n-          
      a[8][iPid] = TMath::Power((np[iPid][1]+np[iPid][0]),2); // 8  n+ + n-           
      a[9][iPid] = (np[iPid][1]+np[iPid][0]) * ((np[iPid][1]+np[iPid][0])-1); // 9  (N(N-1))
    }
    
    b[0]  = a[7][0]*a[7][2];       // 40 N   K
    b[1]  = a[7][1]*a[7][2];       // 41 Pi  K
    b[2]  = a[1][1]*a[1][2];       // 42 pi+ k+
    b[3]  = a[1][1]*a[0][2];       // 43 pi+ k-
    b[4]  = a[0][1]*a[1][2];       // 28 pi- k+  
    b[5]  = a[0][1]*a[0][2];       // 29 pi- k-
    
    b[6]  = a[7][0]*a[7][3];       // 30 N   P
    b[7]  = a[7][2]*a[7][3];       // 31 K   P
    b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
    b[9]  = a[1][2]*a[0][3];       // 33 k+  p-
    b[10] = a[0][2]*a[1][3];       // 34 k-  p+
    b[11] = a[0][2]*a[0][3];       // 35 k-  p-
    
    b[12] = a[7][0]*a[7][1];       // 36 N  Pi
    b[13] = a[7][3]*a[7][1];       // 37 P  Pi
    b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
    b[15] = a[1][3]*a[0][1];       // 39 p+ pi-
    b[16] = a[0][3]*a[1][1];       // 40 p- pi+
    b[17] = a[0][3]*a[0][1];       // 41 p- pi-

    // TList *list_nu = static_cast<TList*>(fOutlistSub->FindObject(Form("f%s_nu",name)));
    Int_t k = 0;
    for (Int_t j = 0; j < 4; j++) {
      for (Int_t i = 0; i < 10; i++) {
	(static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%sNu%02d",iSub, tname.Data(),k))))->Fill(centralityBin,a[i][j]); 
	k++;
      }
    }
    
    for (Int_t j = 0; j < 18; j++) {
      (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%sNu%02d",iSub, tname.Data(),j + 40))))->Fill(centralityBin,b[j]); 
    }
  }
  return;
}




//________________________________________________________________________
Int_t AliEbyENetChargeFluctuationTask::SetupEventCR(AliESDInputHandler *esdHandler, 
						   AliAODInputHandler *aodHandler, 
						   AliMCEvent *mcEvent) {
  if(esdHandler){
    fInputEventHandler = static_cast<AliInputEventHandler*>(esdHandler);
    fESD               = dynamic_cast<AliESDEvent*>(fInputEventHandler->GetEvent());
    if (!fESD) {
      AliError("ESD event handler not available");
      return -1;
    }
  }
  
  else if(aodHandler){
    fInputEventHandler = static_cast<AliInputEventHandler*>(aodHandler);
    fAOD               = dynamic_cast<AliAODEvent*>(fInputEventHandler->GetEvent());
    if (!fAOD) {
      AliError("AOD event handler not available");
      return -1;
    }
  }
  
  //  fPIDResponse = fInputEventHandler->GetPIDResponse();
  fMCEvent     = mcEvent;
  if (fMCEvent)
    fStack     = fMCEvent->Stack();
  
  AliCentrality *centrality = NULL;

  if(esdHandler)
    centrality = fESD->GetCentrality();
  else if(aodHandler)
    centrality = ((AliVAODHeader*)fAOD->GetHeader())->GetCentralityP();

  if (!centrality) {
    AliError("Centrality not available");
    return -1;
  }
  
  if(fIsTen) {
    Int_t centBin = centrality->GetCentralityClass10(fCentralityEstimator.Data());
    if (centBin == 10 || centBin == -1.) { fCentralityBin = -1; }
    else if (centBin >= 0 && centBin < fNCentralityBins) { fCentralityBin = centBin;}
    else {  fCentralityBin = -2; }
  } else {
    Int_t centBin = centrality->GetCentralityClass10(fCentralityEstimator.Data());
    if (centBin == 0) { fCentralityBin = centrality->GetCentralityClass5(fCentralityEstimator.Data()); }
    else if (centBin == 11 || centBin == -1.)           { fCentralityBin = -1; }
    else if (centBin > 0 && centBin < fNCentralityBins) { fCentralityBin = centBin + 1; }
    else {  fCentralityBin = -2; }
  }
  
  if (fCentralityBin >= fCentralityBinMax) fCentralityBin = -2;
  
  Float_t bf =  100./Float_t(fNbwcBin);
  Int_t ai   = Int_t((centrality->GetCentralityPercentile(fCentralityEstimator.Data()))/bf);
  fCentralityPercentile = Double_t(ai);
  
  return 0;
}
//________________________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::TriggeredEvents() {
  Bool_t *aTriggerFired = new Bool_t[fNTriggers];
  for (Int_t ii = 0; ii < fNTriggers; ++ii)
    aTriggerFired[ii] = kFALSE;

  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kMB))          aTriggerFired[0] = kTRUE;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kCentral))     aTriggerFired[1] = kTRUE;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kSemiCentral)) aTriggerFired[2] = kTRUE;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kEMCEJE))      aTriggerFired[3] = kTRUE;
  if ((fInputEventHandler->IsEventSelected() & AliVEvent::kEMCEGA))      aTriggerFired[4] = kTRUE;

  Bool_t isTriggered = kFALSE;

  for (Int_t ii=0; ii<fNTriggers; ++ii) {
    if(aTriggerFired[ii]) {
      if (fNeedQa) (static_cast<TH1F*>(fQaList->FindObject(Form("hTriggerStat"))))->Fill(ii);
    }
  }
  // Fix for only selected events || Histo to keep what is the ratio
  if ((fInputEventHandler->IsEventSelected() & fSelectBit))  isTriggered = kTRUE; 
    delete[] aTriggerFired;
  return isTriggered;
}

//________________________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::RejectedEvent() {
  Int_t *aEventCuts = new Int_t[fHEventStatMax];
    
  for (Int_t ii=0;ii<fHEventStatMax; ++ii)
    aEventCuts[ii] = 0;

  Int_t iStep = 0;
  aEventCuts[iStep] = 0;  //0

  ++iStep;
  if (!TriggeredEvents()) aEventCuts[iStep] = 1; //1

  ++iStep;
  const AliESDVertex* vtxESD = NULL;
  const AliAODVertex* vtxAOD = NULL;
  if (fESD){
    vtxESD = fESD->GetPrimaryVertexTracks();
    if (!vtxESD) aEventCuts[iStep] = 1;   //2
  }
  else if (fAOD){
    vtxAOD = fAOD->GetPrimaryVertex();
    if (!vtxAOD) aEventCuts[iStep] = 1;  //2
  }

  ++iStep;
  if (vtxESD){
    if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvx"))))->Fill(fCentralityPercentile,vtxESD->GetX());
    if(TMath::Abs(vtxESD->GetX()) > fVxMax)  aEventCuts[iStep] = 1;  //3
    else { 
      if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvxA"))))->Fill(fCentralityPercentile,vtxESD->GetX());
    }
  }
  else if(vtxAOD){
    if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvx"))))->Fill(fCentralityPercentile,vtxAOD->GetX());
    if(TMath::Abs(vtxAOD->GetX()) > fVxMax) aEventCuts[iStep] = 1;    //3
    else { 
      if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvxA"))))->Fill(fCentralityPercentile,vtxAOD->GetX());
    }
  }
  else aEventCuts[iStep] = 1; //3
  
  ++iStep;
  if (vtxESD){
   if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvy"))))->Fill(fCentralityPercentile,vtxESD->GetY());
    if(TMath::Abs(vtxESD->GetY()) > fVyMax) aEventCuts[iStep] = 1; //4
    else { 
      if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvyA"))))->Fill(fCentralityPercentile,vtxESD->GetY());
    }  
  }
  else if(vtxAOD){
    if (fNeedQa) if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvy"))))->Fill(fCentralityPercentile,vtxAOD->GetY());
    if(TMath::Abs(vtxAOD->GetY()) > fVyMax) aEventCuts[iStep] = 1; //4
    else {
      if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvyA"))))->Fill(fCentralityPercentile,vtxAOD->GetY());
    }
  }
  else aEventCuts[iStep] = 1; //4


 ++iStep;
  if (vtxESD){
    if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvz"))))->Fill(fCentralityPercentile,vtxESD->GetZ());
    if(TMath::Abs(vtxESD->GetZ()) > fVzMax) aEventCuts[iStep] = 1;  //5
    else {
      if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvzA"))))->Fill(fCentralityPercentile,vtxESD->GetZ());
    }
  }
  else if(vtxAOD){
    if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvz"))))->Fill(fCentralityPercentile,vtxAOD->GetZ());
    if(TMath::Abs(vtxAOD->GetZ()) > fVzMax) aEventCuts[iStep] = 1; //5
    else {
      if (fNeedQa) (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvzA"))))->Fill(fCentralityPercentile,vtxAOD->GetZ());
    }
  }
  else aEventCuts[iStep] = 1; //5


  
  ++iStep;
  if(fCentralityBin == -1.) aEventCuts[iStep] = 1;  //6

  ++iStep;
  if(fCentralityBin == -2.)  aEventCuts[iStep] = 1; //7

  /*for (Int_t ii=0;ii<fHEventStatMax; ++ii)
    printf("%d ", aEventCuts[ii]);
  printf("\n");
  */

  Bool_t isRejected = IsEventStats(aEventCuts);

  

  delete[] aEventCuts;
  return isRejected;
}

//________________________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::IsEventStats(Int_t *aEventCuts) {
  Bool_t isRejected = kFALSE;
  for (Int_t idx = 0; idx < fHEventStatMax ; ++idx) {
    if (aEventCuts[idx])
      isRejected = kTRUE;
    else {
      if (fNeedQa) (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat0"))))->Fill(idx);
    }
  }
  for (Int_t idx = 0; idx < fHEventStatMax; ++idx) {
    if (aEventCuts[idx])
      break;
    if (fNeedQa) (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat1"))))->Fill(idx);
  }
  if (!isRejected) {

    if (fNeedQa) (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityStat"))))->Fill(fCentralityBin+1);
    if (fNeedQa) (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityPercentileAccepted"))))->Fill(fCentralityPercentile+1);
  }
  if (fNeedQa) (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityPercentileAll"))))->Fill(fCentralityPercentile+1);
  return isRejected;
}




void AliEbyENetChargeFluctuationTask::FillCC(Int_t i) {
  if      (i == 0) static_cast<THnSparseF*>(fEffList->FindObject("hmNchContMc"))->Fill(fCurCont);
  else if (i == 1) static_cast<THnSparseF*>(fEffList->FindObject("hmNpiContMc"))->Fill(fCurCont);
  else if (i == 2) static_cast<THnSparseF*>(fEffList->FindObject("hmNkaContMc"))->Fill(fCurCont);
  else if (i == 3) static_cast<THnSparseF*>(fEffList->FindObject("hmNprContMc"))->Fill(fCurCont);
}
void AliEbyENetChargeFluctuationTask::FillRCC(Int_t i) {
  if      (i == 0) static_cast<THnSparseF*>(fEffList->FindObject("hmNchContRec"))->Fill(fCurRec);
  else if (i == 1) static_cast<THnSparseF*>(fEffList->FindObject("hmNpiContRec"))->Fill(fCurRec);
  else if (i == 2) static_cast<THnSparseF*>(fEffList->FindObject("hmNkaContRec"))->Fill(fCurRec);
  else if (i == 3) static_cast<THnSparseF*>(fEffList->FindObject("hmNprContRec"))->Fill(fCurRec);
}
void AliEbyENetChargeFluctuationTask::FillGenCE(Int_t i) {

  // Printf("%d %f %f %f %f %f %f %f %f",i,fCurGen[0],fCurGen[1],fCurGen[2],fCurGen[3],fCurGen[4],fCurGen[5],fCurGen[6],fCurGen[7]);


  if      (i == 0) static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffMc"))->Fill(fCurGen);
  else if (i == 1) static_cast<THnSparseF*>(fEffList->FindObject("hmNpiEffMc"))->Fill(fCurGen);
  else if (i == 2) static_cast<THnSparseF*>(fEffList->FindObject("hmNkaEffMc"))->Fill(fCurGen);
  else if (i == 3) static_cast<THnSparseF*>(fEffList->FindObject("hmNprEffMc"))->Fill(fCurGen);
}
void AliEbyENetChargeFluctuationTask::FillRecCE(Int_t i) {
  if      (i == 0) static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->Fill(fCurRec);
  else if (i == 1) static_cast<THnSparseF*>(fEffList->FindObject("hmNpiEffRec"))->Fill(fCurRec);
  else if (i == 2) static_cast<THnSparseF*>(fEffList->FindObject("hmNkaEffRec"))->Fill(fCurRec);
  else if (i == 3) static_cast<THnSparseF*>(fEffList->FindObject("hmNprEffRec"))->Fill(fCurRec);
}

void AliEbyENetChargeFluctuationTask::FillRecDE(Int_t i) {
  if      (i == 0) static_cast<THnSparseF*>(fDcaList->FindObject("hmNchDca"))->Fill(fCurGenD);
  else if (i == 1) static_cast<THnSparseF*>(fDcaList->FindObject("hmNpiDca"))->Fill(fCurGenD);
  else if (i == 2) static_cast<THnSparseF*>(fDcaList->FindObject("hmNkaDca"))->Fill(fCurGenD);
  else if (i == 3) static_cast<THnSparseF*>(fDcaList->FindObject("hmNprDca"))->Fill(fCurGenD);
}

void AliEbyENetChargeFluctuationTask::FillRecDED(Int_t i) {
  // Printf("%d %f %f %f %f %f %f %f",i,fCurGen[0],fCurRecD[1],fCurRecD[2],fCurRecD[3],fCurRecD[4],fCurRecD[5],fCurRecD[6]);
  if      (i == 0) static_cast<THnSparseF*>(fDcaList->FindObject("hmNchDcaRec"))->Fill(fCurRecD);
  else if (i == 1) static_cast<THnSparseF*>(fDcaList->FindObject("hmNpiDcaRec"))->Fill(fCurRecD);
  else if (i == 2) static_cast<THnSparseF*>(fDcaList->FindObject("hmNkaDcaRec"))->Fill(fCurRecD);
  else if (i == 3) static_cast<THnSparseF*>(fDcaList->FindObject("hmNprDcaRec"))->Fill(fCurRecD);
}

void AliEbyENetChargeFluctuationTask::CreateCE() {

  Int_t     bhepmc[8] = {fCentralityBinMax,fGNBinsSign,2, 2, 2,fGNBinsRap,fGNBinsPhi,fGNBinsPt};
  Double_t mnhepmc[8] = {0.5,fGRngSign[0],-0.5,-0.5,-0.5,fGRngRap[0],fGRngPhi[0],fGRngPt[0]};  
  Double_t mxhepmc[8] = {fCentralityBinMax+0.5,fGRngSign[1],1.5,1.5,1.5,fGRngRap[1],fGRngPhi[1],fGRngPt[1]};  
  
  TString titilemc        = "cent:signMC:findable:recStatus:pidStatus:yMC:phiMC:ptMC";
  TString tiltlelaxmc[8]  = {"Centrality", "sign", "findable","recStatus","recPid","#it{y}_{MC}", 
			     "#varphi_{MC} (rad)","#it{p}_{T,MC} (GeV/#it{c})"};
  
  fEffList->Add(new THnSparseF("hmNpiEffMc",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fEffList->Add(new THnSparseF("hmNkaEffMc",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fEffList->Add(new THnSparseF("hmNprEffMc",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));

  bhepmc[5]  = fGNBinsEta;
  mnhepmc[5] = fGRngEta[0];
  mxhepmc[5] = fGRngEta[1];
  titilemc        = "cent:signMC:findable:recStatus:pidStatus:etaMC:phiMC:ptMC";
  
  fEffList->Add(new THnSparseF("hmNchEffMc",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  
 for (Int_t i = 0; i < 8; i++) { 
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffMc"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
   static_cast<THnSparseF*>(fEffList->FindObject("hmNpiEffMc"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
   static_cast<THnSparseF*>(fEffList->FindObject("hmNkaEffMc"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
   if (i == 5) tiltlelaxmc[5] = "#it{#eta}_{Rec}";
   static_cast<THnSparseF*>(fEffList->FindObject("hmNprEffMc"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());

  }
  
  Int_t    binhnep[5] = {fCentralityBinMax, fGNBinsSign, fGNBinsRap, fGNBinsPhi, fGNBinsPt};
  Double_t minhnep[5] = {0.5,fGRngSign[0],fGRngRap[0],fGRngPhi[0],fGRngPt[0]};
  Double_t maxhnep[5] = {fCentralityBinMax+0.5,fGRngSign[1],fGRngRap[1],fGRngPhi[1],fGRngPt[1]};


  TString titilerec        = "cent:signRec:yRec:phiRec:ptRec";
  TString tiltlelaxrec[5]  = {"Centrality", "sign", "#it{y}_{Rec}", "#varphi_{Rec} (rad)","#it{p}_{T,Rec} (GeV/#it{c})"};

  fEffList->Add(new THnSparseF("hmNpiEffRec",titilerec.Data(),5,binhnep, minhnep, maxhnep));
  fEffList->Add(new THnSparseF("hmNkaEffRec",titilerec.Data(),5,binhnep, minhnep, maxhnep));
  fEffList->Add(new THnSparseF("hmNprEffRec",titilerec.Data(),5,binhnep, minhnep, maxhnep));

  fEffList->Add(new THnSparseF("hmNpiContRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));
  fEffList->Add(new THnSparseF("hmNkaContRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));
  fEffList->Add(new THnSparseF("hmNprContRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));
 
  binhnep[2] = fGNBinsEta;
  minhnep[2] = fGRngEta[0];
  maxhnep[2] = fGRngEta[1];

  titilerec        = "cent:signRec:etaRec:phiRec:ptRec";
    
  fEffList->Add(new THnSparseF("hmNchEffRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));
  fEffList->Add(new THnSparseF("hmNchContRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));

  for (Int_t i = 0; i < 5; i++) { 
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    
    static_cast<THnSparseF*>(fEffList->FindObject("hmNpiContRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNkaContRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
 static_cast<THnSparseF*>(fEffList->FindObject("hmNprContRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    if (i == 2) tiltlelaxrec[2] = "#it{#eta}_{Rec}";
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
static_cast<THnSparseF*>(fEffList->FindObject("hmNchContRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
   
  }  

  //----
  Int_t    binHnCont[6] = {fCentralityBinMax,fGNBinsSign, 5,fGNBinsRap,fGNBinsPhi, fGNBinsPt};  
  Double_t minHnCont[6] = {0.5,fGRngSign[0],-0.5,fGRngRap[0],fGRngPhi[0], fGRngPt[0]};
  Double_t maxHnCont[6] = {fCentralityBinMax+0.5,fGRngSign[1],4.5,fGRngRap[1],fGRngPhi[1], fGRngPt[1]};   

  TString titilecont     = "cent:signMC:contStatus:yMC:phiMC:ptMC";
  TString tiltlelaxcont[6] 
    = {"Centrality","sign","contPart","#it{y}_{MC}","#varphi_{MC} (rad)","#it{p}_{T,MC} (GeV/#it{c})"};
  
  fEffList->Add(new THnSparseF("hmNpiContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont));
  fEffList->Add(new THnSparseF("hmNkaContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont));
  fEffList->Add(new THnSparseF("hmNprContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont));
 
  binHnCont[3] = fGNBinsEta;
  minHnCont[3] = fGRngEta[0];
  maxHnCont[3] = fGRngEta[1];
  titilecont     = "cent:signMC:contStatus:etaMC:phiMC:ptMC";
  fEffList->Add(new THnSparseF("hmNchContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont));

 for (Int_t i = 0; i < 6; i++) {  
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchContMc"))->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNpiContMc"))->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNkaContMc"))->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    if (i == 3) tiltlelaxcont[3] = "#it{#eta}_{Rec}";
    static_cast<THnSparseF*>(fEffList->FindObject("hmNprContMc"))->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
  }


  return;
}

//_____________________________________________________________________________
Double_t * AliEbyENetChargeFluctuationTask::CreateLogAxis(Int_t nbins, Double_t xmin, Double_t xmax) {
  // retun pointer to the array with log axis
  // it is user responsibility to delete the array

  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  Double_t *xbins =  new Double_t[nbins+1];

  xbins[0] = xmin;
  for (Int_t i=1;i<=nbins;i++) {
    xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
  }
  return xbins;
}
void AliEbyENetChargeFluctuationTask::CalEC() {for (Int_t i = 0; i < 4; i++) CalculateCE(i);}
void AliEbyENetChargeFluctuationTask::CalED() {
  if (fIsMC) for (Int_t i = 0; i < 4; i++) CalculateDE(i);
  for (Int_t i = 0; i < 4; i++) CalculateDED(i);
}
//_____________________________________________________________________________
void AliEbyENetChargeFluctuationTask::CalculateCE(Int_t gPid) {
  //  Printf("Step %d",gPid);
  // Loop over ESD/AOD

 Int_t *gLabelsZer = new Int_t[fNTracks];
 Int_t *gLabelsOne = new Int_t[fNTracks];

 for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
   AliVTrack *track = (fESD) ? 
     static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : 
     static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

   // DO track cut.
   if(!AcceptTrackL(track)) continue;
  

   fCurCont[1] = track->Charge() < 0 ? 0 : 1;
   
   Int_t b = 0;

    if (gPid != 0) {
      Int_t a = fHelperPID->GetParticleSpecies(track,kFALSE);
      if(a < 0 || a > 2)  continue;
	b = a + 1;
    }

    if (gPid != b) continue;

    Int_t label  = TMath::Abs(track->GetLabel()); 
    gLabelsZer[idxTrack] = label;    
    
    Double_t rap;
    if (!TrackRapidity(track,rap,b)) continue;
    gLabelsOne[idxTrack] = label;    
    
    AliVParticle* particle = (fESD) ? 
      static_cast<AliVParticle*>(fMCEvent->GetTrack(label)) : 
      static_cast<AliVParticle*>(fArrayMC->At(label));
    if (!particle)
      return;
        
    Bool_t isPhysicalPrimary = (fESD) ? 
      fStack->IsPhysicalPrimary(label):
      (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();

    Bool_t isSecondaryFromWeakDecay = (fESD) ? 
      fStack->IsSecondaryFromWeakDecay(label) : 
      (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();

    Bool_t isSecondaryFromMaterial  = (fESD) ? 
      fStack->IsSecondaryFromMaterial(label)  : 
      (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();
    
    
    if (isPhysicalPrimary) {
      if (b != 0) {
	// -- Check if correctly identified 
	if (particle->PdgCode() == (track->Charge()*GetPDG(b)))
	  fCurCont[2] = 1.;   
	// -- MissIdentification
	else 
	  fCurCont[2] = 2.;
      }
      else
	fCurCont[2] = 1.;   
    } else if(isSecondaryFromWeakDecay)
      fCurCont[2] = 3.;
    else if (isSecondaryFromMaterial)
      fCurCont[2] = 4.;
    else
      fCurCont[2] = -1.;
    
    fCurCont[0] = fCentralityBin+1;
    fCurCont[3] = particle->Y();
    fCurCont[4] = particle->Phi();
    fCurCont[5] = particle->Pt();
    FillCC(b);
    fCurRec[0] = fCentralityBin+1;
    fCurRec[1] = track->Charge();
    fCurRec[2] = rap;
    fCurRec[3] = track->Phi();
    fCurRec[4] = track->Pt();
    FillRCC(b);
    // Printf("%d %d %f %f %f %f %f",gPid,b,fCurRec[0],fCurRec[1],fCurRec[2],fCurRec[3],fCurRec[4]);
 } 

 // Loop Over MC

 Int_t nPart  = (fESD) ? 
   fStack->GetNprimary() : 
   fArrayMC->GetEntriesFast();
 
 for (Int_t idxMC = 0; idxMC < nPart; ++idxMC) {
   AliVParticle* particle = (fESD) ? 
     static_cast<AliVParticle*>(fMCEvent->GetTrack(idxMC)) : 
     static_cast<AliVParticle*>(fArrayMC->At(idxMC));
   
   if (!AcceptTrackLMC(particle, idxMC)) continue;
   
   fCurGen[1] = (particle->PdgCode() < 0) ? -1. : 1.;
      
   Int_t iPid = 0;  
   if (gPid != 0) {
     if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
     else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
     else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
     else    iPid = 0;
   }
  
   if (gPid != iPid) continue;

   Double_t rapmc;
   if (!ParticleRapidity(particle, rapmc, iPid)) continue;

   fCurGen[5] = rapmc;

   fCurGen[2] = (fESD) ? Double_t(IsFindableInTPC(idxMC)) : 1.;
    
   for (Int_t idxRec=0; idxRec < fNTracks; ++idxRec) {
      if (idxMC == gLabelsZer[idxRec]) {
	fCurGen[3] = 1.;
	
	if (idxMC == gLabelsOne[idxRec])
	  fCurGen[4] = 1.;
	
        AliVTrack *track = NULL;
        if(fESD)
          track = fESD->GetTrack(idxRec);
        else if(fAOD)
          track = fAOD->GetTrack(idxRec);
	
        if (track) {
	  Double_t rap;
	  TrackRapidity(track,rap,iPid);
	  fCurRec[0] = fCentralityBin+1;
	  fCurRec[1] = track->Charge();
	  fCurRec[2] = rap;
	  fCurRec[3] = track->Phi();
	  fCurRec[4] = track->Pt();
	  FillRecCE(iPid);
	}      
        break;
      }
    } 
   fCurGen[0] = fCentralityBin+1;
   fCurGen[6] = particle->Phi();
   fCurGen[7] = particle->Pt();

  
   FillGenCE(iPid);
  } 
  

 delete [] gLabelsZer;
 delete [] gLabelsOne;


}

void AliEbyENetChargeFluctuationTask::CreateDED() {

 Int_t      bhepmc[7] = {fCentralityBinMax, fGNBinsSign,   2,   fGNBinsRap, fGNBinsPhi, fGNBinsPt,   140};
  Double_t mnhepmc[7] = {0.5,fGRngSign[0],-0.5, fGRngRap[0],fGRngPhi[0],fGRngPt[0], -3.5};  
  Double_t mxhepmc[7] = {fCentralityBinMax+0.5,fGRngSign[1], 1.5, fGRngRap[1],fGRngPhi[1], fGRngPt[1], 3.5};  
  TString titilemc        = "cent:sign:accepted:y:phi:pt:dcar";

  TString tiltlelaxmc[7]  = {"Centrality", "sign", "Is Accepted","#it{y}","#varphi (rad)","#it{p}_{T} (GeV/#it{c})", "DCAr"};
  
  fDcaList->Add(new THnSparseF("hmNpiDcaRec",titilemc.Data(),7,bhepmc,mnhepmc, mxhepmc));
  fDcaList->Add(new THnSparseF("hmNkaDcaRec",titilemc.Data(),7,bhepmc,mnhepmc, mxhepmc));
  fDcaList->Add(new THnSparseF("hmNprDcaRec",titilemc.Data(),7,bhepmc,mnhepmc, mxhepmc));
  bhepmc[3]  = fGNBinsEta;
  mnhepmc[3] = fGRngEta[0];
  mxhepmc[3] = fGRngEta[1];
  titilemc        = "cent:sign:accepted:eta:phi:pt:dcar";
  fDcaList->Add(new THnSparseF("hmNchDcaRec",titilemc.Data(),7,bhepmc,mnhepmc, mxhepmc));
   
  for (Int_t i = 0; i < 7; i++) { 
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNpiDcaRec"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNkaDcaRec"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNprDcaRec"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    if (i == 4) tiltlelaxmc[3] = "#eta";
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNchDcaRec"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
  }
 
}
void AliEbyENetChargeFluctuationTask::CreateDEM() {

  Int_t     bhepmc[8] = {fCentralityBinMax, fGNBinsSign,  2,   3,   fGNBinsRap, fGNBinsPhi, fGNBinsPt,  140};
  Double_t mnhepmc[8] = {0.5,fGRngSign[0],-0.5, 0.5, fGRngRap[0],fGRngPhi[0],fGRngPt[0],-3.5};  
  Double_t mxhepmc[8] = {fCentralityBinMax+0.5,fGRngSign[1], 1.5, 3.5, fGRngRap[1],fGRngPhi[1],fGRngPt[1], 3.5};  
  TString titilemc    = "cent:sign:cont:accepted:y:phi:pt:dcar";
  TString tiltlelaxmc[8]  = {"Centrality", "sign", "Is Accepted", "1 primary | 2 from WeakDecay | 3 p from Material",
			     "#it{y}","#varphi (rad)","#it{p}_{T} (GeV/#it{c})", "DCAr"};
  
  fDcaList->Add(new THnSparseF("hmNpiDca",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fDcaList->Add(new THnSparseF("hmNkaDca",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fDcaList->Add(new THnSparseF("hmNprDca",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));

  bhepmc[4]  = fGNBinsEta;
  mnhepmc[4] = fGRngEta[0];
  mxhepmc[4] = fGRngEta[1];

  titilemc        = "cent:sign:accepted:cont:eta:phi:pt:dcar";
  fDcaList->Add(new THnSparseF("hmNchDca",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
   
  for (Int_t i = 0; i < 8; i++) { 
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNpiDca"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNkaDca"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNprDca"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    if (i == 4) tiltlelaxmc[4] = "#eta";
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNchDca"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
  }
}

void AliEbyENetChargeFluctuationTask::CalculateDED(Int_t gPid) {

  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliVTrack *track = (fESD) ? 
      static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : 
      static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 
   
    if(!AcceptTrackL(track)) continue; 
     
    fCurRecD[1] = track->Charge() < 0 ? 0. : 1.;      // 1
    
    Int_t b = 0;
    
    if (gPid != 0) {
      Int_t a = fHelperPID->GetParticleSpecies(track,kFALSE);
      if(a < 0 || a > 2)  continue;
      b = a + 1;
    }
    if (gPid != b) continue;
    Double_t rap;
    if (!TrackRapidity(track,rap,b)) continue;
   
    fCurRecD[3] = rap;                                 // 3 
    Float_t dca[2], cov[3]; // 
    if (fESD)
      (dynamic_cast<AliESDtrack*>(track))->GetImpactParameters(dca, cov);
    else  {
      Double_t dcaa[2] = {-999,-999};
      Double_t cova[3] = {-999,-999,-999};
      AliAODTrack* clone =dynamic_cast<AliAODTrack*>(track->Clone("trk_clone"));
      Bool_t propagate = clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),100.,dcaa,cova);
      delete clone;   
      if (!propagate) dca[0] = -999; 
      else dca[0] = Float_t(dcaa[0]);
    }

    if ( TMath::Abs(dca[0]) <= fDcaXy ) fCurRecD[2] = 1.;   // 2
    else fCurRecD[2] = 0.;
    
    fCurRecD[0] = fCentralityBin+1;                // 0
    fCurRecD[4] = track->Phi();                  // 4
    fCurRecD[5] = track->Pt();                   // 5
    fCurRecD[6] = Double_t(dca[0]);              // 6
    //Printf("%d %f %f %f %f %f %f %f",gPid,fCurRecD[0],fCurRecD[1],fCurRecD[2],fCurRecD[3],fCurRecD[4],fCurRecD[5],fCurRecD[6]);
    FillRecDED(gPid);
    
  }
}    

void AliEbyENetChargeFluctuationTask::CalculateDE(Int_t gPid) {

  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliVTrack *track = (fESD) ? 
      static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : 
      static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 
    

    if(!AcceptTrackL(track)) continue;
    fCurGenD[1] = track->Charge() < 0 ? 0. : 1.;            // 1
    
    Int_t b = 0;
    
    if (gPid != 0) {
      Int_t a = fHelperPID->GetParticleSpecies(track,kFALSE);
      if(a < 0 || a > 2)  continue;
      b = a + 1;
    }
    
    if (gPid != b) continue;
    
    Double_t rap;
    if (!TrackRapidity(track,rap,b)) continue;
    
    fCurGenD[4] = rap;                                   // 4
    
    Float_t dca[2], cov[3]; // 
    if (fESD)
      (dynamic_cast<AliESDtrack*>(track))->GetImpactParameters(dca, cov);
    else { 
      Double_t dcaa[2] = {-999,-999};
      Double_t cova[3] = {-999,-999,-999};
      AliAODTrack* clone =dynamic_cast<AliAODTrack*>(track->Clone("trk_clone"));
      Bool_t propagate = clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),100.,dcaa,cova);
      delete clone; 
      if (!propagate) dca[0] = -999; 
      else dca[0] = Float_t(dcaa[0]);
    }
    
    if ( TMath::Abs(dca[0]) <= fDcaXy ) fCurGenD[2] = 1.;   // 2
    else fCurGenD[2] = 0.;
    
    fCurGenD[7] = Double_t(dca[0]);                    // 7
    
    Int_t label  = TMath::Abs(track->GetLabel()); 
    
    AliVParticle* particle = (fESD) ? 
      static_cast<AliVParticle*>(fMCEvent->GetTrack(label)) : 
      static_cast<AliVParticle*>(fArrayMC->At(label));
    if (!particle)
      return;
    
    Bool_t isPhysicalPrimary = (fESD) ? 
      fStack->IsPhysicalPrimary(label):
      (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
    
    Bool_t isSecondaryFromWeakDecay = (fESD) ? 
      fStack->IsSecondaryFromWeakDecay(label) : 
      (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
    
    Bool_t isSecondaryFromMaterial  = (fESD) ? 
      fStack->IsSecondaryFromMaterial(label)  : 
      (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();
    
    
    if (isPhysicalPrimary) fCurGenD[3] = 1.;           // 3
    else if(isSecondaryFromWeakDecay) fCurGenD[3] = 2.;
    else if (isSecondaryFromMaterial) fCurGenD[3] = 3.;
    else fCurGenD[3] = -1.;
    fCurGenD[0] = fCentralityBin+1;                    // 0
    fCurGenD[5] = track->Phi();                      // 5
    fCurGenD[6] = track->Pt();                       // 6
    
    FillRecDE(gPid);
  }
}
