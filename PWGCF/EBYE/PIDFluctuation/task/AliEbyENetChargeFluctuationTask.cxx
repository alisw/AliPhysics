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
//            Version 3:  pp Added with AliVTrack   (2/12/2014) - fixme    //    
//=========================================================================//

#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
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
  fCentralityBin(-1.),    
  fCentralityPercentile(-1.),
// fNp(NULL),
// fMCNp(NULL),
// fRedFactp(NULL),           
   
  fMinTrackLengthMC(80), 
  fAODtrackCutBit(768),   
  fNSubSamples(10),     
  fSubSampleIdx(0),
  fOrder(8),
  fNTriggers(5),
  fHEventStatMax(8),
  fNCentralityBins(11),
  fCentralityBinMax(11),
  fNTracks(0),

  fIsMC(kFALSE),
  fIsRatio(kFALSE),
  fIsAOD(kFALSE),
  fIsSub(kFALSE),
  fIsBS(kFALSE),
  fIsPer(kFALSE),
  fIsEff(kFALSE),
  fDebug(kFALSE),
  fIsQa(kFALSE),
  fIsPhy(kFALSE),
  fIsDca(kFALSE),
  
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
const Int_t   fGNBinsCent = 11 ;
const Float_t fGRngCent[] = {-0.5, 10.5};
const Float_t fGRngEta[]  = {-0.8, 0.8};
const Int_t   fGNBinsEta  = Int_t((fGRngEta[1] - fGRngEta[0])/fGBwRap) +1;

const Float_t fGRngRap[]  = {-0.8, 0.8};
const Int_t   fGNBinsRap  = Int_t((fGRngRap[1] - fGRngRap[0])/fGBwRap) +1;
const Float_t fGRngPhi[]  = {0.0, static_cast<Float_t>(TMath::TwoPi())};
const Int_t   fGNBinsPhi  = 76;

const Float_t fGRngPt[]    = {0.3, 3.0};
const Int_t   fGNBinsPt    = Int_t((fGRngPt[1] - fGRngPt[0])/fGBwPt); 

const Int_t   fGNBinsSign  =  2;
const Float_t fGRngSign[]  = {-0.5, 1.5};


const Char_t* fGEvtNames[] = {"All", "IsTriggered", "HasVertex", "Vx<Vx_{Max}", "Vy<Vy_{Max}", "Vz<Vz_{Max}", "Centrality [0,100]%", "Centrality [<0,>100]%"};
const Char_t* fGCMxNames[] = {"5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100"};
const Char_t* fGTrgNames[] = {"kMB", "kCentral", "kSemiCentral", "kEMCEJE", "kEMCEGA" }; 
const Char_t* fGCntNames[] = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%","50-60%", "60-70%", "70-80%", "80-90%", "90-100%"};

const Char_t* fgkPidName[4] = {"Nch","Npi","Nka","Npr"};
const Char_t* fgkPidShLatex[4] = {"N","#pi","K","p"};
const Char_t* fgkPidLatex[4][2]= {{"N_{-}","N_{+}"}, {"N_{#pi^{-}}","N_{#pi^{+}}"},{"N_{K^{-}}","N_{K^{+}}"}, {"N_{#bar{p}}","N_{p}"}};
const Char_t* fgkPidTitles[4][2]= {{"Negative","Positive"},{"Anti-Pions","Pions"},{"Anti-Kaons","Kaons"}, {"Anti-Protons","Protons"}};


const Char_t* fgkNetHistName[4] = {"","Plus","Minus","Net"};
const Char_t* fgkNetHistLatex[4]      = {"+ + +","+","-","+ - -"};


//---------------------------------------------------------------------------------
void AliEbyENetChargeFluctuationTask::UserCreateOutputObjects() {
  //Bool_t oldStatus = TH1::AddDirectoryStatus();
  //TH1::AddDirectory(kFALSE);
  
  fQaList = new TList();
  fQaList->SetOwner(kTRUE);

  fPhyList = new TList();
  fPhyList->SetOwner(kTRUE);

  fEffList = new TList();
  fEffList->SetOwner(kTRUE);

  fDcaList = new TList();
  fDcaList->SetOwner(kTRUE);
  
    
  fRan = new TRandom3();
  fRan->SetSeed();
  
  fRanIdx = new TRandom3();
  fRanIdx->SetSeed();
  
  Printf(" >>>%d %d %d %d %d %d %d %d %d %d", 
	 fIsAOD, fIsMC, fIsPhy, fIsEff, 
	 fIsDca, fIsQa, 
	 fIsRatio, fIsSub, fIsBS, fIsPer);
  
  CreateBasicQA();
  if (fIsQa)  CreateQA();
  if (fIsPhy) InitPhy();
  if (fIsMC && fIsEff) CreateCE();  
  if (fIsMC && fIsDca) CreateDEM();  
  if (fIsDca) CreateDED();  

  //  if (fIsQa) if (fESDtrackCuts) fQaList->Add(fESDtrackCuts);

  // TH1::AddDirectory(oldStatus);

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
    if (fIsPhy) ExecAA(); 
    if (fIsMC && fIsEff) CalEC();
    if (fIsDca) CalED();
  } 
  else 
    return;
  

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

  fESD = fESDHandler->GetEvent();
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
  Double_t rec[5] = {fCentralityBin,charge,rapp,track->Phi(),track->Pt()};
  
  if (gPid == 0) (static_cast<THnSparseD*>(fQaList->FindObject("fHnNchTrackUnCorr")))->Fill(rec);
  else if (gPid == 1) (static_cast<THnSparseD*>(fQaList->FindObject("fHnNpiTrackUnCorr")))->Fill(rec);
  else if (gPid == 2) (static_cast<THnSparseD*>(fQaList->FindObject("fHnNkaTrackUnCorr")))->Fill(rec);
  else if (gPid == 3) (static_cast<THnSparseD*>(fQaList->FindObject("fHnNprTrackUnCorr")))->Fill(rec);
}

void AliEbyENetChargeFluctuationTask::FillQAThnMc(AliVParticle *particle, Int_t gPid, Double_t rap) {
  Double_t charge = (particle->PdgCode() < 0) ? 0. : 1.;
  Double_t  rapp = (gPid == 0) ? particle->Eta() : rap;
  Double_t rec[5] = {fCentralityBin,charge,rapp,particle->Phi(),particle->Pt()};
  if (gPid == 0) (static_cast<THnSparseD*>(fQaList->FindObject("fHnNchTrackMc")))->Fill(fCurRec);
  else if (gPid == 1) (static_cast<THnSparseD*>(fQaList->FindObject("fHnNpiTrackMc")))->Fill(rec);
  else if (gPid == 2) (static_cast<THnSparseD*>(fQaList->FindObject("fHnNkaTrackMc")))->Fill(rec);
  else if (gPid == 3) (static_cast<THnSparseD*>(fQaList->FindObject("fHnNprTrackMc")))->Fill(rec);
}

//----------------------------------------------------------------------------------
void AliEbyENetChargeFluctuationTask::ExecAA(){
  
  // Printf("Number of Track %d",fNTracks);

  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliVTrack *track = (fESD) ? 
      static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : 
      static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 
    if(!AcceptTrack(track)) continue;

    Int_t icharge = track->Charge() < 0 ? 0 : 1;
    fNp[0][icharge] += 1.; 
    if (fIsQa) FillQAThnRec(track,0,0);
    
    Int_t a = fHelperPID->GetParticleSpecies(track,kTRUE);
    if(a < 0 || a > 2) continue;
    Int_t b = a + 1;

    Double_t rap;
    if (!TrackRapidity(track,rap,b)) continue;
    fNp[b][icharge] += 1.; 
    if (fIsQa) FillQAThnRec(track, b,rap);

  }
  
  //---- - -- - - - - -   -  -- - - - ---- - - - ---
 FillBasicHistos("Phy",kFALSE);
 if (fIsRatio) {
    FillRatioHistos("RatioBin",kFALSE,0);
    if(fIsPer)  FillRatioHistos("RatioPer",kFALSE,1);
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

  //---- - -- - - - - -   -  -- - - - ---- - - - ---
  if (fIsMC) {
    if (fIsAOD) {

      //   Printf("Number of Track %d",fArrayMC->GetEntries());

      for (Int_t idxMC = 0; idxMC < fArrayMC->GetEntries(); idxMC++) {
	AliAODMCParticle *particle = static_cast<AliAODMCParticle*>(fArrayMC->At(idxMC));
	if (!particle) 
	  continue;
	if (!AcceptTrackMC((AliVParticle*)particle, idxMC)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;
	fMCNp[0][icharge]    += 1.;    	   
	if(fIsQa)FillQAThnMc((AliVParticle*)particle,0,0);
	Int_t iPid = 0;  
	if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
	else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
	else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
	else    iPid = 0;

	if(iPid < 1 || iPid > 3) continue;

	Double_t rap;
	if (!ParticleRapidity((AliVParticle*)particle, rap, iPid)) continue;
	fMCNp[iPid][icharge] += 1.;    
	if(fIsQa)FillQAThnMc((AliVParticle*)particle,iPid,rap);   
      }

    } else if (fESD) {
      //---- - -- - - - - -   -  -- - - - ---- - - - --- 
         for (Int_t idxMC = 0; idxMC < fStack->GetNprimary(); ++idxMC) {
	   AliVParticle* particle = fMCEvent->GetTrack(idxMC);
	   if (!particle) 
	     continue;
	   if (!AcceptTrackMC(particle, idxMC)) continue;
	   Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;
	   fMCNp[0][icharge]  += 1.;    	   
	   if(fIsQa)FillQAThnMc(particle, 0, 0);

	   Int_t iPid = 0;  
	   if      (TMath::Abs(particle->PdgCode()) ==  211) iPid = 1; // pion
	   else if (TMath::Abs(particle->PdgCode()) ==  321) iPid = 2; // kaon
	   else if (TMath::Abs(particle->PdgCode()) == 2212) iPid = 3; // proton
	   else  iPid = 0;
	   if(iPid < 1 || iPid > 3) continue;
	   Double_t rap;
	   if (!ParticleRapidity(particle, rap, iPid)) continue;
	   fMCNp[iPid][icharge] += 1.;    
	   if(fIsQa)FillQAThnMc(particle, iPid, rap);    
	 }
    }
    
    //---- - -- - - - - -   -  -- - - - ---- - - - --- 
    FillBasicHistos("MCPhy",kFALSE);
    if (fIsRatio) {
      FillRatioHistos("MCRatioBin",kFALSE,0);
      if(fIsPer)  FillRatioHistos("MCRatioPer",kFALSE,1);
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
    //---- - -- - - - - -   -  -- - - - ---- - - - --- 
  }
  //---- - -- - - - - -   -  -- - - - ---- - - - --- 
  /*
  Printf("%6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f",
	 fNp[0][0],fNp[0][1],
	 fNp[1][0],fNp[1][1],
	 fNp[2][0],fNp[2][1],
	 fNp[3][0],fNp[3][1]);

Printf("%6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f",
	 fMCNp[0][0],fMCNp[0][1],
	 fMCNp[1][0],fMCNp[1][1],
	 fMCNp[2][0],fMCNp[2][1],
	 fMCNp[3][0],fMCNp[3][1]);
  */

}

//___________________________________________________________
void AliEbyENetChargeFluctuationTask::Terminate( Option_t * ){
  Info("AliEbyENetChargeFluctuationTask"," Task Successfully finished");
}

//___________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::AcceptTrack(AliVTrack *track) const {
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
 if (!trackAOD->TestFilterBit(fAODtrackCutBit))
   return kFALSE;
 } else {      // ESDs
   if(!fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))  return kFALSE;
 }


 if(track->Pt() < fPtMin || track->Pt() > fPtMax )  return kFALSE; 
 if (TMath::Abs(track->Eta()) > fEtaMax) return kFALSE; 
 
  
 if (track->Phi() > fPhiMax) return kFALSE;   
 return kTRUE;
}


//___________________________________________________________
Bool_t AliEbyENetChargeFluctuationTask::AcceptTrackMC(AliVParticle *particle, Int_t idxMC) const {
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
  
 Double_t mass = AliPID::ParticleMass(AliPID::kPion);
 if(gCurPid == 1) mass = AliPID::ParticleMass(AliPID::kPion);
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


//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::CreateQA() {
 
  Int_t    bhuc[5] = {fGNBinsCent, fGNBinsSign, fGNBinsRap,fGNBinsPhi,fGNBinsPt };      
  Double_t mnhuc[5] = {fGRngCent[0],fGRngSign[0],fGRngRap[0],fGRngPhi[0],fGRngPt[0]};  
  Double_t mxhuc[5] = {fGRngCent[1],fGRngSign[1],fGRngRap[1],fGRngPhi[1],fGRngPt[1]};  
  const Char_t *ctname = "cent:sign:y:phi:pt";
			     
  fQaList->Add(new THnSparseD("fHnNpiTrackUnCorr",ctname,  5, bhuc, mnhuc, mxhuc));
  fQaList->Add(new THnSparseD("fHnNkaTrackUnCorr", ctname, 5, bhuc, mnhuc, mxhuc));
  fQaList->Add(new THnSparseD("fHnNprTrackUnCorr", ctname, 5, bhuc, mnhuc, mxhuc));

  if (fIsMC) {
    fQaList->Add(new THnSparseD("fHnNpiTrackMc", ctname, 5, bhuc, mnhuc, mxhuc));
    fQaList->Add(new THnSparseD("fHnNkaTrackMc", ctname, 5, bhuc, mnhuc, mxhuc));
    fQaList->Add(new THnSparseD("fHnNprTrackMc", ctname, 5, bhuc, mnhuc, mxhuc));

  }

  ctname = "cent:sign:eta:phi:pt";
  bhuc[2]  = fGNBinsEta;
  mnhuc[2] = fGRngEta[0];
  mxhuc[2] = fGRngEta[1]; 
  
  fQaList->Add(new THnSparseD("fHnNchTrackUnCorr", ctname, 5, bhuc, mnhuc, mxhuc));
  if (fIsMC) {
    fQaList->Add(new THnSparseD("fHnNchTrackMc", ctname, 5, bhuc, mnhuc, mxhuc));
  }
}

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::CreateBasicQA() {
  // -- Basic QA

 
  fQaList->Add(new TH1F("hEventStat0","Event cut statistics 0;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5));
  fQaList->Add(new TH1F("hEventStat1","Event cut statistics 1;Event Cuts;Events", fHEventStatMax,-0.5,fHEventStatMax-0.5));
  
  for ( Int_t ii=0; ii < fHEventStatMax; ii++ ) {
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat0"))))->GetXaxis()->SetBinLabel(ii+1, fGEvtNames[ii]);
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat1"))))->GetXaxis()->SetBinLabel(ii+1, fGEvtNames[ii]);
  }

  (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat0"))))->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", fGCMxNames[fCentralityBinMax-1]));
  (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat1"))))->GetXaxis()->SetBinLabel(fHEventStatMax, Form("Centrality [0-%s]%%", fGCMxNames[fCentralityBinMax-1]));
  
  fQaList->Add(new TH1F("hTriggerStat","Trigger statistics;Trigger;Events", fNTriggers,-0.5,fNTriggers-0.5));

  for ( Int_t ii=0; ii < fNTriggers; ii++ ) {
    (static_cast<TH1F*>(fQaList->FindObject(Form("hTriggerStat"))))->GetXaxis()->SetBinLabel(ii+1, fGTrgNames[ii]);
  }

  
  // -- Initialize trigger statistics histograms
  
  fQaList->Add(new TH1F("hCentralityStat","Centrality statistics;Centrality Bins;Events", 
					   fNCentralityBins,-0.5,fNCentralityBins-0.5));
  
  for ( Int_t ii=0; ii < fNCentralityBins; ii++ ) {
    (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityStat"))))->GetXaxis()->SetBinLabel(ii+1, fGCntNames[ii]);
  }
  
  fQaList->Add(new TH1F("hCentralityPercentileAccepted","Centrality Percentile statistics;Centrality Bins;Events", 
			100,-0.5,99.5));
  
  fQaList->Add(new TH1F("hCentralityPercentileAll","Centrality Percentile statistics;Centrality Bins;Events", 
			100,-0.5,99.5));


  fQaList->Add(new TH2F("fHistQAvx",  "Histo Vx Selected;Centrality;Vx", 100,0,100, 5000, -5., 5.));
  fQaList->Add(new TH2F("fHistQAvy",  "Histo Vy Selected;Centrality;Vy", 100,0,100, 5000, -5., 5.));
  fQaList->Add(new TH2F("fHistQAvz",  "Histo Vz Selected;Centrality;Vz", 100,0,100, 5000, -25., 25.)); 
  fQaList->Add(new TH2F("fHistQAvxA", "Histo Vx;Centrality;Vx", 100,0,100, 5000, -5., 5.));
  fQaList->Add(new TH2F("fHistQAvyA", "Histo Vy;Centrality;Vy", 100,0,100, 5000, -5., 5.));
  fQaList->Add(new TH2F("fHistQAvzA", "Histo Vz;Centrality;Vz", 100,0,100, 5000, -25., 25.));
    
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
void AliEbyENetChargeFluctuationTask::SetAnal(Int_t i){
  switch(i){
  case 0:
    fIsQa  = 1; 
    break;
  case 1:
    fIsDca = 1; 
    break;
  case 2:
    fIsEff = 1; 
    break;
  case 3:
    fIsPhy   = 1; 
    fIsRatio = 1; 
    break;
  case 4:
    fIsPhy   = 1; 
    fIsRatio = 1; fIsPer = 1;
    break;
  case 5:
    fIsPhy = 1; 
    fIsSub = 1; 
    break;
  case 6:
    fIsPhy = 1; 
    fIsSub = 1; fIsPer = 1; 
    break;
  case 7:
    fIsPhy = 1; 
    fIsBS = 1; fIsPer = 1;
    break;
  case 8:
    fIsPhy = 1; 
    fIsSub = 1; fIsBS = 1; 
    break;
  case 9:
    fIsPhy = 1; 
    fIsSub = 1; fIsBS = 1; fIsPer = 1;
    break;
  case 10:
    fIsPhy   = 1;
    fIsRatio = 1; fIsSub = 1; fIsBS  = 1; 
    break;
  case 11:
    fIsPhy   = 1;
    fIsRatio = 1; fIsSub = 1; fIsBS  = 1; fIsPer = 1;
    break;
  case 12:
   fIsPhy   = 1; fIsEff = 1; fIsDca = 1; fIsQa  = 1; 
   fIsRatio = 1; fIsSub = 1; fIsBS  = 1; fIsPer = 1;
   break;
  default:
    cerr<<"Error:  cannot fill histograms!"<<endl;
    fIsPhy   = 0; fIsEff = 0; fIsDca = 0; fIsQa  = 0; 
    fIsRatio = 0; fIsSub = 0; fIsBS  = 0; fIsPer = 0;
    break;
  }

  Printf(" >>> %d %d %d %d %d %d %d %d %d", 
	 i, fIsPhy, fIsEff, fIsDca, fIsQa, 
	 fIsRatio, fIsSub, fIsBS, fIsPer);
  if (fIsEff) {fPtMin = 0.1; fPtMax = 3.;}  

}



//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::InitPhy() {

  Double_t ptRange[2] = {fPtMin,fPtMax};


  TString sTitle("");
  if (fIsRatio) {
    CreateRatioHistos("RatioBin",     Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),0);
    if (fIsPer) CreateRatioHistos("RatioPer",     Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),1);
  }
  CreateBasicHistos("Phy",                      Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]));
  if (fIsSub) {
    CreateGroupHistos("PhyBinSS",        Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fNSubSamples,0);  
    if (fIsPer) CreateGroupHistos("PhyPerSS",        Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fNSubSamples,1);  
  }
  if (fIsBS){
    CreateGroupHistos("PhyBinBS",        Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fNSubSamples,0);  
    if (fIsPer) CreateGroupHistos("PhyPerBS",        Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fNSubSamples,1);  
  }
  if (fIsMC) {
    TString sMCTitle("");
    CreateBasicHistos("MCPhy",                     Form("%s", sTitle.Data()));
    
    if (fIsRatio){
      CreateRatioHistos("MCRatioBin", Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),0);
      if (fIsPer)  CreateRatioHistos("MCRatioPer", Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),1);
    }
    if (fIsSub){
      CreateGroupHistos("MCBinSS",       Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fNSubSamples,0);  
      if (fIsPer)  CreateGroupHistos("MCPerSS",       Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fNSubSamples,1);  
    }
    if (fIsBS){
      CreateGroupHistos("MCBinBS",       Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fNSubSamples,0);  
      if (fIsPer)  CreateGroupHistos("MCPerBS",       Form("%s, #it{p}_{T} [%.1f,%.1f]", sTitle.Data(), ptRange[0], ptRange[1]),fNSubSamples,1);  
    }
  }
}

//________________________________________________________________________
void  AliEbyENetChargeFluctuationTask::CreateBasicHistos(const Char_t *name, const Char_t *title)  {
  TString sName(name);
  TString sTitle(title);
  Double_t etaRange[2] = {fEtaMin,fEtaMax};
  fPhyList->Add(new TList);
  TList *list =  static_cast<TList*>(fPhyList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);
  

  Int_t nBinsCent         =  fGNBinsCent;
  Double_t centBinRange[] = {fGRngCent[0], fGRngCent[1]};

 
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    TString sNetTitle(Form("%s - %s", fgkPidLatex[iPid][1], fgkPidLatex[iPid][0]));

    sTitle = (iPid != 0 ) ? Form("|y| < %.1f", fRapMax) : Form(" |#eta|<%.1f", etaRange[1]);

    list->Add(new TProfile(Form("fProfTot%sPlus%s", fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(100);(%s)",fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   100,-0.5,99.5));

    list->Add(new TProfile(Form("fProfTot%sMinus%s", fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(100);(%s)",fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   100,-0.5,99.5));


    

    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile(Form("fProf%s%sNet%dM", fgkPidName[iPid],name, idx), 
			     Form("(%s)^{%d} : %s;Centrality(100);(%s)^{%d}",sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			     100,-0.5,99.5));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	list->Add(new TProfile(Form("fProf%s%sNetF%02d%02d", fgkPidName[iPid], name, ii, kk),
			       Form("f_{%02d%02d} : %s;Centrality(100);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
			       100,-0.5,99.5));
      }
    }
  
  }  
   
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    TString sNetTitle(Form("%s - %s", fgkPidLatex[iPid][1], fgkPidLatex[iPid][0]));
    sTitle = (iPid != 0 ) ? Form(" |y|<%.1f", fRapMax) : Form(" |#eta| < %.1f", etaRange[1]);

    list->Add(new TProfile(Form("fProfBinTot%sPlus%s", fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(11);(%s)",fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   nBinsCent, centBinRange[0], centBinRange[1]));

    list->Add(new TProfile(Form("fProfBinTot%sMinus%s", fgkPidName[iPid],name), 
			   Form("(%s) : %s;Centrality(11);(%s)",fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
			   nBinsCent, centBinRange[0], centBinRange[1]));



    for (Int_t idx = 1; idx <= fOrder; ++idx) {
      list->Add(new TProfile(Form("fProfBin%s%sNet%dM", fgkPidName[iPid],name, idx), 
			     Form("(%s)^{%d} : %s;Centrality(11);(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
			     nBinsCent, centBinRange[0], centBinRange[1]));
    }
    
    for (Int_t ii = 0; ii <= fOrder; ++ii) {
      for (Int_t kk = 0; kk <= fOrder; ++kk) {
	list->Add(new TProfile(Form("fProfBin%s%sNetF%02d%02d", fgkPidName[iPid], name, ii, kk),
			       Form("f_{%02d%02d} : %s;Centrality(11);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
			       nBinsCent, centBinRange[0], centBinRange[1]));
      }
    }
  
  }  
  
  for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
    list->Add(new TProfile(Form("fProf%sNu%02d",name,iPhy),Form("Physics Variable for index %d | %s ; Centrality;",iPhy,name),100,-0.5,99.5));
  }
  for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
    list->Add(new TProfile(Form("fProfBin%sNu%02d",name,iPhy),Form("Physics Variable for index %d | %s ; Centrality;",iPhy,name),nBinsCent, centBinRange[0], centBinRange[1]));
  }
  
  return;
}



//________________________________________________________________________
void  AliEbyENetChargeFluctuationTask::CreateRatioHistos(const Char_t *name, const Char_t *title, Bool_t isPer)  {
  TString sName(name);
  TString sTitle(title);
  
  fPhyList->Add(new TList);
  TList *list =  static_cast<TList*>(fPhyList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);
  
  Int_t    nRbin  = 15000;
  Double_t mRat[] = {0,1.5};

  Int_t nBinsCent         =  (isPer) ? 100 : fGNBinsCent;
  Double_t centBinRange[2];  
  centBinRange[0]  =  (isPer) ?  0   : fGRngCent[0];
  centBinRange[1]  =  (isPer) ?  100 : fGRngCent[1];

  TString xyz = Form("|y| < %.1f",fRapMax); 

  list->Add(new TH2F(Form("fHistRatioKPi%s",name), 
		     Form("(%s %s) : K/#pi;Centrality(11);K/#pi", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioKpPip%s",name), 
		     Form("(%s %s) : K^{+}/#pi^{+};Centrality(11);K^{+}/#pi^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioKmPip%s",name), 
		     Form("(%s %s) : K^{-}/#pi^{+};Centrality(11);K^{-}/#pi^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioKmPim%s",name), 
		     Form("(%s %s) : K^{-}/#pi^{-};Centrality(11);K^{-}/#pi^{-}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));



 list->Add(new TH2F(Form("fHistRatioPK%s",name), 
		     Form("(%s %s) : P/K;Centrality(11);P/K", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPpKp%s",name), 
		     Form("(%s %s) : P/K^{+};Centrality(11);P/K^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPmKp%s",name), 
		     Form("(%s %s) : #bar{P}/K^{+};Centrality(11);#bar{P}/K^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPmKm%s",name), 
		     Form("(%s %s) : #bar{P}/K^{-};Centrality(11);#bar{P}/K^{-}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));



 list->Add(new TH2F(Form("fHistRatioPPi%s",name), 
		     Form("(%s %s) : P/#pi;Centrality(11);K/#pi", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPpPip%s",name), 
		     Form("(%s %s) : P/#pi^{+};Centrality(11);P/#pi^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPmPip%s",name), 
		     Form("(%s %s) : #bar{P}/#pi^{+};Centrality(11);#bar{P}/#pi^{+}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  
  list->Add(new TH2F(Form("fHistRatioPmPim%s",name), 
		     Form("(%s %s) : #bar{P}/#pi^{-};Centrality(11);#bar{P}/#pi^{-}", xyz.Data(), sTitle.Data()),
		     nBinsCent, centBinRange[0], centBinRange[1], nRbin,mRat[0],mRat[1]));
  

  //------- - - -  -  -   -   - - -   - --- - - --- - - - - - -- --------
  Int_t bin[4] = {2800,2200,1200,600}; 
  Int_t bd[] = {1,2,2,2};
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    Int_t bb = bin[iPid];

    for (Int_t iNet = 0; iNet < 4; ++iNet) {
      Int_t bn     = (iNet == 3) ?  500   : bb/bd[iNet]; 
      Float_t blow = (iNet == 3) ? -250.5 : -0.5;
      Float_t bup  = (iNet == 3) ?  249.5 : bn-0.5;


      list->Add(new TH2F(Form("fHistDist%s%s%s",name, fgkPidName[iPid], fgkNetHistName[iNet]), 
			 Form("(%s %s) : %s Distribution;Centrality(11);%s_{(%s)}", xyz.Data(), sTitle.Data(), 
			      fgkPidShLatex[iPid],fgkPidShLatex[iPid],fgkNetHistLatex[iNet]),
			 nBinsCent, centBinRange[0], centBinRange[1], bn, blow,bup));    
    }
  }
  



return;
}

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::FillBasicHistos(const Char_t *name, Bool_t isMC)  {
  //Double_t **np = (isMC) ? fMCNp : fNp;
  
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

  Float_t centralityBin = fCentralityBin;
  Float_t centralityPer = fCentralityPercentile;
  
  TList *list = static_cast<TList*>(fPhyList->FindObject(Form("f%s",name)));
  
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    Int_t deltaNp = np[iPid][1]-np[iPid][0];  
    Double_t delta = 1.;

    (static_cast<TProfile*>(list->FindObject(Form("fProfBinTot%sPlus%s", fgkPidName[iPid], name))))->Fill(centralityBin, np[iPid][1]);
    (static_cast<TProfile*>(list->FindObject(Form("fProfTot%sPlus%s", fgkPidName[iPid], name))))->Fill(centralityPer, np[iPid][1]);

    (static_cast<TProfile*>(list->FindObject(Form("fProfBinTot%sMinus%s", fgkPidName[iPid], name))))->Fill(centralityBin, np[iPid][0]);
    (static_cast<TProfile*>(list->FindObject(Form("fProfTot%sMinus%s", fgkPidName[iPid], name))))->Fill(centralityPer, np[iPid][0]);


    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
      delta *= deltaNp;

      (static_cast<TProfile*>(list->FindObject(Form("fProfBin%s%sNet%dM", fgkPidName[iPid], name, idxOrder))))->Fill(centralityBin, delta);
      (static_cast<TProfile*>(list->FindObject(Form("fProf%s%sNet%dM", fgkPidName[iPid], name, idxOrder))))->Fill(centralityPer, delta);
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
	(static_cast<TProfile*>(list->FindObject(Form("fProfBin%s%sNetF%02d%02d", fgkPidName[iPid], name, ii, kk))))->Fill(centralityBin, fik);
	(static_cast<TProfile*>(list->FindObject(Form("fProf%s%sNetF%02d%02d", fgkPidName[iPid], name, ii, kk))))->Fill(centralityPer, fik);
      }
    }
  }
 
  //Printf("%6d %20s %6.2f %6d %6d %6d %6d  %6d %6d %6d %6d", idx, name, centralityBin,
  //	 np[0][1],  np[0][0], 
  //	 np[1][1],  np[1][0], 
  ///	 np[2][1],  np[2][0], 
  //	 np[3][1],  np[3][0]);
  //

   Int_t a[6][4]; Int_t b[22];
   for (Int_t iPid = 0; iPid < 4; ++iPid) {
     a[0][iPid] = np[iPid][1]+np[iPid][0];       // 0  n+ + n-
     a[1][iPid] = np[iPid][1];                        // 1  n+
     a[2][iPid] = np[iPid][0];                        // 2  n-
     a[3][iPid] = np[iPid][1]*np[iPid][0];       // 3  n+ . n-
     a[4][iPid] = np[iPid][1]*(np[iPid][1]-1);   // 4  n+ (n+ - 1)
     a[5][iPid] = np[iPid][0]*(np[iPid][0]-1);   // 5  n- (n- - 1)
     
     // Printf("%6d %20s %6.2f %6d %6d %6d ", idx, name, centralityBin,
     //	   a[0][iPid], a[1][iPid], a[2][iPid]);

  }
  
  b[0]  = a[0][0]*a[0][2];       // 24 N   K
  b[1]  = a[0][1]*a[0][2];       // 25 Pi  K
  b[2]  = a[1][1]*a[1][2];       // 26 pi+ k+
  b[3]  = a[1][1]*a[2][2];       // 27 pi+ k-
  b[4]  = a[2][1]*a[1][2];       // 28 pi- k+  
  b[5]  = a[2][1]*a[2][2];       // 29 pi- k-
  
  b[6]  = a[0][0]*a[0][3];       // 30 N   P
  b[7]  = a[0][2]*a[0][3];       // 31 K   P
  b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
  b[9]  = a[1][2]*a[2][3];       // 33 k+  p-
  b[10] = a[2][2]*a[1][3];       // 34 k-  p+
  b[11] = a[2][2]*a[2][3];       // 35 k-  p-
  
  b[12] = a[0][0]*a[0][1];       // 36 N  Pi
  b[13] = a[0][3]*a[0][1];       // 37 P  Pi
  b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
  b[15] = a[1][3]*a[2][1];       // 39 p+ pi-
  b[16] = a[2][3]*a[1][1];       // 40 p- pi+
  b[17] = a[2][3]*a[2][1];       // 41 p- pi-
  
  b[18] = a[0][0]*(a[0][0] - 1); // 42 N ( N - 1 )
  b[19] = a[0][1]*(a[0][1] - 1); // 43 Pi( Pi- 1 )
  b[20] = a[0][2]*(a[0][1] - 1); // 44 K ( K - 1 )
  b[21] = a[0][3]*(a[0][3] - 1); // 45 P ( P - 1 )
  // TList *list_nu = static_cast<TList*>(fPhyList->FindObject(Form("f%s_nu",name)));
  Int_t k = 0;
  for (Int_t j = 0; j < 4; j++) {
    for (Int_t i = 0; i < 6; i++) {
      (static_cast<TProfile*>(list->FindObject(Form("fProfBin%sNu%02d", name,k))))->Fill(centralityBin,a[i][j]); 
      (static_cast<TProfile*>(list->FindObject(Form("fProf%sNu%02d", name,k))))->Fill(centralityPer,a[i][j]); 
      k++;
    }
  }

  for (Int_t j = 0; j < 22; j++) {
    (static_cast<TProfile*>(list->FindObject(Form("fProfBin%sNu%02d", name,j+24))))->Fill(centralityBin,b[j]); 
    (static_cast<TProfile*>(list->FindObject(Form("fProf%sNu%02d", name,j+24))))->Fill(centralityPer,b[j]); 
  }
  
  return;
}


//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::FillRatioHistos(const Char_t *name, Bool_t isMC,Bool_t isPer)  {
   
  //  Double_t **np = (isMC) ? fMCNp : fNp;

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



  Float_t centralityBin = (isPer) ? fCentralityPercentile : fCentralityBin;

  TList *list = static_cast<TList*>(fPhyList->FindObject(Form("f%s",name)));
    
  
  if((Double_t)np[1][1]+(Double_t)np[1][0] != 0 ) {
    Double_t KPi = ((Double_t)np[2][1]+(Double_t)np[2][0])/((Double_t)np[1][1]+(Double_t)np[1][0]);
    Double_t PPi = ((Double_t)np[3][1]+(Double_t)np[3][0])/((Double_t)np[1][1]+(Double_t)np[1][0]);
    
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioKPi%s",name))))->Fill(centralityBin, KPi);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPPi%s",name))))->Fill(centralityBin,   PPi);
  }

  if((Double_t)np[2][1]+(Double_t)np[2][0] != 0 ){
    Double_t PK = ((Double_t)np[3][1]+(Double_t)np[3][0])/((Double_t)np[2][1]+(Double_t)np[2][0]);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPK%s",name))))->Fill(centralityBin, PK);
  }

  if ((Double_t)np[1][1] != 0 ) {
    Double_t KpPip  = ((Double_t)np[2][1])/((Double_t)np[1][1]); 
    Double_t KmPip  = ((Double_t)np[2][0])/((Double_t)np[1][1]); 
    Double_t PpPip  = ((Double_t)np[3][1])/((Double_t)np[1][1]); 
    Double_t PmPip =  ((Double_t)np[3][0])/((Double_t)np[1][1]);

    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioKpPip%s",name))))->Fill(centralityBin, KpPip);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioKmPip%s",name))))->Fill(centralityBin, KmPip);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPpPip%s",name))))->Fill(centralityBin, PpPip);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPmPip%s",name))))->Fill(centralityBin, PmPip);
  }  

  if ((Double_t)np[1][0] != 0) {
    Double_t KmPim = ((Double_t)np[2][0])/((Double_t)np[1][0]);
    Double_t PmPim = ((Double_t)np[3][0])/((Double_t)np[1][0]);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioKmPim%s",name))))->Fill(centralityBin, KmPim);
    (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPmPim%s",name))))->Fill(centralityBin, PmPim);
  }
  
  if ((Double_t)np[2][1] != 0 ) { 
    Double_t PpKp  = ((Double_t)np[3][1])/((Double_t)np[2][1]); 
    Double_t PmKp =  ((Double_t)np[3][0])/((Double_t)np[2][1]); 
   (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPpKp%s",name))))->Fill(centralityBin, PpKp);
   (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPmKp%s",name))))->Fill(centralityBin, PmKp);
  }
  
  if ((Double_t)np[2][0] != 0) {
   Double_t PmKm = ((Double_t)np[3][0])/((Double_t)np[2][0]);
   (static_cast<TH2F*>(list->FindObject(Form("fHistRatioPmKm%s",name))))->Fill(centralityBin, PmKm);
  }
  
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    (static_cast<TH2F*>(list->FindObject(Form("fHistDist%s%s%s",name, fgkPidName[iPid], fgkNetHistName[0]))))->Fill(centralityBin, np[iPid][1]+np[iPid][0]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHistDist%s%s%s",name, fgkPidName[iPid], fgkNetHistName[1]))))->Fill(centralityBin, np[iPid][1]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHistDist%s%s%s",name, fgkPidName[iPid], fgkNetHistName[2]))))->Fill(centralityBin,                  np[iPid][0]); 
    (static_cast<TH2F*>(list->FindObject(Form("fHistDist%s%s%s",name, fgkPidName[iPid], fgkNetHistName[3]))))->Fill(centralityBin, np[iPid][1]-np[iPid][0]); 
  }
  

  return;
}

//________________________________________________________________________
void  AliEbyENetChargeFluctuationTask::CreateGroupHistos(const Char_t *name, const Char_t *title, Int_t nSample, Bool_t isPer)  {

  TString sName(name);
  TString sTitle(title);
  
 

  Float_t etaRange[2] = {fEtaMin, fEtaMax};
  

  //TList *list[4];
  fPhyList->Add(new TList);
  TList *list =  static_cast<TList*>(fPhyList->Last());
  list->SetName(Form("f%s", name));
  list->SetOwner(kTRUE);
  

  TString tname = Form("%s", name);
  Int_t nBinsCent         =  (isPer) ? 100 : fGNBinsCent;
  Double_t centBinRange[2];  
  centBinRange[0]  =  (isPer) ?  0   : fGRngCent[0];
  centBinRange[1]  =  (isPer) ?  100 : fGRngCent[1];
 
  for (Int_t iSub = 0; iSub <= nSample; ++iSub) {
    
    list->Add(new TList);
    TList *listSub = static_cast<TList*>(list->Last());
    listSub->SetName(Form("%s%02d",name, iSub));
    listSub->SetOwner(kTRUE);
    
    for (Int_t iPid = 0; iPid < 4; ++iPid) {
    
      TString sNetTitle(Form("%s - %s", fgkPidLatex[iPid][1], fgkPidLatex[iPid][0]));
      sTitle = (iPid != 0 ) ? Form("|y| < %.1f", fRapMax) : Form(" |#eta|<%.1f", etaRange[1]);
      
      listSub->Add(new TProfile(Form("fProfS%02dTot%sPlus%s", iSub, fgkPidName[iPid],tname.Data()), 
				Form("(%s) : %s;Centrality(11);(%s)",fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
				nBinsCent, centBinRange[0], centBinRange[1]));
      
      listSub->Add(new TProfile(Form("fProfS%02dTot%sMinus%s",iSub, fgkPidName[iPid],tname.Data()), 
				Form("(%s) : %s;Centrality(11);(%s)",fgkPidName[iPid], sTitle.Data(), sNetTitle.Data()),
				nBinsCent, centBinRange[0], centBinRange[1]));
      
      
      
      for (Int_t idx = 1; idx <= fOrder; ++idx) {
	listSub->Add(new TProfile(Form("fProfS%02d%s%sNet%dM",iSub, fgkPidName[iPid],tname.Data(), idx), 
				  Form("(%s)^{%d} : %s;Centrality(11);(%s)^{%d}", sNetTitle.Data(), idx, sTitle.Data(), sNetTitle.Data(), idx),
				  nBinsCent, centBinRange[0], centBinRange[1]));
      }
      
      for (Int_t ii = 0; ii <= fOrder; ++ii) {
	for (Int_t kk = 0; kk <= fOrder; ++kk) {
	  listSub->Add(new TProfile(Form("fProfS%02d%s%sNetF%02d%02d",iSub, fgkPidName[iPid], tname.Data(), ii, kk),
				    Form("f_{%02d%02d} : %s;Centrality(11);f_{%02d%02d}", ii, kk, sTitle.Data(), ii, kk),
				    nBinsCent, centBinRange[0], centBinRange[1]));
	}
      }
                  
    }  // iPid
    //-------------------------------------------
       
    for (Int_t iPhy = 0; iPhy < 46; ++iPhy) { 
      listSub->Add(new TProfile(Form("fProfS%02d%sNu%02d",iSub,tname.Data(),iPhy),
				Form("Physics Variable for index %d | %s | Sub S%02d; Centrality;",iPhy,tname.Data(), iSub),
				nBinsCent, centBinRange[0], centBinRange[1]));
    }
    
  }// isub
  
  return;
}

//________________________________________________________________________
void AliEbyENetChargeFluctuationTask::FillGroupHistos(const Char_t *name,Int_t iSub, Bool_t isMC,Bool_t isPer)  {

  Float_t centralityBin = (isPer) ? fCentralityPercentile : fCentralityBin;
  TString tname = Form("%s", name);
  
  TList *list    = static_cast<TList*>(fPhyList->FindObject(Form("f%s",name)));
  TList *listSub = static_cast<TList*>(list->FindObject(Form("%s%02d",name, iSub)));

  // Double_t **np = (isMC) ? fMCNp : fNp;
     
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
 
  for (Int_t iPid = 0; iPid < 4; ++iPid) {
    Int_t deltaNp = np[iPid][1]-np[iPid][0];  
    Double_t delta = 1.;
    
    (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02dTot%sPlus%s",iSub, fgkPidName[iPid], tname.Data()))))->Fill(centralityBin, np[iPid][1]);
    (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02dTot%sMinus%s",iSub, fgkPidName[iPid], tname.Data()))))->Fill(centralityBin, np[iPid][0]);
    
    for (Int_t idxOrder = 1; idxOrder <= fOrder; ++idxOrder) {
      delta *= deltaNp;
      (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%s%sNet%dM",iSub, fgkPidName[iPid], tname.Data(), idxOrder))))->Fill(centralityBin, delta);
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
       	(static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%s%sNetF%02d%02d",iSub, fgkPidName[iPid], tname.Data(), ii, kk))))->Fill(centralityBin, fik);
      }
    }
  }
 
  //Printf("%6d %20s %6.2f %6d %6d %6d %6d  %6d %6d %6d %6d", idx, name, centralityBin,
  //	 np[0][1],  np[0][0], 
  //	 np[1][1],  np[1][0], 
  ///	 np[2][1],  np[2][0], 
  //	 np[3][1],  np[3][0]);
  //

   Int_t a[6][4]; Int_t b[22];
   for (Int_t iPid = 0; iPid < 4; ++iPid) {
     a[0][iPid] = np[iPid][1]+np[iPid][0];       // 0  n+ + n-
     a[1][iPid] = np[iPid][1];                        // 1  n+
     a[2][iPid] = np[iPid][0];                        // 2  n-
     a[3][iPid] = np[iPid][1]*np[iPid][0];       // 3  n+ . n-
     a[4][iPid] = np[iPid][1]*(np[iPid][1]-1);   // 4  n+ (n+ - 1)
     a[5][iPid] = np[iPid][0]*(np[iPid][0]-1);   // 5  n- (n- - 1)
     
     // Printf("%6d %20s %6.2f %6d %6d %6d ", idx, name, centralityBin,
     //	   a[0][iPid], a[1][iPid], a[2][iPid]);

  }
  
  b[0]  = a[0][0]*a[0][2];       // 24 N   K
  b[1]  = a[0][1]*a[0][2];       // 25 Pi  K
  b[2]  = a[1][1]*a[1][2];       // 26 pi+ k+
  b[3]  = a[1][1]*a[2][2];       // 27 pi+ k-
  b[4]  = a[2][1]*a[1][2];       // 28 pi- k+  
  b[5]  = a[2][1]*a[2][2];       // 29 pi- k-
  
  b[6]  = a[0][0]*a[0][3];       // 30 N   P
  b[7]  = a[0][2]*a[0][3];       // 31 K   P
  b[8]  = a[1][2]*a[1][3];       // 32 k+  p+
  b[9]  = a[1][2]*a[2][3];       // 33 k+  p-
  b[10] = a[2][2]*a[1][3];       // 34 k-  p+
  b[11] = a[2][2]*a[2][3];       // 35 k-  p-
  
  b[12] = a[0][0]*a[0][1];       // 36 N  Pi
  b[13] = a[0][3]*a[0][1];       // 37 P  Pi
  b[14] = a[1][3]*a[1][1];       // 38 p+ pi+
  b[15] = a[1][3]*a[2][1];       // 39 p+ pi-
  b[16] = a[2][3]*a[1][1];       // 40 p- pi+
  b[17] = a[2][3]*a[2][1];       // 41 p- pi-
  
  b[18] = a[0][0]*(a[0][0] - 1); // 42 N ( N - 1 )
  b[19] = a[0][1]*(a[0][1] - 1); // 43 Pi( Pi- 1 )
  b[20] = a[0][2]*(a[0][1] - 1); // 44 K ( K - 1 )
  b[21] = a[0][3]*(a[0][3] - 1); // 45 P ( P - 1 )
  // TList *list_nu = static_cast<TList*>(fOutlistSub->FindObject(Form("f%s_nu",name)));
  Int_t k = 0;
  for (Int_t j = 0; j < 4; j++) {
    for (Int_t i = 0; i < 6; i++) {
      (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%sNu%02d",iSub, tname.Data(),k))))->Fill(centralityBin,a[i][j]); 
       k++;
    }
  }

  for (Int_t j = 0; j < 22; j++) {
    (static_cast<TProfile*>(listSub->FindObject(Form("fProfS%02d%sNu%02d",iSub, tname.Data(),j+24))))->Fill(centralityBin,b[j]); 
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
  
  Int_t centBin = centrality->GetCentralityClass10(fCentralityEstimator.Data());
  if (centBin == 0) { fCentralityBin = centrality->GetCentralityClass5(fCentralityEstimator.Data()); }
  else if (centBin == 11 || centBin == -1.)           { fCentralityBin = -1; }
  else if (centBin > 0 && centBin < fNCentralityBins) { fCentralityBin = centBin + 1; }
  else {  fCentralityBin = -2; }

  if (fCentralityBin >= fCentralityBinMax) fCentralityBin = -2;

  fCentralityPercentile = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
  
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
      isTriggered = kTRUE;
      (static_cast<TH1F*>(fQaList->FindObject(Form("hTriggerStat"))))->Fill(ii);
    }
  }
  
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
    (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvx"))))->Fill(fCentralityPercentile,vtxESD->GetX());
    if(TMath::Abs(vtxESD->GetX()) > fVxMax)  aEventCuts[iStep] = 1;  //3
    else (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvxA"))))->Fill(fCentralityPercentile,vtxESD->GetX());
  }
  else if(vtxAOD){
    (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvx"))))->Fill(fCentralityPercentile,vtxAOD->GetX());
    if(TMath::Abs(vtxAOD->GetX()) > fVxMax) aEventCuts[iStep] = 1;    //3
    else (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvxA"))))->Fill(fCentralityPercentile,vtxAOD->GetX());
  }
  else aEventCuts[iStep] = 1; //3
  
  ++iStep;
  if (vtxESD){
    (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvy"))))->Fill(fCentralityPercentile,vtxESD->GetY());
    if(TMath::Abs(vtxESD->GetY()) > fVyMax) aEventCuts[iStep] = 1; //4
    else (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvyA"))))->Fill(fCentralityPercentile,vtxESD->GetY());
  }
  else if(vtxAOD){
    (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvy"))))->Fill(fCentralityPercentile,vtxAOD->GetY());
    if(TMath::Abs(vtxAOD->GetY()) > fVyMax) aEventCuts[iStep] = 1; //4
    else (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvyA"))))->Fill(fCentralityPercentile,vtxAOD->GetY());
  }
  else aEventCuts[iStep] = 1; //4


 ++iStep;
  if (vtxESD){
    (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvz"))))->Fill(fCentralityPercentile,vtxESD->GetZ());
    if(TMath::Abs(vtxESD->GetZ()) > fVzMax) aEventCuts[iStep] = 1;  //5
    else (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvzA"))))->Fill(fCentralityPercentile,vtxESD->GetZ());
  }
  else if(vtxAOD){
    (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvz"))))->Fill(fCentralityPercentile,vtxAOD->GetZ());
    if(TMath::Abs(vtxAOD->GetZ()) > fVzMax) aEventCuts[iStep] = 1; //5
    else (static_cast<TH2F*>(fQaList->FindObject(Form("fHistQAvzA"))))->Fill(fCentralityPercentile,vtxAOD->GetZ());
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
    else
      (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat0"))))->Fill(idx);
  }
  for (Int_t idx = 0; idx < fHEventStatMax; ++idx) {
    if (aEventCuts[idx])
      break;
    (static_cast<TH1F*>(fQaList->FindObject(Form("hEventStat1"))))->Fill(idx);
  }
  if (!isRejected) {

    (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityStat"))))->Fill(fCentralityBin);
    (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityPercentileAccepted"))))->Fill(fCentralityPercentile);
  }
  (static_cast<TH1F*>(fQaList->FindObject(Form("hCentralityPercentileAll"))))->Fill(fCentralityPercentile);
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

  Int_t    bhepmc[8] = {fGNBinsCent,fGNBinsSign,2, 2, 2,fGNBinsRap,fGNBinsPhi,fGNBinsPt};
  Double_t mnhepmc[8] = {fGRngCent[0],fGRngSign[0],-0.5,-0.5,-0.5,fGRngRap[0],fGRngPhi[0],fGRngPt[0]};  
  Double_t mxhepmc[8] = {fGRngCent[1],fGRngSign[1],1.5,1.5,1.5,fGRngRap[1],fGRngPhi[1],fGRngPt[1]};  
  
  TString titilemc        = "cent:signMC:findable:recStatus:pidStatus:yMC:phiMC:ptMC";
  TString tiltlelaxmc[8]  = {"Centrality", "sign", "findable","recStatus","recPid","#it{y}_{MC}", 
			     "#varphi_{MC} (rad)","#it{p}_{T,MC} (GeV/#it{c})"};
  
  fEffList->Add(new THnSparseF("hmNpiEffMc",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fEffList->Add(new THnSparseF("hmNkaEffMc",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fEffList->Add(new THnSparseF("hmNprEffMc",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fEffList->Add(new THnSparseF("hmNchEffMc",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  
 for (Int_t i = 0; i < 8; i++) { 
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffMc"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
   static_cast<THnSparseF*>(fEffList->FindObject("hmNpiEffMc"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
   static_cast<THnSparseF*>(fEffList->FindObject("hmNkaEffMc"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
   static_cast<THnSparseF*>(fEffList->FindObject("hmNprEffMc"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());

  }
  
  Int_t    binhnep[5] = {fGNBinsCent,fGNBinsSign,fGNBinsRap,fGNBinsPhi,fGNBinsPt};
  Double_t minhnep[5] = {fGRngCent[0],fGRngSign[0],fGRngRap[0],fGRngPhi[0],fGRngPt[0]};
  Double_t maxhnep[5] = {fGRngCent[1],fGRngSign[1],fGRngRap[1],fGRngPhi[1],fGRngPt[1]};


  TString titilerec        = "cent:signRec:yRec:phiRec:ptRec";
  TString tiltlelaxrec[5]  = {"Centrality", "sign", "#it{y}_{Rec}", "#varphi_{Rec} (rad)","#it{p}_{T,Rec} (GeV/#it{c})"};

  fEffList->Add(new THnSparseF("hmNpiEffRec",titilerec.Data(),5,binhnep, minhnep, maxhnep));
  fEffList->Add(new THnSparseF("hmNkaEffRec",titilerec.Data(),5,binhnep, minhnep, maxhnep));
  fEffList->Add(new THnSparseF("hmNprEffRec",titilerec.Data(),5,binhnep, minhnep, maxhnep));
  fEffList->Add(new THnSparseF("hmNchEffRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));

  for (Int_t i = 0; i < 5; i++) { 
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
  }  


  //----
  Int_t    binHnCont[6] = {fGNBinsCent,fGNBinsSign, 5,fGNBinsRap,fGNBinsPhi, fGNBinsPt};  
  Double_t minHnCont[6] = {fGRngCent[0],fGRngSign[0],-0.5,fGRngRap[0],fGRngPhi[0], fGRngPt[0]};
  Double_t maxHnCont[6] = {fGRngCent[1],fGRngSign[1],4.5,fGRngRap[1],fGRngPhi[1], fGRngPt[1]};   

  TString titilecont     = "cent:signMC:contStatus:yMC:phiMC:ptMC";
  TString tiltlelaxcont[6] 
    = {"Centrality","sign","contPart","#it{y}_{MC}","#varphi_{MC} (rad)","#it{p}_{T,MC} (GeV/#it{c})"};
  
  fEffList->Add(new THnSparseF("hmNpiContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont));
  fEffList->Add(new THnSparseF("hmNkaContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont));
  fEffList->Add(new THnSparseF("hmNprContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont));
  fEffList->Add(new THnSparseF("hmNchContMc",titilecont.Data(),6,binHnCont,minHnCont, maxHnCont));

 for (Int_t i = 0; i < 6; i++) {  
    static_cast<THnSparseF*>(fEffList->FindObject("hmNchContMc"))->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNpiContMc"))->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNkaContMc"))->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
    static_cast<THnSparseF*>(fEffList->FindObject("hmNprContMc"))->GetAxis(i)->SetTitle(tiltlelaxcont[i].Data());
  }

 fEffList->Add(new THnSparseF("hmNpiContRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));
 fEffList->Add(new THnSparseF("hmNkaContRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));
 fEffList->Add(new THnSparseF("hmNprContRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));
 fEffList->Add(new THnSparseF("hmNchContRec",titilerec.Data(),5,binhnep,minhnep, maxhnep));

 for (Int_t i = 0; i < 5; i++) {  
   static_cast<THnSparseF*>(fEffList->FindObject("hmNchContRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
   static_cast<THnSparseF*>(fEffList->FindObject("hmNpiContRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
   static_cast<THnSparseF*>(fEffList->FindObject("hmNkaContRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
   static_cast<THnSparseF*>(fEffList->FindObject("hmNprContRec"))->GetAxis(i)->SetTitle(tiltlelaxrec[i].Data());
 }  

 // Magic
 /*
 Double_t *binsPt = 0;
 binsPt = CreateLogAxis(fGNBinsPt,fGRngPt[0],fGRngPt[1]);
 
 static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffMc"))->SetBinEdges(7,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNpiEffMc"))->SetBinEdges(7,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNkaEffMc"))->SetBinEdges(7,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNprEffMc"))->SetBinEdges(7,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->SetBinEdges(4,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->SetBinEdges(4,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->SetBinEdges(4,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNchEffRec"))->SetBinEdges(4,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNchContMc"))->SetBinEdges(5,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNpiContMc"))->SetBinEdges(5,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNkaContMc"))->SetBinEdges(5,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNprContMc"))->SetBinEdges(5,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNchContRec"))->SetBinEdges(4,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNpiContRec"))->SetBinEdges(4,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNkaContRec"))->SetBinEdges(4,binsPt);
 static_cast<THnSparseF*>(fEffList->FindObject("hmNprContRec"))->SetBinEdges(4,binsPt);
 */
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
   if(!AcceptTrack(track)) continue;
   

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
    
    fCurCont[0] = fCentralityBin;
    fCurCont[3] = particle->Y();
    fCurCont[4] = particle->Phi();
    fCurCont[5] = particle->Pt();
    FillCC(b);
    fCurRec[0] = fCentralityBin;
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
   
   if (!AcceptTrackMC(particle, idxMC)) continue;
   
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
	  fCurRec[0] = fCentralityBin;
	  fCurRec[1] = track->Charge();
	  fCurRec[2] = rap;
	  fCurRec[3] = track->Phi();
	  fCurRec[4] = track->Pt();
	  FillRecCE(iPid);
	}      
        break;
      }
    } 
   fCurGen[0] = fCentralityBin;
   fCurGen[6] = particle->Phi();
   fCurGen[7] = particle->Pt();

  
   FillGenCE(iPid);
  } 
  

 delete [] gLabelsZer;
 delete [] gLabelsOne;


}

void AliEbyENetChargeFluctuationTask::CreateDED() {

 Int_t    bhepmc[7] =  {fGNBinsCent, fGNBinsSign,   2,   fGNBinsRap, fGNBinsPhi, fGNBinsPt,800};
  Double_t mnhepmc[7] = {fGRngCent[0],fGRngSign[0],-0.5, fGRngRap[0],fGRngPhi[0],fGRngPt[0],-4};  
  Double_t mxhepmc[7] = {fGRngCent[1],fGRngSign[1], 1.5, fGRngRap[1],fGRngPhi[1], fGRngPt[1], 4.};  
  TString titilemc        = "cent:sign:accepted:y:phi:pt:dcar";

  TString tiltlelaxmc[7]  = {"Centrality", "sign", "Is Accepted","#it{y}","#varphi (rad)","#it{p}_{T} (GeV/#it{c})", "DCAr"};
  
  fDcaList->Add(new THnSparseF("hmNpiDcaRec",titilemc.Data(),7,bhepmc,mnhepmc, mxhepmc));
  fDcaList->Add(new THnSparseF("hmNkaDcaRec",titilemc.Data(),7,bhepmc,mnhepmc, mxhepmc));
  fDcaList->Add(new THnSparseF("hmNprDcaRec",titilemc.Data(),7,bhepmc,mnhepmc, mxhepmc));
  titilemc        = "cent:sign:accepted:eta:phi:pt:dcar";
  fDcaList->Add(new THnSparseF("hmNchDcaRec",titilemc.Data(),7,bhepmc,mnhepmc, mxhepmc));
   
  for (Int_t i = 0; i < 7; i++) { 
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNpiDcaRec"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNkaDcaRec"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNprDcaRec"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    if (i == 4) tiltlelaxmc[3] = "#eta";
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNchDcaRec"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
  }
  /*
  Double_t *binsPt = 0;
  binsPt = CreateLogAxis(fGNBinsPt,fGRngPt[0],fGRngPt[1]);
  static_cast<THnSparseF*>(fDcaList->FindObject("hmNpiDca"))->SetBinEdges(6,binsPt);
  static_cast<THnSparseF*>(fDcaList->FindObject("hmNkaDca"))->SetBinEdges(6,binsPt);
  static_cast<THnSparseF*>(fDcaList->FindObject("hmNprDca"))->SetBinEdges(6,binsPt);
  static_cast<THnSparseF*>(fDcaList->FindObject("hmNchDca"))->SetBinEdges(6,binsPt);
  */

}
void AliEbyENetChargeFluctuationTask::CreateDEM() {

  Int_t    bhepmc[8]  = {fGNBinsCent, fGNBinsSign,  2,     3, fGNBinsRap, fGNBinsPhi,fGNBinsPt,800};
  Double_t mnhepmc[8] = {fGRngCent[0],fGRngSign[0],-0.5, 0.5, fGRngRap[0],fGRngPhi[0],fGRngPt[0],-4};  
  Double_t mxhepmc[8] = {fGRngCent[1],fGRngSign[1],1.5,  3.5, fGRngRap[1],fGRngPhi[1],fGRngPt[1], 4.};  
  TString titilemc    = "cent:sign:cont:accepted:y:phi:pt:dcar";
  TString tiltlelaxmc[8]  = {"Centrality", "sign", "Is Accepted", "1 primary | 2 from WeakDecay | 3 p from Material",
			     "#it{y}","#varphi (rad)","#it{p}_{T} (GeV/#it{c})", "DCAr"};
  
  fDcaList->Add(new THnSparseF("hmNpiDca",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fDcaList->Add(new THnSparseF("hmNkaDca",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
  fDcaList->Add(new THnSparseF("hmNprDca",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));

  titilemc        = "cent:sign:accepted:cont:eta:phi:pt:dcar";
  fDcaList->Add(new THnSparseF("hmNchDca",titilemc.Data(),8,bhepmc,mnhepmc, mxhepmc));
   
  for (Int_t i = 0; i < 8; i++) { 
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNpiDca"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNkaDca"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNprDca"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
    if (i == 4) tiltlelaxmc[4] = "#eta";
    static_cast<THnSparseF*>(fDcaList->FindObject("hmNchDca"))->GetAxis(i)->SetTitle(tiltlelaxmc[i].Data());
  }
  /*
  Double_t *binsPt = 0;
  binsPt = CreateLogAxis(fGNBinsPt,fGRngPt[0],fGRngPt[1]);
  static_cast<THnSparseF*>(fDcaList->FindObject("hmNpiDca"))->SetBinEdges(6,binsPt);
  static_cast<THnSparseF*>(fDcaList->FindObject("hmNkaDca"))->SetBinEdges(6,binsPt);
  static_cast<THnSparseF*>(fDcaList->FindObject("hmNprDca"))->SetBinEdges(6,binsPt);
  static_cast<THnSparseF*>(fDcaList->FindObject("hmNchDca"))->SetBinEdges(6,binsPt);
  */
}

void AliEbyENetChargeFluctuationTask::CalculateDED(Int_t gPid) {

  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    AliVTrack *track = (fESD) ? 
      static_cast<AliVTrack*>(fESD->GetTrack(idxTrack)) : 
      static_cast<AliVTrack*>(fAOD->GetTrack(idxTrack)); 

    if (track->Charge() == 0) continue; // No place for you my friend

    
    if(!AcceptTrack(track)) {fCurRecD[2] = 0.;}       // 2
    else fCurRecD[2] = 1.;   
    
    fCurRecD[1] = track->Charge() < 0 ? 0. : 1.;      // 1
    
    Int_t b = 0;
    
    if (gPid != 0) {
      Int_t a = fHelperPID->GetParticleSpecies(track,kFALSE);
      if(a < 0 || a > 2)  continue;
      b = a + 1;
    }
    if (gPid != b) continue;
    Double_t rap;
    if (!TrackRapidity(track,rap,b)) fCurRecD[2] = 0.; //2
    else fCurRecD[2] = 1.;
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

    fCurRecD[0] = fCentralityBin;                // 0
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

    if (track->Charge() == 0) continue; // No place for you my friend

    if(!AcceptTrack(track)) fCurGenD[2] = 0.;             // 2
    else fCurGenD[2] = 1.;                                // 2
    
    fCurGenD[1] = track->Charge() < 0 ? 0. : 1.;            // 1
    
    Int_t b = 0;
    
    if (gPid != 0) {
      Int_t a = fHelperPID->GetParticleSpecies(track,kFALSE);
      if(a < 0 || a > 2)  continue;
      b = a + 1;
    }
    
    if (gPid != b) continue;
    
    Double_t rap;
    if (!TrackRapidity(track,rap,b)) fCurGenD[2] = 0.;   // 2
    else fCurGenD[2] = 1.;                               // 2
    
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
    
    // Printf(" <<< %f  %f ", dca[0],dca[1]);     

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
    fCurGenD[0] = fCentralityBin;                    // 0
    fCurGenD[5] = track->Phi();                      // 5
    fCurGenD[6] = track->Pt();                       // 6
 
    FillRecDE(gPid);
  }
}
