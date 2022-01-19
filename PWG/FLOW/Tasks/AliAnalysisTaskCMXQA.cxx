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

/* $Id: AliAnalysisTaskCMXQA.cxx  Rihan Haque, 18/09/2019 (ver1) $ */

//-- general include---
#include "TChain.h"
#include "TTree.h"
#include "TGrid.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TMatrixDSym.h"

#include "TMath.h"
#include "stdio.h"
#include "Riostream.h"

//---- manager and handler---
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

//---V0 and ZDC info---
#include "AliAODZDC.h"
#include "AliAODVZERO.h"
#include "AliAODVertex.h"

//---AOD,ESD event--
#include "AliESDEvent.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"

//----for PID-----
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

//----- Vevent and tracks
#include "AliVEventHandler.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliCentrality.h"

//----- must include-------
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliPhysicsSelection.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskCMXQA.h"

using std::cout;
using std::endl;
using std::vector;


ClassImp(AliAnalysisTaskCMXQA)

AliAnalysisTaskCMXQA::AliAnalysisTaskCMXQA(const char *name): AliAnalysisTaskSE(name),
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fV0CutPU(NULL),
  fSPDCutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL),  
  fCentralityMin(0),
  fCentralityMax(90),
  fFilterBit(1),
  fTPCclustMin(70),
  bUseKinkTracks(kFALSE),
  fNSigmaCut(2.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),    
  fdEdxMin(10.0),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  sCentrEstimator("V0M"),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fHistEventCount(NULL)
{

  for(int i=0;i<3;i++){
    for(int j=0;j<4;j++){
    fHistEtaPtwCutChPos[i][j]  = NULL;
    fHistEtaPhiwCutChPos[i][j] = NULL;
    fHistEtaPtwCutChNeg[i][j]  = NULL;
    fHistEtaPhiwCutChNeg[i][j] = NULL;
    fHistEtaVzwCutChPos[i][j]  = NULL;
    fHistEtaVzwCutChNeg[i][j]  = NULL;
    }
  }
  
  //Must be here:
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());
}

//_______________________empty constructor_______________________
AliAnalysisTaskCMXQA::AliAnalysisTaskCMXQA():
  fVevent(NULL),
  fESD(NULL),
  fAOD(NULL),
  fPIDResponse(NULL),
  fMultSelection(NULL),
  fAnalysisUtil(NULL),
  fListHist(NULL),
  fV0CutPU(NULL),
  fSPDCutPU(NULL),
  fMultCutPU(NULL),
  fCenCutLowPU(NULL),
  fCenCutHighPU(NULL), 
  fCentralityMin(0),
  fCentralityMax(90),
  fFilterBit(1),
  fTPCclustMin(70),
  bUseKinkTracks(kFALSE),
  fNSigmaCut(2.0),
  fMinPtCut(0.2),
  fMaxPtCut(5.0),
  fMinEtaCut(-0.8),
  fMaxEtaCut(0.8),
  fTrkChi2Min(0.1),    
  fdEdxMin(10.0),
  fMinVzCut(-10.0),
  fMaxVzCut(10.0),
  sCentrEstimator("V0M"),
  fCentDistBeforCut(NULL),
  fCentDistAfterCut(NULL),
  fHistEventCount(NULL)
{

  for(int i=0;i<3;i++){
    for(int j=0;j<4;j++){
      fHistEtaPtwCutChPos[i][j]  = NULL;
      fHistEtaPhiwCutChPos[i][j] = NULL;
      fHistEtaPtwCutChNeg[i][j]  = NULL;
      fHistEtaPhiwCutChNeg[i][j] = NULL;
      fHistEtaVzwCutChPos[i][j]  = NULL;
      fHistEtaVzwCutChNeg[i][j]  = NULL;
    }
  }
}
  
//__________________ destructor ___________________
AliAnalysisTaskCMXQA::~AliAnalysisTaskCMXQA()
{
  if(fListHist)      delete fListHist;  
  if(fAnalysisUtil)  delete fAnalysisUtil;   // because its 'new' !!

  if(fV0CutPU)       delete fV0CutPU;
  if(fSPDCutPU)      delete fSPDCutPU;
  if(fMultCutPU)     delete fMultCutPU;
  if(fCenCutLowPU)   delete fCenCutLowPU; 
  if(fCenCutHighPU)  delete fCenCutHighPU;   
  
}










//________________ Define Histograms _______________
void AliAnalysisTaskCMXQA::UserCreateOutputObjects()
{
  //std::cout<<"\n UserCreateOutputObject: function begins...\n"<<endl; 

  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
  
  if (!inputHandler) {  printf("\n ...Error: Input handler missing, QUITing Job..!!\n");    return;}

  //PileUp Multi-Vertex
  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);
  
  ///obtain the PID response object if needed:
  fPIDResponse=inputHandler->GetPIDResponse();
    
    
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);

  SetupEventAndTaskConfigInfo();

  fCentDistBeforCut = new TH1F("fCentDistBeforCut","Cent w/o any Cuts; Cent (%); no.Events ",100,0,100);
  fListHist->Add(fCentDistBeforCut);
  fCentDistAfterCut = new TH1F("fCentDistAfterCut","Cent with all Cut; Cent (%); no.Events ",100,0,100);
  fListHist->Add(fCentDistAfterCut);

  TString sPIDnames[4] = {"Chrg","Pion","Kaon","Prot"};

  
  for(int i=0;i<3;i++){
    for(int j=0;j<4;j++){

      const char *part = sPIDnames[j].Data();
      ///eta-pt
      fHistEtaPtwCutChPos[i][j]  = new TH2F(Form("fHistEtaPtwCutChPosCent%dk%s",i,part),Form("(FB%d, wCut);#eta; p_{T};",fFilterBit),32,-0.8,0.8,50,0.2,5.2);
      fListHist->Add(fHistEtaPtwCutChPos[i][j]);
      fHistEtaPtwCutChNeg[i][j]  = new TH2F(Form("fHistEtaPtwCutChNegCent%dk%s",i,part),Form("(FB%d, wCut);#eta; p_{T};)",fFilterBit),32,-0.8,0.8,50,0.2,5.2);
      fListHist->Add(fHistEtaPtwCutChNeg[i][j]);
      ///Eta-phi
      fHistEtaPhiwCutChPos[i][j] = new TH2F(Form("fHistEtaPhiwCutChPosCent%dk%s",i,part),Form("(FB%d, wCut);#phi;#eta;",fFilterBit),100,0,6.284,16,-0.8,0.8);
      fListHist->Add(fHistEtaPhiwCutChPos[i][j]);
      fHistEtaPhiwCutChNeg[i][j] = new TH2F(Form("fHistEtaPhiwCutChPosCent%dk%s",i,part),Form("(FB%d, wCut);#phi;#eta;",fFilterBit),100,0,6.284,16,-0.8,0.8);
      fListHist->Add(fHistEtaPhiwCutChNeg[i][j]);
      ///Eta-Vz:
      fHistEtaVzwCutChPos[i][j] = new TH2F(Form("fHistEtaVzwCutChPosCent%dk%s",i,part),Form("(FB%d, wCut); V_{z} cm; #eta;",fFilterBit),40,-10,10,32,-0.8,0.8);
      fListHist->Add(fHistEtaVzwCutChPos[i][j]);
      fHistEtaVzwCutChNeg[i][j] = new TH2F(Form("fHistEtaVzwCutChNegCent%dk%s",i,part),Form("(FB%d, wCut); V_{z} cm; #eta;",fFilterBit),40,-10,10,32,-0.8,0.8);
      fListHist->Add(fHistEtaVzwCutChNeg[i][j]);
    }
  }




  ///Call this function to setup the pileup removal functions. 
  SetupPileUpRemovalFunctions();

  
  
  PostData(1,fListHist);
 //std::cout<<"\n UserCreateOutputObject: function Ends...\n"<<endl;
}










//____________________________ Call Event by Event ___________________________________
void AliAnalysisTaskCMXQA::UserExec(Option_t*) {
 
  //cout<<"\n Info:UserExec() called ..!!!\n";
  //watch.Start(kTRUE);
  
  Float_t stepCount = 0.5;

  fHistEventCount->Fill(stepCount);   //1
  stepCount++;

  fAOD = dynamic_cast <AliAODEvent*> (InputEvent());
  fESD = dynamic_cast <AliESDEvent*> (InputEvent());
  if(!(fESD || fAOD)) { printf("ERROR: fESD & fAOD not available\n"); return;  }

  
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent) { printf("ERROR: fVevent not available\n");  return;  }
   
  fHistEventCount->Fill(stepCount); //2
  stepCount++;

  //-------------- Vtx cuts ---------------
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t pVtxZ = -999;
  pVtxZ  = pVtx->GetZ();
    
  if(pVtxZ < fMinVzCut || pVtxZ > fMaxVzCut)    return;
  
  fHistEventCount->Fill(stepCount); //3
  stepCount++;

  Float_t centrality = -99.0;
  Float_t centrV0M   = -99.0;
  Float_t centrCL1   = -99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");  // Must never comment this
  if(!fMultSelection) { printf("\n\n **WARNING** \n::UserExec()::\n AliMultSelection object not found. Quit!! \n"); exit(1); }


  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  //centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

 
  centrality = centrV0M;  // This Is Always Default, changes below for other options:
  
  if(sCentrEstimator=="CL0"){
    centrality = fMultSelection->GetMultiplicityPercentile("CL0");;
  }
  else if(sCentrEstimator=="CL1"){
    centrality = fMultSelection->GetMultiplicityPercentile("CL1");;
  }
  else if(sCentrEstimator=="V0C"){
    centrality = fMultSelection->GetMultiplicityPercentile("V0C");
  }
  else if(sCentrEstimator=="V0A"){
    centrality = fMultSelection->GetMultiplicityPercentile("V0A");
  }
  else if(sCentrEstimator=="TRK"){
    centrality = fMultSelection->GetMultiplicityPercentile("TRK");
  }
  else{
    centrality = centrV0M;
  }

  
  fCentDistBeforCut->Fill(centrality);
  
  if(centrality<fCentralityMin || centrality>fCentralityMax){ 
    return;
  }

  fHistEventCount->Fill(stepCount); //4
  stepCount++;

  Int_t ntracks = fAOD->GetNumberOfTracks();
  if(ntracks < 4) return;                        // Minimum 4 tracks per event. 
    
  fHistEventCount->Fill(stepCount); //5
  stepCount++;

  //////----> Get Magnetic field and RunNo.---------
  // Float_t fMagField = fAOD->GetMagneticField();
  // const Int_t QAindex = (fMagField > 0) ? 1 : 0;
  // Int_t runNumber = fAOD->GetRunNumber();
  //------------------------------------------------

  Float_t fMultTPCFull = 0;  // TPC mult estimate
  Float_t fMultGlobal  = 0;  // global track multiplicity
  Float_t fMultwRawFB  = 0;  // Uncorrected Multiplicity
  Float_t fMultCorrFB  = 0;  // Corrected Multiplicity

  
  Int_t   gMultEtaNeg  = 0;
  Int_t   gMultEtaPos  = 0;
  Int_t   gMultEtaAll  = 0;

  Float_t trkPt=0,trkPhi=0,trkEta=0;
  Float_t trkChi2=0,trkdEdx=0;
  Int_t   trkChrg=0, trkTpcNC=0;

  //Prottay's variable:
  Double_t Q2xEtaPos=0., Q2yEtaPos=0.;
  Double_t Q2xEtaNeg=0., Q2yEtaNeg=0.;  


  Int_t icent = 0;

  //if(centrality<=5.0)     icent = 0;
  //else if(centrality<=10) icent = 1;
  //else{
  //icent = (int) centrality/10 + 1;
  //}

  if(centrality<=10) icent = 0;
  else if(centrality<=40) icent = 1;
  else if(centrality<=90) icent = 2;
  else  icent = 10;

  
  if(icent>9) return;  //only fill QA Upto 90% centrality!


  Bool_t kPileupEvent = CheckPileUp2018(fAOD);
  if(kPileupEvent) return;

  Bool_t isItPiontrk1 = kFALSE;
  Bool_t isItKaontrk1 = kFALSE;
  Bool_t isItProttrk1 = kFALSE;

  
    
  ////////-----> Starting 1st track Loop ----->
  
  for(Int_t iTrack = 0; iTrack < ntracks; iTrack++) { 

    AliAODTrack* AODtrack = dynamic_cast <AliAODTrack*> (fVevent->GetTrack(iTrack));
    if(!AODtrack) continue;
    
    if(AODtrack->TestFilterBit(fFilterBit)){

      if(Int_t(AODtrack->GetProdVertex()->GetType()) == AliAODVertex::kKink) continue; /// remove kink Tracks.
      
      trkPt    = AODtrack->Pt();
      trkPhi   = AODtrack->Phi();
      trkEta   = AODtrack->Eta();
      trkChrg  = AODtrack->Charge();
      trkChi2  = AODtrack->Chi2perNDF();
      trkTpcNC = AODtrack->GetTPCNcls();
      trkdEdx  = AODtrack->GetDetPid()->GetTPCsignal();     ///Check After Filter bit is validated!! Otherwise code breaks.!
      
      //Apply track cuts here:
      if((trkPt <= fMaxPtCut) && (trkPt >= fMinPtCut) && (trkEta <= fMaxEtaCut) && (trkEta >= fMinEtaCut) && (trkdEdx >= fdEdxMin) && (trkTpcNC >= fTPCclustMin) && (trkChi2 >= fTrkChi2Min) && (trkChi2 <= 4.0) && TMath::Abs(trkChrg)) {
      ///if((AODtrack->DCA() <= fDCAxyMax) && (AODtrack->ZAtDCA() <= fDCAzMax)) /// DCA cut not applied. Cuts are applied in FB.



	isItPiontrk1 = kFALSE;
	isItKaontrk1 = kFALSE;
	isItProttrk1 = kFALSE; 

	isItPiontrk1 = CheckPIDofParticle(AODtrack,1); // 1=pion
	isItKaontrk1 = CheckPIDofParticle(AODtrack,2); // 2=Kaon
	isItProttrk1 = CheckPIDofParticle(AODtrack,3); // 3=proton
    
	if(trkChrg > 0){
	  
	  fHistEtaPtwCutChPos[icent][0]->Fill(trkEta, trkPt);
	  fHistEtaPhiwCutChPos[icent][0]->Fill(trkPhi,trkEta);
	  fHistEtaVzwCutChPos[icent][0]->Fill(pVtxZ,trkEta);

	  if(isItPiontrk1){
	    fHistEtaPtwCutChPos[icent][1]->Fill(trkEta, trkPt);
	    fHistEtaPhiwCutChPos[icent][1]->Fill(trkPhi,trkEta);
	    fHistEtaVzwCutChPos[icent][1]->Fill(pVtxZ,trkEta);
	  }
	  if(isItKaontrk1){
	    fHistEtaPtwCutChPos[icent][2]->Fill(trkEta, trkPt);
	    fHistEtaPhiwCutChPos[icent][2]->Fill(trkPhi,trkEta);
	    fHistEtaVzwCutChPos[icent][2]->Fill(pVtxZ,trkEta);
	  }
	  if(isItProttrk1){
	    fHistEtaPtwCutChPos[icent][3]->Fill(trkEta, trkPt);
	    fHistEtaPhiwCutChPos[icent][3]->Fill(trkPhi,trkEta);
	    fHistEtaVzwCutChPos[icent][3]->Fill(pVtxZ,trkEta);
	  }
	}
	else if(trkChrg < 0){

	  fHistEtaPtwCutChNeg[icent][0]->Fill(trkEta, trkPt);
	  fHistEtaPhiwCutChNeg[icent][0]->Fill(trkPhi,trkEta);
	  fHistEtaVzwCutChNeg[icent][0]->Fill(pVtxZ,trkEta);
	  
	  if(isItPiontrk1){
	    fHistEtaPtwCutChNeg[icent][1]->Fill(trkEta, trkPt);
	    fHistEtaPhiwCutChNeg[icent][1]->Fill(trkPhi,trkEta);
	    fHistEtaVzwCutChNeg[icent][1]->Fill(pVtxZ,trkEta);
	  }
	  if(isItKaontrk1){
	    fHistEtaPtwCutChNeg[icent][2]->Fill(trkEta, trkPt);
	    fHistEtaPhiwCutChNeg[icent][2]->Fill(trkPhi,trkEta);
	    fHistEtaVzwCutChNeg[icent][2]->Fill(pVtxZ,trkEta);
	  }
	  if(isItProttrk1){
	    fHistEtaPtwCutChNeg[icent][3]->Fill(trkEta, trkPt);
	    fHistEtaPhiwCutChNeg[icent][3]->Fill(trkPhi,trkEta);
	    fHistEtaVzwCutChNeg[icent][3]->Fill(pVtxZ,trkEta);
	  }	  
	}



	
      }//with all trackCuts applied      
    }//-------> if FB is validated    
  }///------> 1st Track loop Ends here.<--------

  
  //// Fill E-by-E quantities:



  
  //Last lines in Event loop
  fCentDistAfterCut->Fill(centrality);
  fHistEventCount->Fill(14.5); //final Event count.
  PostData(1,fListHist);
  
}//---------------- UserExec ----------------------



Bool_t AliAnalysisTaskCMXQA::CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck){
  
  if(pidToCheck==0) return kTRUE;    //// Charge Particles do not need PID check
  
  Bool_t bPIDokay = kFALSE;

  if(!fPIDResponse){
    Printf("\n Could Not access PIDResponse Task, please Add the Task...\n return with kFALSE pid\n");
    return kFALSE;
  }

  /// Rihan todo: To set the low pT cuts for nSigmaTPC from AddTaskMacro!
  /// Although someone barely needs to change it given the purity..
  
  Double_t nSigTPC = 0, nSigTOF = 0, nSigRMS = 0;
  Double_t trkPtPID  = ftrack->Pt();
  Int_t trkChargePID = ftrack->Charge();

  ///Pion => 
  if(pidToCheck==1){ 
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kPion);
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kPion);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);

    if(trkPtPID<=0.5 && TMath::Abs(nSigTPC)<=3.0){
      bPIDokay = kTRUE;
    }
    else if(trkPtPID>0.5 && TMath::Abs(nSigRMS)<=3.0){   // TPC-TOF RMS cut for higher pt. *** Please set high pT range in AddTask**
      bPIDokay = kTRUE;
    }
    return bPIDokay;
  }
  ///Kaon => 
  else if(pidToCheck==2){
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kKaon);	
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kKaon);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
	
    if(trkPtPID<=0.45 && TMath::Abs(nSigTPC)<=3.0){
      bPIDokay = kTRUE;
    }
    else if(trkPtPID>0.45 && TMath::Abs(nSigRMS)<=3.0){  // TPC-TOF RMS cut for higher pt. *** Please set high pT range in AddTask**
      bPIDokay = kTRUE;
    }
    return bPIDokay;
  }
  ///proton => 
  else if(pidToCheck==3){///
    nSigTPC = fPIDResponse->NumberOfSigmasTPC(ftrack, AliPID::kProton);    
    nSigTOF = fPIDResponse->NumberOfSigmasTOF(ftrack, AliPID::kProton);
    nSigRMS = TMath::Sqrt(nSigTPC*nSigTPC + nSigTOF*nSigTOF);
	
    if(trkPtPID<=0.6 && TMath::Abs(nSigTPC)<=3.0){
      bPIDokay = kTRUE;
      if(trkChargePID>0 && trkPtPID<0.4) bPIDokay = kFALSE;  
    }
    else if(trkPtPID>0.6 && TMath::Abs(nSigRMS)<=3.0){
      bPIDokay = kTRUE;
    }
    return bPIDokay;
  }
  else{
    Printf("\n -Ve number not allowed! Choose among: 0,1,2,3,4  (Charge Pion, Kaon, Proton, Lambda)\n return with kFALSE \n");
    return kFALSE;
  }  
  return kFALSE;
}





void AliAnalysisTaskCMXQA::SetupEventAndTaskConfigInfo(){
  fHistEventCount = new TH1F("fHistEventCount","Event counts",15,0,15);
  fHistEventCount->GetXaxis()->SetBinLabel(1,"Called UserExec()");
  fHistEventCount->GetXaxis()->SetBinLabel(2,"Called Exec()");
  fHistEventCount->GetXaxis()->SetBinLabel(3,"AOD Exist");
  fHistEventCount->GetXaxis()->SetBinLabel(4,"Vz < 10");
  fHistEventCount->GetXaxis()->SetBinLabel(5,Form("%2.0f<Cent<%2.0f",fCentralityMin,fCentralityMax));
  fHistEventCount->GetXaxis()->SetBinLabel(6,"noAODtrack > 2 ");
  fHistEventCount->GetXaxis()->SetBinLabel(7,"TPC vs Global");
  fHistEventCount->GetXaxis()->SetBinLabel(8,"TPC128 vs ESD");
  fHistEventCount->GetXaxis()->SetBinLabel(9,"Cent vs TPC");
  fHistEventCount->GetXaxis()->SetBinLabel(10,"mult eta+/- > 2");
  fHistEventCount->GetXaxis()->SetBinLabel(11,"centCL1 < 90");
  fHistEventCount->GetXaxis()->SetBinLabel(15,"Survived Events");
  fListHist->Add(fHistEventCount);
  //fHistEventCount->Fill(1);
}//----------SetupEventAndTaskConfigInfo-----------


Int_t AliAnalysisTaskCMXQA::GetCentralityScaled0to10(Double_t fCent){

 Int_t cIndex = 0;

 if(fCent<5.0) {
   cIndex  = 0; 
 }
 else if(fCent>=5.0 && fCent<10){
   cIndex  = 1;
 }
 else if(fCent>=10.0) {
   cIndex = abs(fCent/10.0)+1;
 }
 return cIndex;
 
}//------------GetCentralityScaled0to10------------






Bool_t AliAnalysisTaskCMXQA::CheckPileUp2018(AliAODEvent *faod){

  Bool_t BisPileup=kFALSE;
  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");

  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(111);
  }

  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = faod->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = faod->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)faod->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
  
  const Int_t nTracks = faod->GetNumberOfTracks();
  Int_t multTrk = 0;
    
  for (Int_t it = 0; it < nTracks; it++) {
    AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);
    if (!aodTrk)  continue;
    
    if (aodTrk->TestFilterBit(32)){
      if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
    }        
  }//AOD track loop
  
  AliAODVZERO* aodV0 = faod->GetVZEROData();
  Float_t  multV0a   = aodV0->GetMTotV0A();
  Float_t  multV0c   = aodV0->GetMTotV0C();
  Float_t  multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On  = multV0aOn + multV0cOn;

  Int_t tpcClsTot = faod->GetNumberOfTPCClusters();
  Float_t nclsDif = Float_t(tpcClsTot) - (60932.9 + 69.2897*multV0Tot - 0.000217837*multV0Tot*multV0Tot);

  if (centrCL0 < fCenCutLowPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (centrCL0 > fCenCutHighPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) {
    BisPileup=kTRUE;
  }     
  if (multV0On < fV0CutPU->Eval(multV0Tot)) {
    BisPileup=kTRUE;
  }
  if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
    BisPileup=kTRUE;
  }
  if (faod->IsIncompleteDAQ()) {
    BisPileup=kTRUE;
  }
  
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;

  Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();
  
  /*
  fHistCentCL0VsV0MBefore->Fill(centrV0M,centrCL0);
  fHistTPCVsESDTrkBefore->Fill(multTrk,multEsd);  
  fHistTPConlyVsCL1Before->Fill(centrCL1,multTrk);
  fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);
  if (!BisPileup) {      
    fHistCentCL0VsV0MAfter->Fill(centrV0M,centrCL0);
    fHistTPCVsESDTrkAfter->Fill(multTrk,multEsd);  
    fHistTPConlyVsCL1After->Fill(centrCL1,multTrk);
    fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);      
  }*/
  
  return BisPileup;  

}


void AliAnalysisTaskCMXQA::SetupPileUpRemovalFunctions(){
  
  ////==========> LHC18q/r PileUp Removal---- Do not Remove them !!! -----------
  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  fV0CutPU  = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);
  
  fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);
  
  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  fCenCutLowPU  = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",  0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
  //-----------------------------------------------------------------------------

}


