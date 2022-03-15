#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TObjArray.h"
#include "TGrid.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskGammaDeltaPID.h"
#endif

void AddTaskLambdaPIDCorr(
  Int_t gFilterBit = 768,
  Double_t fV0PtMin = 0.2,
  Bool_t bDaughterUseTOF = kFALSE,
  Double_t fDautPtMax = 2.0,
  Double_t fDautEtaMin = -0.8, 
  Double_t fDautEtaMax = 0.8,
  Double_t fChi2max = 4.0,
  Int_t gNclustTPC = 70,
  Double_t nSigTPC = 3.0, 
  Double_t nSigTOF = 3.0, 
  Bool_t bSkipPileUp = kFALSE, 
  TString sCentEstimator = "V0M", 
  Float_t fVzMin = -10.0, 
  Float_t fVzMax = 10.0, 
  TString sTrigger = "kINT7", 
  Int_t vnHarmonic = 2, 
  TString sDetForEP = "TPC", 
  TString sMCfilePath="alien:///alice/cern.ch/user/m/mhaque/calib2021/efficiencyBothpol18qnew.root", 
  TString sNUAFilePath = "alien:///alice/cern.ch/user/m/mhaque/calib2021/WgtsNUAChargeAndPion_LHC18qPass3_FB768_AlexPU_DeftMode_Sept2021NoAvgQ.root", 
  TString sDetWgtsFile = "alien:///alice/cern.ch/user/m/mhaque/calib2021/CalibV0GainCorrectionLHC18q_Oct2021.root",
  Bool_t bSkipAnalysis = kFALSE,
  Bool_t bUseZDCSpectatorPlane = kTRUE,
  TString sZDCCorrFile = "alien:///alice/cern.ch/user/s/sqiu/RecenteringResultFinal_2018q.root",
  Double_t fV0MassCut = 0.005, 
  Double_t fV0CosPAmin = 0.995, 
  Double_t fV0RapidityMax = 0.5, 
  Double_t fV0DecLengthMin = 3.0, 
  Double_t fV0DecLengthMax = 100, 
  Double_t fV0DCAToPrimVtx = 1.5, 
  Double_t fV0DCADiffDautMax = 1.0, 
  Double_t fV0DautDCAToPrimVtxMin = 0.02, 
  Bool_t bCheckPIDFlow = kTRUE,
  const char *suffix = "")
{

  printf("===================================================================================\n");
  printf("\n                  Initialising Task: AddTaskLambdaPIDCorr                        \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName; // common outfile filename
  list1OutName += ":Results";         // This directory contains result histograms


  TString TaskName;
  TaskName.Form("gTaskLambdaPIDCorr%d_%d_%s",gFilterBit,gNclustTPC,suffix);

  AliAnalysisTaskGammaDeltaPID *taskLambda = new AliAnalysisTaskGammaDeltaPID(TaskName);

  ///-------> Analysis Object Created, now pass the arguments
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    taskLambda->SelectCollisionCandidates(AliVEvent::kMB);
    printf("\n =========> AddTaskCMW::Info() Trigger = kMB  \n");
  }
  else if(sTrigger=="kSemiCentral" || sTrigger=="SemiCentral" || sTrigger=="semicentral"){
    taskLambda->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    printf("\n =========> AddTaskCMW::Info() Trigger = kSemiCentral \n");
  }
  else if(sTrigger=="kCentral" || sTrigger=="Central" || sTrigger=="central"){
    taskLambda->SelectCollisionCandidates(AliVEvent::kCentral);
    printf("\n =========> AddTaskCMW::Info() Trigger = kCentral \n");
  }
  else if(sTrigger=="kAny" || sTrigger=="kAll"){
    taskLambda->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kSemiCentral | AliVEvent::kCentral);
  }
  else{//if trigger==kINT7 or no trigger provided:
    taskLambda->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskCMW::Info() Trigger = kINT7 \n");
  }
  

  ///Set Event cuts:
  taskLambda->SetVzRangeMin(fVzMin);
  taskLambda->SetVzRangeMax(fVzMax);
  taskLambda->SetFlagSkipPileUpCuts(bSkipPileUp);  
  taskLambda->SetFlagSkipAnalysis(bSkipAnalysis);



  std::cout<<"=========> AddTaskCMW::Info() setting Event Plane Det: "<<sDetForEP<<std::endl;
  taskLambda->SetDetectorforEventPlane(sDetForEP);


  
  if(sCentEstimator=="V0" || sCentEstimator=="V0M"){ 
    taskLambda->SetCentralityEstimator("V0M");    
  }
  else{
    taskLambda->SetCentralityEstimator(sCentEstimator);  //  use the Estimator provided in AddTask.
  }



  
  //Set Track cuts:
  taskLambda->SetParticlePID(0);            // inclusive charge. Note: For this specific Task we check PID of individual particles anyway.
  taskLambda->SetPtRangeMin(0.2);           // Charge Pt min
  taskLambda->SetPtRangeMax(10.0);          // Charge Pt max
  taskLambda->SetFilterBit(gFilterBit);
  taskLambda->SetNSigmaCutTPC(nSigTPC);     // 
  taskLambda->SetNSigmaCutTOF(nSigTOF);
  taskLambda->SetTrackCutChi2Min(0.10);
  taskLambda->SetTrackCutdEdxMin(10.0);  
  taskLambda->SetEtaRangeMin(fDautEtaMin);  // inclusive Eta min (same as V0 Daughter Eta)
  taskLambda->SetEtaRangeMax(fDautEtaMax);  // inclusive Eta max (same as V0 Daughter Eta)
  taskLambda->SetTrackCutChi2Max(fChi2max);
  taskLambda->SetFlagUseKinkTracks(kFALSE);
  taskLambda->SetCumulantHarmonic(vnHarmonic);
  taskLambda->SetTrackCutNclusterMin(gNclustTPC);  

  //Use ZDC Plane
  taskLambda->SetUseZDCSpectatorPlane(bUseZDCSpectatorPlane);

  ///----- Functions For Lambda-X correlation  =>   
  Bool_t bFillLambda=kTRUE;
  taskLambda->SetFlagAnalyseLambda(bFillLambda);
  
  /// Re-use Some AddTask Variables:
  taskLambda->SetV0PtMin(fV0PtMin);                        
  taskLambda->SetV0CPAMin(fV0CosPAmin);                             
  taskLambda->SetV0RapidityMax(fV0RapidityMax);                   
  taskLambda->SetV0DecayLengthMax(fV0DecLengthMax);             
  taskLambda->SetV0DecayLengthMin(fV0DecLengthMin);             
  taskLambda->SetV0DCAToPrimVtxMax(fV0DCAToPrimVtx);           
  //V0 Daughter Cut
  taskLambda->SetDaughtersPIDUseTOF(bDaughterUseTOF);
  taskLambda->SetDaughtersNsigma(nSigTPC);       
  taskLambda->SetDaughtersPtMax(fDautPtMax);                 
  taskLambda->SetDaughtersEtaMax(fDautEtaMax);                        
  taskLambda->SetDaughtersTPCNclsMin(gNclustTPC);    /// Note: same TPC ncluster for Charge and V0 Daughter. Could be changed for different values.
  taskLambda->SetV0DcaBetweenDaughtersMax(fV0DCADiffDautMax);
  taskLambda->SetDaughtersDCAToPrimVtxMin(fV0DautDCAToPrimVtxMin);
  //Lambda Mass Cut
  taskLambda->SetMassMean(1.115683);                             
  taskLambda->SetLambdaMassCut(fV0MassCut);
  taskLambda->IsCheckPIDFlow(bCheckPIDFlow);
  //--------------------------------------------------------------------------

  



  if (!gGrid) TGrid::Connect("alien://");
  //========================= Setup Correction Files ======================> 
  TFile *fMCFile = TFile::Open(sMCfilePath,"READ");
  TList *fListMC=NULL;
  
  if(fMCFile) {
    
    fListMC = dynamic_cast <TList*> (fMCFile->FindObjectAny("fMcEffiHij"));
    std::cout<<" \n ==============> TList found for MC wgts.. GOOD! "<<std::endl;

    if(fListMC) {
      taskLambda->SetListForTrkCorr(fListMC); 
    }
    else{
      printf("\n\n *** AddTask::WARNING \n => MC file Exist, But TList Not Found!!! \n AddTask::Info() ===> NO MC Correction!! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING \n => no MC file!!! \n AddTask::Info() ===> NO MC Correction!! \n\n");
  }
 
 //--------------------------------------------------------------------------
  std::cout<<" NUA file Path "<<sNUAFilePath.Data()<<std::endl;
  TFile* fNUAFile = TFile::Open(sNUAFilePath,"READ");
  TList* fListNUA=NULL;

  //if(fNUAFile->IsOpen()) {
  if(fNUAFile){    
    fListNUA = dynamic_cast <TList*> (fNUAFile->FindObjectAny("fNUA_ChPosChNeg"));
    std::cout<<" \n ==============> TList found for NUA, here is all the histograms : "<<std::endl;
    //fListNUA->ls();
    if(fListNUA) {
      taskLambda->SetListForNUACorr(fListNUA);
    }
    else{
      printf("\n\n *** AddTask::WARNING => NUA file Exist,But TList Not Found!!\n AddTask::Info() ===> NO NUA Correction!! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING => NUA file not Found or Wrong path Set in AddTask Macro!! \n\n");
  }
  //-----------------------------------------------------------------------------
 
  TFile* fV0ZDCWgtsFile = TFile::Open(sDetWgtsFile,"READ");
  TList* fListDetWgts=NULL;

  if(fV0ZDCWgtsFile) {    
    fListDetWgts = dynamic_cast <TList*> (fV0ZDCWgtsFile->FindObjectAny("fWgtsV0ZDC"));
    std::cout<<" \n ==============> TList found for V0 wgts.. GOOD! ";
    // fListDetWgts->ls(); 

    if(fListDetWgts) {
      taskLambda->SetListForV0MCorr(fListDetWgts);
    }
    else{
      printf("\n\n *** AddTask::WARNING => V0/ZDC Weights file Exist, But TList Not Found!!");
      printf("\n May be wrong TList name? No Correction for V0/ZDC !! \n\n");
    }
  }
  else{
    printf("\n\n *** AddTask::WARNING => NO File Found for V0/ZDC Wgts!!\n AddTask::Info() ===> No V0/ZDC Correction!! \n\n");
  }
  //-----------------------------------------------------------------------------
  
  TFile* fZDCRecenterFile = TFile::Open(sZDCCorrFile, "READ");
  TList* fListZDCCorr=NULL;
  
  if(fZDCRecenterFile) {
	fListZDCCorr = dynamic_cast <TList*> (fZDCRecenterFile->FindObjectAny("fOutputRecenter"));
    std::cout<<" \n ==============> TList found for ZDC wgts.. GOOD! "<<std::endl;
	
	if(fListZDCCorr) {
	  taskLambda->SetListForZDCCorr(fListZDCCorr);
    }
    else{
	  printf("\n\n *** AddTask::WARNING => ZDC Recentering file Exist, But TList Not Found!!");
      printf("\n May be wrong TList name? No Correction for ZDC recentering !! \n\n");
	}
  
  } else{
	printf("\n\n *** AddTask::WARNING => NO File Found for ZDC recentering!!\n AddTask::Info() ===> No ZDC Correction!! \n\n");
  }
  //=================================================================================



  

  ///---> Now Pass data and containers to Analysis Object ----
 
  mgr->AddTask(taskLambda);                        // connect the task to the analysis manager
  mgr->ConnectInput(taskLambda, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(taskLambda, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTask was Configured properly... <==================\n\n");


}//Task Ends
