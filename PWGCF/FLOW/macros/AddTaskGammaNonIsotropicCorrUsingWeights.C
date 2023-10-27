
void AddTaskGammaNonIsotropicCorrUsingWeights(Int_t whichData = 2018,TString period = "q",Int_t gFilterBit = 768,Double_t fPtMin=0.2,Double_t fPtMax=5.0,Double_t fEtaMin=-0.8, Double_t fEtaMax=0.8,Double_t fChi2max=4.0,Int_t gNclustTPC=70, Bool_t fUseTPCCrossedRows = kFALSE, Int_t fTPCsharedCut = -1, Int_t fparticle=0,Double_t nSigTPC = 3.0, Double_t nSigTOF = 3.0, Bool_t bSkipPileUp=kFALSE, TString sCentEstimator="V0M", Float_t fVzMin = -10.0, Float_t fVzMax = 10.0, Double_t fdcaxy=0, Double_t fdcaz=0, TString sTrigger="kINT7", Int_t vnHarmonic=2, TString sDetForEP="TPC", TString sDetWgtsFile = "/home/shi/alice/CalibrationFiles/CalibV0GainCorrectionLHC18q_Oct2021.root", TString sZDCCorrFile = "/home/shi/alice/CalibrationFiles/RecenteringResultFinal_2018q_3rdOrderCent_vtx_orbitNo.root", Bool_t bUseMCEfficiencyWeight = kTRUE, TString sMCEfficiencyWeightForTrackFile = "/home/shi/alice/CalibrationFiles/efficiencyWeightAllCharged18qr.root", const char *suffix = "LHC18qAnaTPCEPFB768")
{

  // ====================================================================
  printf("===================================================================================\n");
  printf("           Initialising Task: AliAnalysisTaskGammaDeltaPIDSaveQvec                 \n");
  printf("===================================================================================\n");

  TGrid::Connect("alien://");                                //@Shi
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName; // common outfile filename
  list1OutName += ":Results";         // This directory contains result histograms


  TString TaskName;
  TaskName.Form("gTaskGammaDeltaPID%d_%d_%s", gFilterBit, gNclustTPC, suffix);

  AliAnalysisTaskGammaNonIsotropicCorrUsingWeights *taskGammaPID = new AliAnalysisTaskGammaNonIsotropicCorrUsingWeights(TaskName);

  ///-------> Analysis Object Created, now pass the arguments
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    taskGammaPID->SelectCollisionCandidates(AliVEvent::kMB);
    printf("\n =========> AddTaskCMW::Info() Trigger = kMB  \n");
  }
  else if(sTrigger=="kSemiCentral" || sTrigger=="SemiCentral" || sTrigger=="semicentral"){
    taskGammaPID->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    printf("\n =========> AddTaskCMW::Info() Trigger = kSemiCentral \n");
  }
  else if(sTrigger=="kCentral" || sTrigger=="Central" || sTrigger=="central"){
    taskGammaPID->SelectCollisionCandidates(AliVEvent::kCentral);
    printf("\n =========> AddTaskCMW::Info() Trigger = kCentral \n");
  }
  else if(sTrigger=="kAny" || sTrigger=="kAll"){
    taskGammaPID->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kSemiCentral | AliVEvent::kCentral);
  }
  else if(sTrigger=="kINT7andSemiCentral"){
    taskGammaPID->SelectCollisionCandidates(AliVEvent::kINT7 | AliVEvent::kSemiCentral); 
    printf("\n =========> AddTaskCMW::Info() Trigger = kINT7+kSemiCentral \n");
  }
  else{//if trigger==kINT7 or no trigger provided:
    taskGammaPID->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskCMW::Info() Trigger = kINT7 \n");
  }
  
  

  ///Set Event cuts:
  taskGammaPID->SetVzRangeMin(fVzMin);
  taskGammaPID->SetVzRangeMax(fVzMax);
  taskGammaPID->SetFlagSkipPileUpCuts(bSkipPileUp);  
  
  taskGammaPID->SetwhichData(whichData);
  taskGammaPID->Setperiod(period);

  cout<<"=========> AddTaskCMW::Info() setting Event Plane Det: "<<sDetForEP<<endl;
  taskGammaPID->SetDetectorforEventPlane(sDetForEP);
  
  if(fdcaxy > 0) {
    taskGammaPID->SetDCAXYRangeMax(fdcaxy);
  }
  if(fdcaz > 0) {
    taskGammaPID->SetDCAZRangeMax(fdcaz);
  }
  
  if(sCentEstimator=="V0" || sCentEstimator=="V0M"){ 
    taskGammaPID->SetCentralityEstimator("V0M");    
  }
  else{
    taskGammaPID->SetCentralityEstimator(sCentEstimator);  //  use the Estimator provided in AddTask.
  }



  
  //Set Track cuts:

  taskGammaPID->SetPtRangeMin(fPtMin);
  taskGammaPID->SetPtRangeMax(fPtMax);
  taskGammaPID->SetEtaRangeMin(fEtaMin);
  taskGammaPID->SetEtaRangeMax(fEtaMax);
  taskGammaPID->SetTrackCutChi2Min(0.1);
  taskGammaPID->SetTrackCutdEdxMin(10.0);  
  taskGammaPID->SetFilterBit(gFilterBit);
  taskGammaPID->SetNSigmaCutTPC(nSigTPC);    /// For PID only.Does not apply to Inclusive Charged Tracks
  taskGammaPID->SetNSigmaCutTOF(nSigTOF);
  taskGammaPID->SetParticlePID(fparticle);
  taskGammaPID->SetTrackCutChi2Max(fChi2max);
  taskGammaPID->SetFlagUseKinkTracks(kFALSE);
  taskGammaPID->SetCumulantHarmonic(vnHarmonic);
  taskGammaPID->SetTrackCutNclusterMin(gNclustTPC);  
  taskGammaPID->SetTrackCutUseTPCCrossedRows(fUseTPCCrossedRows);
  taskGammaPID->SetTrackCutTPCsharedMin(fTPCsharedCut);


   ///  -----> Separate AddTask Added For Lambda-X correlation
   ///  AddTaskGammaDeltaPID.C  
 




  //========================= Setup Correction Files ======================> 
  //-----------------------------------------------------------------------------
 
  TFile* fV0ZDCWgtsFile = TFile::Open(sDetWgtsFile,"READ");
  TList* fListDetWgts=NULL;

  if(fV0ZDCWgtsFile) {    
    fListDetWgts = dynamic_cast <TList*> (fV0ZDCWgtsFile->FindObjectAny("fWgtsV0ZDC"));
    std::cout<<" \n ==============> TList found for V0/ZDC wgts.. GOOD! ";
    // fListDetWgts->ls(); 

    if(fListDetWgts) {
      taskGammaPID->SetListForV0MCorr(fListDetWgts);
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
	
	if(fListZDCCorr) {
	  taskGammaPID->SetListForZDCCorr(fListZDCCorr);
    }
    else{
	  printf("\n\n *** AddTask::WARNING => ZDC Recentering file Exist, But TList Not Found!!");
      printf("\n May be wrong TList name? No Correction for ZDC recentering !! \n\n");
	}
  
  } else{
	printf("\n\n *** AddTask::WARNING => NO File Found for ZDC recentering!!\n AddTask::Info() ===> No ZDC Correction!! \n\n");
  }
  //-----------------------------------------------------------------------------
  
  TFile* fMCEfficiencyWeightForTrackFile = TFile::Open(sMCEfficiencyWeightForTrackFile, "READ");
  TList* fListTPCMCEfficiencyWeightForTrack=NULL;
  
  if(bUseMCEfficiencyWeight) {
    if(fMCEfficiencyWeightForTrackFile) {
	  fListTPCMCEfficiencyWeightForTrack = dynamic_cast <TList*> (fMCEfficiencyWeightForTrackFile->FindObjectAny("fOutputRecenter"));
	
	  if(fListTPCMCEfficiencyWeightForTrack) {
	    taskGammaPID->SetListForMCEfficiencyWeightForTrack(fListTPCMCEfficiencyWeightForTrack);
      }
      else{
	    printf("\n\n *** AddTask::WARNING => MC Efficiency Weight For Track file Exist, But TList Not Found!!");
        printf("\n May be wrong TList name? No Correction for MC Efficiency Weight For Track !! \n\n");
	  }
  
    } else{
	  printf("\n\n *** AddTask::WARNING => NO File Found for MC Efficiency Weight For Track file!!\n AddTask::Info() ===> No ZDC Correction!! \n\n");
    }
  }
  
  //=================================================================================


  

  ///---> Now Pass data and containers to Analysis Object ----
 
  mgr->AddTask(taskGammaPID);                        // connect the task to the analysis manager
  mgr->ConnectInput(taskGammaPID, 0, cinput);        // give AOD event to my Task..!!


  AliAnalysisDataContainer  *cOutPut1;
  TString                  sMyOutName;
  sMyOutName.Form("SimpleTask_%s",suffix);
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  mgr->ConnectOutput(taskGammaPID, 1, cOutPut1);
  
 
  printf("\n\n ================> AddTask was Configured properly... <==================\n\n");

  //return taskGammaPID;



}//Task Ends
