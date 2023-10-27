
AliAnalysisTaskCorrPbPbMC *AddTaskCorrPbPbMC(Int_t fCentralityMin=0, Int_t fCentralityMax=90, Double_t fVzMax=10, Int_t fpileupcut = 1, Int_t fCentEstFlag = 0, Int_t fFilterBit=96, Double_t fchi2tpc=2.5, Double_t fchi2its=36, Double_t fnotpccrossedrows=70, Double_t fpidnSigma_Pion=2.0, Double_t fpidnSigma_Kaon=2.0, Double_t fpidnSigma_Proton=2.0, TString OutFileName = "_default", TString sMCfilePath ="alien:///alice/cern.ch/user/s/swati/EfficiencyPbPbLHC20j6a/CentralitywiseEff/EfficiencyHijingPbPb.root")

{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCMWPU                             \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":ResultsMC";         // This directory contains result histograms

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TGrid::Connect("alien://");

  TString TaskMeanpt;
  TaskMeanpt.Form("gTaskMeanpt%d_%d_%s", gCentMin, gCentMax, " ");
  //gROOT->LoadMacro("AliAnalysisTaskResonanceVsMultiplicity.cxx++g");                                                           
  //gInterpreter->ProcessLine(".x AliAnalysisTaskCorrPbPbMC.cxx++g");
  AliAnalysisTaskCorrPbPbMC *task_Mpt = new AliAnalysisTaskCorrPbPbMC(TaskMeanpt);
  task_Mpt->SelectCollisionCandidates(AliVEvent::kINT7); //trigger for analysis                                                               




  /*
  ///-------> Analysis Object Created, now pass the arguments
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    task_Mpt->SelectCollisionCandidates(AliVEvent::kMB);
    printf("\n =========> AddTaskResonancevsMultiplicity::Info() Trigger = kMB  \n");
  }
  else if(sTrigger=="kSemiCentral" || sTrigger=="SemiCentral" || sTrigger=="semicentral"){
    task_Mpt->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kSemiCentral \n");
  }
  else if(sTrigger=="kCentral" || sTrigger=="Central" || sTrigger=="central"){
    task_Mpt->SelectCollisionCandidates(AliVEvent::kCentral);
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kCentral \n");
  }
  else{//if trigger==kINT7 or no trigger provided:
    task_Mpt->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kINT7 \n");
  }
  */
  // ///swati: add event and track cuts
  // task_Mpt->SetCentralityPercentileMin(fCentralityMin);
  // task_Mpt->SetCentralityPercentileMax(fCentralityMax);
  // task_Mpt->SetVzRangeMin(fVzMin);

  
  // ///Event cuts:
  task_Mpt->SetVzRangeMax(fVzMax);
  task_Mpt->SetPileupCutValue(fpileupcut); //default: 1; can vary 2, 3 or 4
  task_Mpt->SetCentralityEstimator(fCentEstFlag); // 0 for V0M, 1 for CL0, 2 for CL1 and 3 for CL2

  // //Track cuts:
  task_Mpt->SetTrackFilterBit(fFilterBit);
  task_Mpt->SetMaxChi2PerTPCClusterRange(fchi2tpc);
  task_Mpt->SetMaxChi2PerITSClusterRange(fchi2its);
  task_Mpt->SetPIDnSigmaCut( fpidnSigma_Pion, fpidnSigma_Kaon, fpidnSigma_Proton);
  task_Mpt->SetMinNoTPCCrossedRows(fnotpccrossedrows);
  
  
  /*
  TString OutTreeName;
  OutTreeName = "fTreeEvent";
  OutTreeName += OutFileName;
  task_Mpt->SetTreeName(OutTreeName);
  */

  TFile *fMCFile = TFile::Open(sMCfilePath,"READ");
  TList *fListMC = NULL;

  if(fMCFile)
    {
      fListMC = dynamic_cast <TList*> (fMCFile->FindObjectAny("fMCEffHijing"));

      if(fListMC)
	{
	  task_Mpt->SetListForTrkCorr(fListMC);
	}
      else
	{
	  printf("\n\n *** AddTask::Warning \n => MC file esists, but TList not found !!! \n AddTask::Info() ===> No MC Correction !! \n\n");
	}
    }
  else
    {
      printf("\n\n *** AddTask::WARNING \n => no MC file!!! \n AddTask::Info() ===> NO MC Correction!! \n\n");
    }

    


  mgr->AddTask(task_Mpt);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_Mpt, 0, cinput);        


  AliAnalysisDataContainer  *cOutPut1;
  AliAnalysisDataContainer  *cOutPut2;
  AliAnalysisDataContainer  *cOutPut3;
  TString                  sMyOutName1;
  TString                  sMyOutName2;
  TString                  sMyOutName3;
  sMyOutName1 += "Analysis_tree";
  sMyOutName1 += OutFileName;
  sMyOutName2 += "Histogram_TrackVariables";
  sMyOutName2 += OutFileName;
  sMyOutName3 += "QAPileupPlots";
  sMyOutName3 += OutFileName;
  
  
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName2,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  cOutPut2 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName3,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  cOutPut3 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName1,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  
  mgr->ConnectOutput(task_Mpt, 1, cOutPut1);
  mgr->ConnectOutput(task_Mpt, 2, cOutPut2);
  mgr->ConnectOutput(task_Mpt, 3, cOutPut3);
  

  return task_Mpt;

}//Task Ends

