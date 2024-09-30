AliAnalysisTaskNetProtonCumulants_pp *AddTaskNetProtonCumulants_pp(Int_t fCentralityMin=0, Int_t fCentralityMax=90, Double_t fVzMax=10, Int_t fpileupcut = 1, Int_t fCentEstFlag = 0, Int_t fFilterBit=96, Double_t fchi2tpc=4, Double_t fchi2its=36, Double_t fnotpccrossedrows=70, Double_t fpidnSigma_Pion=2.0, Double_t fpidnSigma_Kaon=2.0, Double_t fpidnSigma_Proton=2.0, TString OutFileName = "_default", TString sMCfilePath ="alien:///alice/cern.ch/user/s/swati/Efficiency_pp_netprotcum/EfficiencyPythia_pp_LHC18pass1.root", Double_t fetacut=0.5, Int_t lBayesianPIDFlag=1, Double_t lPIDbayesPion=0.9, Double_t lPIDbayesKaon=0.9, Double_t lPIDbayesProton=0.9, Int_t lRapidityCutFlag = 1)
{
  // standard with task
  printf("===================================================================================\n");
  printf("\n                PID: Initialising AliAnalysisTaskCMWPU                             \n");
  printf("===================================================================================\n");


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  TString     outfileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer  *cinput = mgr->GetCommonInputContainer();  // AOD event


  TString list1OutName = outfileName;        // common outfile filename
  list1OutName        += ":ResultsData";         // This directory contains result histograms

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TGrid::Connect("alien://");

  TString TaskMeanpt;
  TaskMeanpt.Form("gTaskMeanpt%d_%d_%s", gCentMin, gCentMax, " ");
  //gROOT->LoadMacro("AliAnalysisTaskResonanceVsMultiplicity.cxx++g");                                                           
  //gInterpreter->ProcessLine(".x AliAnalysisTaskNetProtonCumulants_pp.cxx++g");
  AliAnalysisTaskNetProtonCumulants_pp *task_netprotonCum = new AliAnalysisTaskNetProtonCumulants_pp(TaskMeanpt);
  task_netprotonCum->SelectCollisionCandidates(AliVEvent::kINT7); //trigger for analysis                                                               




  /*
  ///-------> Analysis Object Created, now pass the arguments
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    task_netprotonCum->SelectCollisionCandidates(AliVEvent::kMB);
    printf("\n =========> AddTaskResonancevsMultiplicity::Info() Trigger = kMB  \n");
  }
  else if(sTrigger=="kSemiCentral" || sTrigger=="SemiCentral" || sTrigger=="semicentral"){
    task_netprotonCum->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kSemiCentral \n");
  }
  else if(sTrigger=="kCentral" || sTrigger=="Central" || sTrigger=="central"){
    task_netprotonCum->SelectCollisionCandidates(AliVEvent::kCentral);
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kCentral \n");
  }
  else{//if trigger==kINT7 or no trigger provided:
    task_netprotonCum->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kINT7 \n");
  }
  */
  // ///swati: add event and track cuts
  // task_netprotonCum->SetCentralityPercentileMin(fCentralityMin);
  // task_netprotonCum->SetCentralityPercentileMax(fCentralityMax);
  // task_netprotonCum->SetVzRangeMin(fVzMin);
  
  // ///Event cuts:
  
  task_netprotonCum->SetVzRangeMax(fVzMax);
  task_netprotonCum->SetPileupCutValue(fpileupcut); //default: 1; can vary 2, 3 or 4
  task_netprotonCum->SetCentralityEstimator(fCentEstFlag); // 0 for V0M, 1 for CL0, 2 for CL1 and 3 for CL2
       
  // //Track cuts:
  task_netprotonCum->SetTrackFilterBit(fFilterBit);
  task_netprotonCum->SetMaxChi2PerTPCClusterRange(fchi2tpc);
  task_netprotonCum->SetMaxChi2PerITSClusterRange(fchi2its);
  task_netprotonCum->SetPIDnSigmaCut(fpidnSigma_Pion, fpidnSigma_Kaon, fpidnSigma_Proton);
  task_netprotonCum->SetMinNoTPCCrossedRows(fnotpccrossedrows);
  task_netprotonCum->SetEtaCut(fetacut);
  task_netprotonCum->SetSelectPiKaPrByBayesianPIDFlag(lBayesianPIDFlag);
  task_netprotonCum->SetBayesPIDPionVal(lPIDbayesPion);
  task_netprotonCum->SetBayesPIDKaonVal(lPIDbayesKaon);
  task_netprotonCum->SetBayesPIDProtonVal(lPIDbayesProton);
  task_netprotonCum->SetRapidityCutFlag(lRapidityCutFlag);
  
  
  
  TString OutTreeName;
  OutTreeName = "fTreeEvent";
  OutTreeName += OutFileName;
  task_netprotonCum->SetTreeName(OutTreeName);
  

  TFile *fMCFile = TFile::Open(sMCfilePath,"READ");
  TList *fListMC = NULL;

  if(fMCFile)
    {
      fListMC = dynamic_cast <TList*> (fMCFile->FindObjectAny("fMCEffHijing"));

      if(fListMC)
	{
	  task_netprotonCum->SetListForTrkCorr(fListMC);
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

    

  
  mgr->AddTask(task_netprotonCum);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_netprotonCum, 0, cinput);        


  AliAnalysisDataContainer  *cOutPut1;
  AliAnalysisDataContainer  *cOutPut2;
  AliAnalysisDataContainer  *cOutPut3;
  TString                  sMyOutName1;
  TString                  sMyOutName2;
  TString                  sMyOutName3;
  sMyOutName1 += "SimpleTask_tree";
  sMyOutName1 += OutFileName;
  sMyOutName2 += "Histogram_TrackVariables";
  sMyOutName2 += OutFileName;
  sMyOutName3 += "QAPileupPlots";
  sMyOutName3 += OutFileName;
  
  
  
  cOutPut1 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName2,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  cOutPut2 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName3,TList::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  cOutPut3 = (AliAnalysisDataContainer *) mgr->CreateContainer(sMyOutName1,TTree::Class(),AliAnalysisManager::kOutputContainer,list1OutName.Data());
  
  mgr->ConnectOutput(task_netprotonCum, 1, cOutPut1);
  mgr->ConnectOutput(task_netprotonCum, 2, cOutPut2);
  mgr->ConnectOutput(task_netprotonCum, 3, cOutPut3);
  

  return task_netprotonCum;

}//Task Ends

