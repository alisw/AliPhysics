AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2 *AddTaskDiffPtFluc_V0particles_pTmax5_v2(Int_t lCentralityMin=0,
  Int_t lCentralityMax=90,
  Double_t lVzMax=10,
  Int_t lpileupcut = 1,
  Int_t lCentEstFlag = 0,
  Int_t lFilterBit=96,
  Double_t lchi2tpc=2.5,
  Double_t lchi2its=36,
  Double_t lpidnSigma=2.0,
  Double_t ltpccrossedrows = 70,
  TString OutFileName = "_default",
  TString sMCfilePath = "file:///media/swati/data/Ollitraut_V0/implentation_in_process/latest_data/EfficiencyHijingPbPb.root",
  Double_t letaleftcut = 0.0,
  Double_t letacut=0.4,
  Int_t leffflag=0,
  Int_t leffCorrFlag=0,
  Int_t lElRejectFlag=0,
  Int_t lPIDExclusiveFlag=0,
  Int_t lFillHistTrkQAFlag=1,
  Int_t lFillHistTrkPIDQAFlag=1,
  Int_t lBayesianPIDFlag=1,
  Int_t lMinV0TracksTpcClustersNo=80,
  Int_t lMinV0TracksTpcCrossedRowsNo=80,
  Int_t lMaxV0TracksTpcSharedClustersNo=5,
  Double_t lMaxV0TracksChi2TPCperClstr=2.5,
  Double_t lMaxV0TracksChi2ITSperClstr=36,
  Double_t lRatioTPCcrossedrowsByFindableclusters=0.8,
  Double_t lDcaV0DaughterTracksToPV=0.1,
  Double_t lLambdaPropLifetime=40,
  Double_t lMinLambdaTransDecayRadius=5,
  Double_t lMaxLambdaTransDecayRadius=100,
  Double_t lLambdaDcaV0daughters=0.2,
  Double_t lLambdaDcaV0toPV=0.5,
  Double_t lLambdaCosPAval=0.997,
  Double_t lK0sPropLifetime=10,
  Double_t lMinK0sTransDecayRadius=5,
  Double_t lMaxK0sTransDecayRadius=100,
  Double_t lK0sDcaV0daughters=0.4,
  Double_t lK0sCosPAval=0.999,
  Double_t lArmentousCutVal=0.2,
  Double_t lLambdaDaughtersPIDcut=3.0,
  Double_t lK0sDaughtersPIDcut=3.0,
  Double_t lLambdaMassCut=0.005,
  Double_t lK0sMassCut=0.015,
  Double_t lPIDbayesPion=0.95,
  Double_t lPIDbayesKaon=0.9,
  Double_t lPIDbayesProton=0.9
  )
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

  Int_t gCentMin = lCentralityMin;
  Int_t gCentMax = lCentralityMax;

  TGrid::Connect("alien://");

  TString TaskMeanpt;
  TaskMeanpt.Form("gTaskMeanpt%d_%d_%s", gCentMin, gCentMax, " ");
  //gROOT->LoadMacro("AliAnalysisTaskResonanceVsMultiplicity.cxx++g");                                                           
  //gInterpreter->ProcessLine(".x AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2.cxx++g");
  AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2 *task_v0pT = new AliAnalysisTaskDiffPtFluc_V0particles_pTmax5_v2(TaskMeanpt);
  task_v0pT->SelectCollisionCandidates(AliVEvent::kINT7); //trigger for analysis                                                               




  /*
  ///-------> Analysis Object Created, now pass the arguments
  if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){   // if We want MB Trigger
    task_v0pT->SelectCollisionCandidates(AliVEvent::kMB);
    printf("\n =========> AddTaskResonancevsMultiplicity::Info() Trigger = kMB  \n");
  }
  else if(sTrigger=="kSemiCentral" || sTrigger=="SemiCentral" || sTrigger=="semicentral"){
    task_v0pT->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kSemiCentral \n");
  }
  else if(sTrigger=="kCentral" || sTrigger=="Central" || sTrigger=="central"){
    task_v0pT->SelectCollisionCandidates(AliVEvent::kCentral);
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kCentral \n");
  }
  else{//if trigger==kINT7 or no trigger provided:
    task_v0pT->SelectCollisionCandidates(AliVEvent::kINT7);      // default is kINT7
    printf("\n =========> AddTaskRESONANCEVSMULTIPLICTY::Info() Trigger = kINT7 \n");
  }
  */
  // ///swati: add event and track cuts
  // task_v0pT->SetCentralityPercentileMin(fCentralityMin);
  // task_v0pT->SetCentralityPercentileMax(fCentralityMax);
  // task_v0pT->SetVzRangeMin(fVzMin);
  
  // ///Event cuts:
  
  task_v0pT->SetVzRangeMax(lVzMax);
  task_v0pT->SetPileupCutValue(lpileupcut); //default: 1; can vary 2, 3 or 4
  task_v0pT->SetCentralityEstimator(lCentEstFlag); // 0 for V0M, 1 for CL0, 2 for CL1 and 3 for CL2
       
  // //Track cuts:
  
  task_v0pT->SetTrackFilterBit(lFilterBit);
  task_v0pT->SetMaxChi2PerTPCClusterRange(lchi2tpc);
  task_v0pT->SetMaxChi2PerITSClusterRange(lchi2its);
  task_v0pT->SetPIDnSigmaCut(lpidnSigma);
  task_v0pT->SetEtaLeftCut(letaleftcut);
  task_v0pT->SetEtaCut(letacut);
  task_v0pT->SetMinNCrossedRowsTPC(ltpccrossedrows);
  task_v0pT->SetEfficiencyEffectImposeFlag(leffflag);
  task_v0pT->SetEfficiencyCorrectionFlag(leffCorrFlag);
  task_v0pT->SetRejectElectronFlag(lElRejectFlag);
  task_v0pT->SetExclusivePIDCutFlag(lPIDExclusiveFlag);
  task_v0pT->SetFillTrackQAHistogramsFlag(lFillHistTrkQAFlag);
  task_v0pT->SetFillTrackPIDQAHistogramsFlag(lFillHistTrkPIDQAFlag);
  task_v0pT->SetSelectPiKaPrByBayesianPIDFlag(lBayesianPIDFlag);
  task_v0pT->SetMinV0TracksTpcClustersNo(lMinV0TracksTpcClustersNo);
  task_v0pT->SetMinV0TracksTpcCrossedRowsNo(lMinV0TracksTpcCrossedRowsNo);
  task_v0pT->SetMaxV0TracksTpcSharedClustersNo(lMaxV0TracksTpcSharedClustersNo);
  task_v0pT->SetMaxV0TracksChi2TPCperClstr(lMaxV0TracksChi2TPCperClstr);
  task_v0pT->SetMaxV0TracksChi2ITSperClstr(lMaxV0TracksChi2ITSperClstr);
  task_v0pT->SetRatio_TPCcrossedrowsByFindableclusters(lRatioTPCcrossedrowsByFindableclusters);
  task_v0pT->SetDcaV0DaughterTracksToPV(lDcaV0DaughterTracksToPV);
  task_v0pT->SetLambdaPropLifetime(lLambdaPropLifetime);
  task_v0pT->SetLambdaTransDecayRadius_Min(lMinLambdaTransDecayRadius);
  task_v0pT->SetLambdaTransDecayRadius_Max(lMaxLambdaTransDecayRadius);
  task_v0pT->SetLambdaDCAv0daughters(lLambdaDcaV0daughters);
  task_v0pT->SetLambdaDcaV0toPV(lLambdaDcaV0toPV);
  task_v0pT->SetLambdaCosPA(lLambdaCosPAval);
  task_v0pT->SetK0sPropLifetime(lK0sPropLifetime);
  task_v0pT->SetK0sTransDecayRadius_Min(lMinK0sTransDecayRadius);
  task_v0pT->SetK0sTransDecayRadius_Max(lMaxK0sTransDecayRadius);
  task_v0pT->SetK0sDCAv0daughters(lK0sDcaV0daughters);
  task_v0pT->SetK0sCosPA(lK0sCosPAval);
  task_v0pT->SetK0sArmentousCut(lArmentousCutVal);
  task_v0pT->SetLambdaDaughtersPIDcut(lLambdaDaughtersPIDcut);
  task_v0pT->SetK0sDaughtersPIDcut(lK0sDaughtersPIDcut);
  task_v0pT->SetLambdaFinalMassCut(lLambdaMassCut);
  task_v0pT->SetK0sFinalMassCut(lK0sMassCut);
  task_v0pT->SetBayesPIDPionVal(lPIDbayesPion);
  task_v0pT->SetBayesPIDKaonVal(lPIDbayesKaon);
  task_v0pT->SetBayesPIDProtonVal(lPIDbayesProton);
  
  TString OutTreeName;
  OutTreeName = "fTreeEvent";
  OutTreeName += OutFileName;
  task_v0pT->SetTreeName(OutTreeName);
  

  TFile *fMCFile = TFile::Open(sMCfilePath,"READ");
  TList *fListMC = NULL;

  if(fMCFile)
    {
      fListMC = dynamic_cast <TList*> (fMCFile->FindObjectAny("fMCEffHijing"));

      if(fListMC)
	{
	  task_v0pT->SetListForTrkCorr(fListMC);
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

    

  
  mgr->AddTask(task_v0pT);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_v0pT, 0, cinput);        


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
  
  mgr->ConnectOutput(task_v0pT, 1, cOutPut1);
  mgr->ConnectOutput(task_v0pT, 2, cOutPut2);
  mgr->ConnectOutput(task_v0pT, 3, cOutPut3);
  

  return task_v0pT;

}//Task Ends

