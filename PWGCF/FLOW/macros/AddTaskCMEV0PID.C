#include<TList.h>
#include<TGrid.h>
#include"TSystem.h"

class AliAnalysisTaskCMEV0PID;

void AddTaskCMEV0PID(Int_t gFilterBit = 96, Float_t fPtMin=0.2, Float_t fPtMax=5.0, Float_t fEtaMin=-0.8, Float_t fEtaMax=0.8, Float_t fCentralityMin=0.,Float_t fCentralityMax=90.,Bool_t bUseMC=kFALSE,TString sMCfilePath = "alien:///alice/cern.ch/user/m/mhaque/gain/FB96_Hijing_LHC15o_HI_CorSec.root", const char *suffix = "")
{

  // standard with task
  printf("========================================================================================\n");
  printf("                      Initializing AliAnalysisTaskCMEV0PID                              \n");
  printf("========================================================================================\n");
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString file = AliAnalysisManager::GetCommonFileName();

  TString ContainerFE;
  ContainerFE.Form("FlowEventCont_%s", suffix);

  AliAnalysisDataContainer    *cinput = mgr->GetCommonInputContainer();  //AOD event

  Double_t centrMin[9] = {0, 5,10,20,30,40,50,60,70};
  Double_t centrMax[9] = {5,10,20,30,40,50,60,70,90};


  TString TaskCMEV0PID;

  TString listOutName1 = file;      // file is the common outfile filename
  listOutName1 += ":Results";

  TString fCME_container;

  const Int_t NumberOfCuts = 1;

  Double_t nSigmaCuts[5]  = {2.0, 2.5, 3.0};

  AliAnalysisTaskCMEV0PID *task_CME[10];
  AliAnalysisDataContainer *coutputCont1[10];

  for(int i=0; i<10; i++) {
    task_CME[i] = NULL;
    coutputCont1[i] = NULL;
  }



  for(int i=0;i<NumberOfCuts;i++){
    TaskCMEV0PID.Form("TaskCMEV0PID_Cent%2.0f_%2.0f_%s_%2.1f",fCentralityMin,fCentralityMax, suffix, nSigmaCuts[i]);
    //cout<<"Add taskname = "<<TaskCMEV0PID.Data()<<endl;
    task_CME[i] = new AliAnalysisTaskCMEV0PID(TaskCMEV0PID);
    task_CME[i]->SelectCollisionCandidates(AliVEvent::kINT7);
    task_CME[i]->SetFilterBit(gFilterBit);
    task_CME[i]->SetNSigmaCutTPC(nSigmaCuts[i]);
    task_CME[i]->SetPtRangeMin(fPtMin);
    task_CME[i]->SetPtRangeMax(fPtMax);
    task_CME[i]->SetEtaRangeMin(fEtaMin);
    task_CME[i]->SetEtaRangeMax(fEtaMax);
    task_CME[i]->SetCentralityPercentileMin(fCentralityMin);
    task_CME[i]->SetCentralityPercentileMax(fCentralityMax);
    if(bUseMC){
      task_CME[i]->SetFlagForMCcorrection(kTRUE);
      task_CME[i]->SetFBEfficiencyFilePath(sMCfilePath);
    }
    else{
      task_CME[i]->SetFlagForMCcorrection(kFALSE);
    }

    mgr->AddTask(task_CME[i]);                        // connect the task to the analysis manager
    mgr->ConnectInput(task_CME[i], 0, cinput);        // give AOD event to my Task..!!

    fCME_container.Form("CME_PID_%s_%2.1f", suffix, nSigmaCuts[i]);
    coutputCont1[i] = (AliAnalysisDataContainer *) mgr->CreateContainer(fCME_container,TList::Class(),AliAnalysisManager::kOutputContainer,listOutName1.Data());
    mgr->ConnectOutput(task_CME[i], 1, coutputCont1[i]);
  }
 
  printf("\n ===================> AddTaskCMEV0PID() Configured properly <=====================\n");

}//main ends
