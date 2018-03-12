#include<TList.h>
#include<TGrid.h>
#include "TSystem.h"

class AliAnalysisTaskCMEV0PID;

void AddTaskCMEV0PID(Int_t gFilterBit = 96, Float_t fPtMin=0.2, Float_t fPtMax=5.0, Float_t fEtaMin=-0.8, Float_t fEtaMax=0.8, Float_t fCentralityMin=0.,Float_t fCentralityMax=90.,TString sNuclei="PbPb", TString sTrigger="kINT7", Bool_t bSkipPileUp=kFALSE, Float_t fSlope=3.45, Float_t fConst=100, Int_t gPsiN=2, Bool_t bUseMC=kFALSE, TString sMCfilePath = "alien:///alice/cern.ch/user/m/mhaque/gain/FB96_Hijing_LHC15o_HI_CorSec.root", Bool_t bUseNUA=kFALSE, TString sNUAFilePath = "alien:///alice/cern.ch/user/m/mhaque/gain/NUA15o_pass1_FB96_C15k_CentBin5_AvgEtaFull.root", Bool_t bV0MCorr=kFALSE, TString sV0MFile="alien:///alice/cern.ch/user/m/mhaque/gain/V0GainEq_LHC15o_pass1HI_C15K_RbyR.root", Bool_t bFillNUAPID=kTRUE, const char *suffix = "")
{
  // standard with task
  printf("========================================================================================\n");
  printf("               PID: Initialising AliAnalysisTaskCMEV0PID \n");
  printf("========================================================================================\n");
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString file = AliAnalysisManager::GetCommonFileName();

  TString ContainerFE;
  ContainerFE.Form("FlowEventCont_%s", suffix);

  AliAnalysisDataContainer    *cinput = mgr->GetCommonInputContainer();  //AOD event



  TString TaskCMEV0PID;
  TString fCME_container;

  TString listOutName1 = file;      // file is the common outfile filename
  listOutName1 += ":Results";


  const Int_t NumberOfCuts = 1;  //max 5

  Double_t nSigmaCuts[5]  = {2.0, 2.5, 3.0};


  AliAnalysisTaskCMEV0PID *task_CME[5];
  AliAnalysisDataContainer *coutputCont1[5];

  for(int i=0; i<5; i++) {
    task_CME[i] = NULL;
    coutputCont1[i] = NULL;
  }

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  for(int i=0;i<NumberOfCuts;i++){
    TaskCMEV0PID.Form("TaskCMEV0PID_Cent_%d_%d_%s_%2.1f", gCentMin, gCentMax, suffix, nSigmaCuts[i]);
    //cout<<"Add taskname = "<<TaskCMEV0PID.Data()<<endl;
    task_CME[i] = new AliAnalysisTaskCMEV0PID(TaskCMEV0PID);

    task_CME[i]->SelectCollisionCandidates(AliVEvent::kINT7); //default if kINT7

    if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){
      task_CME[i]->SelectCollisionCandidates(AliVEvent::kMB);
    }

    task_CME[i]->SetFilterBit(gFilterBit);
    task_CME[i]->SetNSigmaCutTPC(nSigmaCuts[i]);
    task_CME[i]->SetPtRangeMin(fPtMin);
    task_CME[i]->SetPtRangeMax(fPtMax);
    task_CME[i]->SetEtaRangeMin(fEtaMin);
    task_CME[i]->SetEtaRangeMax(fEtaMax);
    task_CME[i]->SetCentralityPercentileMin(fCentralityMin);
    task_CME[i]->SetCentralityPercentileMax(fCentralityMax);
    task_CME[i]->SetPileUpCutParam(fSlope,fConst);
    task_CME[i]->SetCollisionSystem(sNuclei);
    task_CME[i]->SetEventPlaneHarmonic(gPsiN);
    task_CME[i]->SetFlagSkipPileUpCuts(bSkipPileUp);

    if(bFillNUAPID){
      task_CME[i]->SetFlagFillNUAforPID(kTRUE);
    }


    if(bUseMC) {
      task_CME[i]->SetFlagForMCcorrection(kTRUE);
      task_CME[i]->SetFBEfficiencyFilePath(sMCfilePath);
    }
    else{
      task_CME[i]->SetFlagForMCcorrection(kFALSE);
    }
    //--------------------------
    if(bUseNUA){      
      TFile* fNUAFile = TFile::Open(sNUAFilePath,"READ");
      if(!fNUAFile) {
	printf("\n\n *** ERROR: NUA wgt file not found! \n  Please check name \n ***(EXIT now)*** \n\n");
	exit(1);
      }
      else { 
	TList* fListNUA = dynamic_cast<TList*>(fNUAFile->FindObjectAny("fNUA_ChPosChNeg"));
	if(fListNUA) {
	  task_CME[i]->SetListForNUACorr(fListNUA);
	}
	else{
	  printf("\n\n *** ERROR: NUA file Exist, But fList name is wrong!!\n please check name \n\n");
	}
      }
    }
    //----------------------------
    if(bV0MCorr){
      TFile* fV0MFile = TFile::Open(sV0MFile,"READ");
      if(!fV0MFile) {
	printf("\n\n !!! ERROR: VOM Gain correction file not found! Please check path \n ***(EXIT now)*** \n\n");
	exit(1);
      } 
      else{
	TList* fListV0MUse = dynamic_cast<TList*>(fV0MFile->FindObjectAny("fV0MChWgts"));
	if(fListV0MUse) {
          task_CME[i]->SetFlagV0MGainCorr(kTRUE);
          task_CME[i]->SetListForV0MCorr(fListV0MUse);
	}
        else{
	  printf("\n\n !!!**** ERROR: TList for VOM Gain Correction found \n Please check name !!!\n\n");
          task_CME[i]->SetListForV0MCorr(NULL);
        }
      }
    }


    
    






    
    mgr->AddTask(task_CME[i]);                        // connect the task to the analysis manager
    mgr->ConnectInput(task_CME[i], 0, cinput);        // give AOD event to my Task..!!

    fCME_container.Form("CME_PID_%s_%2.1f", suffix, nSigmaCuts[i]);
    coutputCont1[i] = (AliAnalysisDataContainer *) mgr->CreateContainer(fCME_container,TList::Class(),AliAnalysisManager::kOutputContainer,listOutName1.Data());
    mgr->ConnectOutput(task_CME[i], 1, coutputCont1[i]);
  }
 



  printf("\n ===================> AddTaskCMEV0PID() Configured properly <=====================\n");

}//main ends
