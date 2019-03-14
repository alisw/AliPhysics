
class AliAnalysisTaskCMEMC;

void AddTaskCMEMC(Int_t gFilterBit = 96, Double_t fPtMin=0.2, Double_t fPtMax=5.0, Double_t fEtaMin=-0.8, Double_t fEtaMax=0.8, Double_t fDCAxyMax=2.4,Double_t fDCAzMax=3.2, Int_t gNclustTPC=70, TString sCentEstimator="V0M", Double_t fCentralityMin=0., Double_t fCentralityMax=90., Float_t fVzMin = -10.0, Float_t fVzMax = 10.0, TString sNuclei="PbPb", TString sTrigger="kINT7", Bool_t bSkipPileUp=kFALSE, Double_t fSlope=3.45, Float_t fConst=100, Int_t gN = 1, Int_t gM = 2, Int_t gPsiN=2, Bool_t bUseMC=kFALSE, TString sMCfilePath = "alien:///alice/cern.ch/user/m/mhaque/gain/FB96_Hijing_LHC15o_HI_CorSec.root", Bool_t bUseNUA=kFALSE, TString sNUAFilePath = "alien:///alice/cern.ch/user/m/mhaque/gain/NUA15o_pass1_FB96_C15k_CentBin5_AvgEtaFull.root", Bool_t bV0MCorr=kFALSE, TString sV0MFile="alien:///alice/cern.ch/user/m/mhaque/gain/V0GainEq_LHC15o_pass1HI_C15K_RbyR.root", Bool_t bFillNUAPID=kTRUE, const char *suffix = "")
{
  // standard with task
  printf("========================================================================================\n");
  printf("               Info: Initialising AliAnalysisTaskCMEMC \n"                               );
  printf("========================================================================================\n");
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  TString file = AliAnalysisManager::GetCommonFileName();

  TString ContainerFE;
  ContainerFE.Form("FlowEventCont_%s", suffix);

  AliAnalysisDataContainer    *cinput = mgr->GetCommonInputContainer();  //AOD event

  TString listOutName1 = file;      // file is the common outfile filename
  listOutName1 += ":Results";

  AliAnalysisDataContainer *coutputCont1;

  Int_t gCentMin = fCentralityMin;
  Int_t gCentMax = fCentralityMax;

  TString TaskCMEV0PID;
  TString fCME_container;

  TaskCMEV0PID.Form("TaskCMEV0PID_Cent_%d_%d_%s", gCentMin, gCentMax, suffix);

  //cout<<"Add taskname = "<<TaskCMEV0PID.Data()<<endl;

  AliAnalysisTaskCMEMC *task_CME  = new AliAnalysisTaskCMEMC(TaskCMEV0PID);

  //task_CME->SelectCollisionCandidates(AliVEvent::kINT7); //default if kINT7
  // if(sTrigger=="kMB" || sTrigger=="kmb" || sTrigger=="MB"){
  //   task_CME->SelectCollisionCandidates(AliVEvent::kMB);
  // }

  //Track cuts:
  task_CME->SetFilterBit(gFilterBit);
  task_CME->SetNSigmaCutTPC(2.0);
  task_CME->SetPtRangeMin(fPtMin);
  task_CME->SetPtRangeMax(fPtMax);
  task_CME->SetEtaRangeMin(fEtaMin);
  task_CME->SetEtaRangeMax(fEtaMax);

  task_CME->SetTrackCutdEdxMin(10.0);
  task_CME->SetTrackCutDCAzMax(fDCAxyMax);
  task_CME->SetTrackCutDCAxyMax(fDCAzMax);
  task_CME->SetTrackCutChi2Min(0.1);
  task_CME->SetTrackCutNclusterMin(gNclustTPC);
  task_CME->SetFlagUseKinkTracks(kFALSE);

 

  //Event cuts:
  task_CME->SetCentralityPercentileMin(fCentralityMin);
  task_CME->SetCentralityPercentileMax(fCentralityMax);
  task_CME->SetPileUpCutParam(fSlope,fConst);
  task_CME->SetCollisionSystem(sNuclei);
  task_CME->SetHarmonicsFor3Particle(gN,gM);  // n and m, Cos(n*phi1 + m*phi2 - (n+m)Psi_{gPsiN})
  task_CME->SetEventPlaneHarmonic(gPsiN);     // gPsiN = order N of Event Plane. Look Here ---->^
  task_CME->SetFlagSkipPileUpCuts(bSkipPileUp);
  task_CME->SetVzRangeMin(fVzMin);
  task_CME->SetVzRangeMax(fVzMax);

  //if(sCentEstimator=="V0" || sCentEstimator=="V0M"){
    //task_CME->SetCentralityEstimator("V0M");    //V0M or V0A or V0C or V0=V0M
  //}
  //else{
    //task_CME->SetCentralityEstimator(sCentEstimator);    //V0M or V0A or V0C or V0=V0M
  //}



  if(bFillNUAPID){
    task_CME->SetFlagFillNUAforPID(kTRUE);
  }

  TFile* fMCFile=NULL;
  TList* fListMC=NULL;

  if(bUseMC) {

    task_CME->SetFlagForMCcorrection(kTRUE);
    fMCFile = TFile::Open(sMCfilePath,"READ");

    if(!fMCFile) {
      Ssiz_t pos = sMCfilePath.Index(".root");
      sMCfilePath.Replace(pos,5,"_copy.root");
      fMCFile = TFile::Open(sMCfilePath,"READ");
      if(!fMCFile) printf("\n\n *** ERROR: MC correction file not found! \n  Please check path/name \n\n");
      //exit(1);
    }

    else { 
      fListMC = dynamic_cast <TList*> (fMCFile->FindObjectAny("fMcEffiHij"));
     
      if(fListMC) {
        //std::cout<<" \n ==============> List found for MC corr, here is all the histograms : ";
        //fListMC->ls();
	task_CME->SetFBEfficiencyList(fListMC);
      }
      else{
	printf("\n\n *** ERROR: MC file Exist, But fList name is wrong!!\n please check name \n\n");
      }
    }
  }
  else{
    task_CME->SetFlagForMCcorrection(kFALSE);
  }
  //--------------------------

  TFile* fNUAFile=NULL;
  TList* fListNUA=NULL;

  if(bUseNUA){      
    fNUAFile = TFile::Open(sNUAFilePath,"READ");

    if(!fNUAFile) {
      Ssiz_t pos = sNUAFilePath.Index(".root");
      sNUAFilePath.Replace(pos,5,"_copy.root");
      fNUAFile = TFile::Open(sNUAFilePath,"READ");
      if(!fNUAFile) printf("\n\n *** ERROR: NUA wgt file not found! \n  Please check path/name \n\n");
      //exit(1);
    }
    else { 
      fListNUA = dynamic_cast <TList*> (fNUAFile->FindObjectAny("fNUA_ChPosChNeg"));
      //std::cout<<" \n ==============> List found for NUA, here is all the histograms : ";
      //fListNUA->ls();

      if(fListNUA) {
	task_CME->SetListForNUACorr(fListNUA);
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
      Ssiz_t pos = sV0MFile.Index(".root");
      sV0MFile.Replace(pos,5,"_copy.root");
      fV0MFile = TFile::Open(sV0MFile,"READ");
      if(!fV0MFile) printf("\n\n !!! ERROR: VOM Gain correction file not found! Please check path \n ***(EXIT now)*** \n\n");
      //exit(1);
    } 
    else{
      TList* fListV0MUse = dynamic_cast<TList*>(fV0MFile->FindObjectAny("fV0MChWgts"));
      if(fListV0MUse) {
	task_CME->SetFlagV0MGainCorr(kTRUE);
	task_CME->SetListForV0MCorr(fListV0MUse);
      }
      else{
	printf("\n\n !!!**** ERROR: TList for VOM Gain Correction found \n Please check name !!!\n\n");
	task_CME->SetListForV0MCorr(NULL);
      }
    }
  }


   
  mgr->AddTask(task_CME);                        // connect the task to the analysis manager
  mgr->ConnectInput(task_CME, 0, cinput);        // give AOD event to my Task..!!

  fCME_container.Form("CME_PID_%s", suffix);
  coutputCont1 = (AliAnalysisDataContainer *) mgr->CreateContainer(fCME_container,TList::Class(),AliAnalysisManager::kOutputContainer,listOutName1.Data());
  mgr->ConnectOutput(task_CME, 1, coutputCont1);
  
 
  printf("\n ===================> AddTaskCMEV0PID() Configured properly Panos <=====================\n");

//return task_CME;



}//main ends
