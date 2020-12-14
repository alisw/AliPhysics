// AddTask for AliAnalysisTaskV0ChCorrelationsys task

AliAnalysisTaskV0ChCorrelationsys* AddTaskV0ChCorrelationsys(float cenMin, float cenMax,  bool effCorr = 0, bool isMC=0,TString fileName_extension = "",TString EffFileNameWithPath = ""){
   AddTaskV0ChCorrelationsys( 
                            cenMin,  cenMax,
                            Form("Cent%d_%d", Int_t(cenMin), Int_t(cenMax)),
                            Form("Cent%d_%d", Int_t(cenMin), Int_t(cenMax)),
                            effCorr, isMC);
}

 AliAnalysisTaskV0ChCorrelationsys* AddTaskV0ChCorrelationsys( float cenMin, float cenMax, TString folderName="myFolder", TString suffixName="waaao", bool effCorr = 0, bool isMC=0 )
{
  // Creates a V0-Ch correlations analysis task and adds it to the analysis manager.

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  TString fileName = AliAnalysisManager::GetCommonFileName();
  //fileName.ReplaceAll(".root","");

    fileName += ":AliAnalysisTaskV0ChCorrelationsys";      // create a subfolder in the file
    fileName += fileName_extension.Data();
    //TString fileName_extension = ""
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskV0ChCorrelationsys", "No analysis manager to connect to.");
    return NULL;
  }
//===================================
if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
//=====================================  
  // create task
  AliAnalysisTaskV0ChCorrelationsys* task = new AliAnalysisTaskV0ChCorrelationsys(folderName,  cenMin,cenMax,effCorr);  

  task->SetAnalysisMC(isMC);
  //------------------------------Mixing part------------------------------
  task->SetMixingTracks(500);
    task->SetPoolSize(100); 
 //--------------------------------Variable--------------------------------
  task->SetVtxCut(9.);
  task->SetNumOfVzBins(9);//   void SetNumOfVzBins(Double_t value){fNumOfVzBins = value;}  
   task->SetCentMin(0);
  task->SetCentMax(10.);
  //-----------------------------Track-------------------------------------
 // task->SetTrackMCPtMin(1);

  task->SetTrackPtMin(1.);
  task->SetTrackPtMax(10.);
  task->SetTrackEta(0.8);
  task->SetFilterBit(768);
  task->SetAssocNcls(70);
  //------------------------------V0--------------------------------------
  //task->SetV0MCPtMin(3);
  task->SetV0PtMin(3.);
  task->SetV0PtMax(16.);
  task->SetV0Eta(0.7);
  task->SetK0sLifeTimeMin(0);
  task->SetK0sLifeTimeMax(20);
  task->SetLambdaLifeTimeMin(0);
  task->SetLambdaLifeTimeMax(25);
  task->SetDCAV0DaughtersMax(0.8);//
  task->SetDCAPostoPrimVertexMink0s(0.1);//
  task->SetDCANegtoPrimVertexMink0s(0.1);//
  task->SetDCAPostoPrimVertexMinLamb(0.1);
  task->SetDCANegtoPrimVertexMinLamb(0.25);//
  task->SetDCAPostoPrimVertexMinALamb(0.25);//
  task->SetDCANegtoPrimVertexMinALamb(0.1);

  task->SetCosPointingAngleMin(0.975);
  task->SetLambdaCPA(0.995);//
  task->Setk0sCPA(0.978);//
  task->Set2DFiducialMin(5);
  task->SetV0DaughterTrackTPCCluster(70.);
  task->SetNCrossedRowsTPCfindable(0.8);
  task->SetPtArmV0AlphaV0(0.2);
  task->SetTrackPileUpCut(kTRUE);
  task->SetV0PileUpCut(kTRUE);
   //-------------------------------------PID--------------------------------
  task->SetV0PIDSigma(3);
  //-------------------------------------------------------------------------
  mgr->AddTask(task);
    
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  //TString outputFileName = AliAnalysisManager::GetCommonFileName();
  //  outputFileName = "XiCh.root";

  // create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput0 
    = mgr->CreateContainer(Form("Output%s", suffixName.Data()),  AliDirList::Class(),   AliAnalysisManager::kOutputContainer, 
			                     Form("%s", fileName.Data()));

AliAnalysisDataContainer *coutput2
    = mgr->CreateContainer(Form("Output2%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//

AliAnalysisDataContainer *coutput3
    = mgr->CreateContainer(Form("Output3%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//
AliAnalysisDataContainer *coutput4
    = mgr->CreateContainer(Form("Output4%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//

AliAnalysisDataContainer *coutput5
    = mgr->CreateContainer(Form("Output5%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//

AliAnalysisDataContainer *coutput6
    = mgr->CreateContainer(Form("Output6%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//

AliAnalysisDataContainer *coutput7
    = mgr->CreateContainer(Form("Output7%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//

  AliAnalysisDataContainer *cinput1 = NULL;

  TList *effList = 0x0;
  
  if(effCorr){

    cinput1 = mgr->CreateContainer(Form("Eff%s", suffixName.Data()),
                                    TList::Class(),
                                    AliAnalysisManager::kInputContainer);
   

//call eff  

    TFile * file = TFile::Open(Form("alien:///alice/cern.ch/user/m/manaam/Efficiency/%s.root",EffFileNameWithPath.Data()));



    if(!file) {
      printf("ERROR: efficiency file is not available!\n",EffFileNameWithPath.Data());
      //return NULL;
    }
    
    effList = (TList*)file->Get("fListCent0_10");

    if(!effList){
      printf("ERROR: no efficiency list fList%s available\n", EffFileNameWithPath.Data());
      //return NULL;
    }
  }


  // connect input/output
  mgr->ConnectInput(task, 0, cinput);

  if(effCorr) mgr->ConnectInput(task, 1, cinput1);
   mgr->ConnectOutput(task, 1, coutput0);
   mgr->ConnectOutput(task, 2, coutput2);
   mgr->ConnectOutput(task, 3, coutput3);
   mgr->ConnectOutput(task, 4, coutput4);
   mgr->ConnectOutput(task, 5, coutput5);
   mgr->ConnectOutput(task, 6, coutput6);
   mgr->ConnectOutput(task, 7, coutput7);


  if(effCorr) cinput1->SetData(effList);

return task;

}


