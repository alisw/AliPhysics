// AddTask for AliAnalysisTaskCascadeChCorrelations task

AliAnalysisTaskV0ChCorrelationpp* AddTaskV0ChCorrelationpp(
                              float cenMin, float cenMax,
                              bool effCorr = 0, bool isMC=0){
   AddTaskV0ChCorrelationpp( 
                            cenMin,  cenMax,
                            Form("Cent%d_%d", Int_t(cenMin), Int_t(cenMax)),
                            Form("Cent%d_%d", Int_t(cenMin), Int_t(cenMax)),
                            effCorr, isMC);
}

 void AddTaskV0ChCorrelationpp(
                              float cenMin, float cenMax,
                              TString folderName="myFolder",
                              TString suffixName="mySuffix",
                              bool effCorr = 0, bool isMC=0       
                              )
{
  // Creates a V0-Ch correlations analysis task and adds it to the analysis manager.

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName.ReplaceAll(".root","");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskV0ChCorrelations", "No analysis manager to connect to.");
    return NULL;
  }
     
  // create task
  AliAnalysisTaskV0ChCorrelationpp* task
  = new AliAnalysisTaskV0ChCorrelationpp(Form("V0ChCorrelation_%s",suffixName.Data()),
                                             cenMin,cenMax,effCorr);  

  task->SetAnalysisMC(isMC);
  //------------------------------Mixing part------------------------------
  task->SetMixingTracks(5000);
  task->SetPoolSize(200); 
 //--------------------------------Variable--------------------------------
  task->SetVtxCut(10.);
  task->SetVtxXMin(10e-5);
  task->SetVtxYMin(10e-5);
  task->SetVtxZMin(10e-5);
  task->SetCentMin(0);
  task->SetCentMax(100.);
  //-----------------------------Track-------------------------------------

  task->SetTrackPtMin(1.);
  task->SetTrackPtMax(8.);
  task->SetTrackEta(0.8);
  task->SetFilterBit(1);
  task->SetAssocNcls(70);
  //------------------------------V0--------------------------------------
  //task->SetV0MCPtMin(3);
  task->SetV0PtMin(3.);
  task->SetV0PtMax(15.);
  task->SetV0Eta(0.8);
  task->SetK0sLifeTimeMin(0);
  task->SetK0sLifeTimeMax(20);
  task->SetLambdaLifeTimeMin(0);
  task->SetLambdaLifeTimeMax(25);
  task->SetDCANegtoPrimVertex(0.06);//
  task->SetDCAPostoPrimVertex(0.06);//
  
  
  task->SetDCAV0DaughtersMax(1);//
  task->Setk0sCPA(0.98);//
  task->SetLambdaCPA(0.998);//
  task->SetCosPointingAngleMin(0.975);
  task->Set2DFiducialMin(0.5);

  task->SetV0DaughterTrackTPCCluster(70.);

  task->SetNCrossedRowsTPCfindable(0.8);

 
  
  task->SetPtArmV0AlphaV0(0.2);
  
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
    = mgr->CreateContainer(Form("Output%s",suffixName.Data()),AliDirList::Class(),AliAnalysisManager::kOutputContainer, 
			                     Form("%s.root", fileName.Data()));
//

  AliAnalysisDataContainer *coutput2
    = mgr->CreateContainer(Form("Output2%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s.root", fileName.Data()));//


AliAnalysisDataContainer *coutput3
    = mgr->CreateContainer(Form("Output3%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s.root", fileName.Data()));//

AliAnalysisDataContainer *coutput4
    = mgr->CreateContainer(Form("Output4%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s.root", fileName.Data()));//



AliAnalysisDataContainer *coutput5
    = mgr->CreateContainer(Form("Output5%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s.root", fileName.Data()));//


AliAnalysisDataContainer *coutput6
    = mgr->CreateContainer(Form("Output6%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s.root", fileName.Data()));//


AliAnalysisDataContainer *coutput7
    = mgr->CreateContainer(Form("Output7%s", suffixName.Data()),AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s.root", fileName.Data()));//



  AliAnalysisDataContainer *cinput1 = NULL;

  TList * effList = 0x0;

  if(effCorr){
    cinput1 = mgr->CreateContainer(Form("Eff%s", suffixName.Data()),
                                    TList::Class(),
                                    AliAnalysisManager::kInputContainer);
    TFile * file = TFile::Open("alien:///alice/cern.ch/user/m/manaam/Efficiencypp/V0Efficiency0_100.root");  
    if(!file) {
      printf("ERROR: efficiency file is not available!\n");
      return;
    }

    effList = (TList*)file->Get(Form("fListCent0_100"));
    if(!effList){
      printf("ERROR: no efficiency list fList%s available\n", suffixName.Data());
      return;
    }
  }
      
  // connect input/output
  mgr->ConnectInput(task, 0, cinput);
  if(effCorr)
   mgr->ConnectInput(task, 1, cinput1);
   mgr->ConnectOutput(task, 1, coutput0);

   mgr->ConnectOutput(task, 2, coutput2);
   mgr->ConnectOutput(task, 3, coutput3);
   mgr->ConnectOutput(task, 4, coutput4);
   mgr->ConnectOutput(task, 5, coutput5);
   mgr->ConnectOutput(task, 6, coutput6);
   mgr->ConnectOutput(task, 7, coutput7);

        
   if(effCorr)
   cinput1->SetData(effList);
  
}

