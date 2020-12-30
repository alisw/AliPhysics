// AddTask for AliAnalysisTaskCascadeChCorrelations task

AliAnalysisTaskV0ChCorrelationppsys* AddTaskV0ChCorrelationppsys(TString taskName = "", float cenMin, float cenMax, bool effCorr = 0, bool isMC=0,TString container_name_extension = "",TString fileName_extension = "",TString EffFileNameWithPath = ""){



  // Creates a V0-Ch correlations analysis task and adds it to the analysis manager.

 //==============================================================================
 
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskV0ChCorrelations", "No analysis manager to connect to.");
    return NULL;
  }

if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

 TString fileName = AliAnalysisManager::GetCommonFileName();


  fileName += ":AliAnalysisTaskV0ChCorrelationppsys";      // create a subfolder in the file
  fileName += fileName_extension.Data();

    
  // create task
  AliAnalysisTaskV0ChCorrelationppsys* task  = new AliAnalysisTaskV0ChCorrelationppsys(taskName, cenMin,cenMax,effCorr); 

  task->SetAnalysisMC(isMC);
  //------------------------------Mixing part------------------------------
  task->SetMixingTracks(5000);
  task->SetPoolSize(200); 
 //--------------------------------Variable--------------------------------
  task->SetVtxCut(10.);
  task->SetNumOfVzBins(10);    
  task->SetVtxXMin(10e-5);
  task->SetVtxYMin(10e-5);
  task->SetVtxZMin(10e-5);
  task->SetCentMin(0);
  task->SetCentMax(100.);
  //-----------------------------Track-------------------------------------

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
  task->SetLambdaLifeTimeMax(30);
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
  
  task->SetTrackPileUpCut(kTRUE);
  task->SetV0PileUpCut(kTRUE);
  task->SetEventPileUpCut(kTRUE);

  //-------------------------------------PID--------------------------------
  task->SetV0PIDSigma(3);
  //-------------------------------------------------------------------------
  mgr->AddTask(task);
    
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  
  // create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput0 
    = mgr->CreateContainer("Output",AliDirList::Class(),AliAnalysisManager::kOutputContainer, 
			                     Form("%s", fileName.Data()));
//

  AliAnalysisDataContainer *coutput2
    = mgr->CreateContainer("Output2",AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//


AliAnalysisDataContainer *coutput3
    = mgr->CreateContainer("Output3", AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//

AliAnalysisDataContainer *coutput4
    = mgr->CreateContainer("Output4",AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//



AliAnalysisDataContainer *coutput5
    = mgr->CreateContainer("Output5",AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//


AliAnalysisDataContainer *coutput6
    = mgr->CreateContainer("Output6",AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//


AliAnalysisDataContainer *coutput7
    = mgr->CreateContainer("Output7",AliDirList::Class(), AliAnalysisManager::kOutputContainer,
                                            Form("%s", fileName.Data()));//



  AliAnalysisDataContainer *cinput1 = NULL;

  TList * effList = 0x0;

  if(effCorr){

        TString eff_container_name = "";
        eff_container_name+=container_name_extension.Data();


    cinput1 = mgr->CreateContainer(Form("%s", eff_container_name.Data()),
                                    TList::Class(),
                                    AliAnalysisManager::kInputContainer);

    TFile * file = TFile::Open(Form("alien:///alice/cern.ch/user/m/manaam/Efficiencypp/%s.root",EffFileNameWithPath.Data())); 

  
if(!file) {
      printf("ERROR: efficiency file is not available!\n",EffFileNameWithPath.Data());
      //return NULL;
    }
    
    effList = (TList*)file->Get("fListCent0_100");

    if(!effList){
      printf("ERROR: no efficiency list fList%s available\n", EffFileNameWithPath.Data());
      //return NULL;
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
      return task;
}

