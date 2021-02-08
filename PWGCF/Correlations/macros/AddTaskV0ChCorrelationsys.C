// AddTask for AliAnalysisTaskV0ChCorrelationsys task

AliAnalysisTaskV0ChCorrelationsys* AddTaskV0ChCorrelationsys(TString taskName = "",float cenMin, float cenMax,  bool effCorr = 0, bool isMC=0,TString container_name_extension = "",TString fileName_extension = "",TString EffFileNameWithPath = ""){
  
 

  // Creates a V0-Ch correlations analysis task and adds it to the analysis manager.

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskV0ChCorrelationsys", "No analysis manager to connect to.");
    return NULL;
  }
//===================================
if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }


TString fileName = AliAnalysisManager::GetCommonFileName();

    fileName += ":AliAnalysisTaskV0ChCorrelationsys";      // create a subfolder in the file
    fileName += fileName_extension.Data();

//=====================================  
  // create task

  AliAnalysisTaskV0ChCorrelationsys* task = new AliAnalysisTaskV0ChCorrelationsys(taskName.Data(),  cenMin,cenMax,effCorr);  

  task->SetAnalysisMC(isMC);
  //------------------------------Mixing part------------------------------
  task->SetMixingTracks(500);
    task->SetPoolSize(100); 
 //--------------------------------Variable--------------------------------
  task->SetVtxCut(9.);
  task->SetNumOfVzBins(9);//    
   task->SetCentMin(0);
  task->SetCentMax(10.);
  //-----------------------------Track-------------------------------------
  task->SetTrackPtMin(1.);
  task->SetTrackPtMax(10.);
  task->SetTrackEta(0.8);
  task->SetFilterBit(768);
  task->SetAssocNcls(70);
  //------------------------------V0--------------------------------------
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
  task->SetV0PIDSigma(3);

  mgr->AddTask(task);
    
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  
  AliAnalysisDataContainer *cinput1 = 0x0;

  TList *effList = 0x0;
  
  if(effCorr){



       TString eff_container_name = "Efficiency";
       eff_container_name+=container_name_extension.Data();

    cinput1 = mgr->CreateContainer(Form("%s", eff_container_name.Data()),
                                    TList::Class(),
                                    AliAnalysisManager::kInputContainer);
//call eff  
     TGrid::Connect("alien://");
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

// your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    if(effCorr){
        cinput1->SetData(effList);
        mgr->ConnectInput(task, 1, cinput1);  
    }


// same for the output
   TString container_name = "MyOutputContainer";
   container_name += container_name_extension.Data();

   mgr->ConnectOutput(task,8,mgr->CreateContainer(container_name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

return task;

}


