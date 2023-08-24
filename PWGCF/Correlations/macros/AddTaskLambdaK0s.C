AliAnalysisTaskLambdaK0s*AddTaskLambdaK0s(TString taskName = "name",float cenMin = 0, float cenMax = 100,  bool effCorr = 0, bool isMC=0,bool eventMixing = 0,TString container_name_extension = "",TString fileName_extension = "",TString EffFileNameWithPath = ""){

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLambdaK0s", "No analysis manager to connect to.");
    return NULL;
  }
//===================================
if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }


TString fileName = AliAnalysisManager::GetCommonFileName();

    fileName += ":AliAnalysisTaskLambdaK0s";      
    fileName += fileName_extension.Data();


 AliAnalysisTaskLambdaK0s* task = new AliAnalysisTaskLambdaK0s(taskName.Data(),  cenMin,cenMax,effCorr); 

  task->SelectCollisionCandidates(AliVEvent::kINT7);
  task->SetPrimVertexCut(8.);  
  task->SetDCAToPrimVtxCut(0.1);
  task->SetDcaV0DaughtersCut(0.8);
  task->SetDCANegtoPrimVertexMinK0s(0.1);
  task->SetDCAPostoPrimVertexMinK0s(0.1);
  task->SetDCANegtoPrimVertexMinLambda(0.25);
  task->SetDCAPostoPrimVertexMinLambda(0.1);
  task->SetDCANegtoPrimVertexMinAntiLambda(0.1);
  task->SetDCAPostoPrimVertexMinAntiLambda(0.25);
  task->SetPLtimeK0s(20.);
  task->SetPLtimeLambda(25.);    
  task->SetCPA(0.95);
  task->SetLambdaCPA(0.98);
  task->SetLambdaDCA(1.);
  task->Set2DFiducial(5);
  task->SetTPCPidNsigmaCut(5);
  task->SetTPCNcls(70);
  task->SetMinCtau(0);
  task->SetMaxCtau(3);
  task->SetMaxEta(0.8);
  task->SetV0DaughterEtaCut(0.8);
  task->SetV0DaughterPtCut(0.15);
  task->SetAssoPtMin(1);
  task->SetAssoPtMax(10);
  task->SetPtTrigMin(5);
  task->SetPtTrigMax(15);
  task->SetAnalysisMC(isMC);
  task->SetEventMixing(eventMixing);
  task->SetAnalysisChXi(kTRUE);
  task->SetAnalysisChV0(kTRUE); 
  task->SetMCpileupgenandrc(1,1);
  task->SetMixingTracks(200);
  task->SetMixingPoolSize(200);

  task->SetV0Eta(0.7);
	task->SetRapidity(kFALSE);

    
   mgr->AddTask(task);
   AliAnalysisDataContainer *cinput1 = 0x0;

  TList *effList = 0x0;
  
  if(effCorr){


       TString eff_container_name = "Efficiency";
       eff_container_name+=container_name_extension.Data();


    cinput1 = mgr->CreateContainer(Form("%s", eff_container_name.Data()),
                                    TList::Class(),
                                    AliAnalysisManager::kInputContainer);
   



   //  TGrid::Connect("alien://");
    TFile * file = TFile::Open(Form("alien:///alice/cern.ch/user/%s.root",EffFileNameWithPath.Data()));
    //TFile * file = TFile::Open(Form("%s.root",EffFileNameWithPath.Data()));

    if(!file) {
      printf("ERROR: efficiency file is not available!\n",EffFileNameWithPath.Data());
    }
    
    effList = (TList*)file->Get("fListCent0_10");

    if(!effList){
      printf("ERROR: no efficiency list fList%s available\n", EffFileNameWithPath.Data());
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

   mgr->ConnectOutput(task,1,mgr->CreateContainer(container_name.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

return task;

}

