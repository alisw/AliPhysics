//=========================================================================//
//                                                                         //
//           Analysis AddTask for  Particle Ratio in FastGen Analysis      //
//              Author: Satyajit Jena || Deepika Rathee                    //
//                      sjena@cern.ch || drathee@cern.ch                   //
//                           v1.0
//                                                                         //
//=========================================================================//

void AddAliEbyEPidRatioTaskOnFly(Double_t ptl=0.1, Double_t pth=3, const Int_t ieta = 20) {
  TString taskname = "FG";
  taskname += "_";
  taskname += "MC";
  taskname += "_";
  taskname += Form("PT_%.1f_%.1f", ptl, pth);
  taskname += "_";
  

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); 

  TString basefilename = AliAnalysisManager::GetCommonFileName();

  AliEbyEPidRatioTaskOnFly *task[ieta];
  AliAnalysisDataContainer *coutt[ieta];

  for(Int_t i = 0; i < ieta ; i ++) {
    Double_t eta = 0.1 + 0.1*i;
    TString taskname1 = taskname;
    taskname1 += Form("ETA_%.2f",eta);
    task[i] = new AliEbyEPidRatioTaskOnFly(taskname1.Data());
    task[i]->SetKinematicsCutsAOD(ptl,pth,eta);
    mgr->AddTask(task[i]);
    

    coutt[i] = mgr->CreateContainer(Form("%s",taskname1.Data()),TList::Class(), 
				   AliAnalysisManager::kOutputContainer,
				   Form("%s:GenPR",basefilename.Data()));
    mgr->ConnectInput(task[i], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[i], 1, coutt[i]);
    
  }
  
  return;
}
