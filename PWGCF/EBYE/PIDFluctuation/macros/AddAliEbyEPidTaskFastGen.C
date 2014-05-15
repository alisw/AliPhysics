//=========================================================================//
//                                                                         //
//           Analysis AddTask for  Particle Ratio in FastGen Analysis      //
//              Author: Satyajit Jena || Deepika Rathee                    //
//                      sjena@cern.ch || drathee@cern.ch                   //
//                           v1.0
//                                                                         //
//=========================================================================//

void AddAliEbyEPidTaskFastGen(Double_t vx, Double_t vy, Double_t vz, Double_t ptl, Double_t pth, const Int_t ieta) {
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

  AliEbyEPidTaskFastGen *task[ieta];
  AliAnalysisDataContainer *cout[ieta];

  for(Int_t i = 0; i < ieta ; i ++) {
    Double_t eta = 0.05 + 0.05*i;
    TString taskname1 = taskname;
    taskname1 += Form("ETA_%.2f",eta);
    task[i] = new AliEbyEPidTaskFastGen(taskname1.Data());
    task[i]->SetVertexDiamond(vx,vy,vz);
    task[i]->SetKinematicsCutsAOD(ptl,pth,eta);
    mgr->AddTask(task[i]);
    

    cout[i] = mgr->CreateContainer(Form("%s",taskname1.Data()),TList::Class(), 
				   AliAnalysisManager::kOutputContainer,
				   Form("%s:CFEbyE_PR",basefilename.Data()));
    mgr->ConnectInput(task[i], 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task[i], 1, cout[i]);
    
  }
  
  return;
}
