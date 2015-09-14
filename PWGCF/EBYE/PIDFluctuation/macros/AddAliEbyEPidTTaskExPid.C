//=========================================================================//
//                                                                         //
//           Analysis AddTask for Particle Ratio Fluctuation Study         //
//              Author: Deepika Rathee  || Satyajit Jenara                 //
//                      drathee@cern.ch || sjena@cern.ch 
//                        Thu Jun 19 11:44:51 CEST 2015
//                                                                         //
//=========================================================================//


AliEbyEPidTTaskExPid *AddAliEbyEPidTTaskExPid(Bool_t isModeAOD = 0,
					 Int_t aodFilterBit = 768, 
					 Int_t cuttype = 9, 
					 Bool_t IisMC  = 1, 
					 Double_t ptl = 0.2, 
					 Double_t pth = 5.0, 
					 Double_t gEta = 0.8, 
					 Double_t dcaxy = 2.4, 
					 Double_t dcaz = 3.2, 
					 TString ctaskname = "D2011") {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskFluctuations", "No analysis manager to connect to.");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFluctuations", "This task requires an input event handler");
    return NULL;
  }
  
  Bool_t isMC = IisMC;

 // Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL);
  if (isMC)
    Info("AddTaskNetParticle", "This task has MC.");
  
   TString taskname = ctaskname;


  AliEbyEPidTTaskExPid *task = new AliEbyEPidTTaskExPid(taskname.Data());
  if (!task) {
    Error("EbyEPidRatio", "Task could not be created.");
    return NULL;
  }
  
  if (isMC) task->SetIsMC();
  
  if (isModeAOD) {
    task->SetIsAOD();                       
    task->SetAODtrackCutBit(aodFilterBit);
  }
  task->SetKinematicsCuts(ptl,pth,gEta);
  if (!isModeAOD) {
    cout << " <<<<<<<<<<<<<<<<<< ADDING ESD TRACK CUTS >>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/EBYE/PIDFluctuation/macros/configureNetChargeTrackCut.C"); 
    AliESDtrackCuts *cuts = configureNetChargeTrackCut(taskname,cuttype,10001006, gEta, dcaxy,dcaz); 
    task->SetAnalysisCutObject(cuts);
  }
  mgr->AddTask(task);
  
  TString basefilename = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutt = mgr->CreateContainer("Event",TTree::Class(),AliAnalysisManager::kOutputContainer, Form("%s",basefilename.Data()));

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutt);
  return task;
}

