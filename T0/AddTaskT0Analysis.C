//_____________________________________________________________________________
AliAnalysisTask  *AddTaskT0Analysis()
{
  //
  // add calibration task
  //
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libT0calib"); 

  cout<<"@@@ AddTaskT0Analysis "<<endl;
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libT0calib"); 


    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskT0Analysis", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskT0Analysis", "This task requires an input event handler");
    return NULL;
  }  

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskT0Analysis", "This task requires an input event handler");
    return NULL;
  }
   TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // setup task
 AliT0CalibAnalysisTask  *task1 = new AliT0CalibAnalysisTask("ObjectsTrain");
   mgr->AddTask(task1);
  

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("T0tree",TList::Class(), AliAnalysisManager::kOutputContainer,"T0AnalysisTree.root");  

  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,1,coutput1);
  return task1;
}
