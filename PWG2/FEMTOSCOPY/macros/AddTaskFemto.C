//=============================================================================
//
// *** AddTaskFemto.C ***
//
// This macro initialize a complete AnalysisTask object for femtoscopy.
//
//=============================================================================

AliAnalysisTaskFemto *AddTaskFemto()
{
// Creates a proton analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFemto", "No analysis manager to connect to.");
    return NULL;
  }  

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler cann also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFemto", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  // C. Create the task, add it to manager.
  //===========================================================================
  gSystem->SetIncludePath("-I$ROOTSYS/include  -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser -I$ALICE_ROOT/include");
  if (dynamic_cast<TProofLite *> gProof) {
    gProof->Exec(".L ConfigFemtoAnalysis.C");
  }
  else if (gProof) {
    gProof->Load("ConfigFemtoAnalysis.C");
  }
  //  gROOT->LoadMacro("ConfigFemtoAnalysis.C++");

  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto");
  mgr->AddTask(taskfemto);

  // D. Configure the analysis task. Extra parameters can be used via optional
  // arguments of the AddTaskXXX() function.
  //===========================================================================
  
  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  AliAnalysisDataContainer *cout_femto  = mgr->CreateContainer("femtolist",  TList::Class(),
							       AliAnalysisManager::kOutputContainer,"Femto.ESD.root");

   mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskfemto, 0, cout_femto);

   // Return task pointer at the end
   return taskfemto;
}
