//=============================================================================
//
// *** AliAnalysisTaskParticleEff.C ***
//
// This macro initialize an AnalysisTask for Particle Efficiency
//
//=============================================================================

AliAnalysisTaskParticleEff *AddTaskEfficiency(TString containerName="femtolist",int method=3, int filterbit=96)
{
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskEfficiency", "No analysis manager to connect to.");
    return NULL;
  }

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEfficiency", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  cout << "Found " <<type << " event handler" << endl;

  // C. Create the task, add it to manager.
  //===========================================================================
//  gSystem->SetIncludePath("-I$ROOTSYS/include  -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser -I$ALICE_ROOT/include");

//   if (TProofMgr::GetListOfManagers()->GetEntries()) {
// //     if (dynamic_cast<TProofLite *> gProof) {
// //       char *macrocommand[10000];
// //       sprintf(macrocommand, ".L %s", configMacroName);
// //       gProof->Exec(macrocommand);
// //     }
// //     else
//     gProof->Load(configMacroName);
//   }
  //  gROOT->LoadMacro("ConfigFemtoAnalysis.C++");

  AliAnalysisTaskParticleEff *taskEfficiency = new AliAnalysisTaskParticleEff("EfficiencyTask",method,filterbit);
  taskEfficiency->SetPidMethod(method);
  taskEfficiency->SetFB(filterbit);
  mgr->AddTask(taskEfficiency);

  // D. Configure the analysis task. Extra parameters can be used via optional
  // arguments of the AddTaskXXX() function.
  //===========================================================================

  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG2FEMTO";
  AliAnalysisDataContainer *cout_femto  = mgr->CreateContainer(containerName,  TList::Class(),
  							       AliAnalysisManager::kOutputContainer,outputfile);


   mgr->ConnectInput(taskEfficiency, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskEfficiency, 1, cout_femto);

   // std::ofstream ofile;
   // ofile.open("test.txt", std::ofstream::app);
   // ofile<<"AddTask:" <<method<<" "<<taskEfficiency->GetPidMethod()<<std::endl;								 
   // ofile.close();
   
   // Return task pointer at the end
   return taskEfficiency;
}
