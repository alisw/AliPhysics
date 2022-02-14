//=============================================================================
//
// *** AddTaskFemto.C ***
// ---train version ---
// This macro initialize a complete AnalysisTask object for femtoscopy.
// from:
// alicepc100/cern/users/erogocha/PbPb2.76/2011/AOD115_0-10_newPID/to_alien_newtag/AddTaskFemto.C
// ---modified to train---
//  KM: March 25, 2013
//=============================================================================

//this line for local: AliAnalysisTaskFemto *AddTaskFemtoKchHBT(const char *configMacroName="ConfigFemtoAnalysis.C", const char *configMacroParameters="" )

AliAnalysisTaskFemto *AddTaskFemtoKchHBT1090(TString configMacroName, const char *containerName="femtolist", const char *configMacroParameters="" )
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
  cout << "Found " <<type << " event handler" << endl;

  // C. Create the task, add it to manager.
  //===========================================================================
//  gSystem->SetIncludePath("-I$ROOTSYS/include  -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser -I$ALICE_PHYSICS/include");
//  gSystem->SetIncludePath("-I$ROOTSYS/include  -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser -I$ALICE_ROOT/include");

  if (TProofMgr::GetListOfManagers()->GetEntries()) {
//     if (dynamic_cast<TProofLite *> gProof) {
//       char *macrocommand[10000];
//       sprintf(macrocommand, ".L %s", configMacroName);
//       gProof->Exec(macrocommand);
//     }
//     else
    gProof->Load(configMacroName);
  }  
  //  gROOT->LoadMacro("ConfigFemtoAnalysis.C++");

  //was befere aliroot 5.04.33: AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto",configMacroName);
  //  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto",configMacroName,kFALSE);
  //March 2013:
  //to check localy before new tag I did symbolic link on my laplot
  //in $ALICE_ROOTALICE_PHYSICS/PWGCF/FEMTOSCOPY/macros/Train/
  //[root@alicethinks Train]# ln -s /scratch/AliWork/PbPb2.76/Train2013/KchHBT KchHBT
  //
  //local------>
  //AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto",configMacroName,configMacroParameters,kFALSE);
  //local------<
  //train::::::::::::::::::::::::::::>
  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto","$ALICE_PHYSICS/"+configMacroName,configMacroParameters,kFALSE);
  //AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto","$ALICE_ROOT/"+configMacroName,configMacroParameters,kFALSE);
  //train::::::::::::::::::::::::::::<
  //10-90% only two triggers: SemiCentral and MB
  //taskfemto->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral);// this a new line for train
  taskfemto->SelectCollisionCandidates(AliVEvent::kINT7);// this a new line for train

  //0-10 % all three triggers
  //taskfemto->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);// this a new line for train
  mgr->AddTask(taskfemto);

  // D. Configure the analysis task. Extra parameters can be used via optional
  // arguments of the AddTaskXXX() function.
  //===========================================================================
  
  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();  
  outputfile += ":PWG2FEMTO";
  AliAnalysisDataContainer *cout_femto  = mgr->CreateContainer("KK07",  TList::Class(),
  							       AliAnalysisManager::kOutputContainer,outputfile);


   mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskfemto, 0, cout_femto);

   // Return task pointer at the end
   return taskfemto;
}
