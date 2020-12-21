//=============================================================================
//
// *** AddTaskFemtoKpm.C ***
//
// This macro initialize a complete AnalysisTask object for femtoscopy
// with MC AOD LEGO Train name: CF_PbPb_MC_AOD
// Konstantin.Mikhaylov@cern.ch May 2017
//Macro parameters:
//"TaskFemtoKpm","PWGCF/FEMTOSCOPY/macros/Train/KpmHBT_MC/ConfigFemtoAnalysis.C","femtolist_KpKm",kFALSE
//=============================================================================

AliAnalysisTaskFemto *AddTaskFemtoKpmLocal/*(const char *configMacroName="ConfigFemtoAnalysis.C", const char *configMacroParameters="")*/
(TString aName,
 TString aConfigMacro,
 TString aConfigParams,
 Bool_t aVerbose)
{
// Creates a K+K- analysis task and adds it to the analysis manager.

  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFemtoKpm", "No analysis manager to connect to.");
    return NULL;
  }

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler cann also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFemtoKpm", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  cout << "AddTaskFemtoKpm->Found " <<type << " event handler" << endl;

  // C. Create the task, add it to manager.
  //===========================================================================
//  gSystem->SetIncludePath("-I$ROOTSYS/include  -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser -I$ALICE_PHYSICS/include");

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

  //TString for taskfemto:
  TString tTaskName = aName;
  TString tConfigMacro = "$ALICE_PHYSICS/"+aConfigMacro;
  TString tConfigParams = aConfigParams;
  Bool_t tVerbose=aVerbose;
  cout<<"------> AddTaskFemtoKpm PARAMETERS set to:----->"<<endl;
  cout<<"AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("
      <<tTaskName<<" , "
      <<tConfigMacro<<" , "
      <<tConfigParams<<" , "
      <<tVerbose<<" )"
      <<endl;
  
  TString tContainerName = tConfigParams;
  tConfigParams = '"' + tConfigParams + '"';

  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto
    (
     tTaskName,
     tConfigMacro,
     //tConfigParams,//this param get an error in local: Symbol femtolist_KpKm is not defined in current scope
     tVerbose
     );
  
  taskfemto->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);// this line for Run1
  // or AliVEvent::kMB | AliVEvent::kCentral

  //taskfemto->SelectCollisionCandidates(AliVEvent::kINT7);// this line for Run2
  mgr->AddTask(taskfemto);

  // D. Configure the analysis task. Extra parameters can be used via optional
  // arguments of the AddTaskXXX() function.
  //===========================================================================

  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWGCFEMTO_KpKm";
  AliAnalysisDataContainer *cout_femto  = mgr->CreateContainer
    //("femtolist_KpKm",
    (tConfigParams,
     TList::Class(),
     AliAnalysisManager::kOutputContainer,outputfile
     );


   mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskfemto, 0, cout_femto);

   // Return task pointer at the end
   return taskfemto;
}
