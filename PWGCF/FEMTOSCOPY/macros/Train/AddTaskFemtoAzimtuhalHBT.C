//=============================================================================
//
// *** AddTaskFemtoAzimtuhalHBT.C ***
//
// This macro initialize a complete AnalysisTask object for femtoscopy.
//
//=============================================================================

AliAnalysisTaskFemto *AddTaskFemtoAzimtuhalHBT(TString configMacroName, const char *containerName="femtolist", const char *configMacroParameters="" )
{
// Creates a proton analysis task and adds it to the analysis manager.
  
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskFemtoAzimtuhalHBT", "No analysis manager to connect to.");
    return NULL;
  }  
	TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
	Bool_t useEtaGap=kFALSE;
	Float_t etaGap=0.;
	Bool_t posTPCAOD=kFALSE;
	TString containername = "EPStat_ttd";

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler cann also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskFemtoAzimtuhalHBT", "This task requires an input event handler");
    return NULL;
  }  
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 // cout << "Found " <<type << " event handler" << endl;

	//gROOT->LoadMacro("AliEPSelectionTask3.cxx+g");
	//gROOT->LoadMacro("AddTaskEventplane.C");
	//AliEPSelectionTask3* epsel = AddTaskEventplane();
	AliEPSelectionTask *eventplaneTask = new AliEPSelectionTask("EventplaneSelection4");
	eventplaneTask->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral | AliVEvent::kCentral);
	if (inputDataType == "AOD"){
		eventplaneTask->SetInput("AOD");
	}
	eventplaneTask->SetTrackType("TPC");
	eventplaneTask->SetUsePtWeight();
	eventplaneTask->SetUsePhiWeight();
	eventplaneTask->SetSaveTrackContribution();
	if(useEtaGap){
		eventplaneTask->SetSubeventsSplitMethod(AliEPSelectionTask3::kEta); 
		eventplaneTask->SetEtaGap(etaGap); 
	}
	if(posTPCAOD){
		eventplaneTask->SetPersonalAODtrackCuts(128,0.,0.8,0.15,20.);
		eventplaneTask->SetSubeventsSplitMethod(AliEPSelectionTask3::kRandom);
	}
	
	
	mgr->AddTask(eventplaneTask);
	
	
	
	
  // C. Create the task, add it to manager.
  //===========================================================================
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

  AliAnalysisTaskFemto *taskfemto = new AliAnalysisTaskFemto("TaskFemto","$ALICE_ROOT/"+configMacroName,configMacroParameters,kFALSE);
  	//taskfemto->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral) ;

	mgr->AddTask(taskfemto);

  // D. Configure the analysis task. Extra parameters can be used via optional
  // arguments of the AddTaskXXX() function.
  //===========================================================================
  
  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();  
 // outputfile += ":PWG2FEMTO";
  AliAnalysisDataContainer *cout_femto  = mgr->CreateContainer(containerName,  TList::Class(),
  							       AliAnalysisManager::kOutputContainer,outputfile);


   mgr->ConnectInput(taskfemto, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskfemto, 0, cout_femto);
   AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(containername, TList::Class(), AliAnalysisManager::kOutputContainer,outputfile);
   
   mgr->ConnectInput(eventplaneTask, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(eventplaneTask,1,coutput1);

   // Return task pointer at the end
   return taskfemto;
}
