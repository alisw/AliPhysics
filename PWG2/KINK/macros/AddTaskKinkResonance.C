AliAnalysisTaskKinkResonance *AddTaskKinkResonance(Short_t lCollidingSystems=0 /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a kink resonance task.
// Get the pointer to the existing analysis manager via the static access method.
//==============================================================================
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 if (!mgr) {
 ::Error("AddTaskKinkResonanceKstar", "No analysis manager to connect to.");
return NULL;
}

// Check the analysis type using the event handlers connected to the analysis manager.
//==============================================================================
if (!mgr->GetInputEventHandler()) {
 ::Error("AddTaskKinkResonanceKstar", "This task requires an input event handler");
 return NULL;
}
TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 if (type != "ESD") {
 ::Error("AddTaskKinkResonanceKstar", "This task needs an ESD input handler");
return NULL;
}
if (!mgr->GetMCtruthEventHandler()) {
 ::Error("AddTaskKinkResonanceKstar", "This task needs an MC handler");
 return NULL;
}

// Setup the analysis object
  
  AliResonanceKink  *kinkResonanceObject=new AliResonanceKink();
  kinkResonanceObject->InitOutputHistograms(60, 0.6, 1.2);
  kinkResonanceObject->SetPDGCodes(kKPlus, kPiPlus, AliResonanceKink::kKstar0); 
  kinkResonanceObject->SetAnalysisType("ESD"); // "ESD" or "MC"
  kinkResonanceObject->SetMaxNsigmaToVertex(4.0);
  kinkResonanceObject->SetMaxDCAxy(3.0);
  kinkResonanceObject->SetMaxDCAzaxis(3.0);
  kinkResonanceObject->SetPtTrackCut(0.25);
  kinkResonanceObject->SetMinTPCclusters(50);
  kinkResonanceObject->SetMaxChi2PerTPCcluster(3.5);
  kinkResonanceObject->SetMaxCov0(2.0);
  kinkResonanceObject->SetMaxCov2(2.0);
  kinkResonanceObject->SetMaxCov5(0.5);
  kinkResonanceObject->SetMaxCov9(0.5);
  kinkResonanceObject->SetMaxCov14(2.0);

// Create and configure the task
AliAnalysisTaskKinkResonance *taskresonanceKstar = new AliAnalysisTaskKinkResonance("TaskResKinkPID");
taskresonanceKstar->SetAnalysisKinkObject(kinkResonanceObject);
mgr->AddTask(taskresonanceKstar);

// Create ONLY the output containers for the data produced by the task.
// Get and connect other common input/output containers via the manager as below
//==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2KINKResonanceKstar";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResKstar",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

mgr->ConnectInput(taskresonanceKstar, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(taskresonanceKstar, 1, coutput1);
return taskresonanceKstar;
} 
