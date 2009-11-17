AliAnalysisTaskKinkResonance *AddTaskKinkResonanceL1520ESD(Short_t lCollidingSystems=0 /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a kink resonance task.
// Get the pointer to the existing analysis manager via the static access method.
//==============================================================================
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 if (!mgr) {
 ::Error("AddTaskKinkResonanceL1520ESD", "No analysis manager to connect to.");
return NULL;
}

// Check the analysis type using the event handlers connected to the analysis manager.
//==============================================================================
if (!mgr->GetInputEventHandler()) {
 ::Error("AddTaskKinkResonanceL1520ESD", "This task requires an input event handler");
 return NULL;
}
TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 if (type != "ESD") {
 ::Error("AddTaskKinkResonanceL1520ESD", "This task needs an ESD input handler");
return NULL;
}
if (!mgr->GetMCtruthEventHandler()) {
 ::Error("AddTaskKinkResonanceL1520ESD", "This task needs an MC handler");
 return NULL;
}

// Setup the analysis object
  
  AliResonanceKink  *kinkResonanceObjectLESD=new AliResonanceKink();
  kinkResonanceObjectLESD->InitOutputHistograms(100,1.4,1.8);
  kinkResonanceObjectLESD->SetPDGCodes(kProton, kKPlus, AliResonanceKink::kLambda1520); 
  kinkResonanceObjectLESD->SetAnalysisType("ESD"); // "ESD" or "MC"
  kinkResonanceObjectLESD->SetMaxNsigmaToVertex(4.0);
  kinkResonanceObjectLESD->SetMaxDCAxy(3.0);
  kinkResonanceObjectLESD->SetMaxDCAzaxis(3.0);
  kinkResonanceObjectLESD->SetPtTrackCut(0.25);
  kinkResonanceObjectLESD->SetMinTPCclusters(50);
  kinkResonanceObjectLESD->SetMaxChi2PerTPCcluster(3.5);
  kinkResonanceObjectLESD->SetMaxCov0(2.0);
  kinkResonanceObjectLESD->SetMaxCov2(2.0);
  kinkResonanceObjectLESD->SetMaxCov5(0.5);
  kinkResonanceObjectLESD->SetMaxCov9(0.5);
  kinkResonanceObjectLESD->SetMaxCov14(2.0);

// Create and configure the task
AliAnalysisTaskKinkResonance *taskresonanceL1520ESD = new AliAnalysisTaskKinkResonance("TaskResL1520ESDKinkPID");
taskresonanceL1520ESD->SetAnalysisKinkObject(kinkResonanceObjectLESD);
mgr->AddTask(taskresonanceL1520ESD);

// Create ONLY the output containers for the data produced by the task.
// Get and connect other common input/output containers via the manager as below
//==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2KINKResonanceL1520ESD";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResL1520ESD",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

mgr->ConnectInput(taskresonanceL1520ESD, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(taskresonanceL1520ESD, 1, coutput1);
return taskresonanceL1520ESD;
} 
