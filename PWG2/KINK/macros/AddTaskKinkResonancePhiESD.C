AliAnalysisTaskKinkResonance *AddTaskKinkResonancePhiESD(Short_t lCollidingSystems=0 /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a kink resonance task.
// Get the pointer to the existing analysis manager via the static access method.
//==============================================================================
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 if (!mgr) {
 ::Error("AddTaskKinkResonancePhiESD", "No analysis manager to connect to.");
return NULL;
}

// Check the analysis type using the event handlers connected to the analysis manager.
//==============================================================================
if (!mgr->GetInputEventHandler()) {
 ::Error("AddTaskKinkResonancePhiESD", "This task requires an input event handler");
 return NULL;
}
TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 if (type != "ESD") {
 ::Error("AddTaskKinkResonancePhiESD", "This task needs an ESD input handler");
return NULL;
}
if (!mgr->GetMCtruthEventHandler()) {
 ::Error("AddTaskKinkResonancePhiESD", "This task needs an MC handler");
 return NULL;
}

// Setup the analysis object
  
  AliResonanceKink  *kinkResonanceObjectPESD=new AliResonanceKink();
  kinkResonanceObjectPESD->InitOutputHistograms(70,0.99,1.088);
  kinkResonanceObjectPESD->SetPDGCodes(kKPlus, kKPlus, AliResonanceKink::kPhi); 
  kinkResonanceObjectPESD->SetAnalysisType("ESD"); // "ESD" or "MC"
  kinkResonanceObjectPESD->SetMaxNsigmaToVertex(4.0);
  kinkResonanceObjectPESD->SetMaxDCAxy(3.0);
  kinkResonanceObjectPESD->SetMaxDCAzaxis(3.0);
  kinkResonanceObjectPESD->SetPtTrackCut(0.25);
  kinkResonanceObjectPESD->SetMinTPCclusters(50);
  kinkResonanceObjectPESD->SetMaxChi2PerTPCcluster(3.5);
  kinkResonanceObjectPESD->SetMaxCov0(2.0);
  kinkResonanceObjectPESD->SetMaxCov2(2.0);
  kinkResonanceObjectPESD->SetMaxCov5(0.5);
  kinkResonanceObjectPESD->SetMaxCov9(0.5);
  kinkResonanceObjectPESD->SetMaxCov14(2.0);

// Create and configure the task
AliAnalysisTaskKinkResonance *taskresonancePhiESD = new AliAnalysisTaskKinkResonance("TaskResPhiESDKinkPID");
taskresonancePhiESD->SetAnalysisKinkObject(kinkResonanceObjectPESD);
mgr->AddTask(taskresonancePhiESD);

// Create ONLY the output containers for the data produced by the task.
// Get and connect other common input/output containers via the manager as below
//==============================================================================
TString outname = "PP";
if (lCollidingSystems) outname = "AA";
if (mgr->GetMCtruthEventHandler()) outname += "-MC-";
outname += "KinkResonanceList.root";
AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResPhiESD",
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
						 	   outname );

mgr->ConnectInput(taskresonancePhiESD, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(taskresonancePhiESD, 1, coutput1);
return taskresonancePhiESD;
} 
