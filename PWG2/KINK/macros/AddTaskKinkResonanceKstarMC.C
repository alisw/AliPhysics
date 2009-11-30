AliAnalysisTaskKinkResonance *AddTaskKinkResonanceKstarMC(Short_t lCollidingSystems=0 /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a kink resonance task.
// Get the pointer to the existing analysis manager via the static access method.
//==============================================================================
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 if (!mgr) {
 ::Error("AddTaskKinkResonanceKstarMC", "No analysis manager to connect to.");
return NULL;
}

// Check the analysis type using the event handlers connected to the analysis manager.
//==============================================================================
if (!mgr->GetInputEventHandler()) {
 ::Error("AddTaskKinkResonanceKstarMC", "This task requires an input event handler");
 return NULL;
}
TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 if (type != "ESD") {
 ::Error("AddTaskKinkResonanceKstarMC", "This task needs an ESD input handler");
return NULL;
}
if (!mgr->GetMCtruthEventHandler()) {
 ::Error("AddTaskKinkResonanceKstarMC", "This task needs an MC handler");
 return NULL;
}

// Setup the analysis object
  
  AliResonanceKink  *kinkResonanceObjectKMC=new AliResonanceKink();
  kinkResonanceObjectKMC->InitOutputHistograms(60, 0.6, 1.2, 36, -0.9, 0.9, 100, 0.0, 10.0);
  kinkResonanceObjectKMC->SetPDGCodes(kKPlus, kPiPlus, AliResonanceKink::kKstar0); 
  kinkResonanceObjectKMC->SetAnalysisType("MC"); // "ESD" or "MC" or "DATA"
  kinkResonanceObjectKMC->SetMaxNsigmaToVertex(4.0);
  kinkResonanceObjectKMC->SetMaxDCAxy(3.0);
  kinkResonanceObjectKMC->SetMaxDCAzaxis(3.0);
  kinkResonanceObjectKMC->SetPtTrackCut(0.25);
  kinkResonanceObjectKMC->SetMinTPCclusters(50);
  kinkResonanceObjectKMC->SetMaxChi2PerTPCcluster(3.5);
  kinkResonanceObjectKMC->SetMaxCov0(2.0);
  kinkResonanceObjectKMC->SetMaxCov2(2.0);
  kinkResonanceObjectKMC->SetMaxCov5(0.5);
  kinkResonanceObjectKMC->SetMaxCov9(0.5);
  kinkResonanceObjectKMC->SetMaxCov14(2.0);
  kinkResonanceObjectKMC->SetMinKinkRadius(120.);
  kinkResonanceObjectKMC->SetMaxKinkRadius(220.);
  kinkResonanceObjectKMC->SetQtLimits(0.05, 0.5);

// Create and configure the task
AliAnalysisTaskKinkResonance *taskresonanceKstarMC = new AliAnalysisTaskKinkResonance("TaskResKstarMCKinkPID");
taskresonanceKstarMC->SetAnalysisKinkObject(kinkResonanceObjectKMC);
mgr->AddTask(taskresonanceKstarMC);

// Create ONLY the output containers for the data produced by the task.
// Get and connect other common input/output containers via the manager as below
//==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2KINKResonanceKstarMC";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResKstarMC",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

mgr->ConnectInput(taskresonanceKstarMC, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(taskresonanceKstarMC, 1, coutput1);
return taskresonanceKstarMC;
} 
