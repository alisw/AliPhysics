AliAnalysisTaskKinkResonance *AddTaskKinkResonancePhiMC(Short_t lCollidingSystems=0 /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a kink resonance task.
// Get the pointer to the existing analysis manager via the static access method.
//==============================================================================
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 if (!mgr) {
 ::Error("AddTaskKinkResonancePhiMC", "No analysis manager to connect to.");
return NULL;
}

// Check the analysis type using the event handlers connected to the analysis manager.
//==============================================================================
if (!mgr->GetInputEventHandler()) {
 ::Error("AddTaskKinkResonancePhiMC", "This task requires an input event handler");
 return NULL;
}
TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 if (type != "ESD") {
 ::Error("AddTaskKinkResonancePhiMC", "This task needs an ESD input handler");
return NULL;
}
if (!mgr->GetMCtruthEventHandler()) {
 ::Error("AddTaskKinkResonancePhiMC", "This task needs an MC handler");
 return NULL;
}

// Setup the analysis object
  
  AliResonanceKink  *kinkResonanceObjectPMC=new AliResonanceKink();
  kinkResonanceObjectPMC->InitOutputHistograms(70,0.99,1.088, 36, -0.9, 0.9, 100, 0.0, 10.0);
  kinkResonanceObjectPMC->SetPDGCodes(kKPlus, kKPlus, AliResonanceKink::kPhi); 
  kinkResonanceObjectPMC->SetAnalysisType("MC"); // "ESD" or "MC" or "DATA"
  kinkResonanceObjectPMC->SetMaxNsigmaToVertex(4.0);
  kinkResonanceObjectPMC->SetMaxDCAxy(3.0);
  kinkResonanceObjectPMC->SetMaxDCAzaxis(3.0);
  kinkResonanceObjectPMC->SetPtTrackCut(0.25);
  kinkResonanceObjectPMC->SetMinTPCclusters(50);
  kinkResonanceObjectPMC->SetMaxChi2PerTPCcluster(3.5);
  kinkResonanceObjectPMC->SetMaxCov0(2.0);
  kinkResonanceObjectPMC->SetMaxCov2(2.0);
  kinkResonanceObjectPMC->SetMaxCov5(0.5);
  kinkResonanceObjectPMC->SetMaxCov9(0.5);
  kinkResonanceObjectPMC->SetMaxCov14(2.0);
  kinkResonanceObjectPMC->SetMinKinkRadius(120.);
  kinkResonanceObjectPMC->SetMaxKinkRadius(220.);
  kinkResonanceObjectPMC->SetQtLimits(0.05, 0.5);
  kinkResonanceObjectPMC->SetUpperAbsEtaCut(0.9);

// Create and configure the task
AliAnalysisTaskKinkResonance *taskresonancePhiMC = new AliAnalysisTaskKinkResonance("TaskResPhiMCKinkPID");
taskresonancePhiMC->SetAnalysisKinkObject(kinkResonanceObjectPMC);
mgr->AddTask(taskresonancePhiMC);

// Create ONLY the output containers for the data produced by the task.
// Get and connect other common input/output containers via the manager as below
//==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2KINKResonancePhiMC";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResPhiMC",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

mgr->ConnectInput(taskresonancePhiMC, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(taskresonancePhiMC, 1, coutput1);
return taskresonancePhiMC;
} 
