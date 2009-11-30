AliAnalysisTaskKinkResonance *AddTaskKinkResonanceL1520MC(Short_t lCollidingSystems=0 /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a kink resonance task.
// Get the pointer to the existing analysis manager via the static access method.
//==============================================================================
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
 if (!mgr) {
 ::Error("AddTaskKinkResonanceL1520MC", "No analysis manager to connect to.");
return NULL;
}

// Check the analysis type using the event handlers connected to the analysis manager.
//==============================================================================
if (!mgr->GetInputEventHandler()) {
 ::Error("AddTaskKinkResonanceL1520MC", "This task requires an input event handler");
 return NULL;
}
TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 if (type != "ESD") {
 ::Error("AddTaskKinkResonanceL1520MC", "This task needs an ESD input handler");
return NULL;
}
if (!mgr->GetMCtruthEventHandler()) {
 ::Error("AddTaskKinkResonanceL1520MC", "This task needs an MC handler");
 return NULL;
}

// Setup the analysis object
  
  AliResonanceKink  *kinkResonanceObjectLMC=new AliResonanceKink();
  kinkResonanceObjectLMC->InitOutputHistograms(100,1.4,1.8, 36, -0.9, 0.9, 100, 0.0, 10.0);
  kinkResonanceObjectLMC->SetPDGCodes(kProton, kKPlus, AliResonanceKink::kLambda1520); 
  kinkResonanceObjectLMC->SetAnalysisType("MC"); // "ESD" or "MC" or "DATA"
  kinkResonanceObjectLMC->SetMaxNsigmaToVertex(4.0);
  kinkResonanceObjectLMC->SetMaxDCAxy(3.0);
  kinkResonanceObjectLMC->SetMaxDCAzaxis(3.0);
  kinkResonanceObjectLMC->SetPtTrackCut(0.25);
  kinkResonanceObjectLMC->SetMinTPCclusters(50);
  kinkResonanceObjectLMC->SetMaxChi2PerTPCcluster(3.5);
  kinkResonanceObjectLMC->SetMaxCov0(2.0);
  kinkResonanceObjectLMC->SetMaxCov2(2.0);
  kinkResonanceObjectLMC->SetMaxCov5(0.5);
  kinkResonanceObjectLMC->SetMaxCov9(0.5);
  kinkResonanceObjectLMC->SetMaxCov14(2.0);
  kinkResonanceObjectLMC->SetMinKinkRadius(120.);
  kinkResonanceObjectLMC->SetMaxKinkRadius(220.);
  kinkResonanceObjectLMC->SetQtLimits(0.05, 0.5);

// Create and configure the task
AliAnalysisTaskKinkResonance *taskresonanceL1520MC = new AliAnalysisTaskKinkResonance("TaskResL1520MCKinkPID");
taskresonanceL1520MC->SetAnalysisKinkObject(kinkResonanceObjectLMC);
mgr->AddTask(taskresonanceL1520MC);

// Create ONLY the output containers for the data produced by the task.
// Get and connect other common input/output containers via the manager as below
//==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2KINKResonanceL1520MC";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResL1520MC",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

mgr->ConnectInput(taskresonanceL1520MC, 0, mgr->GetCommonInputContainer());
mgr->ConnectOutput(taskresonanceL1520MC, 1, coutput1);
return taskresonanceL1520MC;
} 
