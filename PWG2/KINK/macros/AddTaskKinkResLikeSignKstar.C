AliResonanceKinkLikeSign *AddTaskKinkResLikeSignKstar(Short_t lCollidingSystems=0  /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a V0 check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskKinkResonanceLikeSignKstar", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskKinkResonanceLikeSignKstar", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
   if (type != "ESD") {
      ::Error("AddTaskKinkResonanceLikeSignKstar", "This task needs ESD input handler");
      return NULL;
   }   

   // Create and configure the task
	AliResonanceKinkLikeSign *taskkinkreslikesignKstar = new AliResonanceKinkLikeSign("TaskResLikeSignkstar");
        taskkinkreslikesignKstar->SetPDGCodes(kKPlus, kPiPlus);
        taskkinkreslikesignKstar->SetHistoSettings(60, 0.6, 1.2, 100, 0.0, 10.0);
        taskkinkreslikesignKstar->SetEtaLimits(-0.9, 0.9);
        taskkinkreslikesignKstar->SetMaxNsigmaToVertex(4.0);
        taskkinkreslikesignKstar->SetMaxDCAxy(3.0);
        taskkinkreslikesignKstar->SetMaxDCAzaxis(3.0);
        taskkinkreslikesignKstar->SetPtTrackCut(0.25);
        taskkinkreslikesignKstar->SetMinTPCclusters(50);
        taskkinkreslikesignKstar->SetMaxChi2PerTPCcluster(3.5);
        taskkinkreslikesignKstar->SetMaxCov0(2.0);
        taskkinkreslikesignKstar->SetMaxCov2(2.0);
        taskkinkreslikesignKstar->SetMaxCov5(0.5);
        taskkinkreslikesignKstar->SetMaxCov9(0.5);
        taskkinkreslikesignKstar->SetMaxCov14(2.0);
        taskkinkreslikesignKstar->SetMinKinkRadius(120.);
        taskkinkreslikesignKstar->SetMaxKinkRadius(220.);
        taskkinkreslikesignKstar->SetQtLimits(0.05, 0.5);

   mgr->AddTask(taskkinkreslikesignKstar);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2KINKResonanceLikeSignKstar";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResLikeSignkstar",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

   mgr->ConnectInput(taskkinkreslikesignKstar, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskkinkreslikesignKstar, 1, coutput1);
   return taskkinkreslikesignKstar;
}   
