AliResonanceKinkLikeSign *AddTaskKinkResLikeSignL1520(Short_t lCollidingSystems=0  /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a V0 check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskKinkResonanceLikeSignL1520", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskKinkResonanceLikeSignL1520", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
   if (type != "ESD") {
      ::Error("AddTaskKinkResonanceLikeSignL1520", "This task needs ESD input handler");
      return NULL;
   }   

   // Create and configure the task
	AliResonanceKinkLikeSign *taskkinkreslikesignL1520 = new AliResonanceKinkLikeSign("TaskResLikeSignlam");
        taskkinkreslikesignL1520->SetPDGCodes(kProton, kKPlus);
        taskkinkreslikesignL1520->SetHistoSettings(100,1.4,1.8, 100, 0.0, 10.0);
        taskkinkreslikesignL1520->SetEtaLimits(-0.9, 0.9);
        taskkinkreslikesignL1520->SetMaxNsigmaToVertex(4.0);
        taskkinkreslikesignL1520->SetMaxDCAxy(3.0);
        taskkinkreslikesignL1520->SetMaxDCAzaxis(3.0);
        taskkinkreslikesignL1520->SetPtTrackCut(0.25);
        taskkinkreslikesignL1520->SetMinTPCclusters(50);
        taskkinkreslikesignL1520->SetMaxChi2PerTPCcluster(3.5);
        taskkinkreslikesignL1520->SetMaxCov0(2.0);
        taskkinkreslikesignL1520->SetMaxCov2(2.0);
        taskkinkreslikesignL1520->SetMaxCov5(0.5);
        taskkinkreslikesignL1520->SetMaxCov9(0.5);
        taskkinkreslikesignL1520->SetMaxCov14(2.0);
        taskkinkreslikesignL1520->SetMinKinkRadius(120.);
        taskkinkreslikesignL1520->SetMaxKinkRadius(220.);
        taskkinkreslikesignL1520->SetQtLimits(0.05, 0.5);

   mgr->AddTask(taskkinkreslikesignL1520);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2KINKResonanceLikeSignL1520";
   if (lCollidingSystems) outputFileName += "_AA";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResLikeSignL1520",
                                                             TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   mgr->ConnectInput(taskkinkreslikesignL1520, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskkinkreslikesignL1520, 1, coutput1);
   return taskkinkreslikesignL1520;
}   
