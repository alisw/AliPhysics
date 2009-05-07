AliResonanceKinkLikeSign *AddTaskKinkResonanceLikeSign(Short_t lCollidingSystems=0  /*0 = pp, 1 = AA*/)
{
// Creates, configures and attaches to the train a V0 check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskKinkResonanceLikeSign", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskKinkResonanceLikeSign", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
   if (type != "ESD") {
      ::Error("AddTaskKinkResonanceLikeSign", "This task needs ESD input handler");
      return NULL;
   }   

   // Create and configure the task
	AliResonanceKinkLikeSign *taskkinkreslikesign = new AliResonanceKinkLikeSign("TaskResLikeSign");
   mgr->AddTask(taskkinkreslikesign);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   TString outname = "PP";
   if (lCollidingSystems) outname = "AA";
   if (mgr->GetMCtruthEventHandler()) outname += "-MC-";
   outname += "KinkResLikeSignList.root";
	AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("KinkResLikeSign",
								   TList::Class(),
								   AliAnalysisManager::kOutputContainer,
								   outname );
                           
	mgr->ConnectInput(taskkinkreslikesign, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskkinkreslikesign, 1, coutput1);
   return taskkinkreslikesign;
}   
