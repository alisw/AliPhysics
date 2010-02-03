AliAnalysisTaskCheckCascade *AddTaskCheckCascade(Short_t       lCollidingSystems     = 0  /*0 = pp, 1 = AA*/,
						 const TString lMasterJobSessionFlag = "")
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckCascade", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckCascade", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
	AliAnalysisTaskCheckCascade *taskcheckcascade = new AliAnalysisTaskCheckCascade("TaskCheckCascade");
   taskcheckcascade->SetCollidingSystems(lCollidingSystems);
   taskcheckcascade->SetAnalysisType(type);
   
   mgr->AddTask(taskcheckcascade);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)
   /*
   TString lCommonFileName = "sLHC09dxx-CheckCascade";
   if(lMasterJobSessionFlag.Length()){
        lCommonFileName += "-";
        lCommonFileName += lMasterJobSessionFlag.Data();
   }
        lCommonFileName += ".root"; 
   mgr->SetCommonFileName( lCommonFileName.Data() );
   */

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2CheckCascade";
   if (lCollidingSystems) outputFileName += "_AA_";
   else outputFileName += "_PP_";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "MC_";
   if(lMasterJobSessionFlag.Length()) outputFileName += lMasterJobSessionFlag.Data();

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistCasc",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   
   mgr->ConnectInput( taskcheckcascade, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascade, 1, coutput1);
   
   return taskcheckcascade;
}   
