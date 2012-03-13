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
   taskcheckcascade-> SetCollidingSystems           (lCollidingSystems);
   taskcheckcascade-> SetAnalysisType               (type);
   
   taskcheckcascade-> SetTriggerMaskType            ("kMB");
   taskcheckcascade-> SetRelaunchV0CascVertexers    (0);     //NOTE
   taskcheckcascade-> SetQualityCutZprimVtxPos      (kTRUE);
   taskcheckcascade-> SetRejectEventPileUp          (kTRUE);
   taskcheckcascade-> SetQualityCutNoTPConlyPrimVtx (kTRUE);
   taskcheckcascade-> SetQualityCutTPCrefit         (kTRUE);
   taskcheckcascade-> SetQualityCut80TPCcls         (kTRUE);
   taskcheckcascade-> SetAlephParamFor1PadTPCCluster(kTRUE);
        // taskcheckcascade-> SetExtraSelections            (0);
   taskcheckcascade-> SetAngularCorrelationType     ("TrigLeadingTrck-AssoCasc"); // 1.1 - "TrigAnyCasc-AssoAnyPrim", 1.2 - "TrigCascLeading-AssoAnyPrim", 2. - "TrigLeadingTrck-AssoCasc"
   
   mgr->AddTask(taskcheckcascade);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)
   /*
   TString DefaultCommonFileName = AliAnalysisManager::GetCommonFileName();
   
   if(DefaultCommonFileName == "AnalysisResults.root"){
        // Just change the Common File name IF it was not change before
        // -> To avoid screwing-up the analysis train and send the output of the previous task to a non-existing file
        TString lCommonFileName = "sLHC10-CheckCascade";
        if(lMasterJobSessionFlag.Length()){
                lCommonFileName += "-";
                lCommonFileName += lMasterJobSessionFlag.Data();
        }
                lCommonFileName += ".root"; 
        
        mgr->SetCommonFileName( lCommonFileName.Data() );
   }
   */
   
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   outputFileName += ":PWG2CheckCascade";
   if (lCollidingSystems) outputFileName += "_AA_";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   //if(lMasterJobSessionFlag.Length()) outputFileName += lMasterJobSessionFlag.Data();
   
   Printf("AddTaskCheckCascade - Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistCasc",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   
//    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cPaveTextBookKeeping",
// 							     TPaveText::Class(),
// 							     AliAnalysisManager::kOutputContainer,
// 							     outputFileName );
   
   
   mgr->ConnectInput( taskcheckcascade, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascade, 1, coutput1);
   //mgr->ConnectOutput(taskcheckcascade, 2, coutput2);
   
   return taskcheckcascade;
}   
