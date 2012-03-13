AliAnalysisTaskCheckPerformanceCascade *AddTaskCheckPerformanceCascade(Short_t       lCollidingSystems=0  /*0 = pp, 1 = AA*/,
                                                                       const TString lMasterJobSessionFlag = "")
{
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckPerformanceCascade", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckPerformanceCascade", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
        AliAnalysisTaskCheckPerformanceCascade *taskCheckPerfCascade = new AliAnalysisTaskCheckPerformanceCascade("TaskCheckPerfCascade");
   taskCheckPerfCascade->SetCollidingSystems            (lCollidingSystems);
   taskCheckPerfCascade->SetAnalysisType                (type);
   
   taskCheckPerfCascade-> SetTriggerMaskType            ("kMB");
   taskCheckPerfCascade-> SetRelaunchV0CascVertexers    (0);     //NOTE
   taskCheckPerfCascade-> SetQualityCutZprimVtxPos      (kTRUE);
   taskCheckPerfCascade-> SetRejectEventPileUp          (kTRUE);
   taskCheckPerfCascade-> SetQualityCutNoTPConlyPrimVtx (kTRUE);
   taskCheckPerfCascade-> SetQualityCutTPCrefit         (kTRUE);
   taskCheckPerfCascade-> SetQualityCut80TPCcls         (kTRUE);
   taskCheckPerfCascade-> SetAlephParamFor1PadTPCCluster(kTRUE);
        // taskCheckPerfCascade-> SetExtraSelections            (0);
   
   
   mgr->AddTask(taskCheckPerfCascade);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)
   /*
   TString DefaultCommonFileName = AliAnalysisManager::GetCommonFileName();
   
   if(DefaultCommonFileName == "AnalysisResults.root"){
        // Just change the Common File name IF it was not change before
        // -> To avoid screwing-up the analysis train and send the output of the previous task to a non-existing file
        TString lCommonFileName = "sLHC10-CheckPerfCascade";
        if(lMasterJobSessionFlag.Length()){
                lCommonFileName += "-";
                lCommonFileName += lMasterJobSessionFlag.Data();
        }
                lCommonFileName += ".root";

        mgr->SetCommonFileName( lCommonFileName.Data() );
   }
   */

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWG2CheckPerformanceCascade";
   if (lCollidingSystems) outputFileName += "_AA_";
   else outputFileName += "_PP";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   //if(lMasterJobSessionFlag.Length()) outputFileName += lMasterJobSessionFlag.Data();
   
   Printf("AddTaskCheckPerfCascade - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistCascMC",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

   mgr->ConnectInput( taskCheckPerfCascade, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskCheckPerfCascade, 1, coutput1);
   
   return taskCheckPerfCascade;
}   
