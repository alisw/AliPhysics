AliAnalysisTaskCheckCascadePbPb *AddTaskCheckCascadePbPb( Float_t  centrlowlim         = 0.,
                                                          Float_t  centruplim          = 90.,
                                                          TString  centrest            = "V0M",
                                                          Float_t  vtxlim              = 10.,
                                                          Bool_t   kextrasel           = kFALSE,
                                                          Bool_t   krelaunchvertexers  = kFALSE,
                                                          TString  anatype             = "AOD") 
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
	AliAnalysisTaskCheckCascadePbPb *taskcheckcascadepbpb = new AliAnalysisTaskCheckCascadePbPb("TaskCheckCascadePbPb");

  taskcheckcascadepbpb->SetRelaunchV0CascVertexers    (krelaunchvertexers);
  taskcheckcascadepbpb->SetAnalysisType               (anatype);
  taskcheckcascadepbpb->SetQualityCutZprimVtxPos      (kTRUE);             // selects vertices in +-10cm
  taskcheckcascadepbpb->SetQualityCutNoTPConlyPrimVtx (kTRUE);             // retains only events with tracking + SPD vertex
  taskcheckcascadepbpb->SetQualityCutTPCrefit         (kTRUE);             // requires TPC refit flag to be true to select a track
  taskcheckcascadepbpb->SetQualityCut80TPCcls         (kTRUE);             // rejects tracks that have less than 80 clusters in the TPC
  taskcheckcascadepbpb->SetExtraSelections            (kextrasel);         // used to add other selection cuts
  taskcheckcascadepbpb->SetCentralityLowLim           (centrlowlim);       // setting centrality selection vriables
  taskcheckcascadepbpb->SetCentralityUpLim            (centruplim);
  taskcheckcascadepbpb->SetCentralityEst              (centrest);
  taskcheckcascadepbpb->SetVertexRange                (vtxlim);

  mgr->AddTask(taskcheckcascadepbpb);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)
   
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   outputFileName += ":PWGLFCheckCascadePbPb";
   if (mgr->GetMCtruthEventHandler()) outputFileName += "_MC";
   
   Printf("AddTaskCheckCascade - Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistCasc",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   
   
   mgr->ConnectInput( taskcheckcascadepbpb, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascadepbpb, 1, coutput1);
   
   return taskcheckcascadepbpb;
}   
