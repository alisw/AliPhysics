AliAnalysisTaskCheckPerformanceCascadepp276 *AddTaskCheckPerformanceCascadepp276( Int_t    minnTPCcls             = 70,
                                                                                  Float_t  vtxlim                 = 10.,
                                                                                  Bool_t   kextrasel              = kFALSE,
                                                                                  Bool_t   kacccut                = kFALSE,
                                                                                  Bool_t   krelaunchvertexers     = kFALSE,
                                                                                  Bool_t   ksddonselection        = kTRUE,
                                                                                  Float_t  minptondaughtertracks  = 0.,
                                                                                  Float_t  etacutondaughtertracks = 9999999. ) {
// Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckPerformanceCascadepp276", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckPerformanceCascadepp276", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
   AliAnalysisTaskCheckPerformanceCascadepp276 *taskCheckPerfCascadepp276 = new AliAnalysisTaskCheckPerformanceCascadepp276("TaskCheckPerformanceCascadepp276");

   taskCheckPerfCascadepp276->SetAnalysisType               (type);
   taskCheckPerfCascadepp276->SetRelaunchV0CascVertexers    (krelaunchvertexers);     
   taskCheckPerfCascadepp276->SetQualityCutZprimVtxPos      (kTRUE);
   taskCheckPerfCascadepp276->SetRejectEventPileUp          (kFALSE);
   taskCheckPerfCascadepp276->SetQualityCutNoTPConlyPrimVtx (kTRUE);
   taskCheckPerfCascadepp276->SetQualityCutTPCrefit         (kTRUE);
   taskCheckPerfCascadepp276->SetQualityCutnTPCcls          (kTRUE);             
   taskCheckPerfCascadepp276->SetSDDSelection               (ksddonselection);
   taskCheckPerfCascadepp276->SetQualityCutMinnTPCcls       (minnTPCcls);    
   taskCheckPerfCascadepp276->SetExtraSelections            (kextrasel);
   taskCheckPerfCascadepp276->SetApplyAccCut                (kacccut);
   taskCheckPerfCascadepp276->SetVertexRange                (vtxlim);
   taskCheckPerfCascadepp276->SetMinptCutOnDaughterTracks   (minptondaughtertracks); 
   taskCheckPerfCascadepp276->SetEtaCutOnDaughterTracks     (etacutondaughtertracks);
 
   mgr->AddTask(taskCheckPerfCascadepp276);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)

   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWGLFStrangeness.outputCheckPerformanceCascadepp276";
   Printf("AddTaskCheckPerformanceCascadepp276 - Set OutputFileName : \n %s\n", outputFileName.Data() );
   
   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistCascPerf",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );

   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cfcontPIDAsXiM",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cfcontPIDAsXiP",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("cfcontPIDAsOmegaM",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("cfcontPIDAsOmegaP",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("cfcontAsCuts",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );



   mgr->ConnectInput( taskCheckPerfCascadepp276, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 1, coutput1);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 2, coutput2);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 3, coutput3);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 4, coutput4);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 5, coutput5);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 6, coutput6);
   
   return taskCheckPerfCascadepp276;
}   
