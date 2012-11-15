AliAnalysisTaskCheckCascadepp276 *AddTaskCheckCascadepp276( Int_t    minnTPCcls             = 70,
                                                            Float_t  vtxlim                 = 10.,
                                                            Float_t  vtxlimmin              = 0.,
                                                            Bool_t   fwithsdd               = kFALSE,
                                                            Bool_t   kextrasel              = kFALSE,
                                                            Bool_t   krelaunchvertexers     = kFALSE,
                                                            Bool_t   ksddonselection        = kTRUE,
                                                            Float_t  minptondaughtertracks  = 0.,
                                                            Float_t  etacutondaughtertracks = .8,
                                                            Bool_t   standardAnalysis       = kTRUE ) {

   // Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskCheckCascadeipp276", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskCheckCascadepp276", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
   AliAnalysisTaskCheckCascadepp276 *taskcheckcascadepp276 = new AliAnalysisTaskCheckCascadepp276(Form("TaskCheckCascadepp276_vtxlim%2.1f-%2.1f",vtxlim,vtxlimmin));

   taskcheckcascadepp276->SetAnalysisType               (type);
   taskcheckcascadepp276->SetRelaunchV0CascVertexers    (krelaunchvertexers);
   taskcheckcascadepp276->SetSDDSelection               (fwithsdd);          // if TRUE apply the sdd selection
   taskcheckcascadepp276->SetQualityCutZprimVtxPos      (kTRUE);             // selects vertices in +-10cm
   taskcheckcascadepp276->SetQualityCutNoTPConlyPrimVtx (kTRUE);             // retains only events with tracking + SPD vertex
   taskcheckcascadepp276->SetQualityCutTPCrefit         (kTRUE);             // requires TPC refit flag to be true to select a track
   taskcheckcascadepp276->SetQualityCutnTPCcls          (kTRUE);             // rejects tracks that have less than n clusters in the TPC
   taskcheckcascadepp276->SetQualityCutPileup           (kTRUE);
   taskcheckcascadepp276->SetWithSDDOn                  (ksddonselection);   // selects events with SDD on
   taskcheckcascadepp276->SetQualityCutMinnTPCcls       (minnTPCcls);        // minimum number of TPC clusters to accept daughter tracks
   taskcheckcascadepp276->SetExtraSelections            (kextrasel);         // used to add other selection cuts
   taskcheckcascadepp276->SetVertexRange                (vtxlim);
   taskcheckcascadepp276->SetVertexRangeMin             (vtxlimmin);
   taskcheckcascadepp276->SetMinptCutOnDaughterTracks   (minptondaughtertracks);  
   taskcheckcascadepp276->SetEtaCutOnDaughterTracks     (etacutondaughtertracks);

   mgr->AddTask(taskcheckcascadepp276);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)
   
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   if (standardAnalysis) outputFileName += ":PWGLFStrangeness.outputCheckCascadepp276";
   else                  outputFileName += Form(":PWGLFStrangeness.outputCheckCascadepp276_vtxlim%2.1f-%2.1f",vtxlim,vtxlimmin);

   Printf("AddTaskCheckCascade - Set OutputFileName : \n %s\n", outputFileName.Data() );

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("clistCasc",
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer,
							     outputFileName );
   
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cfcontPIDXiM",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cfcontPIDXiP",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   AliAnalysisDataContainer *coutput4 = mgr->CreateContainer("cfcontPIDOmegaM",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   AliAnalysisDataContainer *coutput5 = mgr->CreateContainer("cfcontPIDOmegaP",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   AliAnalysisDataContainer *coutput6 = mgr->CreateContainer("cfcontCuts",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );

   
   mgr->ConnectInput( taskcheckcascadepp276, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascadepp276, 1, coutput1);
   mgr->ConnectOutput(taskcheckcascadepp276, 2, coutput2);
   mgr->ConnectOutput(taskcheckcascadepp276, 3, coutput3);
   mgr->ConnectOutput(taskcheckcascadepp276, 4, coutput4);
   mgr->ConnectOutput(taskcheckcascadepp276, 5, coutput5);
   mgr->ConnectOutput(taskcheckcascadepp276, 6, coutput6);
   
   return taskcheckcascadepp276;
}   
