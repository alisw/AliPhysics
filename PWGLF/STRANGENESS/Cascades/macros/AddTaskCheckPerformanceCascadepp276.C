AliAnalysisTaskCheckPerformanceCascadepp276 *AddTaskCheckPerformanceCascadepp276( Int_t    minnTPCcls             = 70,
                                                                                  Float_t  vtxlim                 = 10.,
                                                                                  Float_t  vtxlimmin              = 0.,
                                                                                  Bool_t   fwithsdd               = kFALSE,
                                                                                  Bool_t   kextrasel              = kFALSE,
                                                                                  Bool_t   kacccut                = kFALSE,
                                                                                  Bool_t   krelaunchvertexers     = kFALSE,
                                                                                  Bool_t   ksddonselection        = kFALSE,
                                                                                  Float_t  minptondaughtertracks  = 0.0,
                                                                                  Float_t  etacutondaughtertracks = 0.8) {
    
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
   //==============================================================================
   TString tasknameperf = "TaskCheckPerformanceCascadepp276";
     tasknameperf += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
   AliAnalysisTaskCheckPerformanceCascadepp276 *taskCheckPerfCascadepp276 = new AliAnalysisTaskCheckPerformanceCascadepp276(tasknameperf);
     taskCheckPerfCascadepp276->SetAnalysisType               (type);                   // "ESD" or "AOD"
     taskCheckPerfCascadepp276->SetRelaunchV0CascVertexers    (krelaunchvertexers);     // choose if reconstruct the vertex of V0 in the cascades
     taskCheckPerfCascadepp276->SetSDDSelection               (fwithsdd);               // choose if apply SDD event selection
     taskCheckPerfCascadepp276->SetQualityCutZprimVtxPos      (kTRUE);                  // choose if apply Z vtx PV position event selection
     taskCheckPerfCascadepp276->SetRejectEventPileUp          (kTRUE);                  // choose if apply no Pileup event selection
     taskCheckPerfCascadepp276->SetQualityCutNoTPConlyPrimVtx (kTRUE);                  // choose if apply no TPC only event selection
     taskCheckPerfCascadepp276->SetQualityCutTPCrefit         (kTRUE);                  // choose if apply TPC refit on daughter tracks
     taskCheckPerfCascadepp276->SetQualityCutnTPCcls          (kTRUE);                  // choose if apply n TPC cluster selection on daughter tracks
     taskCheckPerfCascadepp276->SetWithSDDOn                  (ksddonselection);        // which SDD selection do you want apply? [if kTRUE select SDDon events]
     taskCheckPerfCascadepp276->SetQualityCutMinnTPCcls       (minnTPCcls);             // which value do you want apply for the minTPCcls cut?
     taskCheckPerfCascadepp276->SetExtraSelections            (kextrasel);              // choose if apply the extra selection of cascade reco.
     taskCheckPerfCascadepp276->SetApplyAccCut                (kacccut);                // choose if apply acceptance cut
     taskCheckPerfCascadepp276->SetVertexRange                (vtxlim);                 // which higher value do you want apply for vtx Z cut?
     taskCheckPerfCascadepp276->SetVertexRangeMin             (vtxlimmin);              // which lower value do you want apply for vtx Z cut?
     taskCheckPerfCascadepp276->SetMinptCutOnDaughterTracks   (minptondaughtertracks);  // which value do you want apply for cut on min pt daughter track?
     taskCheckPerfCascadepp276->SetEtaCutOnDaughterTracks     (etacutondaughtertracks); // which value do you want apply for cut on eta daughter track?
    
   mgr->AddTask(taskCheckPerfCascadepp276);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   // Directory name
   TString outputFileNamePerf = Form("%s:PWGLFStrangeness.outputCheckPerformanceCascadepp276", AliAnalysisManager::GetCommonFileName());
   // Objects name
   TString outputnameperf0 = "clistCascPerf";
   TString outputnameperf1 = "cfcontPIDAsXiM";
   TString outputnameperf2 = "cfcontPIDAsXiP";
   TString outputnameperf3 = "cfcontPIDAsOmegaM";
   TString outputnameperf4 = "cfcontPIDAsOmegaP";
   TString outputnameperf5 = "cfcontAsCuts";
     outputnameperf0 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
     outputnameperf1 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
     outputnameperf2 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
     outputnameperf3 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
     outputnameperf4 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
     outputnameperf5 += Form("_minnTPCcls%i_vtxlim%.1f-%.1f_minptdghtrk%.1f_etacutdghtrk%.1f",minnTPCcls,vtxlim,vtxlimmin,minptondaughtertracks,etacutondaughtertracks);
    
    
   //Save objects into the train common file
   AliAnalysisDataContainer *coutputperf1 = mgr->CreateContainer(outputnameperf0,
							                                     TList::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
							                                     outputFileNamePerf );
   AliAnalysisDataContainer *coutputperf2 = mgr->CreateContainer(outputnameperf1,
                                                                 AliCFContainer::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileNamePerf );
   AliAnalysisDataContainer *coutputperf3 = mgr->CreateContainer(outputnameperf2,
                                                                 AliCFContainer::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileNamePerf );
   AliAnalysisDataContainer *coutputperf4 = mgr->CreateContainer(outputnameperf3,
                                                                 AliCFContainer::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileNamePerf );
   AliAnalysisDataContainer *coutputperf5 = mgr->CreateContainer(outputnameperf4,
                                                                 AliCFContainer::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileNamePerf );
   AliAnalysisDataContainer *coutputperf6 = mgr->CreateContainer(outputnameperf5,
                                                                 AliCFContainer::Class(),
                                                                 AliAnalysisManager::kOutputContainer,
                                                                 outputFileNamePerf );

   mgr->ConnectInput( taskCheckPerfCascadepp276, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 1, coutputperf1);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 2, coutputperf2);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 3, coutputperf3);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 4, coutputperf4);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 5, coutputperf5);
   mgr->ConnectOutput(taskCheckPerfCascadepp276, 6, coutputperf6);
   
   return taskCheckPerfCascadepp276;
}   
