AliAnalysisTaskQAMultistrange *AddTaskQAMultistrange( Bool_t   isMC                   = kFALSE,
                                                      Int_t    minnTPCcls             = 70,
                                                      Float_t  centrlowlim            = 0.,
                                                      Float_t  centruplim             = 90.,
                                                      TString  centrest               = "V0M",
                                                      Bool_t   kusecleaning           = kTRUE, 
                                                      Float_t  vtxlim                 = 10.,
                                                      TString  collidingSystem        = "PbPb",
                                                      Bool_t   SDDSelection           = kFALSE,
                                                      Bool_t   withSDD                = kFALSE,
                                                      Float_t  minptondaughtertracks  = 0.,
                                                      Float_t  etacutondaughtertracks = 0.8) {

   // Creates, configures and attaches to the train a cascades check task.
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskQAMultistrange", "No analysis manager to connect to.");
      return NULL;
   }   

   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
      ::Error("AddTaskQAMultistrange", "This task requires an input event handler");
      return NULL;
   }   
   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

   // Create and configure the task
   AliAnalysisTaskQAMultistrange *taskcheckcascade = new AliAnalysisTaskQAMultistrange("TaskCheckCascade");

   taskcheckcascade->SetIsMC                       (isMC);
   taskcheckcascade->SetAnalysisType               (type);
   taskcheckcascade->SetCollidingSystem            (collidingSystem);
   taskcheckcascade->SetSDDselection               (SDDSelection);
   taskcheckcascade->SetQualityCutZprimVtxPos      (kTRUE);             // selects vertices in +-10cm
   taskcheckcascade->SetQualityCutNoTPConlyPrimVtx (kTRUE);             // retains only events with tracking + SPD vertex
   taskcheckcascade->SetQualityCutTPCrefit         (kTRUE);             // requires TPC refit flag to be true to select a track
   taskcheckcascade->SetQualityCutnTPCcls          (kTRUE);             // rejects tracks that have less than n clusters in the TPC
   taskcheckcascade->SetQualityCutMinnTPCcls       (minnTPCcls);        // minimum number of TPC clusters to accept daughter tracks
   taskcheckcascade->SetQualityCutPileup           (kFALSE);
   taskcheckcascade->SetwithSDD                    (withSDD);
   taskcheckcascade->SetCentralityLowLim           (centrlowlim);       // setting centrality selection vriables
   taskcheckcascade->SetCentralityUpLim            (centruplim);
   taskcheckcascade->SetCentralityEst              (centrest);
   taskcheckcascade->SetUseCleaning                (kusecleaning);
   taskcheckcascade->SetVertexRange                (vtxlim);
   taskcheckcascade->SetMinptCutOnDaughterTracks   (minptondaughtertracks);  
   taskcheckcascade->SetEtaCutOnDaughterTracks     (etacutondaughtertracks);

   mgr->AddTask(taskcheckcascade);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   // User file name (if need be)
   
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   outputFileName += ":PWGLFStrangeness.outputCheckCascade";

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cfcontCuts",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cfcontMCgen",
                                                             AliCFContainer::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             outputFileName );
   
   mgr->ConnectInput( taskcheckcascade, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskcheckcascade, 1, coutput1);
   mgr->ConnectOutput(taskcheckcascade, 2, coutput2);  
 
   return taskcheckcascade;
}   

