AliAnalysisTaskESDfilter *AddTaskESDFilter(Bool_t useKineFilter=kTRUE, 
                                           Bool_t writeMuonAOD=kFALSE,
                                           Bool_t writeDimuonAOD=kFALSE,
                                           Bool_t usePhysicsSelection=kFALSE,
                                           Bool_t useCentralityTask=kFALSE, /*obsolete*/
                                           Int_t  tofTimeZeroType=AliESDpid::kTOF_T0,
                                           Bool_t enableTPCOnlyAODTracks=kFALSE,
                                           Bool_t disableCascades=kFALSE,
                                           Bool_t disableKinks=kFALSE)
{
// Creates a filter task and adds it to the analysis manager.

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
      ::Error("AddTaskESDFilter", "No analysis manager to connect to.");
      return NULL;
   }   
   
   // This task requires an ESD input handler and an AOD output handler.
   // Check this using the analysis manager.
   //===============================================================================
   TString type = mgr->GetInputEventHandler()->GetDataType();
   if (!type.Contains("ESD")) {
      ::Error("AddTaskESDFilter", "ESD filtering task needs the manager to have an ESD input handler.");
      return NULL;
   }   
   // Check if AOD output handler exist.
   AliAODHandler *aod_h = (AliAODHandler*)mgr->GetOutputEventHandler();
   if (!aod_h) {
      ::Error("AddTaskESDFilter", "ESD filtering task needs the manager to have an AOD output handler.");
      return NULL;
   }
   // Check if MC handler is connected in case kine filter requested
   AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
   if (!mcH && useKineFilter) {
      ::Error("AddTaskESDFilter", "No MC handler connected while kine filtering requested");
      return NULL;
   }   
   
   // Create the task, add it to the manager and configure it.
   //===========================================================================   
   // Barrel tracks filter
   AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
   esdfilter->SetTimeZeroType(tofTimeZeroType);
   if (disableCascades) esdfilter->DisableCascades();
   if  (disableKinks) esdfilter->DisableKinks();
  
   mgr->AddTask(esdfilter);
  
   // Muons
   Bool_t onlyMuon=kTRUE;
   Bool_t keepAllEvents=kTRUE;
   Int_t mcMode=(useKineFilter ? 2 : 0); // use 1 instead of 2 to get all MC information instead of just ancestors of mu tracks
   AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter",onlyMuon,keepAllEvents,mcMode);
   mgr->AddTask(esdmuonfilter);
   if(usePhysicsSelection){
     esdfilter->SelectCollisionCandidates(AliVEvent::kAny);
     esdmuonfilter->SelectCollisionCandidates(AliVEvent::kAny);
   }  

   // Filtering of MC particles (decays conversions etc)
   // this task has to go AFTER all other filter tasks
   // since it fills the AODMC array with all
   // selected MC Particles, only this way we have the 
   // AODMCparticle information available for following tasks
   AliAnalysisTaskMCParticleFilter *kinefilter = 0;
   if (useKineFilter) {
      kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Kine Filter");
      mgr->AddTask(kinefilter);
   }   

   // Cuts on primary tracks
   AliESDtrackCuts* esdTrackCutsL = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();

   // ITS stand-alone tracks
   AliESDtrackCuts* esdTrackCutsITSsa = new AliESDtrackCuts("ITS stand-alone Track Cuts", "ESD Track Cuts");
   esdTrackCutsITSsa->SetRequireITSStandAlone(kTRUE);

   // Pixel OR necessary for the electrons
   AliESDtrackCuts *itsStrong = new AliESDtrackCuts("ITSorSPD", "pixel requirement for ITS");
   itsStrong->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);


   // PID for the electrons
   AliESDpidCuts *electronID = new AliESDpidCuts("Electrons", "Electron PID cuts");
   electronID->SetTPCnSigmaCut(AliPID::kElectron, 3.);

   // standard cuts with very loose DCA
   AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE); 
   esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
   esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
   esdTrackCutsH->SetDCAToVertex2D(kTRUE);

   // standard cuts with tight DCA cut
   AliESDtrackCuts* esdTrackCutsH2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

   // standard cuts with tight DCA but with requiring the first SDD cluster instead of an SPD cluster
   // tracks selected by this cut are exclusive to those selected by the previous cut
   AliESDtrackCuts* esdTrackCutsH3 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(); 
   esdTrackCutsH3->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
   esdTrackCutsH3->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
 
   // TPC only tracks: Optionally enable the writing of TPConly information
   // constrained to SPD vertex in the filter below
   AliESDtrackCuts* esdTrackCutsTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
   esdTrackCutsTPCOnly->SetMinNClustersTPC(70);

   // Compose the filter
   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   // 1
   trackFilter->AddCuts(esdTrackCutsL);
   // 2
   trackFilter->AddCuts(esdTrackCutsITSsa);
   // 4
   trackFilter->AddCuts(itsStrong);
   itsStrong->SetFilterMask(1);        // AND with Standard track cuts 
   // 8
   trackFilter->AddCuts(electronID);
   electronID->SetFilterMask(4);       // AND with Pixel Cuts
   // 16
   trackFilter->AddCuts(esdTrackCutsH);
   // 32
   trackFilter->AddCuts(esdTrackCutsH2);
   // 64
   trackFilter->AddCuts(esdTrackCutsH3);
   // 128 , 1 << 7
   trackFilter->AddCuts(esdTrackCutsTPCOnly);
   if(enableTPCOnlyAODTracks)esdfilter->SetTPCOnlyFilterMask(128);

   // Filter with cuts on V0s
   AliESDv0Cuts*   esdV0Cuts = new AliESDv0Cuts("Standard V0 Cuts pp", "ESD V0 Cuts");
   esdV0Cuts->SetMinRadius(0.2);
   esdV0Cuts->SetMaxRadius(200);
   esdV0Cuts->SetMinDcaPosToVertex(0.05);
   esdV0Cuts->SetMinDcaNegToVertex(0.05);
   esdV0Cuts->SetMaxDcaV0Daughters(1.5);
   esdV0Cuts->SetMinCosinePointingAngle(0.99);
   AliAnalysisFilter* v0Filter = new AliAnalysisFilter("v0Filter");
   v0Filter->AddCuts(esdV0Cuts);

   esdfilter->SetTrackFilter(trackFilter);
   esdfilter->SetV0Filter(v0Filter);

   // Enable writing of Muon AODs
   esdmuonfilter->SetWriteMuonAOD(writeMuonAOD);
   
   // Enable writing of Dimuon AODs
   esdmuonfilter->SetWriteDimuonAOD(writeDimuonAOD);
 
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================
   mgr->ConnectInput  (esdfilter,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (esdfilter,  0, mgr->GetCommonOutputContainer());
   mgr->ConnectInput  (esdmuonfilter, 0, mgr->GetCommonInputContainer());
   if (useKineFilter) {
      mgr->ConnectInput  (kinefilter,  0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput (kinefilter,  0, mgr->GetCommonOutputContainer());
      AliAnalysisDataContainer *coutputEx = mgr->CreateContainer("cFilterList", TList::Class(),
								   AliAnalysisManager::kOutputContainer,"pyxsec_hists.root");
      mgr->ConnectOutput (kinefilter,  1,coutputEx);
   }   
   return esdfilter;
 }
 
