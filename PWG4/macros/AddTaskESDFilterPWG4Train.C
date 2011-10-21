

AliAnalysisTaskESDfilter *AddTaskESDFilterPWG4Train(Bool_t useKineFilter=kTRUE, 
                                           Bool_t writeMuonAOD=kFALSE,
                                           Bool_t writeDimuonAOD=kFALSE,
					   Bool_t usePhysicsSelection=kFALSE,
					   Bool_t useCentralityTask=kFALSE)
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

   // Make the AOD a little bit lighter and filtering faster
   
   esdfilter->DisableCascades();
   //   esdfilter->DisableV0s();
   esdfilter->DisableKinks();
   //   esdfilter->DisableTracks();
   esdfilter->DisablePmdClusters();
   //   esdfilter->DisableCaloClusters();
   //   esdfilter->DisableCells(); 
   esdfilter->DisableTracklets();
   esdfilter->SetWriteHybridGlobalConstrainedOnly(kTRUE);

   AliAnalysisTaskAODCentralityMaker* ctask = 0;
   if (useCentralityTask) {
       ctask = new AliAnalysisTaskAODCentralityMaker("AODCentralityMaker");
       ctask->SetDeltaAODFileName("AliAODCentrality.root");
       mgr->AddTask(ctask);
   }


   mgr->AddTask(esdfilter);
   // Muons
   //   AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
   //   mgr->AddTask(esdmuonfilter);
   if(usePhysicsSelection){
     esdfilter->SelectCollisionCandidates(AliVEvent::kAny);
     //     esdmuonfilter->SelectCollisionCandidates(AliVEvent::kAny);
   }  

   // Filtering of MC particles (decays conversions etc)
   // this task has to go AFTER all other filter tasks
   // since it fills the AODMC array with all
   // selected MC Particles, only this way we have the 
   // AODMCparticle information available for following tasks
   AliAnalysisTaskMCParticleFilter *kinefilter = 0;
   if (useKineFilter) {
      kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Kine Filter");
      if(usePhysicsSelection)kinefilter->SelectCollisionCandidates(AliVEvent::kAny);
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

   // tighter cuts on primary particles for high pT tracks
   // take the standard cuts, which include already 
   // ITSrefit and use only primaries...

   // ITS cuts for new jet analysis 
   gROOT->LoadMacro("$ALICE_ROOT/PWG4/macros/CreateTrackCutsPWG4.C");
   AliESDtrackCuts* esdTrackCutsHG0 = CreateTrackCutsPWG4(10001005);

   // throw out tracks with too low number of clusters in
   // the first pass (be consistent with TPC only tracks)
   // N.B. the number off crossed rows still acts on the tracks after
   // all iterations if we require tpc standalone, number of clusters
   // and chi2 TPC cuts act on track after the first iteration
   //   esdTrackCutsH0->SetRequireTPCStandAlone(kTRUE);
   //   esdTrackCutsH0->SetMinNClustersTPC(80); // <--- first pass


   // the complement to the one with SPD requirement
   AliESDtrackCuts* esdTrackCutsHG1 = CreateTrackCutsPWG4(10011005);

   // the tracks that must not be taken pass this cut and
   // non HGC1 and HG
   AliESDtrackCuts* esdTrackCutsHG2 = CreateTrackCutsPWG4(10021005);

   


   // standard cuts also used in R_AA analysis
   gROOT->LoadMacro("$ALICE_ROOT/PWG0/dNdPt/macros/CreatedNdPtTrackCuts.C");
   AliESDtrackCuts* esdTrackCutsH2 = CreatedNdPtTrackCuts(200);
   esdTrackCutsH2->SetMaxChi2PerClusterITS(36.);
   esdTrackCutsH2->SetPtRange(0.15,1E10);

   AliESDtrackCuts* esdTrackCutsGCOnly = CreateTrackCutsPWG4(10041005);

   // TPC only tracks
   AliESDtrackCuts* esdTrackCutsTPCCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
   esdTrackCutsTPCCOnly->SetMinNClustersTPC(70);



   // Compose the filter
   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   // 1, 1<<0
   trackFilter->AddCuts(esdTrackCutsL);
   // 2 1<<1
   trackFilter->AddCuts(esdTrackCutsITSsa);
   // 4 1<<2
   trackFilter->AddCuts(itsStrong);
   itsStrong->SetFilterMask(1);        // AND with Standard track cuts 
   // 8 1<<3
   trackFilter->AddCuts(electronID);
   electronID->SetFilterMask(4);       // AND with Pixel Cuts
   // 16 1<<4
   trackFilter->AddCuts(esdTrackCutsHG0);
   // 32 1<<5
   trackFilter->AddCuts(esdTrackCutsHG1);
   // 64 1<<6
   trackFilter->AddCuts(esdTrackCutsHG2);
   // 128 1<<7
   trackFilter->AddCuts(esdTrackCutsHG0); // add once more for tpc only tracks
   // 256 1<<8
   trackFilter->AddCuts(esdTrackCutsGCOnly);
   // 512 1<<9                         
   trackFilter->AddCuts(esdTrackCutsHG1); // add once more for tpc only tracks
   // 1024 1<<10                        
   trackFilter->AddCuts(esdTrackCutsH2); // add r_aa cuts
   // 2048 1<<11                        
   trackFilter->AddCuts(esdTrackCutsTPCCOnly); // add QM TPC only track cuts
   

   esdfilter->SetGlobalConstrainedFilterMask(1<<8|1<<9); // these tracks are written out as global constrained tracks
   esdfilter->SetHybridFilterMaskGlobalConstrainedGlobal((1<<4)); // these normal global tracks will be marked as hybrid

   esdfilter->SetTPCConstrainedFilterMask(1<<11); // these tracks are written out as global constrained tracks

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
   //   esdmuonfilter->SetWriteMuonAOD(writeMuonAOD);
   
   // Enable writing of Dimuon AODs
   //   esdmuonfilter->SetWriteDimuonAOD(writeDimuonAOD);
 
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   mgr->ConnectInput  (esdfilter,  0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (esdfilter,  0, mgr->GetCommonOutputContainer());
   

   //   mgr->ConnectInput  (esdmuonfilter, 0, mgr->GetCommonInputContainer());
   if (useKineFilter) {
      mgr->ConnectInput  (kinefilter,  0, mgr->GetCommonInputContainer());
      mgr->ConnectOutput (kinefilter,  0, mgr->GetCommonOutputContainer());
      AliAnalysisDataContainer *coutputEx = mgr->CreateContainer("cFilterList", TList::Class(),
								   AliAnalysisManager::kOutputContainer,"pyxsec_hists.root");
      mgr->ConnectOutput (kinefilter,  1,coutputEx);
   }   
  if (useCentralityTask) {
       mgr->ConnectInput (ctask, 0, mgr->GetCommonInputContainer());
       mgr->ConnectOutput(ctask, 0, mgr->GetCommonOutputContainer());
   }

   return esdfilter;
 }
 
