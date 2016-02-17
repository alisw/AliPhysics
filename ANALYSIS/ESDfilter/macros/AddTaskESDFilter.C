
Bool_t AddTrackCutsLHC10bcde(AliAnalysisTaskESDfilter* esdFilter);
Bool_t AddTrackCutsLHC10h(AliAnalysisTaskESDfilter* esdFilter);
Bool_t AddTrackCutsLHC11h(AliAnalysisTaskESDfilter* esdFilter);
Bool_t AddTrackCutsLHC15f(AliAnalysisTaskESDfilter* esdFilter);
Bool_t enableTPCOnlyAODTracksLocalFlag=kFALSE;

AliAnalysisTaskESDfilter *AddTaskESDFilter(Bool_t useKineFilter=kTRUE, 
                                           Bool_t writeMuonAOD=kFALSE,
                                           Bool_t writeDimuonAOD=kFALSE, /*obsolete*/
                                           Bool_t usePhysicsSelection=kFALSE,
                                           Bool_t useCentralityTask=kFALSE, /*obsolete*/
                                           Bool_t enableTPCOnlyAODTracks=kFALSE,
                                           Bool_t disableCascades=kFALSE,
                                           Bool_t disableKinks=kFALSE, 
                                           Int_t runFlag = 1500, // The first 2 digits are the year, the second
                                                                 //2 digits are used to distinguish sub-periods (if needed)
                                           Int_t  muonMCMode = 3  ,
                                           Bool_t useV0Filter=kTRUE,
                                           Bool_t muonWithSPDTracklets=kTRUE,
                                           Bool_t isMuonCaloPass=kFALSE)
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
   if (disableCascades) esdfilter->DisableCascades();
   if  (disableKinks) esdfilter->DisableKinks();
  
   if ( isMuonCaloPass )
   {
      esdfilter->SetMuonCaloPass(); // this will call a bunch of DisableXXX methods.
   }
  
   mgr->AddTask(esdfilter);
  
   // Muons
   Bool_t onlyMuon=kTRUE;
   Bool_t keepAllEvents=kTRUE;
   Int_t mcMode= useKineFilter ? muonMCMode : 0;
   AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter",onlyMuon,keepAllEvents,mcMode,muonWithSPDTracklets);
   mgr->AddTask(esdmuonfilter);
   if(usePhysicsSelection)
   {
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

   enableTPCOnlyAODTracksLocalFlag = enableTPCOnlyAODTracks;
  
  if (!isMuonCaloPass)
  {
   if((runFlag/100)==10){
     if((runFlag%100)==0) AddTrackCutsLHC10bcde(esdfilter);
     else AddTrackCutsLHC10h(esdfilter);
   }
   else if ((runFlag/100)==11){
     // default 11h
     AddTrackCutsLHC11h(esdfilter);
   }
   else if ((runFlag/100)==15){
     AddTrackCutsLHC15f(esdfilter);
   }
   else {
     std::cout << "ERROR: illegal runFlag value: " << runFlag  << std::endl;
     return NULL;     
   }
  }
  
   // Filter with cuts on V0s
   if (useV0Filter && !isMuonCaloPass) {
     AliESDv0Cuts*   esdV0Cuts = new AliESDv0Cuts("Standard V0 Cuts pp", "ESD V0 Cuts");
     esdV0Cuts->SetMinRadius(0.2);
     esdV0Cuts->SetMaxRadius(200);
     esdV0Cuts->SetMinDcaPosToVertex(0.05);
     esdV0Cuts->SetMinDcaNegToVertex(0.05);
     esdV0Cuts->SetMaxDcaV0Daughters(1.5);
     esdV0Cuts->SetMinCosinePointingAngle(0.99);
     AliAnalysisFilter* v0Filter = new AliAnalysisFilter("v0Filter");
     v0Filter->AddCuts(esdV0Cuts);

     esdfilter->SetV0Filter(v0Filter);
   }  

   // Enable writing of Muon AODs
   esdmuonfilter->SetWriteMuonAOD(writeMuonAOD);
   
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
 
Bool_t AddTrackCutsLHC10h(AliAnalysisTaskESDfilter* esdfilter){

  Printf("%s%d: Creating Track Cuts for LHC10h",(char*)__FILE__,__LINE__);

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
  electronID->SetTPCnSigmaCut(AliPID::kElectron, 3.5);


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

  
  // tighter cuts on primary particles for high pT tracks
  // take the standard cuts, which include already 
  // ITSrefit and use only primaries...
  
  // ITS cuts for new jet analysis 
  //  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");
  //  AliESDtrackCuts* esdTrackCutsHG0 = CreateTrackCutsPWGJE(10001006);

  AliESDtrackCuts *jetCuts1006 = new AliESDtrackCuts("AliESDtrackCuts"); 

  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  jetCuts1006->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
  jetCuts1006->SetMinNClustersTPC(70);
  jetCuts1006->SetMaxChi2PerClusterTPC(4);
  jetCuts1006->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
  jetCuts1006->SetAcceptKinkDaughters(kFALSE);
  jetCuts1006->SetRequireTPCRefit(kTRUE);
  jetCuts1006->SetMaxFractionSharedTPCClusters(0.4);
  // ITS
  jetCuts1006->SetRequireITSRefit(kTRUE);
  //accept secondaries
  jetCuts1006->SetMaxDCAToVertexXY(2.4);
  jetCuts1006->SetMaxDCAToVertexZ(3.2);
  jetCuts1006->SetDCAToVertex2D(kTRUE);
  //reject fakes
  jetCuts1006->SetMaxChi2PerClusterITS(36);
  jetCuts1006->SetMaxChi2TPCConstrainedGlobal(36);

  jetCuts1006->SetRequireSigmaToVertex(kFALSE);

  jetCuts1006->SetEtaRange(-0.9,0.9);
  jetCuts1006->SetPtRange(0.15, 1E+15.);

  AliESDtrackCuts* esdTrackCutsHG0 = jetCuts1006->Clone("JetCuts10001006");
  esdTrackCutsHG0->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

  // the complement to the one with SPD requirement: tracks with ITS refit but no SPD hit
  //  AliESDtrackCuts* esdTrackCutsHG1 = CreateTrackCutsPWGJE(10011006);
  AliESDtrackCuts* esdTrackCutsHG1 = jetCuts1006->Clone("JetCuts10011006");
  esdTrackCutsHG1->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);

  AliESDtrackCuts* esdTrackCutsHG2 = jetCuts1006->Clone("JetCuts10021006");
  esdTrackCutsHG2->SetMaxChi2PerClusterITS(1E10);

  // all complementary hybrid tracks: no SPD requirement, no ITS refit requirement
  //  AliESDtrackCuts* esdTrackCutsGCOnly = CreateTrackCutsPWGJE(10041006);
  AliESDtrackCuts* esdTrackCutsGCOnly = jetCuts1006->Clone("JetCuts10041006");
  esdTrackCutsGCOnly->SetRequireITSRefit(kFALSE);

  // standard cuts also used in R_AA analysis
  //   "Global track RAA analysis QM2011 + Chi2ITS<36";
  //  AliESDtrackCuts* esdTrackCutsH2 = CreateTrackCutsPWGJE(1000);
  AliESDtrackCuts* esdTrackCutsRAA = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  esdTrackCutsRAA->SetMinNCrossedRowsTPC(120);
  esdTrackCutsRAA->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCutsRAA->SetMaxChi2PerClusterITS(36);
  esdTrackCutsRAA->SetMaxFractionSharedTPCClusters(0.4);
  esdTrackCutsRAA->SetMaxChi2TPCConstrainedGlobal(36);

  esdTrackCutsRAA->SetEtaRange(-0.9,0.9);
  esdTrackCutsRAA->SetPtRange(0.15, 1e10);

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
  trackFilter->AddCuts(esdTrackCutsH);
  // 32 1<<5
  trackFilter->AddCuts(esdTrackCutsH2);
  // 64 1<<6
  trackFilter->AddCuts(esdTrackCutsH3);
  // 128 1<<7
  trackFilter->AddCuts(esdTrackCutsTPCCOnly); // add TPC only track cuts for TPC constrained tracks
  if(enableTPCOnlyAODTracksLocalFlag)esdfilter->SetTPCOnlyFilterMask(128);
  // 256 1<<8
  trackFilter->AddCuts(esdTrackCutsHG0);
  esdfilter->SetHybridFilterMaskGlobalConstrainedGlobal((1<<8)); // these normal global tracks will be marked as hybrid
  // 512 1<<9                         
  trackFilter->AddCuts(esdTrackCutsGCOnly);                      // all complementary hybrids (no SPD req && no ITS refit req && !(1<<8))
  // 1024 1<<10                        
  trackFilter->AddCuts(esdTrackCutsHG1);                         // complementary tracks with ITSrefit & SPD none
  esdfilter->SetGlobalConstrainedFilterMask(1<<9|1<<10);         // these tracks are written out as global constrained tracks
  esdfilter->SetWriteHybridGlobalConstrainedOnly(kTRUE);         // write only the complementary tracks
  // 2048 1<<11
  trackFilter->AddCuts(esdTrackCutsRAA);
  // 4096 1<<12 
  AliESDtrackCuts* esdTrackCutsHG1_tmp = new AliESDtrackCuts(*esdTrackCutsHG1); // avoid double delete
  trackFilter->AddCuts(esdTrackCutsHG1_tmp);
  // 8192 1<<13
  trackFilter->AddCuts(esdTrackCutsHG2);

  esdfilter->SetTrackFilter(trackFilter);
  return kTRUE;
}


Bool_t AddTrackCutsLHC11h(AliAnalysisTaskESDfilter* esdfilter){


  Printf("%s%d: Creating Track Cuts LHC11h",(char*)__FILE__,__LINE__);

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
   electronID->SetTPCnSigmaCut(AliPID::kElectron, 3.5);

   // standard cuts with very loose DCA
   AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
   esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
   esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
   esdTrackCutsH->SetDCAToVertex2D(kTRUE);

   // standard cuts with tight DCA cut
   AliESDtrackCuts* esdTrackCutsH2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

   // standard cuts with tight DCA but with requiring the first SDD cluster instead of an SPD cluster
   // tracks selected by this cut are exclusive to those selected by the previous cut
   AliESDtrackCuts* esdTrackCutsH3 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(); 
   esdTrackCutsH3->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
   esdTrackCutsH3->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
 
   // TPC only tracks: Optionally enable the writing of TPConly information
   // constrained to SPD vertex in the filter below
   AliESDtrackCuts* esdTrackCutsTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
   // The following line is needed for 2010 PbPb reprocessing and pp, but not for 2011 PbPb
   //esdTrackCutsTPCOnly->SetMinNClustersTPC(70);

   // Extra cuts for hybrids
   // first the global tracks we want to take
   AliESDtrackCuts* esdTrackCutsHTG = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
   esdTrackCutsHTG->SetName("Global Hybrid tracks, loose DCA");
   esdTrackCutsHTG->SetMaxDCAToVertexXY(2.4);
   esdTrackCutsHTG->SetMaxDCAToVertexZ(3.2);
   esdTrackCutsHTG->SetDCAToVertex2D(kTRUE);
   esdTrackCutsHTG->SetMaxChi2TPCConstrainedGlobal(36);
   esdTrackCutsHTG->SetMaxFractionSharedTPCClusters(0.4);
   

   // Than the complementary tracks which will be stored as global
   // constraint, complement is done in the ESDFilter task
   AliESDtrackCuts* esdTrackCutsHTGC = new AliESDtrackCuts(*esdTrackCutsHTG);
   esdTrackCutsHTGC->SetName("Global Constraint Hybrid tracks, loose DCA no it requirement");
   esdTrackCutsHTGC->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
   esdTrackCutsHTGC->SetRequireITSRefit(kTRUE);

   // standard cuts with tight DCA cut, using cluster cut instead of crossed rows (a la 2010 default)
   AliESDtrackCuts* esdTrackCutsH2Cluster = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 0);

   // Compose the filter
   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   // 1, 1<<0
   trackFilter->AddCuts(esdTrackCutsL);
   // 2, 1<<1
   trackFilter->AddCuts(esdTrackCutsITSsa);
   // 4, 1<<2
   trackFilter->AddCuts(itsStrong);
   itsStrong->SetFilterMask(1);        // AND with Standard track cuts 
   // 8, 1<<3
   trackFilter->AddCuts(electronID);
   electronID->SetFilterMask(4);       // AND with Pixel Cuts
   // 16, 1<<4
   trackFilter->AddCuts(esdTrackCutsH);
   // 32, 1<<5
   trackFilter->AddCuts(esdTrackCutsH2);
   // 64, 1<<6
   trackFilter->AddCuts(esdTrackCutsH3);
   // 128 , 1 << 7
   trackFilter->AddCuts(esdTrackCutsTPCOnly);
   if(enableTPCOnlyAODTracksLocalFlag)esdfilter->SetTPCOnlyFilterMask(128);
   // 256, 1 << 8 Global Hybrids
   trackFilter->AddCuts(esdTrackCutsHTG);
   esdfilter->SetHybridFilterMaskGlobalConstrainedGlobal((1<<8)); // these normal global tracks will be marked as hybrid    
   // 512, 1<< 9 GlobalConstraint Hybrids
   trackFilter->AddCuts(esdTrackCutsHTGC);
   esdfilter->SetGlobalConstrainedFilterMask(1<<9); // these tracks are written out as global constrained tracks 
   esdfilter->SetWriteHybridGlobalConstrainedOnly(kTRUE); // write only the complement
   // 1024, 1<< 10 // tight DCA cuts
   trackFilter->AddCuts(esdTrackCutsH2Cluster);
   esdfilter->SetTrackFilter(trackFilter);

   return kTRUE;

}

Bool_t AddTrackCutsLHC10bcde(AliAnalysisTaskESDfilter* esdfilter){


  Printf("%s%d: Creating Track Cuts LHC10bcde",(char*)__FILE__,__LINE__);

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
   electronID->SetTPCnSigmaCut(AliPID::kElectron, 3.5);

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
   // The following line is needed for 2010 PbPb reprocessing and pp, but not for 2011 PbPb
   esdTrackCutsTPCOnly->SetMinNClustersTPC(70);

   // Extra cuts for hybrids
   // first the global tracks we want to take
   // take the HTGs from 10h

   //  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");
   //  AliESDtrackCuts* esdTrackCutsHG0 = CreateTrackCutsPWGJE(10001006);

   AliESDtrackCuts *jetCuts1006 = new AliESDtrackCuts("AliESDtrackCuts"); 

   TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
   jetCuts1006->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
   jetCuts1006->SetMinNClustersTPC(70);
   jetCuts1006->SetMaxChi2PerClusterTPC(4);
   jetCuts1006->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
   jetCuts1006->SetAcceptKinkDaughters(kFALSE);
   jetCuts1006->SetRequireTPCRefit(kTRUE);
   jetCuts1006->SetMaxFractionSharedTPCClusters(0.4);
   // ITS
   jetCuts1006->SetRequireITSRefit(kTRUE);
   //accept secondaries
   jetCuts1006->SetMaxDCAToVertexXY(2.4);
   jetCuts1006->SetMaxDCAToVertexZ(3.2);
   jetCuts1006->SetDCAToVertex2D(kTRUE);
   //reject fakes
   jetCuts1006->SetMaxChi2PerClusterITS(36);
   jetCuts1006->SetMaxChi2TPCConstrainedGlobal(36);

   jetCuts1006->SetRequireSigmaToVertex(kFALSE);
   
   jetCuts1006->SetEtaRange(-0.9,0.9);
   jetCuts1006->SetPtRange(0.15, 1E+15.);
   
   AliESDtrackCuts* esdTrackCutsHTG = jetCuts1006->Clone("JetCuts10001006");
   esdTrackCutsHTG->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

   // Than the complementary tracks which will be stored as global
   // constraint, complement is done in the ESDFilter task
   // HGC from 10h

   AliESDtrackCuts* esdTrackCutsHTGC = jetCuts1006->Clone("JetCuts10041006");
   esdTrackCutsHTGC->SetRequireITSRefit(kFALSE);

   // standard cuts with tight DCA cut, using cluster cut instead of crossed rows (a la 2010 default)
   AliESDtrackCuts* esdTrackCutsH2Cluster = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, 0);

   // Compose the filter
   AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
   // 1, 1<<0
   trackFilter->AddCuts(esdTrackCutsL);
   // 2, 1<<1
   trackFilter->AddCuts(esdTrackCutsITSsa);
   // 4, 1<<2
   trackFilter->AddCuts(itsStrong);
   itsStrong->SetFilterMask(1);        // AND with Standard track cuts 
   // 8, 1<<3
   trackFilter->AddCuts(electronID);
   electronID->SetFilterMask(4);       // AND with Pixel Cuts
   // 16, 1<<4
   trackFilter->AddCuts(esdTrackCutsH);
   // 32, 1<<5
   trackFilter->AddCuts(esdTrackCutsH2);
   // 64, 1<<6
   trackFilter->AddCuts(esdTrackCutsH3);
   // 128 , 1 << 7
   trackFilter->AddCuts(esdTrackCutsTPCOnly);
   if(enableTPCOnlyAODTracksLocalFlag)esdfilter->SetTPCOnlyFilterMask(128);
   // 256, 1 << 8 Global Hybrids
   trackFilter->AddCuts(esdTrackCutsHTG);
   esdfilter->SetHybridFilterMaskGlobalConstrainedGlobal((1<<8)); // these normal global tracks will be marked as hybrid    
   // 512, 1<< 9 GlobalConstraint Hybrids
   trackFilter->AddCuts(esdTrackCutsHTGC);
   esdfilter->SetGlobalConstrainedFilterMask(1<<9); // these tracks are written out as global constrained tracks 
   esdfilter->SetWriteHybridGlobalConstrainedOnly(kTRUE); // write only the complement
   // 1024, 1<< 10 // tight DCA cuts
   trackFilter->AddCuts(esdTrackCutsH2Cluster);
   esdfilter->SetTrackFilter(trackFilter);

   return kTRUE;

}

Bool_t AddTrackCutsLHC15f(AliAnalysisTaskESDfilter* esdfilter){
  //
  // filter cuts for RunII pp in 2015
  // basically a duplication of 11h, but with stricter cluster requirement
  //
  Printf("%s%d: Creating Track Cuts for LHC15f",(char*)__FILE__,__LINE__);
  //
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
  electronID->SetTPCnSigmaCut(AliPID::kElectron, 3.5);
  
  // standard cuts with very loose DCA
  AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
  esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
  esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
  esdTrackCutsH->SetDCAToVertex2D(kTRUE);

  // standard cuts with tight DCA cut
  AliESDtrackCuts* esdTrackCutsH2 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  
  // standard cuts with tight DCA but with requiring the first SDD cluster instead of an SPD cluster
  // tracks selected by this cut are exclusive to those selected by the previous cut
  AliESDtrackCuts* esdTrackCutsH3 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(); 
  esdTrackCutsH3->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  esdTrackCutsH3->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);
 
  // TPC only tracks: Optionally enable the writing of TPConly information
  // constrained to SPD vertex in the filter below
  AliESDtrackCuts* esdTrackCutsTPCOnly = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  // The following line is needed for 2010 PbPb reprocessing and pp, but not for 2011 PbPb
  //esdTrackCutsTPCOnly->SetMinNClustersTPC(70);
  
  // Extra cuts for hybrids
  // first the global tracks we want to take
  AliESDtrackCuts* esdTrackCutsHTG = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
  esdTrackCutsHTG->SetName("Global Hybrid tracks, loose DCA");
  esdTrackCutsHTG->SetMaxDCAToVertexXY(2.4);
  esdTrackCutsHTG->SetMaxDCAToVertexZ(3.2);
  esdTrackCutsHTG->SetDCAToVertex2D(kTRUE);
  esdTrackCutsHTG->SetMaxChi2TPCConstrainedGlobal(36);
  esdTrackCutsHTG->SetMaxFractionSharedTPCClusters(0.4);
  
  // Than the complementary tracks which will be stored as global
  // constraint, complement is done in the ESDFilter task
  AliESDtrackCuts* esdTrackCutsHTGC = new AliESDtrackCuts(*esdTrackCutsHTG);
  esdTrackCutsHTGC->SetName("Global Constraint Hybrid tracks, loose DCA no it requirement");
  esdTrackCutsHTGC->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
  esdTrackCutsHTGC->SetRequireITSRefit(kTRUE);

  // standard cuts with tight DCA cut, using cluster cut instead of crossed rows (a la 2010 default)
  AliESDtrackCuts* esdTrackCutsH2Cluster = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 0);
  esdTrackCutsH2Cluster->SetMinNClustersTPC(70); // gain in 2015 is higher than in 2011

  // duplication of 1<<5 = 32 and 1<<6 = 64 with looser requirement 
  // on CrossedRows and CrossedRowsOverFindable in order to go to forward eta (To be used with care!)
  AliESDtrackCuts* esdTrackCutsH2Forward = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  esdTrackCutsH2Forward->SetMinNCrossedRowsTPC(50);
  esdTrackCutsH2Forward->SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);

  AliESDtrackCuts* esdTrackCutsH3Forward = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
  esdTrackCutsH3Forward->SetMinNCrossedRowsTPC(50);
  esdTrackCutsH3Forward->SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);
  esdTrackCutsH3Forward->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
  esdTrackCutsH3Forward->SetClusterRequirementITS(AliESDtrackCuts::kSDD, AliESDtrackCuts::kFirst);


  // Compose the filter
  AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
  // 1, 1<<0
  trackFilter->AddCuts(esdTrackCutsL);
  // 2, 1<<1
  trackFilter->AddCuts(esdTrackCutsITSsa);
  // 4, 1<<2
  trackFilter->AddCuts(itsStrong);
  itsStrong->SetFilterMask(1);        // AND with Standard track cuts 
  // 8, 1<<3
  trackFilter->AddCuts(electronID);
  electronID->SetFilterMask(4);       // AND with Pixel Cuts
   // 16, 1<<4
  trackFilter->AddCuts(esdTrackCutsH);
  // 32, 1<<5
  trackFilter->AddCuts(esdTrackCutsH2);
  // 64, 1<<6
  trackFilter->AddCuts(esdTrackCutsH3);
  // 128 , 1 << 7
  trackFilter->AddCuts(esdTrackCutsTPCOnly);
  if(enableTPCOnlyAODTracksLocalFlag)esdfilter->SetTPCOnlyFilterMask(128);
  // 256, 1 << 8 Global Hybrids
  trackFilter->AddCuts(esdTrackCutsHTG);
  esdfilter->SetHybridFilterMaskGlobalConstrainedGlobal((1<<8)); // these normal global tracks will be marked as hybrid    
  // 512, 1<< 9 GlobalConstraint Hybrids
  trackFilter->AddCuts(esdTrackCutsHTGC);
  esdfilter->SetGlobalConstrainedFilterMask(1<<9); // these tracks are written out as global constrained tracks 
  esdfilter->SetWriteHybridGlobalConstrainedOnly(kTRUE); // write only the complement
  // 1024, 1<< 10 // tight DCA cuts
  trackFilter->AddCuts(esdTrackCutsH2Cluster);
  // 2048, 1<<11 // duplication of 1<<5 with looser CrossedRows requirements for forward eta
  trackFilter->AddCuts(esdTrackCutsH2Forward);
  // 4096, 1<<12 // duplication of 1<<6 with looser CrossedRows requirements for forward eta
  trackFilter->AddCuts(esdTrackCutsH3Forward);

  esdfilter->SetTrackFilter(trackFilter);

  return kTRUE;

}
