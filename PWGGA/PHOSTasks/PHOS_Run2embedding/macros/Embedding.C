Int_t Nevent=0;
const UInt_t trigger = AliVEvent::kINT7 | AliVEvent::kPHI7;
Int_t EventCounter(TChain *chain);

//==================================
void Embedding(const char* dataset="collection.xml")
{
	const Bool_t embpi0   = kTRUE;
	const Bool_t embeta   = kTRUE;
	const Bool_t embgamma = kTRUE;
    
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  
  //load analysis framework
  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice"); //AliAnalysisTaskSE
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include -I$ALICE_ROOT/PHOS -I$ALICE_PHYSICS/include");

  // A task can be compiled dynamically with AClic
  gROOT->LoadMacro("AliPHOSEmbeddingRun2.cxx++g");
 
  // Connect to alien
  TString token = gSystem->Getenv("GRID_TOKEN") ;
  TGrid::Connect("alien://");
  cout << "Pi0Analysis: processing collection " << dataset << endl;
  
  // Create the chain
  TChain* chain = new TChain("esdTree");

  TGridCollection * collection = dynamic_cast<TGridCollection*>(TAlienCollection::Open(dataset));
  
  TAlienResult* result = collection->GetGridResult("",0 ,0);
  TList* rawFileList = result->GetFileInfoList();
  for (Int_t counter=0 ; counter < rawFileList->GetEntries() ; counter++) {
    TFileInfo * fi =  static_cast<TFileInfo*>(rawFileList->At(counter)) ; 
    const char * rawFile = fi->GetCurrentUrl()->GetUrl() ;  
    printf("Processing %s\n", rawFile) ;
    chain->Add(rawFile);
    printf("Chain: %d entries.\n",chain->GetEntries()); 
  }
  
//  chain->Add("../AliESDs.root") ;
   TFileInfo * fi =  static_cast<TFileInfo*>(rawFileList->At(0));
   const char * fn = fi->GetCurrentUrl()->GetUrl() ;

  char runNum[7]={"246982"}; 
  for(Int_t i=0;i<6;i++)runNum[i]=fn[35+i] ;
  runNum[6]=0 ;
  Int_t iRunNum=atoi(runNum) ;

  Nevent = EventCounter(chain);//get the number of event which we need for generating particles in advance.

  printf("Run number=%d \n",iRunNum) ;   

  //Int_t nSimEvents=chain->GetEntries();
  Int_t nSimEvents=Nevent;
  gSystem->Exec("mkdir tmp") ;
  gSystem->ChangeDirectory("tmp");
  gSystem->Exec("mkdir Sim") ;
  gSystem->Exec("cp ../Sim/* Sim") ;

	if(embpi0){
    gSystem->Exec(Form("../dpgsim.sh --run %d --generator kGeneratorCustom --particle 111 --system Pb-Pb --energy 5020. --nevents %d --mode simrecaod",iRunNum,nSimEvents)) ;
    gSystem->Exec("mv AliAOD.root ../AliAODPi0.root") ; //use geometry from OCDB
    gSystem->Exec("mv sim.log ../sim_pi0.log; mv rec.log ../rec_pi0.log; mv aod.log ../aod_pi0.log") ;
    gSystem->Exec("rm *; rm -rf input; rm -rf GRP;") ; //only Sim directory will be left
	}
	if(embeta){
  	gSystem->Exec(Form("../dpgsim.sh --run %d --generator kGeneratorCustom --particle 221 --system Pb-Pb --energy 5020. --nevents %d --mode simrecaod",iRunNum,nSimEvents)) ;
  	gSystem->Exec("mv AliAOD.root ../AliAODEta.root") ; //use geometry from OCDB
  	gSystem->Exec("mv sim.log ../sim_eta.log; mv rec.log ../rec_eta.log; mv aod.log ../aod_eta.log") ;
    gSystem->Exec("rm *; rm -rf input; rm -rf GRP;") ;//only Sim directory will be left
	}
	if(embgamma){
  	gSystem->Exec(Form("../dpgsim.sh --run %d --generator kGeneratorCustom --particle 22 --system Pb-Pb --energy 5020. --nevents %d --mode simrecaod",iRunNum,nSimEvents)) ;
  	gSystem->Exec("mv AliAOD.root ../AliAODGamma.root") ; //use geometry from OCDB
  	gSystem->Exec("mv sim.log ../sim_gamma.log; mv rec.log ../rec_gamma.log; mv aod.log ../aod_gamma.log") ;
    gSystem->Exec("rm *; rm -rf input; rm -rf GRP;") ;//only Sim directory will be left
	}

  gSystem->ChangeDirectory("../");

  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Pi0EmbeddingManager");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);

  // Output
  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName("AliAODout.root");
  mgr->SetOutputEventHandler(aodHandler);

  
  // Debug level
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(2);
  //mgr->SetDebugLevel(AliLog::kDebug);
  //AliLog::SetGlobalLogLevel(AliLog::kDebug);

  // CDB connection
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  AliTaskCDBconnect *taskCDB = AddTaskCDBconnect("raw://");
  if (!taskCDB){
    printf("Can not connect to OCDB\n") ;
    return;
  }
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");

	const Bool_t useMC = kFALSE;

  // Add physics selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  mgr->RegisterExtraFile("event_stat.root");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(useMC);
//      AliOADBPhysicsSelection * oadbDefaultPbPb = CreateOADBphysicsSelection();      
//      physSelTask->GetPhysicsSelection()->SetCustomOADBObjects(oadbDefaultPbPb,0,0);
  mgr->AddStatisticsTask(AliVEvent::kAny);
      
      
  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  //AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();


//  taskCentrality->SetMCInput();


	//for Run2
	gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
	AliMultSelectionTask *taskMult = AddTaskMultSelection();

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"); 
  AliAnalysisTaskPIDResponse *PIDResponse = AddTaskPIDResponse(kFALSE);

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
  AliEPSelectionTask *taskEP = AddTaskEventplane() ;

  AliPHOSEmbeddingRun2 *task1 = new AliPHOSEmbeddingRun2("Embedding");

  //Add local just simlated input

  TChain* chainAOD1=0;
  TChain* chainAOD2=0;
  TChain* chainAOD3=0;

	if(embpi0){
  	chainAOD1 = new TChain("aodTree");
  	chainAOD1->AddFile("AliAODPi0.root") ;
	}
	if(embeta){
	  chainAOD2 = new TChain("aodTree");
 		chainAOD2->AddFile("AliAODEta.root") ;
	}
	if(embgamma){
  	chainAOD3 = new TChain("aodTree");
		chainAOD3->AddFile("AliAODGamma.root") ;  
	}

  task1->SetSignalChains(chainAOD1,chainAOD2,chainAOD3) ;
 	task1->SetEmbeddedFlag(embpi0,embeta,embgamma);

  task1->SelectCollisionCandidates(trigger);
  const TString path = "alien:///alice/cern.ch/user/d/dsekihat/EmbeddingRun2/PrivateOADB/PHOSCalibrations_vAN20161220.root";
  cout << "path = " << path << endl; 

  task1->SetPrivateOADBPath(path);
  task1->SetSignalCalibration(0.985*1.006);

  mgr->AddTask(task1);

  // Connect input/output
  mgr->ConnectInput(task1 , 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task1, 0, mgr->GetCommonOutputContainer());// daiki changed to 0
  //mgr->ConnectOutput(task1, 1, mgr->GetCommonOutputContainer()); //it was set to 1 by Dmitri originally.
 
  
  //-------------Parameters for filtering
  Bool_t useKineFilter=kFALSE ; 
  Bool_t writeMuonAOD=kFALSE ;
  Bool_t usePhysicsSelection=kFALSE;
  Bool_t enableTPCOnlyAODTracks=kFALSE;
  Bool_t disableCascades=kTRUE ;
  Bool_t disableKinks=kTRUE ;
  Int_t runFlag = 1500; // The first 2 digits are the year, the second
                        //2 digits are used to distinguish sub-periods (if needed)
  Int_t  muonMCMode = 3;
  Bool_t useV0Filter=kFALSE;
  Bool_t muonWithSPDTracklets=kFALSE;
  Bool_t isMuonCaloPass=kFALSE;  
  
    if (disableCascades) task1->DisableCascades();
    if  (disableKinks) task1->DisableKinks();   
    if ( isMuonCaloPass )task1->SetMuonCaloPass(); // this will call a bunch of DisableXXX methods.

    // Muons
//   Bool_t onlyMuon=kTRUE;
//   Bool_t keepAllEvents=kTRUE;
//   Int_t mcMode= useKineFilter ? muonMCMode : 0;
//   AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter",onlyMuon,keepAllEvents,mcMode,muonWithSPDTracklets);
//   mgr->AddTask(esdmuonfilter);
//   if(usePhysicsSelection)task1->SelectCollisionCandidates(AliVEvent::kAny);

//    // Filtering of MC particles (decays conversions etc)
//    // this task has to go AFTER all other filter tasks
//    // since it fills the AODMC array with all
//    // selected MC Particles, only this way we have the 
//    // AODMCparticle information available for following tasks
//    AliAnalysisTaskMCParticleFilter *kinefilter = 0;
//    if (useKineFilter) {
//       kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Kine Filter");
//       mgr->AddTask(kinefilter);
//    }   

   Bool_t enableTPCOnlyAODTracksLocalFlag= enableTPCOnlyAODTracks;
  
  if (!isMuonCaloPass)
  {
    if ((runFlag/100)==15){
     AddTrackCutsLHC15f(task1,enableTPCOnlyAODTracksLocalFlag);
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

     task1->SetV0Filter(v0Filter);
   }  
   
//----------------   
   
   
   
   printf("RunNunm===%d \n",iRunNum) ;
  

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
//  gSystem->Exec("aliroot -b -q AnalyzeDiff.C> analyze3.log 2>&1");


} 
Bool_t AddTrackCutsLHC15f(AliAnalysisTaskESDfilter* esdfilter,Bool_t enableTPCOnlyAODTracksLocalFlag){
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
//____________________________________________________________________
Int_t EventCounter(TChain *chain)
{
  gROOT->LoadMacro("AliPHOSEventCounter.cxx++g");

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("EventCounterManager");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);

	const Bool_t useMC = kFALSE;

  // Add physics selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  //mgr->RegisterExtraFile("event_stat.root");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(useMC);
  //mgr->AddStatisticsTask(AliVEvent::kAny);

  AliPHOSEventCounter *task1 = new AliPHOSEventCounter("PHOSCounter");
  task1->SelectCollisionCandidates(trigger);
  mgr->AddTask(task1);

  
  // Connect input/output
  mgr->ConnectInput(task1 , 0, mgr->GetCommonInputContainer());
  //mgr->ConnectOutput(task1, 0, mgr->GetCommonOutputContainer());// daiki changed to 0

  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }

  Int_t N = task1->GetAnalyzedEvent();
  return N;
}

