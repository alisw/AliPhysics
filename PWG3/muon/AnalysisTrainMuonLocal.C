void AnalysisTrainMuonLocal(char* filein = "AliESDs.root", 
                            char* fileout = "AliAODs.root", 
			    Int_t nev=123456789)

// Macro to produce a generic AOD starting from an ESD file. 
// The AOD is filled with two tasks: 
// 1- with the first one (AliAnalysisTaskESDfilter), 
//    all the branches of the AOD are filled apart from the muons. 
// 2- with the second task (AliAnalysisTaskESDMuonFilter) 
//    muons tracks are added to the tracks branch 
// 3- with a third task (AliAnalysisTaskTagCreator) 
//    aod tags are created 
// There is the possibility to apply cuts on the tracks and muon tracks in 
// order to reject them before filling the AOD
// This macro works locally

{
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}/libPWG3muon.so");  // for aliroot based analysis

    // Par files for a par based analysis 
    // SetupPar("STEERBase");
    // SetupPar("ESD");
    // SetupPar("AOD");
    // SetupPar("ANALYSIS");
    // SetupPar("ANALYSISalice");
    // SetupPar("PWG3muon");   

    // Creating ESD Tags on the fly
    // The recreation of the ESD tag file is only needed in order to copy the infos on
    // run/LHC parameters into the AOD tag file. If the ESD tag file is not recreated, the
    // run/LHC info in the AOD tag file will be empty.
    AliESDTagCreator *t = new AliESDTagCreator();
    t->SetStorage(0);
    t->ReadLocalCollection("."); 

    // Input ESD file
    TChain* chain = new TChain("esdTree");  
    chain->Add(filein);
   
    // Define the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Analysis train");
    
    // ESD input handler
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetReadTags();
    
    // AOD output handler
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName(fileout);

    mgr->SetInputEventHandler(esdHandler);
    mgr->SetOutputEventHandler(aodHandler);
    
    // Set of cuts for the ESD filters. 
    // Only tracks surviving the cuts will be copied into the AOD
    // 
    // standard cut for non-muon tracks
    AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
    esdTrackCutsL->SetMinNClustersTPC(50);
    esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
    esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCutsL->SetRequireTPCRefit(kTRUE);
    esdTrackCutsL->SetMinNsigmaToVertex(3);
    esdTrackCutsL->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsL->SetAcceptKingDaughters(kFALSE);
    //
    // hard cuts for non-muon tracks
    AliESDtrackCuts* esdTrackCutsH = new AliESDtrackCuts("AliESDtrackCuts", "Hard");
    esdTrackCutsH->SetMinNClustersTPC(100);
    esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
    esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCutsH->SetRequireTPCRefit(kTRUE);
    esdTrackCutsH->SetMinNsigmaToVertex(2);
    esdTrackCutsH->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsH->SetAcceptKingDaughters(kFALSE);
    esdTrackCutsH->SetPRange(0.,2.);
    //
    //  muon cuts
    AliESDMuonTrackCuts* esdMuonTrackCuts = new AliESDMuonTrackCuts("AliESDMuonTrackCuts", "test");
    esdMuonTrackCuts->SetPRange(0.,12.);
    esdMuonTrackCuts->SetPtRange(0.,2.);   // examples of kinematic cuts that can be applied
    esdMuonTrackCuts->SetHistogramsOn(kTRUE);  // methods to draw control histos
    esdMuonTrackCuts->DefineHistograms();
    esdMuonTrackCuts->DrawHistograms();
    
    // track filter (to reject tracks not surviving the previously defined cuts - 
    // refers to all particles apart from muons)
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsH);
    
    // muon track filter  (to reject muon tracks not surviving the previously defined cuts)
    AliAnalysisFilter* trackMuonFilter = new AliAnalysisFilter("trackMuonFilter");
    trackMuonFilter->AddCuts(esdMuonTrackCuts);

    // ESD filter task to fill standard info in the output generic AOD 
    AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
    //esdfilter->SetTrackFilter(trackFilter); //uncomment to apply cuts on the tracks
    esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);
    
    // ESD filter task filling muon info in the output generic AOD 
    AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
    //esdmuonfilter->SetTrackFilter(trackMuonFilter); //uncomment to apply cuts on the muon tracks
    mgr->AddTask(esdmuonfilter);

    // Tag Creator
    AliAnalysisTaskTagCreator* tagTask = new AliAnalysisTaskTagCreator("AOD Tag Creator");
    mgr->AddTask(tagTask);
    
    // Input container
    AliAnalysisDataContainer *cin_esd = mgr->CreateContainer("cESD",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
    // Output AOD container. 
    AliAnalysisDataContainer *cout_aod = mgr->CreateContainer("cAOD", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");
    // Tag container							      
    AliAnalysisDataContainer *cout_tags = mgr->CreateContainer("cTag",TTree::Class(), 
                                                               AliAnalysisManager::kOutputContainer, "AOD.tag.root");
    
    // Connect containers to tasks slots
    mgr->ConnectInput  (esdfilter,  0, cin_esd  );
    mgr->ConnectOutput (esdfilter,  0, cout_aod );

    mgr->ConnectInput  (esdmuonfilter,  0, cin_esd);
    mgr->ConnectOutput (esdmuonfilter,  0, cout_aod );
 
    mgr->ConnectInput  (tagTask, 0, cin_esd);
    mgr->ConnectOutput (tagTask, 1, cout_tags);

    //
    // Run the analysis
    //    
    if (mgr->InitAnalysis()) {
        mgr->PrintStatus();
        mgr->StartAnalysis("local",chain,nev);
    }   
}

//______________________________________________________________________________
void SetupPar(char* pararchivename)
{
    if (pararchivename) {
	char processline[1024];
	sprintf(processline,".! tar xvzf %s.par",pararchivename);
	gROOT->ProcessLine(processline);
	TString ocwd = gSystem->WorkingDirectory();
	gSystem->ChangeDirectory(pararchivename);
	
	// check for BUILD.sh and execute
	if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
	    printf("*******************************\n");
	    printf("*** Building PAR archive    ***\n");
	    printf("*******************************\n");
	    
	    if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
		Error("runProcess","Cannot Build the PAR Archive! - Abort!");
		return -1;
	    }
	}
	// check for SETUP.C and execute
	if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
	    printf("*******************************\n");
	    printf("*** Setup PAR archive       ***\n");
	    printf("*******************************\n");
	    gROOT->Macro("PROOF-INF/SETUP.C");
	}
	
	gSystem->ChangeDirectory(ocwd.Data());
   printf("Current dir: %s\n", ocwd.Data());
    } 
}
