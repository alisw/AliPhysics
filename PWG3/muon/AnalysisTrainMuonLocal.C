void AnalysisTrainMuonLocal(char* filein = "AliESDs.root", char* fileout = "AliAODs.root")

// Macro to produce a generic AOD starting from an ESD file. 
// The AOD is filled with two tasks: 
// 1- with the first one (AliAnalysisTaskESDfilter), 
//    all the branches of the AOD are filled apart from the muons. 
// 2- with the second task (AliAnalysisTaskESDMuonFilter) 
//    muons tracks are added to the tracks branch 
// There is the possibility to apply cuts on the muon tracks in order 
// to reject muons before filling the AOD
// This macro works locally
// R. Arnaldi 5/5/08

{
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");

    // If analysis is .par based:

    // Common packages
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
    // Analysis-specific packages
    SetupPar("PWG3muon");      
  
    // Input ESD file
    TChain* chain = new TChain("esdTree");  
    chain->Add(filein);
   
    // Define the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "Analysis train");
    
    // ESD input handler
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetInactiveBranches("FMD CaloCluster");
    
    // AOD output handler
    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName(fileout);

    mgr->SetInputEventHandler(esdHandler);
    mgr->SetOutputEventHandler(aodHandler);
    
    // Set of cuts for the ESD filter
    // 
    // standard cut
    AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
    esdTrackCutsL->SetMinNClustersTPC(50);
    esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
    esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCutsL->SetRequireTPCRefit(kTRUE);
    esdTrackCutsL->SetMinNsigmaToVertex(3);
    esdTrackCutsL->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsL->SetAcceptKingDaughters(kFALSE);
    //
    // hard cuts
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
    esdMuonTrackCuts->SetPRange(0.,20.);
    //esdMuonTrackCuts->SetPtRange(0.,0.5);   // examples of kinematic cuts that can be applied
    esdMuonTrackCuts->SetHistogramsOn(kTRUE);  // methods to draw control histos
    esdMuonTrackCuts->DefineHistograms();
    esdMuonTrackCuts->DrawHistograms();
    
    // track filter (to reject tracks not surviving the cuts - refers to all particles apart from muons)
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsH);
    
    // muon track filter  (to reject muon tracks not surviving the cuts)
    AliAnalysisFilter* trackMuonFilter = new AliAnalysisFilter("trackMuonFilter");
    trackMuonFilter->AddCuts(esdMuonTrackCuts);

    // ESD filter task putting standard info in the output generic AOD 
    AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
    //esdfilter->SetTrackFilter(trackFilter);
    esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);
    
    // ESD filter task putting muon info in the output generic AOD 
    AliAnalysisTaskESDMuonFilter *esdmuonfilter = new AliAnalysisTaskESDMuonFilter("ESD Muon Filter");
    esdmuonfilter->SetTrackFilter(trackMuonFilter);
    mgr->AddTask(esdmuonfilter);

    // Containers for input/output
    AliAnalysisDataContainer *cin_esd = mgr->CreateContainer("cESD",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
    // Output AOD container. 
    AliAnalysisDataContainer *cout_aod = mgr->CreateContainer("cAOD", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");
							      
    // Connect containers to tasks slots
    mgr->ConnectInput  (esdfilter,  0, cin_esd  );
    mgr->ConnectOutput (esdfilter,  0, cout_aod );

    mgr->ConnectInput  (esdmuonfilter,  0, cin_esd);
    mgr->ConnectOutput (esdmuonfilter,  0, cout_aod );
 
    //
    // Run the analysis
    //    
    if (mgr->InitAnalysis()) {
        mgr->PrintStatus();
        mgr->StartAnalysis("local",chain);
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
