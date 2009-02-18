void runPWG2AODCreation(const char *esdsource="ESD1503X_v1.txt", int nfiles=1)
{
  producePWG2AOD(esdsource, nfiles);
}

void producePWG2AOD(const char *esdsource, int nfiles)
{
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libXMLIO.so");
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("PWG0base");
    SetupPar("PWG2AOD");
    // Make the analysis manager
    //
    // Chain from CAF
    gROOT->LoadMacro("CreateESDChain.C");
    TChain* chain = CreateESDChain(esdsource, nfiles);

    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("aod.root");
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetInactiveBranches("FMD CaloCluster");
    
    AliMCEventHandler* mcHandler = new AliMCEventHandler();

    AliAnalysisManager *mgr  = new AliAnalysisManager("esd to pwg2 aod", "testing aod analysis");
    mgr->SetInputEventHandler(esdHandler);
    mgr->SetOutputEventHandler(aodHandler);
    mgr->SetMCtruthEventHandler(mcHandler);
    mgr->SetDebugLevel(10);
    AliLog::EnableDebug(kTRUE);
    AliLog::SetGlobalLogLevel(2);


    // Set of cuts
    // 
    // standard
    AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
    esdTrackCutsL->SetMinNClustersTPC(50);
    esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
    esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCutsL->SetRequireTPCRefit(kTRUE);
    esdTrackCutsL->SetMinNsigmaToVertex(3);
    esdTrackCutsL->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsL->SetAcceptKingDaughters(kFALSE);
    //
    // hard
    AliESDtrackCuts* esdTrackCutsH = new AliESDtrackCuts("AliESDtrackCuts", "Hard");
    esdTrackCutsH->SetMinNClustersTPC(95);
    esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
    esdTrackCutsH->SetMaxChi2PerClusterITS(3.0);
    esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCutsH->SetRequireTPCRefit(kTRUE);
    esdTrackCutsH->SetRequireITSRefit(kTRUE);
    esdTrackCutsH->SetMinNsigmaToVertex(2);
    esdTrackCutsH->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsH->SetAcceptKingDaughters(kFALSE);
    esdTrackCutsH->SetPtRange(0.1,1.5);
    esdTrackCutsH->SetEtaRange(-1.0,1.0);
    //
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsL);
    trackFilter->AddCuts(esdTrackCutsH);
    //
    AliAnalysisTaskPWG2ESDfilter *esdfilter = new AliAnalysisTaskPWG2ESDfilter("PWG2 ESD Filter");
    esdfilter->SetTrackFilter(trackFilter);
    esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);

    // Create containers for input/output
    // Top ESD container
    AliAnalysisDataContainer *cin_esd = mgr->GetCommonInputContainer();

    // Output AOD container
    AliAnalysisDataContainer *cout_aod = mgr->GetCommonOutputContainer();

    mgr->ConnectInput  (esdfilter,  0, cin_esd  );
    mgr->ConnectOutput (esdfilter,  0, cout_aod );
    //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
    delete mgr;
}

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
