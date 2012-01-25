void runESDQA()
{
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libXMLIO.so");

    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libJETAN.so");
    gSystem->Load("libPWG4JetTasks.so");
#if 0
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");    
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
#endif

    TGrid::Connect("alien://");
    AliTagAnalysis *analysis = new AliTagAnalysis();
    //TChain *chain = analysis->GetChainFromCollection("pythia_10tev_10files.xml","esdTree");
    TChain *chain = analysis->GetChainFromCollection("pythia_10tev_ph2_10files.xml","esdTree");
    
    // Make the analysis manager

    //
    // Chain from CAF

    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetInactiveBranches("FMD CaloCluster");
    
    AliAnalysisManager *mgr  = new AliAnalysisManager("esd to aod to histos", "testing aod analysis");
    mgr->SetInputEventHandler(esdHandler);

    mgr->SetDebugLevel(10);
    AliLog::EnableDebug(kTRUE);
    AliLog::SetGlobalLogLevel(2);


    // Set of cuts
    // 
    // standard
    /*
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
    esdTrackCutsH->SetMinNClustersTPC(100);
    esdTrackCutsH->SetMaxChi2PerClusterTPC(2.0);
    esdTrackCutsH->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
    esdTrackCutsH->SetRequireTPCRefit(kTRUE);
    esdTrackCutsH->SetMinNsigmaToVertex(2);
    esdTrackCutsH->SetRequireSigmaToVertex(kTRUE);
    esdTrackCutsH->SetAcceptKingDaughters(kFALSE);
    //
    AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
    trackFilter->AddCuts(esdTrackCutsL);
    trackFilter->AddCuts(esdTrackCutsH);
    //
    AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
    esdfilter->SetTrackFilter(trackFilter);
    esdfilter->SetDebugLevel(10);
    mgr->AddTask(esdfilter);
*/

    
    // Take it from the library no need to compile directly
    // Standalone does not need ANALYISalice/AOD/JETAN/PWG4JetsTasks
    //    gProof->Load("AliAnaESDSpectraQA.cxx++g");


    AliAnaESDSpectraQA *pwg4QA = new AliAnaESDSpectraQA("AliAnaSpectraQA");
    mgr->AddTask(pwg4QA);

    // Create containers for input/output
    // Top ESD container
    AliAnalysisDataContainer *cin_esd = mgr->GetCommonInputContainer();

    // Histos
    //AliAnalysisDataContainer *cout_hist = mgr->CreateContainer("qa_hists", TObjArray::Class(),
    AliAnalysisDataContainer *cout_hist = mgr->CreateContainer("qa_hists", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "PWG4QAHists.root");

    mgr->ConnectInput (pwg4QA,  0, cin_esd );
    mgr->ConnectOutput (pwg4QA,  0, cout_hist );
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
