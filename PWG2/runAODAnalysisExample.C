void produceAOD()
{
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libXMLIO.so");
    SetupPar("STEERBase");
//    gSystem->Load("STEERBase/libSTEERBase.so");
    SetupPar("ESD");
//    gSystem->Load("ESD/libESD.so");
    SetupPar("AOD");
//    gSystem->Load("AOD/libAOD.so");
    SetupPar("ANALYSIS");
//    gSystem->Load("ANALYSIS/libANALYSIS.so");
    SetupPar("PWG0base");
    // Make the analysis manager
    //
    // Chain from CAF
    gROOT->LoadMacro("CreateESDChain.C");
//     TChain* chain = CreateESDChain("ESD1503X_v1.txt", 10);
    TChain* chain = CreateESDChain("/mnt/data/alice/pbpb_therminator/central.2/", 885);

    AliAODHandler* aodHandler   = new AliAODHandler();
    aodHandler->SetOutputFileName("aod.root");
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetInactiveBranches("FMD CaloCluster");
    
    AliMCEventHandler* mcHandler = new AliMCEventHandler();

    AliAnalysisManager *mgr  = new AliAnalysisManager("esd to aod to histos", "testing aod analysis");
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

void runPWG2()
{
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libXMLIO.so");
    SetupPar("STEERBase");
//    gSystem->Load("STEERBase/libSTEERBase.so");
    SetupPar("ESD");
//    gSystem->Load("ESD/libESD.so");
    SetupPar("AOD");
//    gSystem->Load("AOD/libAOD.so");
    SetupPar("ANALYSIS");
//    gSystem->Load("ANALYSIS/libANALYSIS.so");
    SetupPar("PWG0base");
//    gSystem->Load("PWG0base/libPWG0base.so");
//    SetupPar("TASKFILTER");
//    gSystem->Load("TASKFILTER/libTASKFILTER.so");
//    gSystem->Load("JETAN/libJETAN.so");
  SetupPar("PWG2spectra");
  gROOT->ProcessLine(".L $ALICE_ROOT/PWG2/AliAnalysisTaskProtons.cxx++g");
   //
//    if (gApplication) gApplication->InitializeGraphics();
    // Create the chain
    //
    // Chain from aod.root
     TChain* chain = new TChain("aodTree");
     chain->Add("aod.root");
 //   TChain* chain = CreateESDChain("ESD12001.txt", 2);  
    // Chain from local files
//    gROOT->LoadMacro("CreateLocalChain.C");
//    TChain* chain = CreateLocalChain();  
    // Chain from ALIEN files
//    TGrid::Connect("alien://"); 
//    gROOT->LoadMacro("CreateChain.C");
//    TChain* chain = CreateChain("global.xml");  

    //
    // Make the analysis manager
    //
    AliAODInputHandler *aodHandler = new AliAODInputHandler();

    AliAnalysisManager *mgr  = new AliAnalysisManager("esd to aod to histos", "testing aod analysis");
    mgr->SetInputEventHandler(aodHandler);
    mgr->SetDebugLevel(10);
    AliLog::EnableDebug(kTRUE);
    AliLog::SetGlobalLogLevel(2);


   //
  //____________________________________________//
  // 1st Proton task
  AliAnalysisTaskProtons *taskproton = new AliAnalysisTaskProtons("TaskProtons");
  taskproton->SetType("AOD");
  TFile *f = TFile::Open("$ALICE_ROOT/PWG2/data/PriorProbabilities.root ");
  TF1 *fitElectrons = (TF1 *)f->Get("fitElectrons");
  TF1 *fitMuons = (TF1 *)f->Get("fitMuons");
  TF1 *fitPions = (TF1 *)f->Get("fitPions");
  TF1 *fitKaons = (TF1 *)f->Get("fitKaons");
  TF1 *fitProtons = (TF1 *)f->Get("fitProtons");
  taskproton->SetPriorProbabilityFunctions(fitElectrons,
				      fitMuons,
				      fitPions,
				      fitKaons,
				      fitProtons);
  mgr->AddTask(taskproton);
    // Create containers for input/output
    // Top AOD container
    AliAnalysisDataContainer *cin_aod = mgr->GetCommonInputContainer();

    // Output histogram container
    AliAnalysisDataContainer *cout_hist = mgr->CreateContainer("protonhistos", TList::Class(),AliAnalysisManager::kOutputContainer,"protonhistos.root");

    mgr->ConnectInput(taskproton,0,cin_aod);
    mgr->ConnectOutput(taskproton,0,cout_hist);
   //
    // Run the analysis
    //    
    mgr->InitAnalysis();
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
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
