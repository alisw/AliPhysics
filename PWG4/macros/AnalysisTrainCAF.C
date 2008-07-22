//______________________________________________________________________________
void AnalysisTrainCAF()
{
// Example of running analysis train in CAF. To run in debug mode:
//  - export ROOTSYS=debug  on your local client
//  - un-comment gProof->ClearPackages()
//  - un-comme lines with debugging info

// WHY AOD is not a exchange container when running from ESD->AOD

    Bool_t debug         = kTRUE;
    Bool_t useMC         = kTRUE;
    Bool_t readTR        = kFALSE;
    
    Int_t iAODanalysis   = 0;
    Int_t iAODhandler    = 1;
    Int_t iESDfilter     = 1;  // Only active if iAODanalysis=0
    Int_t iJETAN         = 1;
    Int_t iJETANMC       = 1;
    Int_t iPWG4UE      = 1;

    if (iAODanalysis) {
       useMC = kFALSE;
       readTR = kFALSE;
       iESDfilter = 0;
       iMUONfilter = 0;
    }    
    if (iJETAN) iESDfilter=1;
    if (iESDfilter) iAODhandler=1;

    // Dataset from CAF
    TString dataset = "/PWG4/arian/jetjet15-50";
    printf("==================================================================\n");
    printf("===========    RUNNING ANALYSIS TRAIN IN CAF MODE    =============\n");
    printf("==================================================================\n");
    if (iAODanalysis) printf("=  AOD analysis on dataset: %s\n", dataset.Data());
    else              printf("=  ESD analysis on dataset: %s\n", dataset.Data());
    if (iESDfilter)   printf("=  ESD filter                                                    =\n");
    if (iJETAN)       printf("=  Jet analysis                                                  =\n");
    if (iJETANMC)       printf("=  Jet analysis from Kinematics                                  =\n");
    if (iPWG4UE)  printf("=  PWG4 UE                                                   =\n");
    printf("==================================================================\n");
    if (useMC) printf(":: use MC    TRUE\n");
    else       printf(":: use MC    FALSE\n");
    if (readTR) printf(":: read TR   TRUE\n");
    else        printf(":: read TR   FALSE\n");
    if (debug) printf(":: debugging TRUE\n");
    else       printf(":: debugging FALSE\n");
    
    // Load common libraries
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");


    // Reset user processes if CAF if not responding anymore
    //TProof::Reset("lxb6046"); 
    // One may enable a different ROOT version on CAF

    const char* proofNode = "lxb6046";
    //    const char* proofNode = "kleinb@localhost";

    TProof::Mgr(proofNode)->ShowROOTVersions();
    //    TProof::Mgr(proofNode)->SetROOTVersion("vHEAD-r24503_dbg");

    // Connect to proof
    TProof::Open(proofNode); // may be username@lxb6046 if user not the same as on local

    // Clear packages if changing ROOT version on CAF or local
    //    gProof->ClearPackages();
    // Enable proof debugging if needed
    //    gProof->SetLogLevel(5);
    // To debug the train in PROOF mode, type in a root session:
    // root[0] TProof::Mgr("lxb6064")->GetSessionLogs()->Display("*",0,10000);
    // Common packages
    // --- Enable the STEERBase Package
    gProof->UploadPackage("pars/STEERBase.par");
    gProof->EnablePackage("STEERBase");
    // --- Enable the ESD Package
    gProof->UploadPackage("pars/ESD.par");
    gProof->EnablePackage("ESD");
    // --- Enable the AOD Package
    gProof->UploadPackage("pars/AOD.par");
    gProof->EnablePackage("AOD");
    // --- Enable the ANALYSIS Package
    gProof->UploadPackage("pars/ANALYSIS.par");
    gProof->EnablePackage("ANALYSIS");
    // --- Enable the ANALYSISalice Package
    gProof->UploadPackage("pars/ANALYSISalice.par");
    gProof->EnablePackage("ANALYSISalice");

    AliPDG::AddParticlesToPdgDataBase();

    // --- Enable the JETAN Package
    if (iJETAN||iJETANMC) {
       gProof->UploadPackage("pars/JETAN.par");
       gProof->EnablePackage("JETAN");
    }   
    // --- Enable particle correlation analysis
    if (iPWG4UE) {
       gProof->UploadPackage("pars/PWG4JetTasks.par");
       gProof->EnablePackage("PWG4JetTasks");
    }   

    // Create the chain
    //

    // Make the analysis manager
    AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", "A test setup for the analysis train");
    // Top container for input 
    AliAnalysisDataContainer *cinput = mgr->CreateContainer("cInput",TChain::Class(), 
							     AliAnalysisManager::kInputContainer);
    if (iAODanalysis) {
    // AOD input handler
       AliAODInputHandler *aodH = new AliAODInputHandler();
       mgr->SetInputEventHandler(aodH);
    } else {   
    // ESD input handler
       AliESDInputHandler *esdHandler = new AliESDInputHandler();
       mgr->SetInputEventHandler(esdHandler);
//       esdHandler->SetInactiveBranches("FMD CaloCluster");
    }
    // Monte Carlo handler
    if (useMC && !iAODanalysis) {
       AliMCEventHandler* mcHandler = new AliMCEventHandler();
       mgr->SetMCtruthEventHandler(mcHandler);
       mcHandler->SetReadTR(readTR); 
    }   
    
    // This container is managed by the AOD handler
    AliAnalysisDataContainer *cout_aod = 0;
    if (iAODhandler) {
       // AOD output handler
       AliAODHandler* aodHandler   = new AliAODHandler();
       mgr->SetOutputEventHandler(aodHandler);
       aodHandler->SetOutputFileName("AliAODs.root");
//       aodHandler->SetCreateNonStandardAOD();
       cout_aod = mgr->CreateContainer("cAOD", TTree::Class(),
							      AliAnalysisManager::kOutputContainer, "default");
       cout_aod->SetSpecialOutput();
    }   

    // Debugging if needed
    if (debug) mgr->SetDebugLevel(3);
//    AliLog::EnableDebug(kTRUE);
//    AliLog::SetGlobalLogLevel(2);


    if (iESDfilter && !iAODanalysis) {
       // Set of cuts plugged into the ESD filter
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
       AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
       trackFilter->AddCuts(esdTrackCutsL);
       //
       // ESD filter task putting standard info to output AOD (with cuts)
       AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
       esdfilter->SetTrackFilter(trackFilter);
       esdfilter->SetDebugLevel(10);
       mgr->AddTask(esdfilter);
       // Connect to data containers
       mgr->ConnectInput  (esdfilter,  0, cinput  );
       mgr->ConnectOutput (esdfilter,  0, cout_aod );
    }   
    // Jet analysis
    AliAnalysisDataContainer *c_aodjet = 0;
    if (iJETAN && !iAODanalysis) {
       AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
       jetana->SetDebugLevel(10);
       jetana->SetConfigFile("${ALICE_ROOT}/JETAN/ConfigJetAnalysis.C");
       mgr->AddTask(jetana);
       // Output histograms list for jet analysis                       
       AliAnalysisDataContainer *cout_jet = mgr->CreateContainer("jethist", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "jethist.root");
       // Dummy AOD output container for jet analysis (no client yet)
       c_aodjet = mgr->CreateContainer("cAODjet", TTree::Class(),
                           AliAnalysisManager::kExchangeContainer);
       // Connect to data containers
       mgr->ConnectInput  (jetana,     0, cinput  );
       mgr->ConnectOutput (jetana,     0, c_aodjet );
       // mgr->ConnectOutput (jetana,     0, cout_aod );
       mgr->ConnectOutput (jetana,     1, cout_jet );
    }   
    // Jet analysisMC
    AliAnalysisDataContainer *c_aodjetMC = 0;
    if (iJETANMC && useMC) {
       AliAnalysisTaskJets *jetanaMC = new AliAnalysisTaskJets("JetAnalysisMC");
       jetanaMC->SetDebugLevel(10);
       jetanaMC->SetConfigFile("${ALICE_ROOT}/JETAN/ConfigJetAnalysisMC.C");
       jetanaMC->SetNonStdBranch("jetsMC");
       mgr->AddTask(jetanaMC);
       // Output histograms list for jet analysis                       
       AliAnalysisDataContainer *cout_jetMC = mgr->CreateContainer("jethistMC", TList::Class(),
							      AliAnalysisManager::kOutputContainer, "jethistMC.root");
       // Dummy AOD output container for jet analysis (no client yet)
       c_aodjetMC = mgr->CreateContainer("cAODjetMC", TTree::Class(),
                           AliAnalysisManager::kExchangeContainer);
       // Connect to data containers
       mgr->ConnectInput  (jetanaMC,     0, cinput  );
       mgr->ConnectOutput (jetanaMC,     0, c_aodjetMC );
       // mgr->ConnectOutput (jetanaMC,     0, cout_aod );
       mgr->ConnectOutput (jetanaMC,     1, cout_jetMC );
    }   
    
    // Particle correlation analysis
    if (iPWG4UE) {
      AliAnalysisTaskUE* ueana = new  AliAnalysisTaskUE("Undelying Event");
      

      // default parameters use a switch via iPWGUE
      // or a config file
      Int_t anaType =1; 
      Int_t regType =1;
      Double_t jetEtaCut=0.2;
      Double_t trackPtCut=0.5; 
      Double_t trackEtaCut= 0.9; 
      Double_t rad=0.7; 
      Double_t deltaPhiCut = 2.616;

      ueana->SetDebugLevel(10); 
      ueana->SetPtRangeInHist(25, 0., 250.);
      ueana->SetAnaTopology(anaType);      
      ueana->SetRegionType(regType);        
      ueana->SetJet1EtaCut(jetEtaCut);     
      ueana->SetTrackPtCut(trackPtCut); 
      ueana->SetPtSumOrdering(2);
      ueana->SetConeRadius(rad);   
      ueana->SetTrackEtaCut(trackEtaCut);
      ueana->SetJet2DeltaPhiCut(deltaPhiCut);
      mgr->AddTask(ueana);


      AliAnalysisDataContainer *coutput1_UE = mgr->CreateContainer("histosUE", TList::Class(),AliAnalysisManager::kOutputContainer, "histosUE.root");

      //      mgr->ConnectInput  (ueana,  0, cinput);    
      mgr->ConnectInput  (ueana,  0, c_aodjet);    
      mgr->ConnectOutput (ueana,     0, coutput1_UE );
    }   

    // Run the analysis
    //    
    if (mgr->InitAnalysis()) {
       mgr->PrintStatus();
       mgr->StartAnalysis("proof",dataset.Data(), 20000);
    }   
}
