void runAODFilterMC()
{
  // PROOF example

      gSystem->Load("libTree.so");
      gSystem->Load("libPhysics.so");
      gSystem->Load("libGeom.so");
      gSystem->Load("libVMC.so");

      bool bKineFilter = true;
      //      TProof::Mgr("alicecaf")->SetROOTVersion("v5-21-01-alice_dbg");
      TProof::Open("alicecaf");
      //      gProof->SetParallel(1);

      char *dataset = "/COMMON/COMMON/LHC08c11_10TeV_0.5T";
      //      char *dataset = "/PWG4/kleinb/LHC08q_jetjet100";
      // gProof->ClearPackages();
      gProof->UploadPackage("${ALICE_ROOT}/STEERBase.par");
      gProof->EnablePackage("STEERBase");
      gProof->UploadPackage("${ALICE_ROOT}/ESD.par");
      gProof->EnablePackage("ESD");
      gProof->UploadPackage("${ALICE_ROOT}/AOD.par");
      gProof->EnablePackage("AOD");
      gProof->UploadPackage("${ALICE_ROOT}/ANALYSIS.par");
      gProof->EnablePackage("ANALYSIS");
      gProof->UploadPackage("${ALICE_ROOT}/ANALYSISalice.par");
      gProof->EnablePackage("ANALYSISalice");
      
      //
      if (gApplication) gApplication->InitializeGraphics();
      // Create the chain
      //


	/////////////////////////////////////////////////////////////////////////////////// 
	// Create the analysis manager
	//
	// Input 
      AliESDInputHandler* inpHandler = new AliESDInputHandler();
      // Output
      AliAODHandler* aodHandler = new AliAODHandler();
      aodHandler->SetOutputFileName("aod_ckb2.root");
      // MC Truth
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      AliAnalysisManager *mgr  = new AliAnalysisManager("Filter Manager", "Filter Manager");
     if(bKineFilter){
       mgr->SetMCtruthEventHandler(mcHandler);
     }

      mgr->SetInputEventHandler  (inpHandler);
      mgr->SetOutputEventHandler (aodHandler);
      aodHandler->Dump();

      mgr->SetDebugLevel(10);

      // Filtering of MC particles (decays conversions etc)
      // this task is also needed to set the MCEventHandler
      // to the AODHandler, this will not be needed when
      // AODHandler goes to ANALYSISalice
      AliAnalysisTaskMCParticleFilter *kinefilter = new AliAnalysisTaskMCParticleFilter("Particle Filter");
      if(bKineFilter)mgr->AddTask(kinefilter);
      
      
      // 
      AliESDtrackCuts* esdTrackCutsL = new AliESDtrackCuts("AliESDtrackCuts", "Loose");
      esdTrackCutsL->SetMinNClustersTPC(50);
      esdTrackCutsL->SetMaxChi2PerClusterTPC(3.5);
      esdTrackCutsL->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
      esdTrackCutsL->SetRequireTPCRefit(kTRUE);
      esdTrackCutsL->SetDCAToVertexZ(3.0);
      esdTrackCutsL->SetDCAToVertexXY(3.0);
      esdTrackCutsL->SetDCAToVertex2D(kTRUE);
      esdTrackCutsL->SetRequireSigmaToVertex(kFALSE);
      esdTrackCutsL->SetAcceptKingDaughters(kFALSE);
      
      AliAnalysisFilter* trackFilter = new AliAnalysisFilter("trackFilter");
      trackFilter->AddCuts(esdTrackCutsL);
      
      AliAnalysisTaskESDfilter *esdfilter = new AliAnalysisTaskESDfilter("ESD Filter");
      esdfilter->SetTrackFilter(trackFilter);
      
      mgr->AddTask(esdfilter);
    
    
      //
      // Create containers for input/output
      AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
      AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
      
      coutput1->SetSpecialOutput();

      if(bKineFilter){
	mgr->ConnectInput  (kinefilter,     0, cinput1  );
	mgr->ConnectOutput (kinefilter,     0, coutput1 );
      }

	mgr->ConnectInput  (esdfilter,     0, cinput1  );
	mgr->ConnectOutput (esdfilter,     0, coutput1 );
      
      //
      // Run the analysis
      //    
      mgr->InitAnalysis();
      mgr->PrintStatus();
      mgr->StartAnalysis("proof",dataset,10000);

}
