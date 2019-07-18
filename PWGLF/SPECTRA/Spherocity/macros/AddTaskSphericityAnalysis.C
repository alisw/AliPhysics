AliAnalysisSphericityTask* AddTaskSphericityAnalysis(
          Bool_t AnalysisMC = kFALSE,             //if TRUE change StoreMcIn as well
					Int_t typerun =1,                       //0 for pp and 1 for Pb-Pb or pPb
					UInt_t kTriggerInt =AliVEvent::kMB,     //for pPb kINT7, for pp or PbPb kMB
					Float_t minCent = 0.,                   //minimum centrality
					Float_t maxCent = 80.,                  //maximum centrality
		  		Bool_t ispileuprej = kTRUE,             //def. kFALSE
          Bool_t Strange = kTRUE                  //selection of event with analysed strange particles
					)
{
  // Get the pointer to the existing analysis manager via the static
  //access method
  //=========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHighPtDeDx", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskHighPtDeDx", "This task requires an input event handler");
    return NULL;
  }  
  ////////////////////////////////////// analysis cuts //////////////////////////////////////////////

  gROOT->LoadMacro("$(ALICE_PHYSICS)/PWGJE/macros/CreateTrackCutsPWGJE.C");

  AliAnalysisFilter* trackFilterGolden = new AliAnalysisFilter("trackFilter");
  AliESDtrackCuts* esdTrackCutsGolden = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE,1);
  trackFilterGolden->AddCuts(esdTrackCutsGolden);

  //Hybrid track cuts,  https://twiki.cern.ch/twiki/bin/viewauth/ALICE/HybridTracks
  //                    https://twiki.cern.ch/twiki/pub/ALICE/PWGPPAODTrackCuts/FilterBitProposal.pdf    
  AliAnalysisFilter* trackFilterHybrid1 = new AliAnalysisFilter("trackFilterHybrid");       //abs (eta range)= 0.9; pT (0.15;1e10)
  AliESDtrackCuts*  trackCutsHybrid1 = CreateTrackCutsPWGJE(10001006);
  trackFilterHybrid1->AddCuts(trackCutsHybrid1);
  //second part, these tracks need to be constrained
  AliAnalysisFilter* trackFilterHybrid2 = new AliAnalysisFilter("trackFilterHybrid2");
  AliESDtrackCuts*  trackCutsHybrid2 = CreateTrackCutsPWGJE(10041006);
  trackFilterHybrid2->AddCuts(trackCutsHybrid2);

  //////////////////////////////////////////////////////////////////////////////////////////////////////
   
  Double_t fV0Sels[7];
  Double_t fCascSels[15];

  fV0Sels[0] =    33.   ;  // max allowed chi2
  fV0Sels[1] =    0.02  ;  // min allowed impact parameter for the 1st daughter     (0.05) 
  fV0Sels[2] =    0.02  ;  // min allowed impact parameter for the 2nd daughter     (0.05)
  fV0Sels[3] =    3.0   ;  // max allowed DCA between the daughter tracks           (1.5)
  fV0Sels[4] =    0.95  ;  // min allowed cosine of V0's pointing angle             (0.98)
  fV0Sels[5] =    0.1   ;  // min radius of the fiducial volume                     (0.2)
  fV0Sels[6] =    100.  ;  // max radius of the fiducial volume                     (100.0)

  fCascSels[0] =  33.   ;  // max allowed chi2 v0
  fCascSels[1] =  0.01  ;  // min allowed impact parameter for the 1st daughter v0  (0.02) 
  fCascSels[2] =  0.01  ;  // min allowed impact parameter for the 2nd daughter v0  (0.02)
  fCascSels[3] =  4.0   ;  // max allowed DCA between the v0 daughter tracks        (2.0)
  fCascSels[4] =  0.95  ;  // min allowed cosine of Casc's pointing angle v0        (0.97)
  fCascSels[5] =  0.5   ;  // min v0 radius of the fiducial volume                  (1.0)
  fCascSels[6] =  100.  ;  // max v0 radius of the fiducial volume                  (100.0)
  fCascSels[7] =  33.   ;  // max allowed chi2 xi
  fCascSels[8] =  0.02  ;  // min allowed V0 impact parameter                       (0.05)
  fCascSels[9] =  0.0012;  // "window" around the Lambda mass                       (0.008)
  fCascSels[10] = 0.01  ;  // min allowed bachelor's impact parameter               (0.03)
  fCascSels[11] = 4.0   ;  // max allowed DCA between the V0 and the bachelor       (2.0)
  fCascSels[12] = 0.96  ;  // min allowed cosine of the cascade pointing angle      (0.97)
  fCascSels[13] = 0.01  ;  // min radius of the fiducial volume                     (0.04)
  fCascSels[14] = 100.  ;  // max radius of the fiducial volume                     (100)


  //////////////////////////////////////   "ESD" or "AOD", pp or PbPb, MC, centrality   ///////////////////////////////////
  AliAnalysisSphericityTask* taskESA = new AliAnalysisSphericityTask("taskESA");
  TString type = mgr->GetInputEventHandler()->GetDataType();                // can be "ESD" or "AOD"
  taskESA->SetAnalysisType(type);
  taskESA->SetAnalysisMC(AnalysisMC);
  taskESA->SetMinCent(minCent);                         
  taskESA->SetMaxCent(maxCent);                         
  taskESA->SetCentralityEstimator("V0M");

  if(typerun==1)
  {
    taskESA->SetAnalysisPbPb(kTRUE);
  }
  else
    taskESA->SetAnalysisPbPb(kFALSE);

  ////////////////////////////////////// PID - Strangeness ////////////////////////////////////////////////////////////////

  taskESA->SetMinDaughterTpcClusters(70);                 // select TPC clusters>80
  taskESA->SetQualityCutTPCrefit(kFALSE);                  // Boolean : kTRUE = ask for TPCrefit for the 3 daughter tracks
  taskESA->SetPIDMode(AliAnalysisSphericityTask::kSigma,3,100,1.5);
  taskESA->SetExtraSelections(kFALSE);                     //set basic cuts of background
  taskESA->SetExtraSelectionsCut(kTRUE);                  //set cuts on invariant mass
  taskESA->SetRerunVertexers(kFALSE);
  taskESA->SetV0Cuts(fV0Sels);
  taskESA->SetCascadeCuts(fCascSels);
  taskESA->SetMaxV0Rapidity(0.8);                         // select |y|<0.75
  taskESA->SetInvMassCutKaon(0.008);                      // 0.008
  taskESA->SetInvMassCutLambda(0.008);                    // 0.008
  taskESA->SetInvMassCutXi(0.008);                        // 0.008
  taskESA->SetInvMassCutOmega(0.010);                     // 0.008

  ////////////////////////////////////// Event Shape Analysis //////////////////////////////////////////////////
  
  taskESA->SetUseHybridESA(kTRUE);                        //def. kTRUE - higher multiplicity cuts
  taskESA->SetTrackFilterESAHyb1(trackFilterHybrid1);
  taskESA->SetTrackFilterESAHyb2(trackFilterHybrid2);
  taskESA->SetTrackFilterESA(trackFilterGolden);
  taskESA->SetMinMultForESA(3);                           //set condition in AliTransverseEventShape::GetSpherocity evente tracks controll if condition false nTrack<MinMult  return -0.5 to hso and hst
  taskESA->SetStepSizeESA(0.1);
  taskESA->SetIsEtaAbsESA(kTRUE);                        
  taskESA->SetTrackEtaMinESA(0.0);
  taskESA->SetTrackEtaMaxESA(0.8);
  taskESA->SetTrackPtMinESA(0.15);
  taskESA->SetTrackPtMaxESA(5);                           //pT maximum cut               
  taskESA->SetDebugLevel(0);
  taskESA->SetEtaCut(0.8);
  taskESA->SetVtxCut(10.0);
  taskESA->SetTrigger(kTriggerInt);                       //trigger for data - ftrigBit
  taskESA->SetPileUpRej(ispileuprej);
  taskESA->SetStoreMcIn(AnalysisMC);
  taskESA->SetStrangeness(Strange);                    
  
  mgr->AddTask(taskESA);

  ////////////////////////////////////////////////////     OUTPUT   ///////////////////////////////////////////////////////////

  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
  //=======================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cout_histdedx;
  cout_histdedx=0;
  cout_histdedx = mgr->CreateContainer("OutputAna", TList::Class(), AliAnalysisManager::kOutputContainer, outputFileName);
  mgr->ConnectInput (taskESA, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskESA, 1, cout_histdedx);

  // Return task pointer at the end
  return taskESA;

}
