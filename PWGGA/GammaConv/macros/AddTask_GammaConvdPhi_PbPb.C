AliAnalysisTask *AddTask_GammaConvdPhi_PbPb(TString v0Cut = "109000200209297002322000000",
					 TString pionCut = "01522045009000",
					 TString photoncut = "",
					 Bool_t pbpb = kTRUE) {

  // standard with task
  printf("========================================================================================\n");
  printf("dPhiAnalysis: Initialising AliAnalysisTaskdPhi\n");
  printf("========================================================================================\n");
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTask_slindal_dPhi", "No analysis manager found.");
    return 0;
  }

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  AliAnalysisTaskdPhi *task = new AliAnalysisTaskdPhi((TString("slindalTask_dPhi")+"_" + v0Cut));

  ///Axes for histrograms
  Double_t cptbins[13] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 25, 50, 100};
  task->GetAxiscPt().Set(12, cptbins);

  Double_t tptbins[10] = {2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15, 50, 100};
  task->GetAxistPt().Set(9, tptbins);

  if(pbpb) {
    Double_t centBins[7] = {0, 5, 10, 30, 50, 60, 90};
    task->GetAxisCent().Set(6, centBins);
  } else {
    Double_t centBins[2] = {-9999, 9999};
    task->GetAxisCent().Set(1, centBins);
  }

  Double_t zbins[6] = { -10, -5, -1.5, 1.5, 5, 10};
  task->GetAxisZ().Set(5, zbins);

  Double_t mbins[17] = {0.07, 0.09, 0.1, 0.11, 0.12, 0.125, 0.1275, 0.13, 0.14, 0.1425, 0.145, 0.15, 0.16, 0.18, 0.2, 0.24, 0.26};
  task->GetAxisPiMass().Set(16, mbins);


  //AliConversionTrackCuts * cuts = new AliConversionTrackCuts();
  AliESDtrackCuts * trackCuts = new AliESDtrackCuts();

  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  trackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
  trackCuts->SetMaxChi2PerClusterTPC(4);
  trackCuts->SetRequireTPCStandAlone(kTRUE);
  trackCuts->SetAcceptKinkDaughters(kFALSE);
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetMaxFractionSharedTPCClusters(0.4);

  trackCuts->SetMaxDCAToVertexXY(2.4);
  trackCuts->SetMaxDCAToVertexZ(3.2);
  trackCuts->SetDCAToVertex2D(kTRUE);

  trackCuts->SetMaxChi2PerClusterITS(36);
  trackCuts->SetMaxChi2TPCConstrainedGlobal(36);

  trackCuts->SetRequireSigmaToVertex(kFALSE);

  trackCuts->SetEtaRange(-0.9, 0.9);
  trackCuts->SetPtRange(cptbins[0], 1000000.0);

  trackCuts->SetRequireITSRefit(kFALSE);

  //cuts->SetESDCuts(esdTrackCuts);
  task->SetTrackCuts(trackCuts);
  
 
  ///Pion cuts
  AliConversionMesonCuts * picuts = new AliConversionMesonCuts("dphi_pioncuts");
  picuts->InitializeCutsFromCutString(pionCut);
  task->SetMesonFilter(picuts);

  if(photoncut.Length() > 0) {
    ///V0 analysis cuts (applied before pion analysis)
    AliConversionCuts * gcuts = new AliConversionCuts();
    gcuts->InitializeCutsFromCutString(photoncut);
    task->SetPhotonFilter(gcuts);
  }

  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString baseString("slindal_dPhi_");
  baseString += v0Cut + "_" + pionCut;

  //dumm output container
  AliAnalysisDataContainer *coutput0 =
	mgr->CreateContainer(baseString+"tree",
						 TTree::Class(),
						 AliAnalysisManager::kExchangeContainer,
						 "slindal_default");
  
  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
	mgr->CreateContainer(baseString+"me", TList::Class(),
						 AliAnalysisManager::kOutputContainer,baseString+".root");

  //define output containers, please use 'username'_'somename'
  
  
  //========= Add PID Reponse to ANALYSIS manager ====
  Bool_t isMC = kFALSE;
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }
   
  //=========  Set Cutnumber for V0Reader ================================
  //  TString cutnumber = "100000000008400100150000000"; 
  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  AliV0ReaderV1 * fV0ReaderV1 = (AliV0ReaderV1*)mgr->GetTask("V0ReaderV1");
  if( !fV0reader ){
    AliV0ReaderV1 *fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
      
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

    // Set AnalysisCut Number
    AliConversionCuts *fCuts=NULL;
    if(v0Cut!=""){
      fCuts= new AliConversionCuts(v0Cut.Data(),v0Cut.Data());
      if(pbpb) {
	fCuts->SetPreSelectionCutFlag(kFALSE);
      } else {
	fCuts->SetPreSelectionCutFlag(kTRUE);
      }
      if(fCuts->InitializeCutsFromCutString(v0Cut)){
	fV0ReaderV1->SetConversionCuts(fCuts);
	//fCuts->SetFillCutHistograms("",kTRUE);
      }


    }
      
    fV0ReaderV1->Init();
    mgr->AddTask(fV0ReaderV1);
    mgr->ConnectInput(fV0ReaderV1,0,cinput);
    task->SetV0Reader(fV0ReaderV1);
    task->SaveReaderHists();
  } else {
    ///V0 analysis cuts (applied before pion analysis)
    AliConversionCuts * v0cuts = new AliConversionCuts();
    v0cuts->InitializeCutsFromCutString(v0Cut);
    task->SetV0Reader(fV0ReaderV1);
    task->SetV0Filter(v0cuts);
  }

  mgr->AddTask(task);
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);
  
  return task;
}


