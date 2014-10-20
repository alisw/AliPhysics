AliAnalysisTask *AddTask_GammaConvdPhi_pp(TString v0Cut = "0000011002093663003800000000",
					 TString pionCut = "01631031009000",
					 Bool_t pbpb = kFALSE) {

	////////////////////CURRENTLY NOT WORKING ///////////////////////////////
	
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
  Double_t cptbins[9] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10};
  task->GetAxiscPt().Set(8, cptbins);

  Double_t tptbins[3] = {3.0, 5.0, 10.0};
  task->GetAxistPt().Set(2, tptbins);

  if(pbpb) {
    Double_t centBins[6] = {0, 5, 10, 30, 50, 60, 90};
    task->GetAxisCent().Set(5, centBins);
  } else {
    Double_t centBins[2] = {-9999, 9999};
    task->GetAxisCent().Set(1, centBins);
  }

  Double_t zbins[6] = { -10, -5, -1.5, 1.5, 5, 10};
  task->GetAxisZ().Set(5, zbins);

  Double_t mbins[17] = {0.07, 0.09, 0.1, 0.11, 0.12, 0.125, 0.1275, 0.13, 0.14, 0.1425, 0.145, 0.15, 0.16, 0.18, 0.2, 0.24, 0.26};
  task->GetAxisPiMass().Set(16, mbins);


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

  AliConversionTrackCuts* tc = new AliConversionTrackCuts();
  tc->SetEsdTrackCuts(trackCuts);
  task->SetTrackFilter(tc);


  ///open side nclusters graph
  AliESDtrackCuts * trackCuts2 = new AliESDtrackCuts();
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","65.+30./20.*x");
  trackCuts2->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
  trackCuts2->SetMaxChi2PerClusterTPC(4);
  trackCuts2->SetRequireTPCStandAlone(kTRUE);
  trackCuts2->SetAcceptKinkDaughters(kFALSE);
  trackCuts2->SetRequireTPCRefit(kTRUE);
  trackCuts2->SetMaxFractionSharedTPCClusters(0.4);

  trackCuts2->SetMaxDCAToVertexXY(2.4);
  trackCuts2->SetMaxDCAToVertexZ(3.2);
  trackCuts2->SetDCAToVertex2D(kTRUE);

  trackCuts2->SetMaxChi2PerClusterITS(36);
  trackCuts2->SetMaxChi2TPCConstrainedGlobal(36);

  trackCuts2->SetRequireSigmaToVertex(kFALSE);

  trackCuts2->SetEtaRange(-0.9, 0.9);
  trackCuts2->SetPtRange(cptbins[0], 1000000.0);

  trackCuts2->SetRequireITSRefit(kFALSE);
  //task->AddTrackFilter(trackCuts2, kTRUE);

  AliConversionTrackCuts tc2 = new AliConversionTrackCuts();
  tc2->SetEsdTrackCuts(trackCuts2);
  task->AddTrackFilter(tc2, kTRUE);


  ///Tight side nclusters graphs
  AliESDtrackCuts * trackCuts3 = new AliESDtrackCuts();
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","75.+30./20.*x");
  trackCuts3->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
  trackCuts3->SetMaxChi2PerClusterTPC(4);
  trackCuts3->SetRequireTPCStandAlone(kTRUE);
  trackCuts3->SetAcceptKinkDaughters(kFALSE);
  trackCuts3->SetRequireTPCRefit(kTRUE);
  trackCuts3->SetMaxFractionSharedTPCClusters(0.4);

  trackCuts3->SetMaxDCAToVertexXY(2.4);
  trackCuts3->SetMaxDCAToVertexZ(3.2);
  trackCuts3->SetDCAToVertex2D(kTRUE);

  trackCuts3->SetMaxChi2PerClusterITS(36);
  trackCuts3->SetMaxChi2TPCConstrainedGlobal(36);

  trackCuts3->SetRequireSigmaToVertex(kFALSE);

  trackCuts3->SetEtaRange(-0.9, 0.9);
  trackCuts3->SetPtRange(cptbins[0], 1000000.0);

  trackCuts3->SetRequireITSRefit(kFALSE);
  //task->AddTrackFilter(trackCuts3, kFALSE);

  AliConversionTrackCuts tc3 = new AliConversionTrackCuts();
  tc3->SetEsdTrackCuts(trackCuts3);
  task->AddTrackFilter(tc3, kFALSE);



  ////Tight xyz track cuts
  AliESDtrackCuts * trackCuts4 = new AliESDtrackCuts();
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  trackCuts4->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
  trackCuts4->SetMaxChi2PerClusterTPC(4);
  trackCuts4->SetRequireTPCStandAlone(kTRUE);
  trackCuts4->SetAcceptKinkDaughters(kFALSE);
  trackCuts4->SetRequireTPCRefit(kTRUE);
  trackCuts4->SetMaxFractionSharedTPCClusters(0.4);

  trackCuts4->SetMaxDCAToVertexXY(2.2);
  trackCuts4->SetMaxDCAToVertexZ(3.0);
  trackCuts4->SetDCAToVertex2D(kTRUE);

  trackCuts4->SetMaxChi2PerClusterITS(36);
  trackCuts4->SetMaxChi2TPCConstrainedGlobal(36);

  trackCuts4->SetRequireSigmaToVertex(kFALSE);

  trackCuts4->SetEtaRange(-0.9, 0.9);
  trackCuts4->SetPtRange(cptbins[0], 1000000.0);

  trackCuts4->SetRequireITSRefit(kFALSE);
  //task->AddTrackFilter(trackCuts4, kFALSE);

  AliConversionTrackCuts tc4 = new AliConversionTrackCuts();
  tc4->SetEsdTrackCuts(trackCuts4);
  task->AddTrackFilter(tc4, kFALSE);


  ////lose xyz track cuts
  AliESDtrackCuts * trackCuts5 = new AliESDtrackCuts();
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  trackCuts5->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
  trackCuts5->SetMaxChi2PerClusterTPC(4);
  trackCuts5->SetRequireTPCStandAlone(kTRUE);
  trackCuts5->SetAcceptKinkDaughters(kFALSE);
  trackCuts5->SetRequireTPCRefit(kTRUE);
  trackCuts5->SetMaxFractionSharedTPCClusters(0.4);

  trackCuts5->SetMaxDCAToVertexXY(2.6);
  trackCuts5->SetMaxDCAToVertexZ(3.6);
  trackCuts5->SetDCAToVertex2D(kTRUE);

  trackCuts5->SetMaxChi2PerClusterITS(36);
  trackCuts5->SetMaxChi2TPCConstrainedGlobal(36);

  trackCuts5->SetRequireSigmaToVertex(kFALSE);

  trackCuts5->SetEtaRange(-0.9, 0.9);
  trackCuts5->SetPtRange(cptbins[0], 1000000.0);

  trackCuts5->SetRequireITSRefit(kFALSE);
  //task->AddTrackFilter(trackCuts5);

  AliConversionTrackCuts tc5 = new AliConversionTrackCuts();
  tc5->SetEsdTrackCuts(trackCuts5);
  task->AddTrackFilter(tc5, kTRUE);


  ///open side shared tpt clusters
  AliESDtrackCuts * trackCuts6 = new AliESDtrackCuts();
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  trackCuts6->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
  trackCuts6->SetMaxChi2PerClusterTPC(4);
  trackCuts6->SetRequireTPCStandAlone(kTRUE);
  trackCuts6->SetAcceptKinkDaughters(kFALSE);
  trackCuts6->SetRequireTPCRefit(kTRUE);
  trackCuts6->SetMaxFractionSharedTPCClusters(0.44);

  trackCuts6->SetMaxDCAToVertexXY(2.4);
  trackCuts6->SetMaxDCAToVertexZ(3.2);
  trackCuts6->SetDCAToVertex2D(kTRUE);

  trackCuts6->SetMaxChi2PerClusterITS(36);
  trackCuts6->SetMaxChi2TPCConstrainedGlobal(36);

  trackCuts6->SetRequireSigmaToVertex(kFALSE);

  trackCuts6->SetEtaRange(-0.9, 0.9);
  trackCuts6->SetPtRange(cptbins[0], 1000000.0);

  trackCuts6->SetRequireITSRefit(kFALSE);
  //task->AddTrackFilter(trackCuts6, kTRUE);

  AliConversionTrackCuts tc6 = new AliConversionTrackCuts();
  tc6->SetEsdTrackCuts(trackCuts6);
  task->AddTrackFilter(tc6, kTRUE);


  ///tight side shared tpt clusters
  AliESDtrackCuts * trackCuts7 = new AliESDtrackCuts();
  TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
  trackCuts7->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
  trackCuts7->SetMaxChi2PerClusterTPC(4);
  trackCuts7->SetRequireTPCStandAlone(kTRUE);
  trackCuts7->SetAcceptKinkDaughters(kFALSE);
  trackCuts7->SetRequireTPCRefit(kTRUE);
  trackCuts7->SetMaxFractionSharedTPCClusters(0.36);

  trackCuts7->SetMaxDCAToVertexXY(2.4);
  trackCuts7->SetMaxDCAToVertexZ(3.2);
  trackCuts7->SetDCAToVertex2D(kTRUE);

  trackCuts7->SetMaxChi2PerClusterITS(36);
  trackCuts7->SetMaxChi2TPCConstrainedGlobal(36);

  trackCuts7->SetRequireSigmaToVertex(kFALSE);

  trackCuts7->SetEtaRange(-0.9, 0.9);
  trackCuts7->SetPtRange(cptbins[0], 1000000.0);

  trackCuts7->SetRequireITSRefit(kFALSE);
  //task->AddTrackFilter(trackCuts7, kFALSE);
  
  AliConversionTrackCuts tc7 = new AliConversionTrackCuts();
  tc7->SetEsdTrackCuts(trackCuts7);
  task->AddTrackFilter(tc7, kFALSE);


  ///Pion cuts
  AliConversionMesonCuts * picuts = new AliConversionMesonCuts("dphi_pioncuts");
  picuts->InitializeCutsFromCutString(pionCut);
  task->SetMesonFilter(picuts);


  //================================================
  //              data containers
  //================================================
  //            find input container
  //below the trunk version
  AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
  TString baseString("slindal_dPhi_");
  baseString += v0Cut + pionCut;

  //dumm output container
  AliAnalysisDataContainer *coutput0 =
	mgr->CreateContainer(baseString+"tree",
						 TTree::Class(),
						 AliAnalysisManager::kExchangeContainer,
						 "slindal_default");
  
  //define output containers, please use 'username'_'somename'
  AliAnalysisDataContainer *coutput1 = 
	mgr->CreateContainer(baseString+"me", TList::Class(),
						 AliAnalysisManager::kOutputContainer,baseString+"me.root");

  //define output containers, please use 'username'_'somename'
  
  
  //========= Add PID Reponse to ANALYSIS manager ====
  Bool_t isMC = kFALSE;
  if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AddTaskPIDResponse(isMC);
  }
   
  //=========  Set Cutnumber for V0Reader ================================
  TString cutnumber = "000000000208400000220000000"; 

  //========= Add V0 Reader to  ANALYSIS manager if not yet existent =====
  AliV0ReaderV1 * fV0ReaderV1 = (AliV0ReaderV1*)mgr->GetTask("V0ReaderV1");
  if( !fV0ReaderV1 ){
    fV0ReaderV1 = new AliV0ReaderV1("V0ReaderV1");
      
    fV0ReaderV1->SetUseOwnXYZCalculation(kTRUE);
    fV0ReaderV1->SetCreateAODs(kFALSE);// AOD Output
    fV0ReaderV1->SetUseAODConversionPhoton(kTRUE);

    if (!mgr) {
      Error("AddTask_V0ReaderV1", "No analysis manager found.");
      return;
    }

    // Set AnalysisCut Number
    AliConversionCuts *fCuts=NULL;
    if(cutnumber!=""){
      fCuts= new AliConversionCuts(cutnumber.Data(),cutnumber.Data());
      fCuts->SetPreSelectionCutFlag(kTRUE);
      }
      if(fCuts->InitializeCutsFromCutString(cutnumber)){
	fV0ReaderV1->SetConversionCuts(fCuts);
      }
      
      fV0ReaderV1->Init();
      mgr->AddTask(fV0ReaderV1);
      mgr->ConnectInput(fV0ReaderV1,0,cinput);
      task->SaveReaderHists();
  }
 
  task->SetV0Reader(fV0ReaderV1);

  AliConversionCuts * v0cuts = new AliConversionCuts();
  v0cuts->InitializeCutsFromCutString(v0Cut);
  task->SetV0Filter(v0cuts);

  Int_t numberOfCuts = 7;
  TString *tcutarray  =new TString[7];
  TString *wcutarray = new TString[7];
  
  //"0000011002093663003800000000",
                // "0000011002093663003800000000"; // -4 5 sigma electron line
  tcutarray[0] = "0000011002092663003800000000"; // -3 5 sigma
  wcutarray[0] = "0000011002091663003800000000"; //dedex -5 -5 sigma

             // "0000011002093663003800000000"; // 663: 2 sigma low pt (0.25)  0.5 sigma: high pt (3.5)
  tcutarray[1] = "0000011002093863003800000000"; // 863  2                      1 
  wcutarray[1] = "0000011002093563003800000000"; // 563  2                      -10

  tcutarray[2] = "0000011002093663003800000000";
  wcutarray[2] = "0000011002093663003800000000";

  tcutarray[3] = "0000011002693663003800000000"; // single pt 6: 0.04 tight
  wcutarray[3] = "0000011002493663003800000000"; // single pt 4: 0.075 loose

                            //                 //tpc cluster cut 9 : 0.6 default
  tcutarray[4] = "0000011002063663003800000000"; //6 : > 0.7 tight 
  wcutarray[4] = "0000011002083663003800000000"; //8 : > 0.35 loose

                                   //          // qt max 3: 0.05 default
  tcutarray[5] = "0000011002093663004800000000"; //4 : < 0.03 tight
  wcutarray[5] = "0000011002093663002800000000"; //2 : < 0.07 loose 

                                    //         // chi2 8 default <20
  tcutarray[6] = "0000011002093663003900000000"; //9 : < 15 tight 
  wcutarray[6] = "0000011002093663003200000000"; //2 : < 30 loose


  ///Add the tight cuts
  for(Int_t i = 0; i < numberOfCuts; i++) {
    AliConversionCuts * v0cuts = new AliConversionCuts("aaa","aaa");
    v0cuts->SetPreSelectionCutFlag(kTRUE);
    if(v0cuts->InitializeCutsFromCutString(tcutarray[i].Data())){
      task->AddV0Filter(v0cuts, kFALSE);
      v0cuts->SetFillCutHistograms("aaa",kFALSE);
    }
  }

  ///Add the tight cuts
  for(Int_t i = 0; i < numberOfCuts; i++) {
    AliConversionCuts * v0cuts = new AliConversionCuts("aaa","aaa");
    v0cuts->SetPreSelectionCutFlag(kTRUE);
    if(v0cuts->InitializeCutsFromCutString(wcutarray[i].Data())){
      task->AddV0Filter(v0cuts, kTRUE);
      v0cuts->SetFillCutHistograms("aaa",kFALSE);
    }
  }

  mgr->AddTask(task);
  task->SelectCollisionCandidates(AliVEvent::kMB);
  mgr->ConnectInput  (task,  0, cinput );
  mgr->ConnectOutput (task,  0, coutput0);
  mgr->ConnectOutput (task,  1, coutput1);
  
  return task;
}


