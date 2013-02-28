void AddTask_dNdPtPbPb2011()
{
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  if (!mgr) {
    Error("AddTask_dNdPtPbPb2011", "No analysis manager found.");
    return 0;
  }
  
  
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  
  cout << "HasMC: " << hasMC << endl;
  
  // Switch off all AliInfo (too much output!!!)
  AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->SetDebugLevel(0);
  
  //
  // Create physics trigger selection class
  //
  
  
  // AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();
  
  //
  // Create event cuts
  //
  Float_t zvWindow = 30. ;
  
  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kTRUE);
  
  //
  // Create geom. acceptance cuts
  //
  Float_t etaWindow = 1. ;
  Float_t ptMin = 0.10 ;
  
  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);
  
  // cut out trouble sector(s)
  Float_t phiMin = TMath::Pi() / 9. * 6.;
  Float_t phiMax = TMath::Pi() / 9. * 7.;
  // additional cut to remove one sector
  AlidNdPtAcceptanceCuts *recCuts = new AlidNdPtAcceptanceCuts("AlidNdPtRecAcceptanceCuts","Geom. acceptance cuts for Recostructed tracks");  
  recCuts->SetExcludeEtaPhiRange(-1.,0.,phiMin,phiMax);
  
  
  //
  // Create standard esd track cuts
  //
  //   Int_t cutMode = 200;
  
  //gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");  
  //   gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/ChargedHadrons/dNdPt/macros/CreatedNdPtTrackCuts.C");
  //   AliESDtrackCuts* esdTrackCuts = CreatedNdPtTrackCuts(cutMode);
  //   if (!esdTrackCuts) {
 //     printf("ERROR: esdTrackCuts could not be created\n");
 //     return;
 //   } else {
   //esdTrackCuts->SetHistogramsOn(kTRUE);
   //     esdTrackCuts->SetHistogramsOn(kTRUE);
   //   }
   //   esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
   //   esdTrackCuts->SetMaxChi2PerClusterITS(36.);
   //   esdTrackCuts->SetMaxDCAToVertexXY(3.0); // is in Michaels cutmode 2222
   
   Float_t minNCrossedRowsTPC = 120;
   Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8;
   Float_t maxFractionSharedTPCCluster = 0.4;
   Double_t maxchi2perTPCcl=4.;
   Double_t maxdcazITSTPC=2.0;
   Double_t maxdaczTPC=3.0;
   Double_t maxdcaxyTPC=3.0;
   
   AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
   esdTrackCuts->SetRequireTPCRefit(kTRUE);
   esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
   esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
   esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
   esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);
   esdTrackCuts->SetMaxDCAToVertexZ(maxdaczTPC);
   esdTrackCuts->SetMaxDCAToVertexXY(maxdcaxyTPC);
   esdTrackCuts->SetRequireITSRefit(kTRUE);
   esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
   esdTrackCuts->SetMaxChi2PerClusterITS(36.);
   esdTrackCuts->SetDCAToVertex2D(kFALSE);
   esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
   esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
   esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
   esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
   esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);
   esdTrackCuts->SetHistogramsOn(kTRUE);
   
   //
   // Create task
   //
   AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask");
   task->SetUseMCInfo(hasMC);
   
   // trigger selection: MB
//    task->SelectCollisionCandidates(AliVEvent::kINT7); 
   
   //
   // set analysis options from the Helper here
   //
   AlidNdPtHelper::OutputObject outputObject = AlidNdPtHelper::kAnalysisPbPb;
   AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPCITS ;
   AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart ;
   
   
   //
   // Create cut analysis object
   //
   if(outputObject==AlidNdPtHelper::kAnalysisPbPb){
     
     //gROOT->LoadMacro("dNdPtPbPb/AlidNdPtAnalysisPbPb.cxx+");
     AlidNdPtAnalysisPbPb2011 *fdNdPtAnalysisPbPb = new AlidNdPtAnalysisPbPb2011("dNdPtAnalysisPbPb2011","dN/dPt Analysis");
     fdNdPtAnalysisPbPb->SetEventCuts(evtCuts);
     fdNdPtAnalysisPbPb->SetAcceptanceCuts(accCuts);
     fdNdPtAnalysisPbPb->SetTrackCuts(esdTrackCuts);
     fdNdPtAnalysisPbPb->SetAnalysisMode(analysisMode);
     fdNdPtAnalysisPbPb->SetParticleMode(particleMode); 
     //      fdNdPtAnalysisPbPb->SetCentralityEstimator("ZNA");
     fdNdPtAnalysisPbPb->SetCentralityEstimator("V0A");
     fdNdPtAnalysisPbPb->SetTriggerMask(AliVEvent::kMB);
//      fdNdPtAnalysisPbPb->SetTriggerMask(AliVEvent::kCINT5);
     //fdNdPtAnalysisPbPb->SetTriggerMask(AliVEvent::kEMC1);
     
     // cut to remove tpc sector
     //fdNdPtAnalysisPbPb->SetRecAcceptanceCuts(recCuts);
     
     // change binning
     // change binning
//      Double_t centralitybins[12] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
     Double_t centralitybins[6] = {0., 20., 40., 60., 80., 100.};  
     
     Double_t ptbins[85] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0, 300.0, 400.0, 500.0};
     Double_t ptcorrbins[85] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0, 300.0, 400.0, 500.0};
     Double_t multbins[48] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,19.5, 20.5, 30.5, 40.5 , 50.5 , 60.5 , 70.5 , 80.5 , 90.5 , 100.5,200.5, 300.5, 400.5, 500.5, 600.5, 700.5, 800.5, 900.5, 1000.5, 2000.5, 3000.5, 4000.5, 5000.5, 6000.5, 7000.5, 8000.5, 9000.5, 10000.5 };
//      Double_t etabins[31] = {-1.465409,-1.365409,-1.265409,-1.165409,-1.065409,-0.965409,-0.865409,-0.765409,-0.665409,-0.565409,-0.465409,-0.365409,-0.265409,-0.165409,-0.065409,0.034591,0.134591,0.234591,0.334591,0.434591,0.534591,0.634591,0.734591,0.834591,0.934591,1.034591,1.134591,1.234591,1.334591,1.434591,1.534591};
     Double_t zvbins[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
     
     
     //Double_t* binsPt = new Double_t[85];
     //for (int i=0; i<85; i++) {binsPt[i] = ptbins[i];}
     
     //fdNdPtAnalysisPbPb->SetBinsPt(85, ptbins);
     //fdNdPtAnalysisPbPb->SetBinsPtCorr(85, ptbins);  
     //fdNdPtAnalysisPbPb->SetBinsCentrality(4, centralitybins);
     
     if(ptbins)
     {
       Int_t nptbins = sizeof(ptbins) / sizeof(Double_t);
       Printf("Setting %i ptbins", nptbins);
       fdNdPtAnalysisPbPb->SetBinsPt(nptbins, ptbins);
     }
     
     if(ptbins)
     {
       Int_t nptcorrbins = sizeof(ptcorrbins) / sizeof(Double_t);
       Printf("Setting %i ptcorrbins", nptcorrbins);
       fdNdPtAnalysisPbPb->SetBinsPtCorr(nptcorrbins, ptcorrbins);
     }
     
     if(centralitybins)
     {
       Int_t ncentralitybins = sizeof(centralitybins) / sizeof(Double_t);
       Printf("Setting %i centralitybins", ncentralitybins);
       fdNdPtAnalysisPbPb->SetBinsCentrality(ncentralitybins, centralitybins);
     }
     
//      if(etabins)
//      {
//        Int_t netabins = sizeof(etabins) / sizeof(Double_t);
//        Printf("Setting %i etabins", netabins);
//        fdNdPtAnalysisPbPb->SetBinsEta(netabins, etabins);
//      }
     
     if(hasMC) 
     {
       fdNdPtAnalysisPbPb->SetUseMCInfo(kTRUE);
       fdNdPtAnalysisPbPb->SetHistogramsOn(kTRUE);
     }
     else 
     { 
       // online trigger
       // fdNdPtAnalysisPbPb->SetPhysicsTriggerSelection(physTrigSel); 
     }
     
     
     task->AddAnalysisObject( fdNdPtAnalysisPbPb );
   }
   
   // Add task
   mgr->AddTask(task);
   
   // Create containers for input
   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   
   TString outputFileName = AliAnalysisManager::GetCommonFileName();
   
   AliAnalysisDataContainer *coutput  = mgr->CreateContainer("dNdPtPbPb", 
							     TList::Class(),
							     AliAnalysisManager::kOutputContainer, 	
							     Form("%s:dNdPtHistos", mgr->GetCommonFileName()));
							     
							     mgr->ConnectInput(task, 0, cinput);
							     mgr->ConnectOutput(task, 1, coutput);
							     
							     
}

