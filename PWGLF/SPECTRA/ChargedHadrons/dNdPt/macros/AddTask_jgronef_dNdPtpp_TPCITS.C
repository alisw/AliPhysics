/// \file AddTask_jgronef_dNdPtpp_TPCITS.C
/// \brief Basic \b AddTask for pp
///
/// Copied from Michaels train directory on Feb 18, 2015do. To test running over pp-pass4.
///
/// The object can be opend with: 
/// \code AlidNdPtAnalysis *obj = jgronef_dNdPtpp_TPCITS->FindObject("dNdPtAnalysis")
/// THnSparse* fRecTrackHist2 = obj->GetRecTrackHist();
/// fRecTrackHist2->GetNdimensions() \endcode
///
/// Functions often used are:
/// \li GetRecTrackHist() - For pT, Eta and Phi distribution
///
/// \author Anton, Jacek, Michael, Philipp...
/// \author Julius Gronefeld <j.gronefeld@cern.ch>, GSI-Darmstadt
/// \date Feb 18, 2014
/// \see AlidNdPtAnalysis.h

void AddTask_jgronef_dNdPtpp_TPCITS()
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {Error("AddTask_dNdPtAnalysis_TPCITS", "No analysis manager found.");return 0;}

//   Switch off all AliInfo (too much output!!!)
//   AliLog::SetGlobalLogLevel(AliLog::kError);
//   mgr->SetDebugLevel(0);

/// Create physics trigger selection class

  //AliPhysicsSelection *physTrigSel =  new AliPhysicsSelection();
  //physTrigSel->AddBackgroundIdentification(new AliBackgroundSelection());


/// <b> Create event cuts:</b> \li Vertex z-window #pm 10 cm \li Trigger Requirement set to \b kFALSE if no physics selection exists.

  Float_t zvWindow = 10. ;

  AlidNdPtEventCuts *evtCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  evtCuts->SetZvRange(-zvWindow,zvWindow);
  evtCuts->SetMeanXYZv(0.0,0.0,0.0);
  evtCuts->SetSigmaMeanXYZv(1.0,1.0,10.0);
  evtCuts->SetTriggerRequired(kFALSE);


/// <b> Create geom. acceptance cuts:</b> \li Eta window #pm 1 \li Minimum pT = 0.1 Gev/c

  Float_t etaWindow = 1.0 ;
  Float_t ptMin = 0.10;

  AlidNdPtAcceptanceCuts *accCuts = new AlidNdPtAcceptanceCuts("AlidNdPtAcceptanceCuts","Geom. acceptance cuts");
  accCuts->SetEtaRange(-etaWindow,etaWindow);
  accCuts->SetPtRange(ptMin,1.e10);


/// <b> Create standard esd track cuts:</b> \li Stop loading external macro. \li Trackcuts according to "CutMode=222"

  
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");

  //Int_t    minclsTPC=70;
  Float_t minNCrossedRowsTPC = 120; 
  Float_t minRatioCrossedRowsOverFindableClustersTPC = 0.8; 
  Float_t maxFractionSharedTPCCluster = 0.4;
  Double_t maxchi2perTPCcl=4.;
  Double_t maxdcazITSTPC=2.0;
  //
  // TPC
  //
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  
  //esdTrackCuts->SetMinNClustersTPC(minclsTPC);
  esdTrackCuts->SetMinNCrossedRowsTPC(minNCrossedRowsTPC);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minRatioCrossedRowsOverFindableClustersTPC);
  esdTrackCuts->SetMaxChi2PerClusterTPC(maxchi2perTPCcl);
  esdTrackCuts->SetMaxFractionSharedTPCClusters(maxFractionSharedTPCCluster);
  //
  // ITS
  //
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMaxChi2PerClusterITS(36.);
  //
  // primary selection
  //
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMaxDCAToVertexZ(maxdcazITSTPC);
  // DCArphi parametrization (LHC10c pass2)
  // 7*(0.0026+0.0050/pt^1.01)
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  // tpcc cut
  esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36.);    
  
  
  
  Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

/// <b>Create task</b>
/// \code AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask_TPCITS"); task->SetUseMCInfo(hasMC); \endcode
  
  AlidNdPtTask *task = new AlidNdPtTask("AlidNdPtTask_TPCITS");  
  task->SetUseMCInfo(hasMC);
  

  // trigger  
  task->SelectCollisionCandidates(AliVEvent::kMB); 


/// Set analysis options from the helper

  //AliTriggerAnalysis::Trigger trigger = AliTriggerAnalysis::kMB1;
  AlidNdPtHelper::AnalysisMode analysisMode = AlidNdPtHelper::kTPCITS;
  AlidNdPtHelper::ParticleMode particleMode = AlidNdPtHelper::kAllPart;


/// <b>Create analysis object</b>


  AlidNdPtAnalysis *fdNdPtAnalysis = new AlidNdPtAnalysis("dNdPtAnalysis","dN/dPt Analysis with TPC-ITS tracking");
  fdNdPtAnalysis->SetEventCuts(evtCuts);
  fdNdPtAnalysis->SetAcceptanceCuts(accCuts);
  fdNdPtAnalysis->SetTrackCuts(esdTrackCuts);
  //fdNdPtAnalysis->SetBackgroundCuts(backCuts);
  fdNdPtAnalysis->SetAnalysisMode(analysisMode); 
  fdNdPtAnalysis->SetParticleMode(particleMode); 
  //fdNdPtAnalysis->SetTrigger(trigger);
  fdNdPtAnalysis->SetTriggerMask(AliVEvent::kMB);
  //fdNdPtAnalysis->SetTriggerMask(AliVEvent::kEMC1);  
  if(hasMC) 
  {
    //physTrigSel->SetAnalyzeMC();
    //fdNdPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 
    fdNdPtAnalysis->SetUseMCInfo(kTRUE);
    fdNdPtAnalysis->SetHistogramsOn(kTRUE);
    //fdNdPtAnalysis->SetHistogramsOn(kFALSE);
  }
  else { // online trigger
//     physTrigSel->SetUseBXNumbers();
//     physTrigSel->SetComputeBG();
//     fdNdPtAnalysis->SetPhysicsTriggerSelection(physTrigSel); 
  }  
  
/// \li Change binning
  Int_t multNbins = 252;  
  Double_t binsMult[253];
  for (int i=0; i<=multNbins; i++) { binsMult[i] = -0.5 + i; }
  binsMult[252] = 1000.;
  
  const Int_t ptNbins = 81;
  Double_t bins[82] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 180.0, 200.0};
  Double_t* binsPt = new Double_t[82];
  for (int i=0; i<82; i++) {binsPt[i] = bins[i];}
  fdNdPtAnalysis->SetBinsPt(ptNbins, binsPt);
  fdNdPtAnalysis->SetBinsPtCorr(ptNbins, binsPt);  
  fdNdPtAnalysis->SetBinsMult(multNbins, binsMult);
   
  task->AddAnalysisObject( fdNdPtAnalysis );
  mgr->AddTask(task);
  
/// <b>Create containers for input and output </b>
/// \code AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer(); \endcode
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
///  Adjust the Data Container to the Grid:
/// \code AliAnalysisDataContainer *coutput = mgr->CreateContainer("dNdPtpp-Standard",TList::Class(), AliAnalysisManager::kOutputContainer, "dNdPt_jgronef");  \endcode
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("dNdPtpp-Standard",
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   Form("%s:dNdPt_jgronef", mgr->GetCommonFileName())); 

  
  
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput);
  
//   return task; 

}

