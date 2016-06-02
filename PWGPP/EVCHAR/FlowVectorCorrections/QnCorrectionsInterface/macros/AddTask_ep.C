// event plane calibration task
// author: Jaap Onderwaater, jacobus.onderwaater@cern.ch
//         2014/Dec/06


//////////////////////////////////////////////////////////////
//
//    Event plane correction options:
//
//    
//    Calibration method: 0: Calibration of Q/sqrt(M)
//                        1: Calibration of Q/M
//                        2: Calibration of Q/|Q|
//                        3: Calibration of Q
//
//    Calibration step  : 0: Raw
//                        1: Equalization
//                        2: Recentering
//                        3: Twist
//                        4: Scaling
//
//    Equalization method : 0: M/<M>
//                          1: 1+(M-<M>)/sigma(M)
//
//    Channel list      : Array of channel numbers that are included in Q-vector calculation
//
//    Twist and Scaling method: to be implemented
//                              0: Double harmonic track wise (advisable for TPC)
//                              1: Double harmonic Q wise
//                              2: Correlations
//
//    Event plane detector name : Name to give your event plane (has to be unique)
//
//    Correlation detector names: The detectors used to perform the twist and scaling with correlations 
//
//    
///////////////////////////////////////////////////////////////


void DefineHistograms(AliQnCorrectionsManager* QnManager, AliQnCorrectionsHistos* histos, TString histClass);

void AddVZERO(AliQnCorrectionsManager* QnManager);
void AddTPC(AliQnCorrectionsManager* QnManager, Bool_t isESD, Bool_t UseTPConlyTracks);
void AddTZERO(AliQnCorrectionsManager* QnManager);
void AddFMD(AliAnalysisManager *mgr, AliQnCorrectionsManager* QnManager, Bool_t isESD);
void AddFMDTaskForESDanalysis(AliAnalysisManager *mgr);
void AddRawFMD(AliQnCorrectionsManager* QnManager);
void AddZDC(AliQnCorrectionsManager* QnManager);
void AddSPD(AliQnCorrectionsManager* QnManager);


AliAnalysisDataContainer* AddTask_ep() {
  /* temporal flag to use multiplicity instead of centrality and to inhibit detectors for 2015 dataset */
  Bool_t bUseMultiplicity = kFALSE;
  Bool_t b2015DataSet = kFALSE;


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    mgr = new AliAnalysisManager("AnalysisManagerQcorrections");

    Error("AddTask_ep", "No analysis manager found.");
    return 0;
  }

  Bool_t isESD=mgr->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();

  AliQnCorrectionsManager* QnManager = new AliQnCorrectionsManager();
  AliAnalysisTaskFlowVectorCorrections* taskEP = new AliAnalysisTaskFlowVectorCorrections("FlowVectorCorrections");
  taskEP->SetRunListPath("$(ALICE_PHYSICS)/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/runsLHC10h.txt");

  AliQnCorrectionsCuts *eventCuts = new AliQnCorrectionsCuts();
  eventCuts->AddCut(AliQnCorrectionsVarManager::kVtxZ,-10.0,10.0);
  if (bUseMultiplicity) {
    eventCuts->AddCut(AliQnCorrectionsVarManager::kVZEROMultPercentile,0.0,90.0);
  }
  else {
    eventCuts->AddCut(AliQnCorrectionsVarManager::kCentVZERO,0.0,90.0);
  }



  taskEP->SetEventCuts(eventCuts);
  if (!b2015DataSet) {
    taskEP->SelectCollisionCandidates(AliVEvent::kMB);  // Events passing trigger and physics selection for analysis
  }
  else
    taskEP->SelectCollisionCandidates(AliVEvent::kMB|AliVEvent::kINT7);  // Events passing trigger and physics selection for analysis

  if(kFALSE){  // attach calibration file directly to correction framework here in the AddTask
    TGrid::Connect("alien://");
    TFile* inputfile = TFile::Open("alien:///alice/cern.ch/user/j/jonderwa/CalibrationFiles/LHC10h.root");
    Bool_t foundInput = QnManager->SetCalibrationFile(inputfile); 
  }
  else {      // access calibration file from job site in the analysis task
    taskEP->SetCalibrationFilePath("alien:///alice/cern.ch/user/j/jonderwa/CalibrationFiles/LHC10h.root");
  }

  Bool_t UseTPConlyTracks=kFALSE;
  taskEP->GetFillEvent()->SetUseTPCStandaloneTracks(UseTPConlyTracks);  // Use of TPC standalone tracks or Global tracks (only for ESD analysis)

  if (b2015DataSet) {
    AddTPC(QnManager, isESD, UseTPConlyTracks);
    AddSPD(QnManager);
    AddVZERO(QnManager);
    AddTZERO(QnManager);
    AddRawFMD(QnManager);
    AddZDC(QnManager);
  }
  else {
    AddVZERO(QnManager);
    AddTPC(QnManager, isESD, UseTPConlyTracks);
    if(isESD) AddTZERO(QnManager);
    AddZDC(QnManager);
    AddFMD(mgr,QnManager, isESD);
    //AddRawFMD(QnManager);
    //AddSPD(QnManager);
  }


  //QnManager->SetFillTreeQnVectors();
  QnManager->SetFillHistogramsQA();
  QnManager->SetFillHistogramsQnCorrections();

  taskEP->FillExchangeContainerWithQvectors(kTRUE);

  taskEP->SetEventPlaneManager(QnManager);
  taskEP->DefineInOutput();

  AliQnCorrectionsHistos* hists = taskEP->GetHistograms();
  TString histClass = "";
  histClass += "Event_NoCuts;";
  histClass += "Event_Analysis;";
  DefineHistograms(QnManager, hists, histClass);

  mgr->AddTask(taskEP);

  //create output container
  AliAnalysisDataContainer *cOutputHist;
  if(QnManager->ShouldFillHistogramsQnCorrections()){
    cOutputHist = mgr->CreateContainer("CalibrationHistos",
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        "CalibrationHistograms.root");
  }
  
  AliAnalysisDataContainer *cOutputQvec;
  if(QnManager->ShouldFillTreeQnVectors()){
     cOutputQvec = mgr->CreateContainer("CalibratedQvector",
        TTree::Class(),
        AliAnalysisManager::kOutputContainer,
        "QvectorsTree.root");
  }
  
  AliAnalysisDataContainer *cOutputQvecList;
  if(taskEP->IsFillExchangeContainerWithQvectors()){
    cOutputQvecList =
    mgr->CreateContainer("CalibratedQvectorList",
        TList::Class(),
        AliAnalysisManager::kExchangeContainer,
        "QvectorsList.root");
  }
  
  AliAnalysisDataContainer *cOutputHistQA;
  if(QnManager->ShouldFillHistogramsQA()){
    cOutputHistQA = mgr->CreateContainer("CalibrationQA",
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        "CalibrationQA.root");
  }

  AliAnalysisDataContainer *cOutputQnEventQA;
  if(taskEP->IsFillEventQA()){
    cOutputQnEventQA = mgr->CreateContainer("QnEventQA",
        TList::Class(),
        AliAnalysisManager::kOutputContainer,
        "QnEventQA.root");
  }


  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  mgr->ConnectInput(taskEP,  0, mgr->GetCommonInputContainer());
  if(QnManager->ShouldFillHistogramsQnCorrections())    mgr->ConnectOutput(taskEP, taskEP->OutputSlotHistQn()        , cOutputHist );
  if(QnManager->ShouldFillTreeQnVectors())              mgr->ConnectOutput(taskEP, taskEP->OutputSlotTree()          , cOutputQvec );
  if(QnManager->ShouldFillHistogramsQA())               mgr->ConnectOutput(taskEP, taskEP->OutputSlotHistQA()        , cOutputHistQA );
  if(taskEP->IsFillExchangeContainerWithQvectors())     mgr->ConnectOutput(taskEP, taskEP->OutputSlotGetListQnVectors() , cOutputQvecList );
  if(taskEP->IsFillEventQA())     mgr->ConnectOutput(taskEP, taskEP->OutputSlotEventQA()       , cOutputQnEventQA );

  return cOutputQvecList;
}




void AddVZERO(AliQnCorrectionsManager* QnManager){


  Short_t VZEROchannels[4][64];
  for(Int_t iv0=0; iv0<4; iv0++) for(Int_t ich=0; ich<64; ich++) VZEROchannels[iv0][ich] = 0;

  for(Int_t ich=32; ich<64; ich++) VZEROchannels[0][ich] = 1;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<32; ich++) VZEROchannels[1][ich] = 1;


  TArrayS* channelGroups=new TArrayS(64);
  for(Int_t ich=0; ich<64; ich++) channelGroups->SetAt((Int_t) ich/8, ich);

  //-----------------------------------------------------------
  // Binning for Q-vector calibration
  //
  const Int_t nVZEROdim=2;
  AliQnCorrectionsAxes * VZERObinning = new AliQnCorrectionsAxes(nVZEROdim);  //  declare binning object with number of calibration dimensions
  TAxis VZERObinningAxis[nVZEROdim+1];

  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }}; 
  VZERObinningAxis[0] = VZERObinning->MakeAxis(Ctbinning);

  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};  // array of pairs, where the 1st element of each pair is the lower edge of a coarse bin, and the 2nd element is the number of fine bins inside the coarse bin. The 2nd element of the first` pair is the number of coarse bins plus one (i.e. total number of pairs).
  VZERObinningAxis[1] = VZERObinning->MakeAxis(VtxZbinning);

  Double_t VarIdMap[nVZEROdim] = { AliQnCorrectionsVarManager::kCentVZERO , AliQnCorrectionsVarManager::kVtxZ};

  for(Int_t idim=0; idim<nVZEROdim; idim++) VZERObinning->SetAxis(idim, VarIdMap[idim], VZERObinningAxis[idim], AliQnCorrectionsVarManager::VarName(VarIdMap[idim]));
  
  //-----------------------------------------------------------
  // Binning for channel equalization


  AliQnCorrectionsAxes * VZEROchannelbinning = new AliQnCorrectionsAxes(nVZEROdim+1); // +1 is for channel axis
  for(Int_t idim=0; idim<nVZEROdim; idim++) VZEROchannelbinning->SetAxis(idim, VarIdMap[idim], VZERObinningAxis[idim], AliQnCorrectionsVarManager::VarName(VarIdMap[idim]));
  VZEROchannelbinning->SetNchannels(64);

  ////////// end of binning


  AliQnCorrectionsConfiguration * VZEROAconf = new AliQnCorrectionsConfiguration();
  VZEROAconf->SetQnNormalization(1);
  VZEROAconf->SetQnHarmonicsRange(1,4);
  VZEROAconf->SetHarmonicForAlignment(2);
  VZEROAconf->SetDataVectorEqualizationMethod(0);
  VZEROAconf->SetTwistAndRescalingMethod(2);
  VZEROAconf->SetDataVectorIdList(new TArrayS(64, VZEROchannels[0]));
  VZEROAconf->SetDataVectorEqualizationGroups(channelGroups);
  VZEROAconf->SetQnConfigurationName("VZEROA");
  VZEROAconf->SetReferenceQnForAlignment("TPC");
  VZEROAconf->SetReferenceQnForTwistAndRescaling("VZEROC","TPC");
  VZEROAconf->SetQnCorrectionsCommonAxes(VZERObinning);
  VZEROAconf->SetDataVectorEqualizationAxes(VZEROchannelbinning);
  VZEROAconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run

  VZEROAconf->SetQnCorrectionDataVectorEqualization();
  VZEROAconf->SetQnCorrectionRecentering();
  VZEROAconf->SetQnCorrectionAlignment();

  //VZEROAconf->SetCorrectionSteps(AliQnCorrectionsSteps::kDataVectorEqualization, AliQnCorrectionsSteps::kRecentering, AliQnCorrectionsSteps::kAlignment);

  QnManager->AddQnConfiguration(VZEROAconf, AliQnCorrectionsVarManager::kVZERO);

  
  AliQnCorrectionsConfiguration * VZEROCconf = new AliQnCorrectionsConfiguration();
  VZEROCconf->SetQnNormalization(1);
  VZEROCconf->SetQnHarmonicsRange(1,4);
  VZEROCconf->SetHarmonicForAlignment(2);
  VZEROCconf->SetDataVectorEqualizationMethod(0);
  VZEROCconf->SetTwistAndRescalingMethod(2);
  VZEROCconf->SetDataVectorIdList(new TArrayS(64, VZEROchannels[1]));
  VZEROCconf->SetDataVectorEqualizationGroups(channelGroups);
  VZEROCconf->SetQnConfigurationName("VZEROC");
  VZEROCconf->SetReferenceQnForAlignment("TPC");
  VZEROCconf->SetReferenceQnForTwistAndRescaling("VZEROA","TPC");
  VZEROCconf->SetQnCorrectionsCommonAxes(VZERObinning);
  VZEROCconf->SetDataVectorEqualizationAxes(VZEROchannelbinning);
  VZEROCconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run

  VZEROCconf->SetQnCorrectionDataVectorEqualization();
  VZEROCconf->SetQnCorrectionRecentering();
  VZEROCconf->SetQnCorrectionAlignment();

  //VZEROCconf->SetCorrectionSteps(AliQnCorrectionsSteps::kDataVectorEqualization, AliQnCorrectionsSteps::kRecentering, AliQnCorrectionsSteps::kAlignment);

  QnManager->AddQnConfiguration(VZEROCconf, AliQnCorrectionsVarManager::kVZERO);



}

void AddTPC(AliQnCorrectionsManager* QnManager, Bool_t isESD, Bool_t UseTPConlyTracks){

  /////////////// Add TPC subdetectors ///////////////////

  //-----------------------------------------------------------
  // Binning for Q-vector calibration
  //
  const Int_t nTPCdim=2;
  AliQnCorrectionsAxes * TPCbinning = new AliQnCorrectionsAxes(nTPCdim);  //  declare binning object with number of calibration dimensions
  TAxis TPCbinningAxis[nTPCdim];

  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};  // array of pairs, where the 1st element of each pair is the lower edge of a coarse bin, and the 2nd element is the number of fine bins inside the coarse bin. The 2nd element of the first` pair is the number of coarse bins plus one (i.e. total number of pairs).
  TPCbinningAxis[0] = TPCbinning->MakeAxis(VtxZbinning);

  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }}; 
  TPCbinningAxis[1] = TPCbinning->MakeAxis(Ctbinning);

  Double_t VarIdMap[nTPCdim] = { AliQnCorrectionsVarManager::kVtxZ , AliQnCorrectionsVarManager::kCentVZERO};

  for(Int_t idim=0; idim<nTPCdim; idim++) TPCbinning->SetAxis(idim, VarIdMap[idim], TPCbinningAxis[idim], AliQnCorrectionsVarManager::VarName(VarIdMap[idim]));
  

  ////////// end of binning


  AliQnCorrectionsConfiguration * TPCconf = new AliQnCorrectionsConfiguration();
  //TPCconf->SetCorrectionSteps( AliQnCorrectionsSteps::kRecentering, AliQnCorrectionsSteps::kRecentering, AliQnCorrectionsSteps::kRescaling);
  TPCconf->SetQnNormalization(1);
  TPCconf->SetQnHarmonicsRange(1,4);
  TPCconf->SetTwistAndRescalingMethod(0);
  TPCconf->SetQnConfigurationName("TPC");
  //TPCconf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
  TPCconf->SetQnCorrectionsCommonAxes(TPCbinning);
  TPCconf->SetTracking();   // this is a tracking detector
  TPCconf->UseCalibrationDirectoryNameAllEvents();


  AliQnCorrectionsCuts *trackTPC = new AliQnCorrectionsCuts();
  if(!isESD){
    trackTPC->AddFlag(AliQnCorrectionsVarManager::kFilterBitMask768,0.5,1.5);
  //trackTPC->AddFlag(QnCorrectionsVarManager::kFilterBit,QnCorrectionsVarManager::kFilterBit0,kTRUE);
    trackTPC->AddCut(AliQnCorrectionsVarManager::kEta,-0.8,0.8);
    trackTPC->AddCut(AliQnCorrectionsVarManager::kPt,0.2,5.);
  }
  else{
    if(UseTPConlyTracks){
      trackTPC->AddCut(AliQnCorrectionsVarManager::kDcaXY,-3.0,3.0);
      trackTPC->AddCut(AliQnCorrectionsVarManager::kDcaZ ,-3.0,3.0);
      trackTPC->AddCut(AliQnCorrectionsVarManager::kEta,-0.8,0.8);
      //trackTPC->AddCut(AliQnCorrectionsVarManager::kEta,-0.3,0.3, kTRUE);  // exclusion region
      trackTPC->AddCut(AliQnCorrectionsVarManager::kPt,0.2,5.);
      trackTPC->AddCut(AliQnCorrectionsVarManager::kTPCnclsIter1,70.0,161.0);
      trackTPC->AddCut(AliQnCorrectionsVarManager::kTPCchi2,0.2,4.0);
    }
    else{
      trackTPC->AddCut(AliQnCorrectionsVarManager::kDcaXY,-0.3,0.3);
      trackTPC->AddCut(AliQnCorrectionsVarManager::kDcaZ ,-0.3,0.3);
      trackTPC->AddCut(AliQnCorrectionsVarManager::kEta,-0.8,0.8);
      //trackTPC->AddCut(AliQnCorrectionsVarManager::kEta,-0.3,0.3, kTRUE);  // exclusion region
      trackTPC->AddCut(AliQnCorrectionsVarManager::kPt,0.2,5.);
      trackTPC->AddCut(AliQnCorrectionsVarManager::kTPCncls,70.0,161.0);
      trackTPC->AddCut(AliQnCorrectionsVarManager::kTPCchi2,0.2,4.0);
    }
  }

  TPCconf->SetDataVectorCuts(trackTPC);

  TPCconf->SetQnCorrectionRecentering();

  QnManager->AddQnConfiguration(TPCconf, AliQnCorrectionsVarManager::kTPC);


} 


void AddSPD(AliQnCorrectionsManager* QnManager){

  /////////////// Add SPD subdetectors ///////////////////

  //-----------------------------------------------------------
  // Binning for Q-vector calibration
  //
  const Int_t nSPDdim=2;
  AliQnCorrectionsAxes * SPDbinning = new AliQnCorrectionsAxes(nSPDdim);  //  declare binning object with number of calibration dimensions
  TAxis SPDbinningAxis[nSPDdim];

  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};  // array of pairs, where the 1st element of each pair is the lower edge of a coarse bin, and the 2nd element is the number of fine bins inside the coarse bin. The 2nd element of the first` pair is the number of coarse bins plus one (i.e. total number of pairs).
  SPDbinningAxis[0] = SPDbinning->MakeAxis(VtxZbinning);

  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }}; 
  SPDbinningAxis[1] = SPDbinning->MakeAxis(Ctbinning);

  Double_t VarIdMap[nSPDdim] = { AliQnCorrectionsVarManager::kVtxZ , AliQnCorrectionsVarManager::kVZEROMultPercentile};

  for(Int_t idim=0; idim<nSPDdim; idim++) SPDbinning->SetAxis(idim, VarIdMap[idim], SPDbinningAxis[idim], AliQnCorrectionsVarManager::VarName(VarIdMap[idim]));
  

  ////////// end of binning



  AliQnCorrectionsConfiguration * SPDconf = new AliQnCorrectionsConfiguration();
  //SPDconf->SetCorrectionSteps( AliQnCorrectionsSteps::kRecentering, AliQnCorrectionsSteps::kRecentering, AliQnCorrectionsSteps::kRescaling);
  SPDconf->SetQnNormalization(1);
  SPDconf->SetQnHarmonicsRange(1,4);
  SPDconf->SetQnConfigurationName("SPD");
  SPDconf->SetQnCorrectionsCommonAxes(SPDbinning);
  SPDconf->SetTracking();   // this is a tracking detector
  SPDconf->UseCalibrationDirectoryNameAllEvents();

  SPDconf->SetQnCorrectionRecentering();

  QnManager->AddQnConfiguration(SPDconf, AliQnCorrectionsVarManager::kSPD);


} 



void AddTZERO(AliQnCorrectionsManager* QnManager){

  /////////////// Add TZERO subdetectors ///////////////////

  Short_t TZEROchannels[2][24];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<24; ich++) TZEROchannels[iv0][ich] = 0;

  for(Int_t ich=12; ich<24; ich++) TZEROchannels[0][ich] = 1;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<12; ich++) TZEROchannels[1][ich] = 1;


  //-----------------------------------------------------------
  // Binning for Q-vector calibration
  //
  const Int_t nTZEROdim=2;
  AliQnCorrectionsAxes * TZERObinning = new AliQnCorrectionsAxes(nTZEROdim);  //  declare binning object with number of calibration dimensions
  TAxis TZERObinningAxis[nTZEROdim+1];

  Double_t VtxZbinning[][2] = { { -10.0, 4} , {-7.0, 1}, {7.0, 8}, {10.0, 1}};  // array of pairs, where the 1st element of each pair is the lower edge of a coarse bin, and the 2nd element is the number of fine bins inside the coarse bin. The 2nd element of the first` pair is the number of coarse bins plus one (i.e. total number of pairs).
  TZERObinningAxis[0] = TZERObinning->MakeAxis(VtxZbinning);

  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }}; 
  TZERObinningAxis[1] = TZERObinning->MakeAxis(Ctbinning);

  Double_t VarIdMap[nTZEROdim] = {AliQnCorrectionsVarManager::kVtxZ , AliQnCorrectionsVarManager::kCentVZERO };

  for(Int_t idim=0; idim<nTZEROdim; idim++) TZERObinning->SetAxis(idim, VarIdMap[idim], TZERObinningAxis[idim], AliQnCorrectionsVarManager::VarName(VarIdMap[idim]));
  
  //-----------------------------------------------------------
  // Binning for channel equalization


  AliQnCorrectionsAxes * TZEROchannelbinning = new AliQnCorrectionsAxes(nTZEROdim+1); // +1 is for channel axis
  for(Int_t idim=0; idim<nTZEROdim; idim++) TZEROchannelbinning->SetAxis(idim, VarIdMap[idim], TZERObinningAxis[idim], AliQnCorrectionsVarManager::VarName(VarIdMap[idim]));
  TZEROchannelbinning->SetNchannels(24);

  ////////// end of binning


  AliQnCorrectionsConfiguration * TZEROAconf = new AliQnCorrectionsConfiguration();
  TZEROAconf->SetQnNormalization(1);
  TZEROAconf->SetQnHarmonicsRange(1,4);
  TZEROAconf->SetHarmonicForAlignment(2);
  TZEROAconf->SetDataVectorEqualizationMethod(0);
  TZEROAconf->SetTwistAndRescalingMethod(2);
  TZEROAconf->SetDataVectorIdList(new TArrayS(24, TZEROchannels[0]));
  TZEROAconf->SetQnConfigurationName("TZEROA");
  TZEROAconf->SetReferenceQnForAlignment("TPC");
  TZEROAconf->SetReferenceQnForTwistAndRescaling("TZEROC","TPC");
  TZEROAconf->SetQnCorrectionsCommonAxes(TZERObinning);
  TZEROAconf->SetDataVectorEqualizationAxes(TZEROchannelbinning);
  TZEROAconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run

  TZEROAconf->SetQnCorrectionDataVectorEqualization();
  TZEROAconf->SetQnCorrectionRecentering();
  TZEROAconf->SetQnCorrectionAlignment();

  //TZEROAconf->SetCorrectionSteps(AliQnCorrectionsSteps::kDataVectorEqualization, AliQnCorrectionsSteps::kRecentering, AliQnCorrectionsSteps::kAlignment);

  QnManager->AddQnConfiguration(TZEROAconf, AliQnCorrectionsVarManager::kTZERO);

  
  AliQnCorrectionsConfiguration * TZEROCconf = new AliQnCorrectionsConfiguration();
  TZEROCconf->SetQnNormalization(1);
  TZEROCconf->SetQnHarmonicsRange(1,4);
  TZEROCconf->SetHarmonicForAlignment(2);
  TZEROCconf->SetDataVectorEqualizationMethod(0);
  TZEROCconf->SetTwistAndRescalingMethod(2);
  TZEROCconf->SetDataVectorIdList(new TArrayS(24, TZEROchannels[1]));
  TZEROCconf->SetQnConfigurationName("TZEROC");
  TZEROCconf->SetReferenceQnForAlignment("TPC");
  TZEROCconf->SetReferenceQnForTwistAndRescaling("TZEROA","TPC");
  TZEROCconf->SetQnCorrectionsCommonAxes(TZERObinning);
  TZEROCconf->SetDataVectorEqualizationAxes(TZEROchannelbinning);
  TZEROCconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run

  TZEROCconf->SetQnCorrectionDataVectorEqualization();
  TZEROCconf->SetQnCorrectionRecentering();
  TZEROCconf->SetQnCorrectionAlignment();

  //TZEROCconf->SetCorrectionSteps(AliQnCorrectionsSteps::kDataVectorEqualization, AliQnCorrectionsSteps::kRecentering, AliQnCorrectionsSteps::kAlignment);

  QnManager->AddQnConfiguration(TZEROCconf, AliQnCorrectionsVarManager::kTZERO);



}


void AddZDC(AliQnCorrectionsManager* QnManager){
  /////////////// Add ZDC subdetectors ///////////////////

  AliQnCorrectionsVarManager::SetDefaultVarNames();

  Short_t ZDCchannels[2][10];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<10; ich++) ZDCchannels[iv0][ich] = 0;

  for(Int_t ich=6; ich<10; ich++) ZDCchannels[0][ich] = 1;
  for(Int_t ich=1; ich<5; ich++)  ZDCchannels[1][ich] = 1;

  //-----------------------------------------------------------
  // Binning for Q-vector calibration
  //
  const Int_t nZDCdim=3;
  AliQnCorrectionsAxes * ZDCbinning = new AliQnCorrectionsAxes(nZDCdim);  //  declare binning object with number of calibration dimensions
  TAxis ZDCbinningAxis[nZDCdim+1];

  // vertex X bins
  Double_t VtXbinning[][2] = {{ -0.06, 2}, {0.04, 5 }}; 
  ZDCbinningAxis[0] = ZDCbinning->MakeAxis(VtXbinning);

  // vertex Y bins
  Double_t VtYbinning[][2] = {{ 0.12, 2}, {0.24, 5 }}; 
  ZDCbinningAxis[1] = ZDCbinning->MakeAxis(VtYbinning);

  // centrality bins
  Double_t Ctbinning[][2] = {{ 0.0, 2}, {100.0, 100 }}; 
  ZDCbinningAxis[2] = ZDCbinning->MakeAxis(Ctbinning);

  Double_t VarIdMap[nZDCdim] = {AliQnCorrectionsVarManager::kVtxX , AliQnCorrectionsVarManager::kVtxY ,AliQnCorrectionsVarManager::kCentVZERO };

  for(Int_t idim=0; idim<nZDCdim; idim++) ZDCbinning->SetAxis(idim, VarIdMap[idim], ZDCbinningAxis[idim], AliQnCorrectionsVarManager::VarName(VarIdMap[idim]));
  
  Short_t ZDCchannels[2][10];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<10; ich++) ZDCchannels[iv0][ich] = 0;

  for(Int_t ich=6; ich<10; ich++) ZDCchannels[0][ich] = 1;
  for(Int_t ich=1; ich<5; ich++)  ZDCchannels[1][ich] = 1;


  AliQnCorrectionsConfiguration * ZDCAconf = new AliQnCorrectionsConfiguration();
  ZDCAconf->SetQnNormalization(1);
  ZDCAconf->SetQnHarmonicsRange(1,3);
  ZDCAconf->SetDataVectorIdList(new TArrayS(10, ZDCchannels[0]));
  ZDCAconf->SetQnConfigurationName("ZDCA");
  ZDCAconf->SetQnCorrectionsCommonAxes(ZDCbinning);
  ZDCAconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run

  ZDCAconf->SetQnCorrectionRecentering();

  QnManager->AddQnConfiguration(ZDCAconf, AliQnCorrectionsVarManager::kZDC);
  

  AliQnCorrectionsConfiguration * ZDCCconf = new AliQnCorrectionsConfiguration();
  ZDCCconf->SetQnNormalization(1);
  ZDCCconf->SetQnHarmonicsRange(1,3);
  ZDCCconf->SetDataVectorIdList(new TArrayS(10, ZDCchannels[1]));
  ZDCCconf->SetQnConfigurationName("ZDCC");
  ZDCCconf->SetQnCorrectionsCommonAxes(ZDCbinning);
  ZDCCconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run

  ZDCCconf->SetQnCorrectionRecentering();


  QnManager->AddQnConfiguration(ZDCCconf, AliQnCorrectionsVarManager::kZDC);


}

void AddFMDTaskForESDanalysis(AliAnalysisManager *mgr){

  gSystem->Load("libPWGLFforward2");  // for FMD

  // Create the FMD task and add it to the manager
  //===========================================================================


  //--- AOD output handler -----------------------------------------
  AliAODHandler* ret = new AliAODHandler();
  //ret->SetFillAOD(kTRUE);
  ret->SetOutputFileName("AliAOD.pass2.root");
  mgr->SetOutputEventHandler(ret);

  //gROOT->LoadClass("AliAODForwardMult", "libPWGLFforward2");
  gSystem->Load("libESDfilter.so");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/ESDfilter/macros/AddTaskESDFilter.C");
  AliAnalysisTaskESDfilter *esdfilter = AddTaskESDFilter(kTRUE, kFALSE, kFALSE, kFALSE, kFALSE, kTRUE);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (esdfilter,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (esdfilter,  0, mgr->GetCommonOutputContainer());

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/AddTaskForwardMult.C");

  Bool_t   mc  = false; // false: real data, true: simulated data
  ULong_t run = 0; // 0: get from data???
  UShort_t sys = 0; // 0: get from data, 1: pp, 2: AA 
  UShort_t sNN = 0; // 0: get from data, otherwise center of mass energy (per nucleon pair)
  Short_t  fld = 0; // 0: get from data, otherwise L3 field in kG
  //AliAnalysisTask *taskFmd  = AddTaskForwardMult(mc, sys, sNN, fld);
  const Char_t* config = "$ALICE_PHYSICS/PWGPP/EVCHAR/FlowVectorCorrections/QnCorrectionsInterface/macros/ForwardAODConfig2.C";
  AliAnalysisTask *taskFmd  = AddTaskForwardMult(mc, run, sys, sNN, fld, config);



  // --- Make the output container and connect it --------------------
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
        AliAnalysisManager::kExchangeContainer,
        "Forward");
  //AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer *output = 
    mgr->CreateContainer("ForwardResultsP", TList::Class(), 
        AliAnalysisManager::kParamContainer, 
        "ForwardResultsP");
  //AliAnalysisManager::GetCommonFileName());

  mgr->ConnectInput(taskFmd, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskFmd, 1, histOut);

}

void AddFMD(AliAnalysisManager *mgr, AliQnCorrectionsManager* QnManager, Bool_t isESD){

  if(isESD) AddFMDTaskForESDanalysis(mgr);

  const Int_t gkFMDstep=0;
  Short_t FMDchannels[2][4000];
  for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<4000; ich++) FMDchannels[iv0][ich] = 0;

  for(Int_t ich=2000; ich<4000; ich++) FMDchannels[0][ich] = 1;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<2000; ich++) FMDchannels[1][ich] = 1;

  //-----------------------------------------------------------
  // Binning for Q-vector calibration

  const Int_t FMDdim              = 2;
  AliQnCorrectionsAxes * FMDbinning = new AliQnCorrectionsAxes(FMDdim);


  // vertex Z bins
  const Int_t nVtxZwidths = 3;   
  Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
  Int_t VtxZbins[nVtxZwidths] = {1,8,1};

  FMDbinning->SetAxis(0, AliQnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges,AliQnCorrectionsVarManager::VarName(AliQnCorrectionsVarManager::kVtxZ));
  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {100};

  FMDbinning->SetAxis(1, AliQnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges,AliQnCorrectionsVarManager::VarName(AliQnCorrectionsVarManager::kCentVZERO));

  //-----------------------------------------------------------
  // Binning for channel equalization
  AliQnCorrectionsAxes * FMDchannelbinning = new AliQnCorrectionsAxes(FMDdim+1); // +1 is for channel axis
  FMDchannelbinning->SetAxis(0, AliQnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges,AliQnCorrectionsVarManager::VarName(AliQnCorrectionsVarManager::kVtxZ));
  FMDchannelbinning->SetAxis(1, AliQnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges,AliQnCorrectionsVarManager::VarName(AliQnCorrectionsVarManager::kCentVZERO));
  FMDchannelbinning->SetNchannels(4000);

  ////////// end of binning


  AliQnCorrectionsConfiguration * FMDAconf = new AliQnCorrectionsConfiguration();
  FMDAconf->SetQnNormalization(1);
  FMDAconf->SetQnHarmonicsRange(1,4);
  FMDAconf->SetHarmonicForAlignment(2);
  //FMDAconf->SetDataVectorEqualizationMethod(1);
  FMDAconf->SetDataVectorIdList(new TArrayS(4000, FMDchannels[0]));
  FMDAconf->SetQnConfigurationName("FMDA");
  FMDAconf->SetReferenceQnForAlignment("TPC");
  FMDAconf->SetQnCorrectionsCommonAxes(FMDbinning);
  FMDAconf->SetDataVectorEqualizationAxes(FMDchannelbinning);
  FMDAconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run

  FMDAconf->SetQnCorrectionRecentering();
  FMDAconf->SetQnCorrectionAlignment();
  QnManager->AddQnConfiguration(FMDAconf, AliQnCorrectionsVarManager::kFMD);

  AliQnCorrectionsConfiguration * FMDCconf = new AliQnCorrectionsConfiguration();
  FMDCconf->SetQnNormalization(1);
  FMDCconf->SetQnHarmonicsRange(1,4);
  FMDCconf->SetHarmonicForAlignment(2);
  //FMDCconf->SetDataVectorEqualizationMethod(1);
  FMDCconf->SetDataVectorIdList(new TArrayS(4000, FMDchannels[1]));
  FMDCconf->SetQnConfigurationName("FMDC");
  FMDCconf->SetReferenceQnForAlignment("TPC");
  FMDCconf->SetQnCorrectionsCommonAxes(FMDbinning);
  FMDCconf->SetDataVectorEqualizationAxes(FMDchannelbinning);
  FMDCconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run

  FMDCconf->SetQnCorrectionRecentering();
  FMDCconf->SetQnCorrectionAlignment();
  QnManager->AddQnConfiguration(FMDCconf, AliQnCorrectionsVarManager::kFMD);





}



void AddRawFMD(AliQnCorrectionsManager* QnManager){


  const Int_t FMDdim              = 2;
  AliQnCorrectionsAxes * FMDbinning = new AliQnCorrectionsAxes(FMDdim);


  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {100};

  FMDbinning->SetAxis(0, AliQnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges,AliQnCorrectionsVarManager::VarName(AliQnCorrectionsVarManager::kVZEROMultPercentile));


  // vertex Z bins
  const Int_t nVtxZwidths = 3;   
  Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
  Int_t VtxZbins[nVtxZwidths] = {1,8,1};

  FMDbinning->SetAxis(1, AliQnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges,AliQnCorrectionsVarManager::VarName(AliQnCorrectionsVarManager::kVtxZ));

  //-----------------------------------------------------------
  // Binning for channel equalization
  AliQnCorrectionsAxes * FMDchannelbinning = new AliQnCorrectionsAxes(FMDdim+1); // +1 is for channel axis
  FMDchannelbinning->SetAxis(0, AliQnCorrectionsVarManager::kVZEROMultPercentile, nCtwidths, Ctbins, Ctedges,AliQnCorrectionsVarManager::VarName(AliQnCorrectionsVarManager::kVZEROMultPercentile));
  FMDchannelbinning->SetAxis(1, AliQnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges,AliQnCorrectionsVarManager::VarName(AliQnCorrectionsVarManager::kVtxZ));
  FMDchannelbinning->SetNchannels(52000);

  ////////// end of binning


  AliQnCorrectionsCuts *cutFMDA = new AliQnCorrectionsCuts();
  cutFMDA->AddCut(AliQnCorrectionsVarManager::kFMDEta,0.0,6.0);

  AliQnCorrectionsCuts *cutFMDC = new AliQnCorrectionsCuts();
  cutFMDC->AddCut(AliQnCorrectionsVarManager::kFMDEta,-6.0,0.0);

  AliQnCorrectionsConfiguration * FMDAconf = new AliQnCorrectionsConfiguration();
  FMDAconf->SetQnNormalization(1);
  FMDAconf->SetQnHarmonicsRange(1,4);
  //FMDAconf->SetDataVectorEqualizationMethod(1);
  //FMDAconf->SetDataVectorIdList(new TArrayS(4000, FMDchannels[0]));
  FMDAconf->SetQnConfigurationName("FMDAraw");
  FMDAconf->SetQnCorrectionsCommonAxes(FMDbinning);
  //FMDAconf->SetDataVectorEqualizationAxes(FMDchannelbinning);
  FMDAconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run
  FMDAconf->SetDataVectorCuts(cutFMDA);

  FMDAconf->SetQnCorrectionRecentering();
  QnManager->AddQnConfiguration(FMDAconf, AliQnCorrectionsVarManager::kFMDraw);

  AliQnCorrectionsConfiguration * FMDCconf = new AliQnCorrectionsConfiguration();
  FMDCconf->SetQnNormalization(1);
  FMDCconf->SetQnHarmonicsRange(1,4);
  //FMDCconf->SetDataVectorEqualizationMethod(1);
  //FMDCconf->SetDataVectorIdList(new TArrayS(4000, FMDchannels[0]));
  FMDCconf->SetQnConfigurationName("FMDCraw");
  FMDCconf->SetQnCorrectionsCommonAxes(FMDbinning);
  //FMDCconf->SetDataVectorEqualizationAxes(FMDchannelbinning);
  FMDCconf->UseCalibrationDirectoryNameAllEvents(kFALSE);   // in the current analysis task, label is run number, so with this setting corrections are made run-by-run
  FMDCconf->SetDataVectorCuts(cutFMDC);

  FMDCconf->SetQnCorrectionRecentering();
  QnManager->AddQnConfiguration(FMDCconf, AliQnCorrectionsVarManager::kFMDraw);





}





//__________________________________________________________________
void DefineHistograms(AliQnCorrectionsManager* QnManager, AliQnCorrectionsHistos* histos, TString histClass) {
  //
  // define the histograms
  //
  //#define VAR AliQnCorrectionsVarManager

  histClass+= "TrackQA_NoCuts;";
  for(Int_t iconf=0; iconf<QnManager->GetNumberOfQnConfigurations(); iconf++){
  AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) QnManager->GetQnConfiguration(iconf);
   if(!QnConf) continue;
   if(QnConf->IsTracking()) if(QnConf->QnConfigurationName().Contains("TPC")) histClass+= "TrackQA_"+QnConf->QnConfigurationName()+";";
   if(QnConf->IsTracking()) if(QnConf->QnConfigurationName().Contains("SPD")) histClass+= "TrackletQA_"+QnConf->QnConfigurationName()+";";
  }

  const Char_t* histClasses = histClass.Data();

  cout << "Defining histograms ..." << endl;
  cout << "histogram classes: " << histClass<< endl;

  //fHistosFile=new TFile(output,"RECREATE");

  TString classesStr(histClasses);
  TObjArray* arr=classesStr.Tokenize(";");

  const Int_t kNRunBins = 3000;
  Double_t runHistRange[2] = {137000.,140000.};

  for(Int_t iclass=0; iclass<arr->GetEntries(); ++iclass) {
    TString classStr = arr->At(iclass)->GetName();
    cout << "hist class: " << classStr.Data() << endl;

    // Event wise histograms
    if(classStr.Contains("Event")) {
      histos->AddHistClass(classStr.Data());
      histos->AddHistogram(classStr.Data(),"RunNo","Run numbers;Run", kFALSE, kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo);
      histos->AddHistogram(classStr.Data(),"BC","Bunch crossing;BC", kFALSE,3000,0.,3000.,AliQnCorrectionsVarManager::kBC);
      histos->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag;;", kFALSE,
          2,-0.5,1.5,AliQnCorrectionsVarManager::kIsPhysicsSelection, 0,0.0,0.0,AliQnCorrectionsVarManager::kNothing, 0,0.0,0.0,AliQnCorrectionsVarManager::kNothing, "off;on");

      histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,0.0,0.0,AliQnCorrectionsVarManager::kVtxZ);
      //histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-15.,15.,AliQnCorrectionsVarManager::kVtxZ);
      histos->AddHistogram(classStr.Data(),"VtxX","Vtx X;vtx X (cm)", kFALSE,300,-1.,1.,AliQnCorrectionsVarManager::kVtxX);
      histos->AddHistogram(classStr.Data(),"VtxY","Vtx Y;vtx Y (cm)", kFALSE,300,-1.,1.,AliQnCorrectionsVarManager::kVtxY);


      histos->AddHistogram(classStr.Data(),"CentVZEROvsMultPVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kVZEROMultPercentile, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentVZEROvsCentSPD","Centrality(VZERO);centrality VZERO (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentSPD","Centrality(TPC);centrality TPC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentVZERO","Centrality(TPC);centrality TPC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentZDC","Centrality(TPC);centrality TPC (percents);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentVZERO","Centrality(ZDC);centrality ZDC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentSPD","Centrality(ZDC);centrality ZDC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);

      histos->AddHistogram(classStr.Data(),"MultVZEROvsCentVZERO","Multiplicity;multiplicity VZERO;VZERO centrality", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManager::kVZEROTotalMult, 100,0.,100., AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"MultSPDvsCentPSD","Multiplicity;SPD tracklets;SPD centrality", kFALSE,
          100, 0.0, 3000.0, AliQnCorrectionsVarManager::kSPDntracklets, 100,0.,100., AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"MultTPCvsCentSPD","Multiplicity;TPC selected tracks;TPC centrality", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManager::kNtracksSelected, 100,0.,100., AliQnCorrectionsVarManager::kCentTPC);
      histos->AddHistogram(classStr.Data(),"MultZDCvsCentZDC","Multiplicity;multiplicity ZDC;ZDC centrality", kFALSE,
          100, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy, 100,0.,100., AliQnCorrectionsVarManager::kCentZDC);


      histos->AddHistogram(classStr.Data(),"MultTPCvsMultVZERO","Multiplicity;tracks TPC;multiplicity VZERO", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManager::kNtracksSelected, 100, 0.0, 32000.0, AliQnCorrectionsVarManager::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultSPD","Multiplicity;tracklets SPD;tracks TPC", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManager::kNtracksSelected, 100, 0.0, 3000.0, AliQnCorrectionsVarManager::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultVZERO","Multiplicity;tracklets SPD;multiplicity VZERO", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManager::kVZEROTotalMult, 100, 0.0, 3000.0, AliQnCorrectionsVarManager::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultZDC","Multiplicity;tracks TPC;energy ZDC", kFALSE,
          100, 0.0, 3500.0, AliQnCorrectionsVarManager::kNtracksSelected, 100, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultVZEROvsMultZDC","Multiplicity;multiplicity VZERO;energy ZDC", kFALSE,
          100, 0.0, 32000.0, AliQnCorrectionsVarManager::kVZEROTotalMult, 100, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultZDC","Multiplicity;tracklets SPD;energy ZDC", kFALSE,
          100, 0.0, 3000.0, AliQnCorrectionsVarManager::kSPDntracklets, 100, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy);



      histos->AddHistogram(classStr.Data(),"MultVZERO","Multiplicity;multiplicity VZERO", kFALSE,
          320, 0.0, 25000.0, AliQnCorrectionsVarManager::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultVZEROA","Multiplicity;multiplicity VZEROA", kFALSE,
          250, 0.0, 0.0, AliQnCorrectionsVarManager::kVZEROATotalMult);//10000.0
      histos->AddHistogram(classStr.Data(),"MultVZEROC","Multiplicity;multiplicity VZEROC", kFALSE,
          250, 0.0, 0.0, AliQnCorrectionsVarManager::kVZEROCTotalMult);//15000.0
      histos->AddHistogram(classStr.Data(),"MultZDC","Multiplicity;multiplicity ZDC", kFALSE,
          200, 0.0, 300000.0, AliQnCorrectionsVarManager::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCA","Multiplicity;multiplicity ZDCA", kFALSE,
          200, 0.0, 150000.0, AliQnCorrectionsVarManager::kZDCATotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCC","Multiplicity;multiplicity ZDCC", kFALSE,
          200, 0.0, 150000.0, AliQnCorrectionsVarManager::kZDCCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultFMD1","Multiplicity;multiplicity FMD1", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD1TotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2I","Multiplicity;multiplicity FMD2I", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD2ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2O","Multiplicity;multiplicity FMD2O", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD2OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3I","Multiplicity;multiplicity FMD3I", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD3ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3O","Multiplicity;multiplicity FMD3O", kFALSE,
          300, 0.0, 10000.0, AliQnCorrectionsVarManager::kFMD3OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROA","Multiplicity;multiplicity TZEROA", kFALSE,
          300, 0.0, 3000.0, AliQnCorrectionsVarManager::kTZEROATotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROC","Multiplicity;multiplicity TZEROC", kFALSE,
          300, 0.0, 3000.0, AliQnCorrectionsVarManager::kTZEROCTotalMult);




      histos->AddHistogram(classStr.Data(),"MultPercentVZERO","Multiplicity percentile (VZERO);multiplicity VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kVZEROMultPercentile);
      histos->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC);centrality TPC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC);

      histos->AddHistogram(classStr.Data(),"CentQuality","Centrality quality;centrality quality", kFALSE,
          100, -50.5, 49.5, AliQnCorrectionsVarManager::kCentQuality);
      histos->AddHistogram(classStr.Data(),"CentVZERO_Run_prof","<Centrality(VZERO)> vs run;Run; centrality VZERO (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD_Run_prof","<Centrality(SPD)> vs run;Run; centrality SPD (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC_Run_prof","<Centrality(TPC)> vs run;Run; centrality TPC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC_Run_prof","<Centrality(ZDC)> vs run;Run; centrality ZDC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC);


      histos->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event;# pairs", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManager::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event;# pairs", kFALSE,
          1000,0.,10000.,AliQnCorrectionsVarManager::kNV0selected);
      histos->AddHistogram(classStr.Data(),"NPairs","Number of candidates per event;# pairs", kFALSE,
          5000,0.,5000.,AliQnCorrectionsVarManager::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NPairsSelected", "Number of selected pairs per event; #pairs", kFALSE,
          5000,0.,5000.,AliQnCorrectionsVarManager::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event;# tracks", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManager::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event;# tracks", kFALSE,
          1000,0.,30000.,AliQnCorrectionsVarManager::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets; tracklets", kFALSE,
          3000, -0.5, 2999.5, AliQnCorrectionsVarManager::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"SPDnSingleClusters", "SPD #single clusters; tracklets", kFALSE,
          3000, -0.5, 2999.5, AliQnCorrectionsVarManager::kSPDnSingleClusters);

      histos->AddHistogram(classStr.Data(),"NV0total_Run_prof", "<Number of total V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0selected_Run_prof", "<Number of selected V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNV0selected);
      histos->AddHistogram(classStr.Data(),"Ndielectrons_Run_prof", "<Number of dielectrons> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NpairsSelected_Run_prof", "<Number of selected pairs> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected_Run_prof", "<Number of selected tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 100, 0., 10000., AliQnCorrectionsVarManager::kSPDntracklets);

      histos->AddHistogram(classStr.Data(),"VtxZ_CentVZERO","Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManager::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentSPD","Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManager::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentTPC","Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)", kFALSE,
          300,-15.,15.,AliQnCorrectionsVarManager::kVtxZ, 100, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC);
      continue;
    }  // end if className contains "Event"    


    // Offline trigger histograms
    if(classStr.Contains("OfflineTriggers")) {
      histos->AddHistClass(classStr.Data());

      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += Form("%s",AliQnCorrectionsVarManager::fOfflineTriggerNames[i]); triggerNames+=";";}

      histos->AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTrigger, 2, -0.5, 1.5, AliQnCorrectionsVarManager::kOfflineTriggerFired, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data(), "off;on");
      histos->AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentTPC_Triggers2", "Offline triggers fired vs centrality TPC; ; centrality TPC;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentTPC, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentSPD_Triggers2", "Offline triggers fired vs centrality SPD; ; centrality SPD;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentSPD, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentZDC_Triggers2", "Offline triggers fired vs centrality ZDC; ; centrality ZDC;", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentZDC, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);", kFALSE,
          64, -0.5, 63.5, AliQnCorrectionsVarManager::kOfflineTriggerFired2, 200, -20.0, 20.0, AliQnCorrectionsVarManager::kVtxZ, 0, 0.0, 0.0, AliQnCorrectionsVarManager::kNothing, triggerNames.Data());
      continue;
    }

    // Track histograms
    if(classStr.Contains("Tracks")) {
      histos->AddHistClass(classStr.Data());
      for(Int_t ih=0; ih<6; ++ih) {
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManager::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kSinNPhi+ih);
      }
    }


    // Track histograms
    if(classStr.Contains("TrackQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE,
          1000, 0.0, 50.0, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -1.5, 1.5, AliQnCorrectionsVarManager::kEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          1000, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE,
          1000, -10.0, 10.0, AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE,
          1000, -10.0, 10.0, AliQnCorrectionsVarManager::kDcaZ);
      histos->AddHistogram(classStr.Data(), "TPCncls", "TPCncls; TPCncls", kFALSE,
          160, 0.0, 160.0, AliQnCorrectionsVarManager::kTPCncls);

      // run dependence
      histos->AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, 0.0, 50.0, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, -1.5, 1.5, AliQnCorrectionsVarManager::kEta);      
      histos->AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi);      
      histos->AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, -10.0, 10.0, AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], AliQnCorrectionsVarManager::kRunNo, 1000, -10.0, 10.0, AliQnCorrectionsVarManager::kDcaZ);

      // correlations between parameters
      histos->AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE,
          300, -1.5, +1.5, AliQnCorrectionsVarManager::kEta, 100, 0.0, 10.0, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt", "p_{T} vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kFALSE,
          300, -0.01, 6.3, AliQnCorrectionsVarManager::kPhi, 100, 0.0, 2.2, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kTRUE,
          300, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi, 100, 0.0, 10.0, AliQnCorrectionsVarManager::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -1.0, +1.0, AliQnCorrectionsVarManager::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi);
      histos->AddHistogram(classStr.Data(), "TPCncls_Eta_Phi_prof", "<TPC ncls> vs #varphi vs #eta; #eta; #varphi (rad.);TPC ncls", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManager::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManager::kTPCncls);
      histos->AddHistogram(classStr.Data(), "DCAxy_Eta_Phi_prof", "<DCAxy> vs #varphi vs #eta; #eta; #varphi (rad.);DCAxy (cm)", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManager::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Eta_Phi_prof", "<DCAz> vs #varphi vs #eta; #eta; #varphi (rad.);DCAz (cm)", kTRUE,
          200, -1.0, +1.0, AliQnCorrectionsVarManager::kEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kPhi, 10, 0.0, 200., AliQnCorrectionsVarManager::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
          100, 0.0, 10.0, AliQnCorrectionsVarManager::kPt, 500, -2.0, 2.0, AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE,
          100, 0.0, 10.0, AliQnCorrectionsVarManager::kPt, 500, -2.0, 2.0, AliQnCorrectionsVarManager::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE,
          100, -1.0, 1.0, AliQnCorrectionsVarManager::kEta, 500, -2.0, 2.0, AliQnCorrectionsVarManager::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE,
          100, -1.0, 1.0, AliQnCorrectionsVarManager::kEta, 500, -2.0, 2.0, AliQnCorrectionsVarManager::kDcaZ);

      for(Int_t ih=0; ih<6; ++ih) {
        //histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManager::kCosNPhi+ih);
        //histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_Pt_Eta",ih+1), Form("<cos%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, AliQnCorrectionsVarManager::kEta, 30, 0.0, 3.0, AliQnCorrectionsVarManager::kPt, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_Pt_Eta",ih+1), Form("<sin%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, AliQnCorrectionsVarManager::kEta, 30, 0.0, 3.0, AliQnCorrectionsVarManager::kPt, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO_VtxZ",ih+1), Form("<cos%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, AliQnCorrectionsVarManager::kVtxZ, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1., 1., AliQnCorrectionsVarManager::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO_VtxZ",ih+1), Form("<sin%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, AliQnCorrectionsVarManager::kVtxZ, 20, 0.0, 100.0, AliQnCorrectionsVarManager::kCentVZERO, 500, -1.0, 1.0, AliQnCorrectionsVarManager::kSinNPhi+ih);
      }
    }

    // Tracklet histograms
    if(classStr.Contains("TrackletQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -3.0, 3.0, AliQnCorrectionsVarManager::kSPDtrackletEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          300, -0.01, 6.3, AliQnCorrectionsVarManager::kSPDtrackletPhi);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -3.0, +3.0, AliQnCorrectionsVarManager::kSPDtrackletEta, 100, 0.0, 6.3, AliQnCorrectionsVarManager::kSPDtrackletPhi);
    }

  }

  cout << " done" << endl;
}







