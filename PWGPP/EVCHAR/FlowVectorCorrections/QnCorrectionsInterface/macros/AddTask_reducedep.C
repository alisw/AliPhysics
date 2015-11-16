// event plane calibration task
// author: Jaap Onderaater, jacobus.onderwaater@cern.ch
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

void DefineHistograms(QnCorrectionsManager* EPmanager, QnCorrectionsHistos* histos, TString histClass);

void AddVZERO(QnCorrectionsManager* EPmanager);
void AddTPC(QnCorrectionsManager* EPmanager);
void AddTZERO(QnCorrectionsManager* EPmanager);
void AddFMD(QnCorrectionsManager* EPmanager);
void AddZDC(QnCorrectionsManager* EPmanager);




AliAnalysisTaskEventPlaneCalibration* AddTask_ep() {
  //gSystem->Load("libESD");
  //gSystem->Load("libPWGPPevcharEp.so");

  //TChain* chain = CreateESDChain("/u/jonderw/files.txt", 1, 0);
  //AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");
  //
  //AliVEventHandler* esdH = new AliESDInputHandler;
  //mgr->SetInputEventHandler(esdH);

  //AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  //if (!mgr) {
  //  Error("AddTask_jonderw_ep", "No analysis manager found.");
  //  return 0;
  //}

 //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
 //AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

 // Apply the event selection
 //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
 //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
  
  //Do we have an MC handler?
  //Bool_t hasMC=(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);
  //Bool_t isAOD=mgr->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  


////  AliAnalysisDataContainer *coutput1 =
////    mgr->CreateContainer("diele_defaultTree",
////                         TTree::Class(),
////                         AliAnalysisManager::kExchangeContainer,
////                         "diele_default");
  //
  //AliAnalysisDataContainer *cOutputHist1 =
  //  mgr->CreateContainer("qaHistos",
  //                       TList::Class(),
  //                       AliAnalysisManager::kOutputContainer,
  //                       "dst_Histos.root");

  //AliAnalysisDataContainer *cOutputHist2 =
  //  mgr->CreateContainer("dstTree",
  //                       TTree::Class(),
  //                       AliAnalysisManager::kOutputContainer,
  //                       "dstTree.root");
 
 
  //AliAnalysisDataContainer *cReducedEvent =
  //  mgr->CreateContainer("ReducedEvent",
  //                       AliReducedEvent::Class(),
  //                       AliAnalysisManager::kExchangeContainer,
  //                       "reducedEvent");


 
  QnCorrectionsManager* EPmanager = new QnCorrectionsManager();
  AliAnalysisTaskEventPlaneCalibration* taskEP = new AliAnalysisTaskEventPlaneCalibration("EventPlaneCalibration");

  QnCorrectionsCuts *eventCuts = new QnCorrectionsCuts();
  eventCuts->AddCut(QnCorrectionsVarManager::kVtxZ,-10.0,10.0);
  eventCuts->AddCut(QnCorrectionsVarManager::kCentVZERO,0.0,90.0);
  eventCuts->AddCut(QnCorrectionsVarManager::kIsPhysicsSelection,0.5,1.5);

  taskEP->SetEventCuts(eventCuts);
  //taskEP->SelectCollisionCandidates(AliVEvent::kMB);

  AddVZERO(EPmanager);
  AddTPC(EPmanager);
  //AddTZERO(EPmanager);
  //AddZDC(EPmanager);
  //AddFMD(EPmanager);

  taskEP->SetEventPlaneManager(EPmanager);
  taskEP->DefineInOutput();


  QnCorrectionsHistos* hists = taskEP->GetHistograms();
  TString histClass = "";
  histClass += "Event_NoCuts;";
  histClass += "Event_Analysis;";
  DefineHistograms(EPmanager, hists, histClass);


  //create output container
  //  AliAnalysisDataContainer *cOutputHist =
  //  mgr->CreateContainer("CalibrationHistos",
  //                       TList::Class(),
  //                       AliAnalysisManager::kOutputContainer,
  //                       "CalibrationHistograms.root");
  //  AliAnalysisDataContainer *cOutputQvec =
  //  mgr->CreateContainer("CalibratedQvector",
  //                       TTree::Class(),
  //                       AliAnalysisManager::kOutputContainer,
  //                       "Qvectors.root");
  //  AliAnalysisDataContainer *cOutputHistQA =
  //  mgr->CreateContainer("CalibrationQA",
  //                       TList::Class(),
  //                       AliAnalysisManager::kOutputContainer,
  //                       "CalibrationQA.root");


  //AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  //cout<<cinput1<<endl;

  //mgr->ConnectInput(taskDst,  0, cinput1);
  ////mgr->ConnectOutput(taskDst, 0, coutput1 );
  //mgr->ConnectOutput(taskDst, 1, cOutputHist1);
  //mgr->ConnectOutput(taskDst, 2, cOutputHist2);
  //mgr->ConnectOutput(taskDst, 3, cReducedEvent);

  //mgr->AddTask(taskEP);
  //mgr->ConnectInput(taskEP,  0, cReducedEvent);
  //mgr->ConnectInput(taskEP,  0, cinput1);
  //mgr->ConnectInput(taskEP,  0, mgr->GetCommonInputContainer());
  //mgr->ConnectOutput(taskEP, 1, cOutputHist );
  //mgr->ConnectOutput(taskEP, 2, cOutputQvec );
  //mgr->ConnectOutput(taskEP, 3, cOutputHistQA );


  return taskEP;
}




void AddVZERO(QnCorrectionsManager* EPmanager){

  const Char_t* gkVzeroEqualizationPath    = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
  const Char_t* gkVzeroRecenteringPath     = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
  const Char_t* gkVzeroDiagonalizationPath = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
  const Int_t gkVZEROstep=0;
  Short_t VZEROchannels[4][64];
  for(Int_t iv0=0; iv0<4; iv0++) for(Int_t ich=0; ich<64; ich++) VZEROchannels[iv0][ich] = 0;
  
  for(Int_t ich=32; ich<64; ich++) VZEROchannels[0][ich] = 1;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<32; ich++) VZEROchannels[1][ich] = 1;

  //-----------------------------------------------------------
  // Binning for Q-vector calibration

  const Int_t VZEROdim              = 2;
  QnCorrectionsAxes * VZERObinning = new QnCorrectionsAxes(VZEROdim);
   

  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {100};

  VZERObinning->AddAxis(0, QnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

    
  // vertex Z bins
  const Int_t nVtxZwidths = 3;   
  Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
  Int_t VtxZbins[nVtxZwidths] = {1,8,1};

  VZERObinning->AddAxis(1, QnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);
    
  //-----------------------------------------------------------
  // Binning for channel equalization
  QnCorrectionsAxes * VZEROchannelbinning = new QnCorrectionsAxes(VZEROdim+1); // +1 is for channel axis
  VZEROchannelbinning->AddAxis(0, QnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);
  VZEROchannelbinning->AddAxis(1, QnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);
  VZEROchannelbinning->SetNchannels(64);
  
  ////////// end of binning
  
 
  QnCorrectionsConfiguration * VZEROAconf = new QnCorrectionsConfiguration();
  VZEROAconf->SetCalibrationStep(gkVZEROstep);
  VZEROAconf->SetQnNormalization(1);
  VZEROAconf->SetDataVectorEqualizationMethod(1);
  VZEROAconf->SetTwistAndRescalingMethod(2);
  VZEROAconf->SetDataVectorIdList(new TArrayS(64, VZEROchannels[0]));
  VZEROAconf->SetCalibrationDetectorName("VZEROA");
  VZEROAconf->SetEqualizationDetectorName("VZEROA");
  VZEROAconf->SetQnConfigurationName("VZEROA");
  VZEROAconf->SetReferenceQnForTwistAndRescaling("TPC","VZEROC");
  VZEROAconf->SetEqualizationHistPath(gkVzeroEqualizationPath);
  VZEROAconf->SetRecenteringHistPath(gkVzeroRecenteringPath);
  VZEROAconf->SetCorrelationHistPath(gkVzeroDiagonalizationPath);
  VZEROAconf->SetQnCorrectionsCommonAxes(VZERObinning);
  VZEROAconf->SetDataVectorEqualizationAxes(VZEROchannelbinning);

  EPmanager->AddQnConfiguration(VZEROAconf, QnCorrectionsVarManager::kVZERO);



  QnCorrectionsConfiguration * VZEROCconf = new QnCorrectionsConfiguration();
  VZEROCconf->SetCalibrationStep(gkVZEROstep);
  VZEROCconf->SetQnNormalization(1);
  VZEROCconf->SetDataVectorEqualizationMethod(1);
  VZEROCconf->SetTwistAndRescalingMethod(2);
  VZEROCconf->SetDataVectorIdList(new TArrayS(64, VZEROchannels[1]));
  VZEROCconf->SetCalibrationDetectorName("VZEROC");
  VZEROCconf->SetEqualizationDetectorName("VZEROC");
  VZEROCconf->SetQnConfigurationName("VZEROC");
  VZEROCconf->SetReferenceQnForTwistAndRescaling("TPC","VZEROA");
  VZEROCconf->SetEqualizationHistPath(gkVzeroEqualizationPath);
  VZEROCconf->SetRecenteringHistPath(gkVzeroRecenteringPath);
  VZEROCconf->SetCorrelationHistPath(gkVzeroDiagonalizationPath);
  VZEROCconf->SetQnCorrectionsCommonAxes(VZERObinning);
  VZEROCconf->SetDataVectorEqualizationAxes(VZEROchannelbinning);
  EPmanager->AddQnConfiguration(VZEROCconf, QnCorrectionsVarManager::kVZERO);
  
}
 
  void AddTPC(QnCorrectionsManager* EPmanager){

     /////////////// Add TPC subdetectors ///////////////////
        
const Int_t gkTPCstep=0;
const Char_t* gkTpcRecenteringPath       = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
const Char_t* gkTpcDiagonalizationPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";



    //-----------------------------------------------------------
    // Binning for Q-vector calibration
    //
    const Int_t TPCdim              = 2;    // number of calibration parameters
    QnCorrectionsAxes * TPCbinning = new QnCorrectionsAxes(TPCdim);  //  declare binning object
    
    // centrality bins
    const Int_t nCtwidths = 1;   
    Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
    Int_t Ctbins[nCtwidths] = {100};

    TPCbinning->AddAxis(0, QnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

      
    // vertex Z bins
    const Int_t nVtxZwidths = 3;   
    Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
    Int_t VtxZbins[nVtxZwidths] = {1,8,1};

    TPCbinning->AddAxis(1, QnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);

    ////////// end of binning
    


    QnCorrectionsConfiguration * TPCconf = new QnCorrectionsConfiguration();
    TPCconf->SetCalibrationStep(gkTPCstep);
    TPCconf->SetQnNormalization(1);
    TPCconf->SetTwistAndRescalingMethod(0);
    TPCconf->SetCalibrationDetectorName("TPC");
    TPCconf->SetQnConfigurationName("TPC");
    TPCconf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPCconf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPCconf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    //TPCconf->CreateQvectorHistograms(TPCdim, TPCvar, binLimits);
    TPCconf->SetQnCorrectionsCommonAxes(TPCbinning);
    TPCconf->SetTracking();   // this is a tracking detector
   
    QnCorrectionsCuts *trackTPC = new QnCorrectionsCuts();
    trackTPC->AddCut(QnCorrectionsVarManager::kDcaXY,-1.0,1.0);
    trackTPC->AddCut(QnCorrectionsVarManager::kDcaZ,-3.0,3.0);
    trackTPC->AddCut(QnCorrectionsVarManager::kEta,-0.8,0.8);
    //trackTPC->AddCut(QnCorrectionsVarManager::kEta,-0.3,0.3, kTRUE);  // exclusion region
    trackTPC->AddCut(QnCorrectionsVarManager::kPt,0.2,5.);
    trackTPC->AddCut(QnCorrectionsVarManager::kTPCnclsIter1,70.0,161.0);
    trackTPC->AddCut(QnCorrectionsVarManager::kTPCchi2Iter1,0.2,4.0);
    TPCconf->SetDataVectorCuts(trackTPC);

    //AliFlowCuts *trackTPC = new AliFlowCuts();
    //trackTPC = trackTPC->GetStandardTPCStandaloneCuts2010();
    //TPCconf->SetDataVectorCuts(trackTPC);

    EPmanager->AddQnConfiguration(TPCconf, QnCorrectionsVarManager::kTPC);


    return;


    QnCorrectionsConfiguration * TPCpConf = new QnCorrectionsConfiguration();
    TPCpConf->SetCalibrationStep(gkTPCstep);
    TPCpConf->SetQnNormalization(1);
    TPCpConf->SetTwistAndRescalingMethod(0);
    TPCpConf->SetCalibrationDetectorName("TPC");
    TPCpConf->SetQnConfigurationName("TPCp");
    TPCpConf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPCpConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPCpConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    TPCpConf->SetQnCorrectionsCommonAxes(TPCbinning);
   
    QnCorrectionsCuts *trackTPCp = new QnCorrectionsCuts();
    trackTPCp->CopyCuts(trackTPC);
    trackTPCp->AddCut(QnCorrectionsVarManager::kCharge, 0.5, 1.5);
    TPCpConf->SetDataVectorCuts(trackTPCp);

    QnCorrectionsConfiguration * TPCnConf = new QnCorrectionsConfiguration();
    TPCnConf->SetCalibrationStep(gkTPCstep);
    TPCnConf->SetQnNormalization(1);
    TPCnConf->SetTwistAndRescalingMethod(0);
    TPCnConf->SetCalibrationDetectorName("TPC");
    TPCnConf->SetQnConfigurationName("TPCn");
    TPCnConf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPCnConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPCnConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    TPCnConf->SetQnCorrectionsCommonAxes(TPCbinning);
   
    QnCorrectionsCuts *trackTPCn = new QnCorrectionsCuts();
    trackTPCn->CopyCuts(trackTPC);
    trackTPCn->AddCut(QnCorrectionsVarManager::kCharge, -1.5, -0.5);
    TPCnConf->SetDataVectorCuts(trackTPCn);

    QnCorrectionsConfiguration * TPClConf = new QnCorrectionsConfiguration();
    TPClConf->SetCalibrationStep(gkTPCstep);
    TPClConf->SetQnNormalization(1);
    TPClConf->SetTwistAndRescalingMethod(0);
    TPClConf->SetCalibrationDetectorName("TPC");
    TPClConf->SetQnConfigurationName("TPCl");
    TPClConf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPClConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPClConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    TPClConf->SetQnCorrectionsCommonAxes(TPCbinning);
   
    QnCorrectionsCuts *trackTPCl = new QnCorrectionsCuts();
    trackTPCl->CopyCuts(trackTPC);
    trackTPCl->AddCut(QnCorrectionsVarManager::kEta, -0.8, -0.4);
    TPClConf->SetDataVectorCuts(trackTPCl);


    QnCorrectionsConfiguration * TPCrConf = new QnCorrectionsConfiguration();
    TPCrConf->SetCalibrationStep(gkTPCstep);
    TPCrConf->SetQnNormalization(1);
    TPCrConf->SetTwistAndRescalingMethod(0);
    TPCrConf->SetCalibrationDetectorName("TPC");
    TPCrConf->SetQnConfigurationName("TPCr");
    TPCrConf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPCrConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPCrConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    TPCrConf->SetQnCorrectionsCommonAxes(TPCbinning);
   
    QnCorrectionsCuts *trackTPCr = new QnCorrectionsCuts();
    trackTPCr->CopyCuts(trackTPC);
    trackTPCr->AddCut(QnCorrectionsVarManager::kEta, 0.4, 0.8);
    TPCrConf->SetDataVectorCuts(trackTPCr);


    QnCorrectionsConfiguration * TPClpConf = new QnCorrectionsConfiguration();
    TPClpConf->SetCalibrationStep(gkTPCstep);
    TPClpConf->SetQnNormalization(1);
    TPClpConf->SetTwistAndRescalingMethod(0);
    TPClpConf->SetCalibrationDetectorName("TPC");
    TPClpConf->SetQnConfigurationName("TPClp");
    TPClpConf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPClpConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPClpConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    TPClpConf->SetQnCorrectionsCommonAxes(TPCbinning);
   
    QnCorrectionsCuts *trackTPClp = new QnCorrectionsCuts();
    trackTPClp->CopyCuts(trackTPC);
    trackTPClp->AddCut(QnCorrectionsVarManager::kCharge, 0.5, 1.5);
    trackTPClp->AddCut(QnCorrectionsVarManager::kEta, -0.8, -0.4);
    TPClpConf->SetDataVectorCuts(trackTPClp);


    QnCorrectionsConfiguration * TPClnConf = new QnCorrectionsConfiguration();
    TPClnConf->SetCalibrationStep(gkTPCstep);
    TPClnConf->SetQnNormalization(1);
    TPClnConf->SetTwistAndRescalingMethod(0);
    TPClnConf->SetCalibrationDetectorName("TPC");
    TPClnConf->SetQnConfigurationName("TPCln");
    TPClnConf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPClnConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPClnConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    TPClnConf->SetQnCorrectionsCommonAxes(TPCbinning);
   
    QnCorrectionsCuts *trackTPCln = new QnCorrectionsCuts();
    trackTPCln->CopyCuts(trackTPC);
    trackTPCln->AddCut(QnCorrectionsVarManager::kCharge, -1.5, -0.5);
    trackTPCln->AddCut(QnCorrectionsVarManager::kEta, -0.8, -0.4);
    TPClnConf->SetDataVectorCuts(trackTPCln);



    QnCorrectionsConfiguration * TPCrpConf = new QnCorrectionsConfiguration();
    TPCrpConf->SetCalibrationStep(gkTPCstep);
    TPCrpConf->SetQnNormalization(1);
    TPCrpConf->SetTwistAndRescalingMethod(0);
    TPCrpConf->SetCalibrationDetectorName("TPC");
    TPCrpConf->SetQnConfigurationName("TPCrp");
    TPCrpConf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPCrpConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPCrpConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    TPCrpConf->SetQnCorrectionsCommonAxes(TPCbinning);
   
    QnCorrectionsCuts *trackTPCrp = new QnCorrectionsCuts();
    trackTPCrp->CopyCuts(trackTPC);
    trackTPCrp->AddCut(QnCorrectionsVarManager::kCharge, 0.5, 1.5);
    trackTPCrp->AddCut(QnCorrectionsVarManager::kEta, 0.4, 0.8);
    TPCrpConf->SetDataVectorCuts(trackTPCrp);


    QnCorrectionsConfiguration * TPCrnConf = new QnCorrectionsConfiguration();
    TPCrnConf->SetCalibrationStep(gkTPCstep);
    TPCrnConf->SetQnNormalization(1);
    TPCrnConf->SetTwistAndRescalingMethod(0);
    TPCrnConf->SetCalibrationDetectorName("TPC");
    TPCrnConf->SetQnConfigurationName("TPCrn");
    TPCrnConf->SetReferenceQnForTwistAndRescaling("VZEROC","VZEROA");
    TPCrnConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPCrnConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    TPCrnConf->SetQnCorrectionsCommonAxes(TPCbinning);
   
    QnCorrectionsCuts *trackTPCrn = new QnCorrectionsCuts();
    trackTPCrn->CopyCuts(trackTPC);
    trackTPCrn->AddCut(QnCorrectionsVarManager::kCharge, -1.5, -0.5);
    trackTPCrn->AddCut(QnCorrectionsVarManager::kEta, 0.4, 0.8);
    TPCrnConf->SetDataVectorCuts(trackTPCrn);


    //EPmanager->AddQnConfiguration(TPClConf, QnCorrectionsVarManager::kTPC);
    EPmanager->AddQnConfiguration(TPCrConf, QnCorrectionsVarManager::kTPC);
    //EPmanager->AddQnConfiguration(TPClpConf, QnCorrectionsVarManager::kTPC);
    //EPmanager->AddQnConfiguration(TPClnConf, QnCorrectionsVarManager::kTPC);
    EPmanager->AddQnConfiguration(TPCrpConf, QnCorrectionsVarManager::kTPC);
    EPmanager->AddQnConfiguration(TPCrnConf, QnCorrectionsVarManager::kTPC);





    //QnCorrectionsConfiguration * TPCposConf = new QnCorrectionsConfiguration();
    //TPCposConf->SetCalibrationStep(gkTPCstep);
    //TPCposConf->SetQnNormalization(1);
    //TPCposConf->SetTwistAndRescalingMethod(0);
    //TPCposConf->SetCalibrationDetectorName("TPC");
    //TPCposConf->SetQnConfigurationName("TPCpos");
    //TPCposConf->SetReferenceQnForTwistAndRescaling("TPCneg","VZEROA");
    //TPCposConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    //TPCposConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
   
    //QnCorrectionsCuts *trackTPCpos = new QnCorrectionsCuts();
    //trackTPCpos->CopyCuts(trackTPC);
    //trackTPCpos->AddCut(kCharge, 0.5, 1.5);
    //TPCposConf->SetDataVectorCuts(trackTPCpos);

    //EPmanager->AddQnConfiguration(TPCposConf);

    //QnCorrectionsConfiguration * TPCnegConf = new QnCorrectionsConfiguration();
    //TPCnegConf->SetCalibrationStep(gkTPCstep);
    //TPCnegConf->SetQnNormalization(1);
    //TPCnegConf->SetTwistAndRescalingMethod(0);
    //TPCnegConf->SetCalibrationDetectorName("TPC");
    //TPCnegConf->SetQnConfigurationName("TPCneg");
    //TPCnegConf->SetReferenceQnForTwistAndRescaling("TPCpos","VZEROA");
    //TPCnegConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    //TPCnegConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
   
    //QnCorrectionsCuts *trackTPCneg = new QnCorrectionsCuts();
    //trackTPCneg->CopyCuts(trackTPC);
    //trackTPCneg->AddCut(kCharge, -0.5, -1.5);
    //TPCnegConf->SetDataVectorCuts(trackTPCneg);

    //EPmanager->AddQnConfiguration(TPCnegConf);


    

    //QnCorrectionsConfiguration * TPCrightConf = new QnCorrectionsConfiguration();
    //TPCrightConf->SetCalibrationStep(gkTPCstep);
    //TPCrightConf->SetQnNormalization(1);
    //TPCrightConf->SetTwistAndRescalingMethod(0);
    //TPCrightConf->SetCalibrationDetectorName("TPC");
    //TPCrightConf->SetQnConfigurationName("TPCright");
    //TPCrightConf->SetReferenceQnForTwistAndRescaling("TPCleft","VZEROC");
    //TPCrightConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    //TPCrightConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    //TPCrightConf->SetQnCorrectionsCommonAxes(TPCbinning);
    //TPCrightconf->CreateQvectorHistograms(TPCdim, TPCvar, binLimits);
   
    //QnCorrectionsCuts *trackTPCright = new QnCorrectionsCuts();
    //trackTPCright->CopyCuts(trackTPC);
    //trackTPCright->AddCut(QnCorrectionsVarManager::kEta, 0.0, 1.0);
    //TPCrightConf->SetDataVectorCuts(trackTPCright);

    //EPmanager->AddQnConfiguration(TPCrightConf, QnCorrectionsVarManager::kTPC);

    //QnCorrectionsConfiguration * TPCleftConf = new QnCorrectionsConfiguration();
    //TPCleftConf->SetCalibrationStep(gkTPCstep);
    //TPCleftConf->SetQnNormalization(1);
    //TPCleftConf->SetTwistAndRescalingMethod(0);
    //TPCleftConf->SetCalibrationDetectorName("TPC");
    //TPCleftConf->SetQnConfigurationName("TPCleft");
    //TPCleftConf->SetReferenceQnForTwistAndRescaling("TPCright","VZEROA");
    //TPCleftConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    //TPCleftConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    //TPCleftConf->SetQnCorrectionsCommonAxes(TPCbinning);
    //TPCleftConf->CreateQvectorHistograms(TPCdim, TPCvar, binLimits);
   
    //QnCorrectionsCuts *trackTPCleft = new QnCorrectionsCuts();
    //trackTPCleft->CopyCuts(trackTPC);
    //trackTPCleft->AddCut(QnCorrectionsVarManager::kEta, -1.0, 0);
    //TPCleftConf->SetDataVectorCuts(trackTPCleft);

    //EPmanager->AddQnConfiguration(TPCleftConf, QnCorrectionsVarManager::kTPC);




  } 
   void AddTZERO(QnCorrectionsManager* EPmanager){
    
    /////////////// Add TZERO subdetectors ///////////////////
const Int_t gkTZEROstep=0;
const Char_t* gkTzeroEqualizationPath      = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
const Char_t* gkTzeroRecenteringPath       = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
const Char_t* gkTzeroDiagonalizationPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";

    

    Short_t TZEROchannels[2][24];
    for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<24; ich++) TZEROchannels[iv0][ich] = 0;
    
    for(Int_t ich=12; ich<24; ich++) TZEROchannels[0][ich] = 1;
    for(Int_t ich=0; ich<12; ich++)  TZEROchannels[1][ich] = 1;
    TZEROchannels[1][8] = 0; // Channel 8 in 2010 PbPb is not ok


  //-----------------------------------------------------------
  // Binning for Q-vector calibration

  const Int_t TZEROdim              = 2;
  QnCorrectionsAxes * TZERObinning = new QnCorrectionsAxes(TZEROdim);
   

  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {100};

  TZERObinning->AddAxis(0, QnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

    
  // vertex Z bins
  const Int_t nVtxZwidths = 3;   
  Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
  Int_t VtxZbins[nVtxZwidths] = {1,8,1};

  TZERObinning->AddAxis(1, QnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);
    
  //-----------------------------------------------------------
  // Binning for channel equalization
  QnCorrectionsAxes * TZEROchannelbinning = new QnCorrectionsAxes(TZEROdim+1); // +1 is for channel axis
  TZEROchannelbinning->AddAxis(0, QnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);
  TZEROchannelbinning->AddAxis(1, QnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);
  TZEROchannelbinning->SetNchannels(24);
  
  ////////// end of binning
  
 
      QnCorrectionsConfiguration * TZEROAconf = new QnCorrectionsConfiguration();
      TZEROAconf->SetCalibrationStep(gkTZEROstep);
      TZEROAconf->SetQnNormalization(1);
      TZEROAconf->SetDataVectorEqualizationMethod(0);
      TZEROAconf->SetTwistAndRescalingMethod(2);
      TZEROAconf->SetDataVectorIdList(new TArrayS(24, TZEROchannels[0]));
      TZEROAconf->SetCalibrationDetectorName("TZEROA");
      TZEROAconf->SetEqualizationDetectorName("TZEROA");
      TZEROAconf->SetQnConfigurationName("TZEROA");
      TZEROAconf->SetReferenceQnForTwistAndRescaling("TPC","TZEROC");
      TZEROAconf->SetEqualizationHistPath(gkTzeroEqualizationPath);
      TZEROAconf->SetRecenteringHistPath(gkTzeroRecenteringPath);
      TZEROAconf->SetCorrelationHistPath(gkTzeroDiagonalizationPath);
      TZEROAconf->SetQnCorrectionsCommonAxes(TZERObinning);
      TZEROAconf->SetDataVectorEqualizationAxes(TZEROchannelbinning);
      EPmanager->AddQnConfiguration(TZEROAconf, QnCorrectionsVarManager::kTZERO);

      QnCorrectionsConfiguration * TZEROCconf = new QnCorrectionsConfiguration();
      TZEROCconf->SetCalibrationStep(gkTZEROstep);
      TZEROCconf->SetQnNormalization(1);
      TZEROCconf->SetDataVectorEqualizationMethod(0);
      TZEROCconf->SetTwistAndRescalingMethod(2);
      TZEROCconf->SetDataVectorIdList(new TArrayS(24, TZEROchannels[1]));
      TZEROCconf->SetCalibrationDetectorName("TZEROC");
      TZEROCconf->SetEqualizationDetectorName("TZEROC");
      TZEROCconf->SetQnConfigurationName("TZEROC");
      TZEROCconf->SetReferenceQnForTwistAndRescaling("TPC","TZEROA");
      TZEROCconf->SetEqualizationHistPath(gkTzeroEqualizationPath);
      TZEROCconf->SetRecenteringHistPath(gkTzeroRecenteringPath);
      TZEROCconf->SetCorrelationHistPath(gkTzeroDiagonalizationPath);
      TZEROCconf->SetQnCorrectionsCommonAxes(TZERObinning);
      TZEROCconf->SetDataVectorEqualizationAxes(TZEROchannelbinning);
      EPmanager->AddQnConfiguration(TZEROCconf, QnCorrectionsVarManager::kTZERO);

  }
  
 void AddZDC(QnCorrectionsManager* EPmanager){
    /////////////// Add ZDC subdetectors ///////////////////

const Int_t gkZDCstep=0;
const Char_t* gkZdcEqualizationPath      = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
const Char_t* gkZdcRecenteringPath       = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
const Char_t* gkZdcDiagonalizationPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";


  //-----------------------------------------------------------
  // Binning for Q-vector calibration

  const Int_t ZDCdim              = 3;
  QnCorrectionsAxes * ZDCbinning = new QnCorrectionsAxes(ZDCdim);
   

  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {100};

  ZDCbinning->AddAxis(0, QnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

    
  // vertex Z bins
  const Int_t nVtxXwidths = 1;   
  Double_t VtxXedges[nVtxXwidths+1] = {-0.1,0.05};
  Int_t VtxXbins[nVtxXwidths] = {5};
  ZDCbinning->AddAxis(1, QnCorrectionsVarManager::kVtxX, nVtxXwidths, VtxXbins, VtxXedges);
  const Int_t nVtxYwidths = 1;   
  Double_t VtxYedges[nVtxYwidths+1] = {0.11,0.25};
  Int_t VtxYbins[nVtxYwidths] = {5};
  ZDCbinning->AddAxis(2, QnCorrectionsVarManager::kVtxY, nVtxYwidths, VtxYbins, VtxYedges);
    
  //-----------------------------------------------------------
  // Binning for channel equalization
  QnCorrectionsAxes * ZDCchannelbinning = new QnCorrectionsAxes(ZDCdim+1); // +1 is for channel axis
  ZDCchannelbinning->AddAxis(0, QnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);
  ZDCchannelbinning->AddAxis(1, QnCorrectionsVarManager::kVtxX, nVtxXwidths, VtxXbins, VtxXedges);
  ZDCchannelbinning->AddAxis(2, QnCorrectionsVarManager::kVtxY, nVtxYwidths, VtxYbins, VtxYedges);
  ZDCchannelbinning->SetNchannels(10);
  
  ////////// end of binning

    Short_t ZDCchannels[2][10];
    for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<10; ich++) ZDCchannels[iv0][ich] = 0;
    
    for(Int_t ich=6; ich<10; ich++) ZDCchannels[0][ich] = 1;
    for(Int_t ich=1; ich<5; ich++) ZDCchannels[1][ich] = 1;

      QnCorrectionsConfiguration * ZDCAconf = new QnCorrectionsConfiguration();
      ZDCAconf->SetCalibrationStep(gkZDCstep);
      ZDCAconf->SetQnNormalization(1);
      //ZDCAconf->SetDataVectorEqualizationMethod(2);
      //ZDCAconf->SetTwistAndRescalingMethod(2);
      ZDCAconf->SetDataVectorIdList(new TArrayS(10, ZDCchannels[0]));
      ZDCAconf->SetCalibrationDetectorName("ZDCA");
      //ZDCAconf->SetEqualizationDetectorName("ZDCA");
      ZDCAconf->SetQnConfigurationName("ZDCA");
      ZDCAconf->SetReferenceQnForTwistAndRescaling("ZDCC","TPC");
      //ZDCAconf->SetEqualizationHistPath(gkZdcEqualizationPath);
      ZDCAconf->SetRecenteringHistPath(gkZdcRecenteringPath);
      ZDCAconf->SetCorrelationHistPath(gkZdcDiagonalizationPath);
      ZDCAconf->SetQnCorrectionsCommonAxes(ZDCbinning);
      //ZDCAconf->SetDataVectorEqualizationAxes(ZDCchannelbinning);
      EPmanager->AddQnConfiguration(ZDCAconf, QnCorrectionsVarManager::kZDC);

      QnCorrectionsConfiguration * ZDCCconf = new QnCorrectionsConfiguration();
      ZDCCconf->SetCalibrationStep(gkZDCstep);
      ZDCCconf->SetQnNormalization(1);
      //ZDCCconf->SetDataVectorEqualizationMethod(2);
      //ZDCCconf->SetTwistAndRescalingMethod(2);
      ZDCCconf->SetDataVectorIdList(new TArrayS(10, ZDCchannels[1]));
      ZDCCconf->SetCalibrationDetectorName("ZDCC");
      //ZDCCconf->SetEqualizationDetectorName("ZDCC");
      ZDCCconf->SetQnConfigurationName("ZDCC");
      ZDCCconf->SetReferenceQnForTwistAndRescaling("ZDCA","TPC");
      //ZDCCconf->SetEqualizationHistPath(gkZdcEqualizationPath);
      ZDCCconf->SetRecenteringHistPath(gkZdcRecenteringPath);
      ZDCCconf->SetCorrelationHistPath(gkZdcDiagonalizationPath);
      ZDCCconf->SetQnCorrectionsCommonAxes(ZDCbinning);
      //ZDCCconf->SetDataVectorEqualizationAxes(ZDCchannelbinning);
      EPmanager->AddQnConfiguration(ZDCCconf, QnCorrectionsVarManager::kZDC);

  }








void AddFMD(QnCorrectionsManager* EPmanager){

  const Char_t* gkFmdEqualizationPath      = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
  const Char_t* gkFmdRecenteringPath       = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
  const Char_t* gkFmdDiagonalizationPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20150215_12212b/%d/CalibrationHistogramsMerged.root";
  
  
  const Int_t gkFMDstep=0;


  //-----------------------------------------------------------
  // Binning for Q-vector calibration

  const Int_t FMDdim              = 2;
  QnCorrectionsAxes * FMDbinning = new QnCorrectionsAxes(FMDdim);
   

  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {100};

  FMDbinning->AddAxis(0, QnCorrectionsVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

    
  // vertex Z bins
  const Int_t nVtxZwidths = 3;   
  Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
  Int_t VtxZbins[nVtxZwidths] = {1,8,1};

  FMDbinning->AddAxis(1, QnCorrectionsVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);

  QnCorrectionsConfiguration * FMDAcor = new QnCorrectionsConfiguration();
  FMDAcor->SetCalibrationStep(gkFMDstep);
  FMDAcor->SetQnNormalization(1);
  //FMDAcor->SetDataVectorEqualizationMethod(1);
  //FMDAcor->SetTwistAndRescalingMethod(2);
  FMDAcor->SetCalibrationDetectorName("FMD");
  //FMDAcor->SetEqualizationDetectorName("FMD");
  FMDAcor->SetQnConfigurationName("FMDAcor");
  FMDAcor->SetReferenceQnForTwistAndRescaling("TPC","FMDCcor");
  //FMDAcor->SetEqualizationHistPath(gkFmdEqualizationPath);
  FMDAcor->SetRecenteringHistPath(gkFmdRecenteringPath);
  //FMDAcor->SetCorrelationHistPath(gkFmdDiagonalizationPath);
  FMDAcor->SetQnCorrectionsCommonAxes(FMDbinning);
  EPmanager->AddQnConfiguration(FMDAcor, QnCorrectionsVarManager::kFMD);

  QnCorrectionsConfiguration * FMDCcor = new QnCorrectionsConfiguration();
  FMDCcor->SetCalibrationStep(gkFMDstep);
  FMDCcor->SetQnNormalization(1);
  //FMDCcor->SetDataVectorEqualizationMethod(1);
  FMDCcor->SetTwistAndRescalingMethod(2);
  FMDCcor->SetCalibrationDetectorName("FMD");
  //FMDCcor->SetEqualizationDetectorName("FMD");
  FMDCcor->SetQnConfigurationName("FMDCcor");
  FMDCcor->SetReferenceQnForTwistAndRescaling("TPC","VZEROA");
  //FMDCcor->SetEqualizationHistPath(gkFmdEqualizationPath);
  FMDCcor->SetRecenteringHistPath(gkFmdRecenteringPath);
  //FMDCcor->SetCorrelationHistPath(gkFmdDiagonalizationPath);
  FMDCcor->SetQnCorrectionsCommonAxes(FMDbinning);
  EPmanager->AddQnConfiguration(FMDCcor, QnCorrectionsVarManager::kFMD);
 



//    /////////////// Add FMD subdetectors ///////////////////
//
//    
//      QnCorrectionsConfiguration * FMD1conf = new QnCorrectionsConfiguration();
//      FMD1conf->SetDetectorId(QnCorrectionsVarManager::kFMD1);
//      FMD1conf->SetCalibrationStep(gkFMDstep);
//      FMD1conf->SetQnNormalization(1);
//      FMD1conf->SetDataVectorEqualizationMethod(2);
//      FMD1conf->SetCalibrationDetectorName("FMD1");
//      FMD1conf->SetEqualizationDetectorName("FMD1");
//      FMD1conf->SetQnConfigurationName("FMD1");
//      FMD1conf->SetReferenceQnForTwistAndRescaling("FMD2I","TPC");
//      FMD1conf->SetEqualizationHistPath(gkFmdEqualizationPath);
//      FMD1conf->SetRecenteringHistPath(gkFmdRecenteringPath);
//      FMD1conf->SetCorrelationHistPath(gkFmdDiagonalizationPath);
//      EPmanager->AddQnConfiguration(FMD1conf);
//
//      QnCorrectionsConfiguration * FMD2Iconf = new QnCorrectionsConfiguration();
//      FMD2Iconf->SetDetectorId(QnCorrectionsVarManager::kFMD2I);
//      FMD2Iconf->SetCalibrationStep(gkFMDstep);
//      FMD2Iconf->SetQnNormalization(1);
//      FMD2Iconf->SetDataVectorEqualizationMethod(2);
//      FMD2Iconf->SetCalibrationDetectorName("FMD2I");
//      FMD2Iconf->SetEqualizationDetectorName("FMD2I");
//      FMD2Iconf->SetQnConfigurationName("FMD2I");
//      FMD2Iconf->SetReferenceQnForTwistAndRescaling("FMD1","TPC");
//      FMD2Iconf->SetEqualizationHistPath(gkFmdEqualizationPath);
//      FMD2Iconf->SetRecenteringHistPath(gkFmdRecenteringPath);
//      FMD2Iconf->SetCorrelationHistPath(gkFmdDiagonalizationPath);
//      EPmanager->AddQnConfiguration(FMD2Iconf);
//
//
//
//
//
//    ////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////
//    ////////////////////////////////////////////////////////

  }








// Helper macros for creating chains
// from: CreateESDChain.C,v 1.10 jgrosseo Exp

TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
{
  // creates chain of files in a given directory or file containing a list.
  // In case of directory the structure is expected as:
  // <aDataDir>/<dir0>/AliESDs.root
  // <aDataDir>/<dir1>/AliESDs.root
  // ...
  
  if (!aDataDir)
    return 0;
  
  Long_t id, size, flags, modtime;
  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
    {
      printf("%s not found.\n", aDataDir);
      return 0;
    }
  
  TChain* chain = new TChain("esdTree");
  TChain* chaingAlice = 0;
  
  if (flags & 2)
    {
      TString execDir(gSystem->pwd());
      TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
      TList* dirList            = baseDir->GetListOfFiles();
      Int_t nDirs               = dirList->GetEntries();
      gSystem->cd(execDir);
      
      Int_t count = 0;
      
      for (Int_t iDir=0; iDir<nDirs; ++iDir)
	{
	  TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
	  if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
	    continue;
	  
	  if (offset > 0)
	    {
	      --offset;
	      continue;
	    }
	  
	  if (count++ == aRuns)
	    break;
	  
	  TString presentDirName(aDataDir);
	  presentDirName += "/";
	  presentDirName += presentDir->GetName();	  
	  chain->Add(presentDirName + "/AliESDs.root/esdTree");
	  //  cerr<<presentDirName<<endl;
	}
      
    }
  else
    {
      // Open the input stream
      ifstream in;
      in.open(aDataDir);
      
      Int_t count = 0;
      
      // Read the input list of files and add them to the chain
      TString esdfile;
      while(in.good()) {
	in >> esdfile;
	if (!esdfile.Contains("root")) continue; // protection
	
	if (offset > 0)
	  {
	    --offset;
	    continue;
	  }
	
	if (count++ == aRuns)
	  break;
	
	// add esd file
	chain->Add(esdfile);
      }
      
      in.close();
    }
  
  return chain;
}




//__________________________________________________________________
void DefineHistograms(QnCorrectionsManager* EPmanager, QnCorrectionsHistos* histos, TString histClass) {
  //
  // define the histograms
  //
  #define VAR QnCorrectionsVarManager

  for(Int_t iconf=0; iconf<EPmanager->NumberOfQnConfigurations(); iconf++){
  QnCorrectionsConfiguration* QnConf = (QnCorrectionsConfiguration*) EPmanager->GetQnConfiguration(iconf);
   if(!QnConf) continue;
   if(QnConf->IsTracking()) histClass+= "TrackQA_"+QnConf->QnConfigurationName()+";";
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
      histos->AddHistogram(classStr.Data(),"RunNo","Run numbers;Run", kFALSE, kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo);
      histos->AddHistogram(classStr.Data(),"BC","Bunch crossing;BC", kFALSE,3000,0.,3000.,VAR::kBC);
      histos->AddHistogram(classStr.Data(),"IsPhysicsSelection","Physics selection flag;;", kFALSE,
          2,-0.5,1.5,VAR::kIsPhysicsSelection, 0,0.0,0.0,VAR::kNothing, 0,0.0,0.0,VAR::kNothing, "off;on");

      histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,0.0,0.0,VAR::kVtxZ);
      //histos->AddHistogram(classStr.Data(),"VtxZ","Vtx Z;vtx Z (cm)", kFALSE,300,-15.,15.,VAR::kVtxZ);
      histos->AddHistogram(classStr.Data(),"VtxX","Vtx X;vtx X (cm)", kFALSE,300,-1.,1.,VAR::kVtxX);
      histos->AddHistogram(classStr.Data(),"VtxY","Vtx Y;vtx Y (cm)", kFALSE,300,-1.,1.,VAR::kVtxY);


      histos->AddHistogram(classStr.Data(),"CentVZEROvsCentSPD","Centrality(VZERO);centrality VZERO (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentVZERO, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentSPD","Centrality(TPC);centrality TPC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentVZERO","Centrality(TPC);centrality TPC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentTPCvsCentZDC","Centrality(TPC);centrality TPC (percents);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC, 100, 0.0, 100.0, VAR::kCentZDC);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentVZERO","Centrality(ZDC);centrality ZDC (percents);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentZDCvsCentSPD","Centrality(ZDC);centrality ZDC (percents);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC, 100, 0.0, 100.0, VAR::kCentSPD);

      histos->AddHistogram(classStr.Data(),"MultVZEROvsCentVZERO","Multiplicity;multiplicity VZERO;VZERO centrality", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100,0.,100., VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"MultSPDvsCentPSD","Multiplicity;SPD tracklets;SPD centrality", kFALSE,
          100, 0.0, 3000.0, VAR::kSPDntracklets, 100,0.,100., VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"MultTPCvsCentSPD","Multiplicity;TPC selected tracks;TPC centrality", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100,0.,100., VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"MultZDCvsCentZDC","Multiplicity;multiplicity ZDC;ZDC centrality", kFALSE,
          100, 0.0, 300000.0, VAR::kZDCTotalEnergy, 100,0.,100., VAR::kCentZDC);


      histos->AddHistogram(classStr.Data(),"MultTPCvsMultVZERO","Multiplicity;tracks TPC;multiplicity VZERO", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 32000.0, VAR::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultSPD","Multiplicity;tracklets SPD;tracks TPC", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 3000.0, VAR::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultVZERO","Multiplicity;tracklets SPD;multiplicity VZERO", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100, 0.0, 3000.0, VAR::kSPDntracklets);
      histos->AddHistogram(classStr.Data(),"MultTPCvsMultZDC","Multiplicity;tracks TPC;energy ZDC", kFALSE,
          100, 0.0, 3500.0, VAR::kNtracksSelected, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultVZEROvsMultZDC","Multiplicity;multiplicity VZERO;energy ZDC", kFALSE,
          100, 0.0, 32000.0, VAR::kVZEROTotalMult, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultSPDvsMultZDC","Multiplicity;tracklets SPD;energy ZDC", kFALSE,
          100, 0.0, 3000.0, VAR::kSPDntracklets, 100, 0.0, 300000.0, VAR::kZDCTotalEnergy);



      histos->AddHistogram(classStr.Data(),"MultVZERO","Multiplicity;multiplicity VZERO", kFALSE,
          320, 0.0, 25000.0, VAR::kVZEROTotalMult);
      histos->AddHistogram(classStr.Data(),"MultVZEROA","Multiplicity;multiplicity VZEROA", kFALSE,
          250, 0.0, 0.0, VAR::kVZEROATotalMult);//10000.0
      histos->AddHistogram(classStr.Data(),"MultVZEROC","Multiplicity;multiplicity VZEROC", kFALSE,
          250, 0.0, 0.0, VAR::kVZEROCTotalMult);//15000.0
      histos->AddHistogram(classStr.Data(),"MultZDC","Multiplicity;multiplicity ZDC", kFALSE,
          200, 0.0, 300000.0, VAR::kZDCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCA","Multiplicity;multiplicity ZDCA", kFALSE,
          200, 0.0, 150000.0, VAR::kZDCATotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultZDCC","Multiplicity;multiplicity ZDCC", kFALSE,
          200, 0.0, 150000.0, VAR::kZDCCTotalEnergy);
      histos->AddHistogram(classStr.Data(),"MultFMD1","Multiplicity;multiplicity FMD1", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD1TotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2I","Multiplicity;multiplicity FMD2I", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD2ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD2O","Multiplicity;multiplicity FMD2O", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD2OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3I","Multiplicity;multiplicity FMD3I", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD3ITotalMult);
      histos->AddHistogram(classStr.Data(),"MultFMD3O","Multiplicity;multiplicity FMD3O", kFALSE,
          300, 0.0, 10000.0, VAR::kFMD3OTotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROA","Multiplicity;multiplicity TZEROA", kFALSE,
          300, 0.0, 3000.0, VAR::kTZEROATotalMult);
      histos->AddHistogram(classStr.Data(),"MultTZEROC","Multiplicity;multiplicity TZEROC", kFALSE,
          300, 0.0, 3000.0, VAR::kTZEROCTotalMult);




      histos->AddHistogram(classStr.Data(),"CentVZERO","Centrality(VZERO);centrality VZERO (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD","Centrality(SPD);centrality SPD (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC","Centrality(TPC);centrality TPC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC","Centrality(ZDC);centrality ZDC (percents)", kFALSE,
          100, 0.0, 100.0, VAR::kCentZDC);

      histos->AddHistogram(classStr.Data(),"CentQuality","Centrality quality;centrality quality", kFALSE,
          100, -50.5, 49.5, VAR::kCentQuality);
      histos->AddHistogram(classStr.Data(),"CentVZERO_Run_prof","<Centrality(VZERO)> vs run;Run; centrality VZERO (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"CentSPD_Run_prof","<Centrality(SPD)> vs run;Run; centrality SPD (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"CentTPC_Run_prof","<Centrality(TPC)> vs run;Run; centrality TPC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentTPC);
      histos->AddHistogram(classStr.Data(),"CentZDC_Run_prof","<Centrality(ZDC)> vs run;Run; centrality ZDC (%)", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0.0, 100.0, VAR::kCentZDC);


      histos->AddHistogram(classStr.Data(),"NV0sTotal","Number of V0 candidates per event;# pairs", kFALSE,
          1000,0.,30000.,VAR::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0sSelected","Number of selected V0 candidates per event;# pairs", kFALSE,
          1000,0.,10000.,VAR::kNV0selected);
      histos->AddHistogram(classStr.Data(),"NPairs","Number of candidates per event;# pairs", kFALSE,
          5000,0.,5000.,VAR::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NPairsSelected", "Number of selected pairs per event; #pairs", kFALSE,
          5000,0.,5000.,VAR::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal","Number of total tracks per event;# tracks", kFALSE,
          1000,0.,30000.,VAR::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected","Number of selected tracks per event;# tracks", kFALSE,
          1000,0.,30000.,VAR::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets", "SPD #tracklets in |#eta|<1.0; tracklets", kFALSE,
          3000, -0.5, 2999.5, VAR::kSPDntracklets);

      histos->AddHistogram(classStr.Data(),"NV0total_Run_prof", "<Number of total V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNV0total);
      histos->AddHistogram(classStr.Data(),"NV0selected_Run_prof", "<Number of selected V0s> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNV0selected);
      histos->AddHistogram(classStr.Data(),"Ndielectrons_Run_prof", "<Number of dielectrons> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNdielectrons);
      histos->AddHistogram(classStr.Data(),"NpairsSelected_Run_prof", "<Number of selected pairs> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNpairsSelected);
      histos->AddHistogram(classStr.Data(),"NTracksTotal_Run_prof", "<Number of tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNtracksTotal);
      histos->AddHistogram(classStr.Data(),"NTracksSelected_Run_prof", "<Number of selected tracks> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kNtracksSelected);
      histos->AddHistogram(classStr.Data(),"SPDntracklets_Run_prof", "<SPD ntracklets> per run; Run; #tracks", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 100, 0., 10000., VAR::kSPDntracklets);

      histos->AddHistogram(classStr.Data(),"VtxZ_CentVZERO","Centrality(VZERO) vs vtx. Z;vtx Z (cm); centrality VZERO (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentVZERO);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentSPD","Centrality(SPD) vs vtx. Z;vtx Z (cm); centrality SPD (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentSPD);
      histos->AddHistogram(classStr.Data(),"VtxZ_CentTPC","Centrality(TPC) vs vtx. Z;vtx Z (cm); centrality TPC (%)", kFALSE,
          300,-15.,15.,VAR::kVtxZ, 100, 0.0, 100.0, VAR::kCentTPC);
      continue;
    }  // end if className contains "Event"    


    // Offline trigger histograms
    if(classStr.Contains("OfflineTriggers")) {
      histos->AddHistClass(classStr.Data());

      TString triggerNames = "";
      for(Int_t i=0; i<64; ++i) {triggerNames += Form("%s",VAR::fOfflineTriggerNames[i]); triggerNames+=";";}

      histos->AddHistogram(classStr.Data(), "Triggers", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTrigger, 2, -0.5, 1.5, VAR::kOfflineTriggerFired, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data(), "off;on");
      histos->AddHistogram(classStr.Data(), "Triggers2", "Offline triggers fired; ; ;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 0, 0.0, 0.0, VAR::kNothing, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentVZERO_Triggers2", "Offline triggers fired vs centrality VZERO; ; centrality VZERO;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentVZERO, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentTPC_Triggers2", "Offline triggers fired vs centrality TPC; ; centrality TPC;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentTPC, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentSPD_Triggers2", "Offline triggers fired vs centrality SPD; ; centrality SPD;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentSPD, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "CentZDC_Triggers2", "Offline triggers fired vs centrality ZDC; ; centrality ZDC;", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 20, 0.0, 100.0, VAR::kCentZDC, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      histos->AddHistogram(classStr.Data(), "VtxZ_Triggers2", "Offline triggers fired vs vtxZ; ; vtx Z (cm.);", kFALSE,
          64, -0.5, 63.5, VAR::kOfflineTriggerFired2, 200, -20.0, 20.0, VAR::kVtxZ, 0, 0.0, 0.0, VAR::kNothing, triggerNames.Data());
      continue;
    }

    // Track histograms
    if(classStr.Contains("Tracks")) {
      histos->AddHistClass(classStr.Data());
      for(Int_t ih=0; ih<6; ++ih) {
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
      }
    }


    // Track histograms
    if(classStr.Contains("TrackQA")) {
      histos->AddHistClass(classStr.Data());

      histos->AddHistogram(classStr.Data(), "Pt", "p_{T} distribution; p_{T} (GeV/c^{2});", kFALSE,
          1000, 0.0, 50.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta", "#eta illumination; #eta;", kFALSE,
          1000, -1.5, 1.5, VAR::kEta);
      histos->AddHistogram(classStr.Data(), "Phi", "#varphi illumination; #varphi;", kFALSE,
          1000, 0.0, 6.3, VAR::kPhi);
      histos->AddHistogram(classStr.Data(), "DCAxy", "DCAxy; DCAxy (cm.)", kFALSE,
          1000, -10.0, 10.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz", "DCAz; DCAz (cm.)", kFALSE,
          1000, -10.0, 10.0, VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "TPCncls", "TPCncls; TPCncls", kFALSE,
          160, 0.0, 160.0, VAR::kTPCncls);

      // run dependence
      histos->AddHistogram(classStr.Data(), "Pt_Run", "<p_{T}> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, 0.0, 50.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Run", "<#eta> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -1.5, 1.5, VAR::kEta);      
      histos->AddHistogram(classStr.Data(), "Phi_Run", "<#varphi> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, 0.0, 6.3, VAR::kPhi);      
      histos->AddHistogram(classStr.Data(), "DCAxy_Run", "<DCAxy> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -10.0, 10.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Run", "<DCAz> vs run; run;", kTRUE,
          kNRunBins, runHistRange[0], runHistRange[1], VAR::kRunNo, 1000, -10.0, 10.0, VAR::kDcaZ);

      // correlations between parameters
      histos->AddHistogram(classStr.Data(), "Eta_Pt_prof", "<p_{T}> vs #eta; #eta; p_{T} (GeV/c);", kTRUE,
          300, -1.5, +1.5, VAR::kEta, 100, 0.0, 10.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt", "p_{T} vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kFALSE,
          300, -0.01, 6.3, VAR::kPhi, 100, 0.0, 2.2, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Phi_Pt_prof", "<p_{T}> vs #varphi; #varphi (rad.); p_{T} (GeV/c)", kTRUE,
          300, 0.0, 6.3, VAR::kPhi, 100, 0.0, 10.0, VAR::kPt);
      histos->AddHistogram(classStr.Data(), "Eta_Phi", "#varphi vs #eta; #eta; #varphi (rad.);", kFALSE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi);
      histos->AddHistogram(classStr.Data(), "TPCncls_Eta_Phi_prof", "<TPC ncls> vs #varphi vs #eta; #eta; #varphi (rad.);TPC ncls", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kTPCncls);
      histos->AddHistogram(classStr.Data(), "DCAxy_Eta_Phi_prof", "<DCAxy> vs #varphi vs #eta; #eta; #varphi (rad.);DCAxy (cm)", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "DCAz_Eta_Phi_prof", "<DCAz> vs #varphi vs #eta; #eta; #varphi (rad.);DCAz (cm)", kTRUE,
          200, -1.0, +1.0, VAR::kEta, 100, 0.0, 6.3, VAR::kPhi, 10, 0.0, 200., VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Pt_DCAxy", "DCAxy vs p_{T}; p_{T} (GeV/c); DCA_{xy} (cm)", kFALSE,
          100, 0.0, 10.0, VAR::kPt, 500, -2.0, 2.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Pt_DCAz", "DCAz vs p_{T}; p_{T} (GeV/c); DCA_{z} (cm)", kFALSE,
          100, 0.0, 10.0, VAR::kPt, 500, -2.0, 2.0, VAR::kDcaZ);
      histos->AddHistogram(classStr.Data(), "Eta_DCAxy", "DCAxy vs #eta; #eta; DCA_{xy} (cm)", kFALSE,
          100, -1.0, 1.0, VAR::kEta, 500, -2.0, 2.0, VAR::kDcaXY);
      histos->AddHistogram(classStr.Data(), "Eta_DCAz", "DCAz vs #eta; #eta; DCA_{z} (cm)", kFALSE,
          100, -1.0, 1.0, VAR::kEta, 500, -2.0, 2.0, VAR::kDcaZ);

      for(Int_t ih=0; ih<6; ++ih) {
        //histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO",ih+1), Form("<cos%d #varphi> vs (CentVZERO); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        //histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO",ih+1), Form("<sin%d #varphi> vs (CentVZERO); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
        //                   20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_Pt_Eta",ih+1), Form("<cos%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <cos%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, VAR::kEta, 30, 0.0, 3.0, VAR::kPt, 500, -1.0, 1.0, VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_Pt_Eta",ih+1), Form("<sin%d #varphi> vs (#eta,p_{T}); #eta; p_{T} (GeV/c); <sin%d #varphi>", ih+1, ih+1), kTRUE,
            20, -1.0, 1.0, VAR::kEta, 30, 0.0, 3.0, VAR::kPt, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Cos%dPhi_CentVZERO_VtxZ",ih+1), Form("<cos%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <cos%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, VAR::kVtxZ, 20, 0.0, 100.0, VAR::kCentVZERO, 500, -1., 1., VAR::kCosNPhi+ih);
        histos->AddHistogram(classStr.Data(), Form("Sin%dPhi_CentVZERO_VtxZ",ih+1), Form("<sin%d #varphi> vs (CentVZERO,VtxZ); Z (cm.); centrality VZERO; <sin%d #varphi>", ih+1, ih+1), kTRUE,
            30, -15.0, 15.0, VAR::kVtxZ, 20, 0.0, 100.0, VAR::kCentVZERO, 500, -1.0, 1.0, VAR::kSinNPhi+ih);
      }
    }


  }

  cout << " done" << endl;
}







