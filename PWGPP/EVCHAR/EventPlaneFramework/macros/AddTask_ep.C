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



#include "AliTrigger.h"
//#include<iostream>
//#include<TSystem.h>
//#include<TROOT.h>
//#include"AliAnalysisManager.h"
//#include<AliEventPlaneManager.h>
//#include<AliEventPlaneCuts.h>
//#include<AliEventPlaneVarManager.h>
//#include<AliAnalysisTaskEventPlaneCalibration.h>


void AddVZERO(AliEventPlaneManager* EPmanager);
void AddTPC(AliEventPlaneManager* EPmanager);
void AddTZERO(AliEventPlaneManager* EPmanager);
//void AddFMD(AliEventPlaneManager* EPmanager);
void AddZDC(AliEventPlaneManager* EPmanager);


void AddTask_ep() {
  //gSystem->Load("libESD");
  //gSystem->Load("libPWGPPevcharEp.so");

  //TChain* chain = CreateESDChain("/u/jonderw/files.txt", 1, 0);
  //AliAnalysisManager *mgr = new AliAnalysisManager("FlowAnalysisManager");
  //

  // Load common libraries
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  // gSystem->Load("libNetx.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  
  gSystem->Load("libPWGPPevcharEP.so");





  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    //mgr = new AliAnalysisManager("EventPlaneCorrections");

  //Error("AddTask_jonderw_ep", "No analysis manager found.");
  //return 0;
  }
  //AliVEventHandler* esdH = new AliESDInputHandler;
  //mgr->SetInputEventHandler(esdH);

 //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
 //AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

 //// Apply the event selection
 //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
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



 
  AliEventPlaneManager* EPmanager = new AliEventPlaneManager();
  AliAnalysisTaskEventPlaneCalibration* taskEP = new AliAnalysisTaskEventPlaneCalibration("EventPlaneCalibration");


  AliEventPlaneCuts *eventCuts = new AliEventPlaneCuts();
  eventCuts->AddCut(AliEventPlaneVarManager::kVtxZ,-7.0,7.0);
  eventCuts->AddCut(AliEventPlaneVarManager::kCentVZERO,0.0,90.0);
  //eventCuts->AddCut(AliEventPlaneVarManager::kIsPhysicsSelection,0.5,1.5);

  
  taskEP->SetEventCuts(eventCuts);
  taskEP->SetTrigger(AliTrigger::kMB);                 // Trigger stream to be used for calibration and QA histograms
  taskEP->SelectCollisionCandidates(AliTrigger::kMB);  // Events passing trigger and physics selection for analysis

  AddVZERO(EPmanager);
  AddTPC(EPmanager);
  //AddTZERO(EPmanager);
  AddZDC(EPmanager);
  //AddFMD(EPmanager);

  taskEP->SetEventPlaneManager(EPmanager);

  //create output container
    AliAnalysisDataContainer *cOutputHist =
    mgr->CreateContainer("CalibrationHistos",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "CalibrationHistograms.root");
    AliAnalysisDataContainer *cOutputQvec =
    mgr->CreateContainer("CalibratedQvector",
                         TTree::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "Qvectors.root");
    AliAnalysisDataContainer *cOutputHistQA =
    mgr->CreateContainer("CalibrationQA",
                         TList::Class(),
                         AliAnalysisManager::kOutputContainer,
                         "CalibrationQA.root");


  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  //mgr->ConnectInput(taskDst,  0, cinput1);
  ////mgr->ConnectOutput(taskDst, 0, coutput1 );
  //mgr->ConnectOutput(taskDst, 1, cOutputHist1);
  //mgr->ConnectOutput(taskDst, 2, cOutputHist2);
  //mgr->ConnectOutput(taskDst, 3, cReducedEvent);

  mgr->AddTask(taskEP);
  //mgr->ConnectInput(taskEP,  0, cReducedEvent);
  //mgr->ConnectInput(taskEP,  0, cinput1);
  mgr->ConnectInput(taskEP,  0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(taskEP, 1, cOutputHist );
  mgr->ConnectOutput(taskEP, 2, cOutputQvec );
  mgr->ConnectOutput(taskEP, 3, cOutputHistQA );


  return;
}




void AddVZERO(AliEventPlaneManager* EPmanager){

  //const Char_t* gkVzeroEqualizationPath = "";
  const Char_t* gkVzeroEqualizationPath = "/alien/alice/cern.ch/user/j/jonderwa/AnalysisTesting5/output/000%d/CalibrationHistograms.root";
  const Char_t* gkVzeroRecenteringPath  = "";
  const Char_t* gkVzeroDiagonalizationPath = "";
  //const Char_t* gkVzeroEqualizationPath = "/hera/alice/jonderw/PbPb2010/EventPlane_20141210_000-1-1/";
  //const Char_t* gkVzeroRecenteringPath  = "/hera/alice/jonderw/PbPb2010/EventPlane_20141210_000-1-1/";
  //const Char_t* gkVzeroDiagonalizationPath = "/hera/alice/jonderw/PbPb2010/EventPlane_20131010_3222-1/";
  const Int_t gkVZEROstep=1;
  Short_t VZEROchannels[4][64];
  for(Int_t iv0=0; iv0<4; iv0++) for(Int_t ich=0; ich<64; ich++) VZEROchannels[iv0][ich] = 0;
  
  for(Int_t ich=32; ich<64; ich++) VZEROchannels[0][ich] = 1;  // channel list: value 1 if channel should be used
  for(Int_t ich=0; ich<32; ich++) VZEROchannels[1][ich] = 1;

  //-----------------------------------------------------------
  // Binning for Q-vector calibration

  const Int_t VZEROdim              = 2;
  AliEventPlaneBinning * VZERObinning = new AliEventPlaneBinning(VZEROdim);
   

  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {100};

  VZERObinning->AddBinAxis(0, AliEventPlaneVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

    
  // vertex Z bins
  const Int_t nVtxZwidths = 3;   
  Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
  Int_t VtxZbins[nVtxZwidths] = {1,8,1};

  VZERObinning->AddBinAxis(1, AliEventPlaneVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);
    
  //-----------------------------------------------------------
  // Binning for channel equalization
  AliEventPlaneBinning * VZEROchannelbinning = new AliEventPlaneBinning(VZEROdim+1); // +1 is for channel axis
  VZEROchannelbinning->AddBinAxis(0, AliEventPlaneVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);
  VZEROchannelbinning->AddBinAxis(1, AliEventPlaneVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);
  VZEROchannelbinning->SetNchannels(64);
  
  ////////// end of binning
  
 
  AliEventPlaneConfiguration * VZEROAconf = new AliEventPlaneConfiguration();
  VZEROAconf->SetCalibrationStep(gkVZEROstep);
  VZEROAconf->SetCalibrationMethod(1);
  VZEROAconf->SetEqualizationMethod(1);
  VZEROAconf->SetTwistAndScalingMethod(2);
  VZEROAconf->SetChannelList(new TArrayS(64, VZEROchannels[0]));
  VZEROAconf->SetCalibrationDetectorName("VZEROA");
  VZEROAconf->SetEqualizationDetectorName("VZEROA");
  VZEROAconf->SetEventPlaneDetectorName("VZEROA");
  VZEROAconf->SetCorrelationDetectorNames("TPC","VZEROC");
  VZEROAconf->SetEqualizationHistPath(gkVzeroEqualizationPath);
  VZEROAconf->SetRecenteringHistPath(gkVzeroRecenteringPath);
  VZEROAconf->SetCorrelationHistPath(gkVzeroDiagonalizationPath);
  VZEROAconf->SetCalibrationBinning(VZERObinning);
  VZEROAconf->SetEqualizationBinning(VZEROchannelbinning);

  EPmanager->AddEventPlaneConfiguration(VZEROAconf, AliEventPlaneManager::kVZERO);



  AliEventPlaneConfiguration * VZEROCconf = new AliEventPlaneConfiguration();
  VZEROCconf->SetCalibrationStep(gkVZEROstep);
  VZEROCconf->SetCalibrationMethod(1);
  VZEROCconf->SetEqualizationMethod(1);
  VZEROCconf->SetTwistAndScalingMethod(2);
  VZEROCconf->SetChannelList(new TArrayS(64, VZEROchannels[1]));
  VZEROCconf->SetCalibrationDetectorName("VZEROC");
  VZEROCconf->SetEqualizationDetectorName("VZEROC");
  VZEROCconf->SetEventPlaneDetectorName("VZEROC");
  VZEROCconf->SetCorrelationDetectorNames("TPC","VZEROA");
  VZEROCconf->SetEqualizationHistPath(gkVzeroEqualizationPath);
  VZEROCconf->SetRecenteringHistPath(gkVzeroRecenteringPath);
  VZEROCconf->SetCorrelationHistPath(gkVzeroDiagonalizationPath);
  VZEROCconf->SetCalibrationBinning(VZERObinning);
  VZEROCconf->SetEqualizationBinning(VZEROchannelbinning);

  EPmanager->AddEventPlaneConfiguration(VZEROCconf, AliEventPlaneManager::kVZERO);
  
}
 
  void AddTPC(AliEventPlaneManager* EPmanager){

     /////////////// Add TPC subdetectors ///////////////////
        
const Int_t gkTPCstep=0;
//const Char_t* gkTpcRecenteringPath   = "";
const Char_t* gkTpcRecenteringPath   = "/alien/alice/cern.ch/user/j/jonderwa/AnalysisTesting5/output/000%d/CalibrationHistograms.root";
//const Char_t* gkTpcRecenteringPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20131120_00000/";
const Char_t* gkTpcDiagonalizationPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20141207_00-1-1-1/";



    //-----------------------------------------------------------
    // Binning for Q-vector calibration
    //
    const Int_t TPCdim              = 2;    // number of calibration parameters
    AliEventPlaneBinning * TPCbinning = new AliEventPlaneBinning(TPCdim);  //  declare binning object
    
    // centrality bins
    const Int_t nCtwidths = 1;   
    Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
    Int_t Ctbins[nCtwidths] = {100};

    TPCbinning->AddBinAxis(0, AliEventPlaneVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

      
    // vertex Z bins
    const Int_t nVtxZwidths = 3;   
    Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
    Int_t VtxZbins[nVtxZwidths] = {1,8,1};

    TPCbinning->AddBinAxis(1, AliEventPlaneVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);

    ////////// end of binning
    


    AliEventPlaneConfiguration * TPCconf = new AliEventPlaneConfiguration();
    TPCconf->SetCalibrationStep(gkTPCstep);
    TPCconf->SetCalibrationMethod(1);
    TPCconf->SetTwistAndScalingMethod(0);
    TPCconf->SetCalibrationDetectorName("TPC");
    TPCconf->SetEventPlaneDetectorName("TPC");
    TPCconf->SetCorrelationDetectorNames("VZEROC","VZEROA");
    TPCconf->SetRecenteringHistPath(gkTpcRecenteringPath);
    TPCconf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    //TPCconf->CreateQvectorHistograms(TPCdim, TPCvar, binLimits);
    TPCconf->SetCalibrationBinning(TPCbinning);
   
    AliEventPlaneCuts *trackTPC = new AliEventPlaneCuts();
    trackTPC->AddCut(AliEventPlaneVarManager::kDcaXY,-1.0,1.0);
    trackTPC->AddCut(AliEventPlaneVarManager::kDcaZ,-3.0,3.0);
    trackTPC->AddCut(AliEventPlaneVarManager::kEta,-0.8,0.8);
    trackTPC->AddCut(AliEventPlaneVarManager::kEta,-0.3,0.3, kTRUE);  // exclusion region
    trackTPC->AddCut(AliEventPlaneVarManager::kPt,0.2,2.);
    trackTPC->AddCut(AliEventPlaneVarManager::kTPCnclsIter1,70.0,161.0);
    trackTPC->AddCut(AliEventPlaneVarManager::kTPCchi2Iter1,0.2,4.0);
    TPCconf->SetTrackCuts(trackTPC);

    //AliFlowCuts *trackTPC = new AliFlowCuts();
    //trackTPC = trackTPC->GetStandardTPCStandaloneCuts2010();
    //TPCconf->SetTrackCuts(trackTPC);

    EPmanager->AddEventPlaneConfiguration(TPCconf, AliEventPlaneManager::kTPC);



    //AliEventPlaneConfiguration * TPCposConf = new AliEventPlaneConfiguration();
    //TPCposConf->SetCalibrationStep(gkTPCstep);
    //TPCposConf->SetCalibrationMethod(1);
    //TPCposConf->SetTwistAndScalingMethod(0);
    //TPCposConf->SetCalibrationDetectorName("TPC");
    //TPCposConf->SetEventPlaneDetectorName("TPCpos");
    //TPCposConf->SetCorrelationDetectorNames("TPCneg","VZEROA");
    //TPCposConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    //TPCposConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
   
    //AliEventPlaneCuts *trackTPCpos = new AliEventPlaneCuts();
    //trackTPCpos->CopyCuts(trackTPC);
    //trackTPCpos->AddCut(kCharge, 0.5, 1.5);
    //TPCposConf->SetTrackCuts(trackTPCpos);

    //EPmanager->AddEventPlaneConfiguration(TPCposConf);

    //AliEventPlaneConfiguration * TPCnegConf = new AliEventPlaneConfiguration();
    //TPCnegConf->SetCalibrationStep(gkTPCstep);
    //TPCnegConf->SetCalibrationMethod(1);
    //TPCnegConf->SetTwistAndScalingMethod(0);
    //TPCnegConf->SetCalibrationDetectorName("TPC");
    //TPCnegConf->SetEventPlaneDetectorName("TPCneg");
    //TPCnegConf->SetCorrelationDetectorNames("TPCpos","VZEROA");
    //TPCnegConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    //TPCnegConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
   
    //AliEventPlaneCuts *trackTPCneg = new AliEventPlaneCuts();
    //trackTPCneg->CopyCuts(trackTPC);
    //trackTPCneg->AddCut(kCharge, -0.5, -1.5);
    //TPCnegConf->SetTrackCuts(trackTPCneg);

    //EPmanager->AddEventPlaneConfiguration(TPCnegConf);


    

    //AliEventPlaneConfiguration * TPCrightConf = new AliEventPlaneConfiguration();
    //TPCrightConf->SetCalibrationStep(gkTPCstep);
    //TPCrightConf->SetCalibrationMethod(1);
    //TPCrightConf->SetTwistAndScalingMethod(0);
    //TPCrightConf->SetCalibrationDetectorName("TPC");
    //TPCrightConf->SetEventPlaneDetectorName("TPCright");
    //TPCrightConf->SetCorrelationDetectorNames("TPCleft","VZEROC");
    //TPCrightConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    //TPCrightConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    //TPCrightConf->SetCalibrationBinning(TPCbinning);
    //TPCrightconf->CreateQvectorHistograms(TPCdim, TPCvar, binLimits);
   
    //AliEventPlaneCuts *trackTPCright = new AliEventPlaneCuts();
    //trackTPCright->CopyCuts(trackTPC);
    //trackTPCright->AddCut(AliEventPlaneVarManager::kEta, 0.0, 1.0);
    //TPCrightConf->SetTrackCuts(trackTPCright);

    //EPmanager->AddEventPlaneConfiguration(TPCrightConf, AliEventPlaneManager::kTPC);

    //AliEventPlaneConfiguration * TPCleftConf = new AliEventPlaneConfiguration();
    //TPCleftConf->SetCalibrationStep(gkTPCstep);
    //TPCleftConf->SetCalibrationMethod(1);
    //TPCleftConf->SetTwistAndScalingMethod(0);
    //TPCleftConf->SetCalibrationDetectorName("TPC");
    //TPCleftConf->SetEventPlaneDetectorName("TPCleft");
    //TPCleftConf->SetCorrelationDetectorNames("TPCright","VZEROA");
    //TPCleftConf->SetRecenteringHistPath(gkTpcRecenteringPath);
    //TPCleftConf->SetCorrelationHistPath(gkTpcDiagonalizationPath);
    //TPCleftConf->SetCalibrationBinning(TPCbinning);
    //TPCleftConf->CreateQvectorHistograms(TPCdim, TPCvar, binLimits);
   
    //AliEventPlaneCuts *trackTPCleft = new AliEventPlaneCuts();
    //trackTPCleft->CopyCuts(trackTPC);
    //trackTPCleft->AddCut(AliEventPlaneVarManager::kEta, -1.0, 0);
    //TPCleftConf->SetTrackCuts(trackTPCleft);

    //EPmanager->AddEventPlaneConfiguration(TPCleftConf, AliEventPlaneManager::kTPC);




  } 
   void AddTZERO(AliEventPlaneManager* EPmanager){
    
    /////////////// Add TZERO subdetectors ///////////////////
const Int_t gkTZEROstep=0;
const Char_t* gkTzeroEqualizationPath      = "/hera/alice/jonderw/PbPb2010/EventPlane_20141210_000-1-1/";
const Char_t* gkTzeroRecenteringPath       = "/hera/alice/jonderw/PbPb2010/EventPlane_20141210_000-1-1/";
const Char_t* gkTzeroDiagonalizationPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20141210_000-1-1/";

    

    Short_t TZEROchannels[2][24];
    for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<24; ich++) TZEROchannels[iv0][ich] = 0;
    
    for(Int_t ich=12; ich<24; ich++) TZEROchannels[0][ich] = 1;
    for(Int_t ich=0; ich<12; ich++) TZEROchannels[1][ich] = 1;


  //-----------------------------------------------------------
  // Binning for Q-vector calibration

  const Int_t TZEROdim              = 2;
  AliEventPlaneBinning * TZERObinning = new AliEventPlaneBinning(TZEROdim);
   

  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {100};

  TZERObinning->AddBinAxis(0, AliEventPlaneVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

    
  // vertex Z bins
  const Int_t nVtxZwidths = 3;   
  Double_t VtxZedges[nVtxZwidths+1] = {-10.0, -7.0, 7.0, 10.0};
  Int_t VtxZbins[nVtxZwidths] = {1,8,1};

  TZERObinning->AddBinAxis(1, AliEventPlaneVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);
    
  //-----------------------------------------------------------
  // Binning for channel equalization
  AliEventPlaneBinning * TZEROchannelbinning = new AliEventPlaneBinning(TZEROdim+1); // +1 is for channel axis
  TZEROchannelbinning->AddBinAxis(0, AliEventPlaneVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);
  TZEROchannelbinning->AddBinAxis(1, AliEventPlaneVarManager::kVtxZ, nVtxZwidths, VtxZbins, VtxZedges);
  TZEROchannelbinning->SetNchannels(24);
  
  ////////// end of binning
  
 
      AliEventPlaneConfiguration * TZEROAconf = new AliEventPlaneConfiguration();
      TZEROAconf->SetCalibrationStep(gkTZEROstep);
      TZEROAconf->SetCalibrationMethod(1);
      TZEROAconf->SetEqualizationMethod(2);
      TZEROAconf->SetTwistAndScalingMethod(2);
      TZEROAconf->SetChannelList(new TArrayS(24, TZEROchannels[0]));
      TZEROAconf->SetCalibrationDetectorName("TZEROA");
      TZEROAconf->SetEqualizationDetectorName("TZEROA");
      TZEROAconf->SetEventPlaneDetectorName("TZEROA");
      TZEROAconf->SetCorrelationDetectorNames("TPC","TZEROC");
      TZEROAconf->SetEqualizationHistPath(gkTzeroEqualizationPath);
      TZEROAconf->SetRecenteringHistPath(gkTzeroRecenteringPath);
      TZEROAconf->SetCorrelationHistPath(gkTzeroDiagonalizationPath);
      TZEROAconf->SetCalibrationBinning(TZERObinning);
      TZEROAconf->SetEqualizationBinning(TZEROchannelbinning);
      EPmanager->AddEventPlaneConfiguration(TZEROAconf, AliEventPlaneManager::kTZERO);

      AliEventPlaneConfiguration * TZEROCconf = new AliEventPlaneConfiguration();
      TZEROCconf->SetCalibrationStep(gkTZEROstep);
      TZEROCconf->SetCalibrationMethod(1);
      TZEROCconf->SetEqualizationMethod(2);
      TZEROCconf->SetTwistAndScalingMethod(2);
      TZEROCconf->SetChannelList(new TArrayS(24, TZEROchannels[1]));
      TZEROCconf->SetCalibrationDetectorName("TZEROC");
      TZEROCconf->SetEqualizationDetectorName("TZEROC");
      TZEROCconf->SetEventPlaneDetectorName("TZEROC");
      TZEROCconf->SetCorrelationDetectorNames("TPC","TZEROA");
      TZEROCconf->SetEqualizationHistPath(gkTzeroEqualizationPath);
      TZEROCconf->SetRecenteringHistPath(gkTzeroRecenteringPath);
      TZEROCconf->SetCorrelationHistPath(gkTzeroDiagonalizationPath);
      TZEROCconf->SetCalibrationBinning(TZERObinning);
      TZEROCconf->SetEqualizationBinning(TZEROchannelbinning);
      EPmanager->AddEventPlaneConfiguration(TZEROCconf, AliEventPlaneManager::kTZERO);

  }
  
 void AddZDC(AliEventPlaneManager* EPmanager){
    /////////////// Add ZDC subdetectors ///////////////////

const Int_t gkZDCstep=0;
const Char_t* gkZdcEqualizationPath   = "";
const Char_t* gkZdcRecenteringPath   = "";
const Char_t* gkZdcDiagonalizationPath   = "";


    Short_t ZDCchannels[2][10];
    for(Int_t iv0=0; iv0<2; iv0++) for(Int_t ich=0; ich<10; ich++) ZDCchannels[iv0][ich] = 0;
    
    for(Int_t ich=6; ich<10; ich++) ZDCchannels[0][ich] = 1;
    for(Int_t ich=1; ich<5; ich++)  ZDCchannels[1][ich] = 1;


  //-----------------------------------------------------------
  // Binning for Q-vector calibration

  const Int_t ZDCdim              = 3;
  AliEventPlaneBinning * ZDCbinning = new AliEventPlaneBinning(ZDCdim);
   

  // centrality bins
  const Int_t nCtwidths = 1;   
  Double_t Ctedges[nCtwidths+1] = {0.0, 100.0};
  Int_t Ctbins[nCtwidths] = {50};

  ZDCbinning->AddBinAxis(0, AliEventPlaneVarManager::kCentVZERO, nCtwidths, Ctbins, Ctedges);

    
  // vertex Z bins
  const Int_t nVtxXwidths = 1;   
  Double_t VtxXedges[nVtxXwidths+1] = {-0.06, 0.04};
  Int_t VtxXbins[nVtxXwidths] = {5};

  ZDCbinning->AddBinAxis(1, AliEventPlaneVarManager::kVtxX, nVtxXwidths, VtxXbins, VtxXedges);


    
  // vertex Z bins
  const Int_t nVtxYwidths = 1;   
  Double_t VtxYedges[nVtxYwidths+1] = {0.12,0.24};
  Int_t VtxYbins[nVtxYwidths] = {5};

  ZDCbinning->AddBinAxis(2, AliEventPlaneVarManager::kVtxY, nVtxYwidths, VtxYbins, VtxYedges);



      AliEventPlaneConfiguration * ZDCAconf = new AliEventPlaneConfiguration();
      ZDCAconf->SetCalibrationStep(gkZDCstep);
      ZDCAconf->SetCalibrationMethod(1);
      //ZDCAconf->SetEqualizationMethod(2);
      ZDCAconf->SetTwistAndScalingMethod(2);
      ZDCAconf->SetChannelList(new TArrayS(10, ZDCchannels[0]));
      ZDCAconf->SetCalibrationDetectorName("ZDCA");
      ZDCAconf->SetEqualizationDetectorName("ZDCA");
      ZDCAconf->SetEventPlaneDetectorName("ZDCA");
      ZDCAconf->SetCorrelationDetectorNames("ZDCC","TPC");
      //ZDCAconf->SetEqualizationHistPath(gkZdcEqualizationPath);
      ZDCAconf->SetRecenteringHistPath(gkZdcRecenteringPath);
      ZDCAconf->SetCorrelationHistPath(gkZdcDiagonalizationPath);
      ZDCAconf->SetCalibrationBinning(ZDCbinning);
      EPmanager->AddEventPlaneConfiguration(ZDCAconf, AliEventPlaneManager::kZDC);

      AliEventPlaneConfiguration * ZDCCconf = new AliEventPlaneConfiguration();
      ZDCCconf->SetCalibrationStep(gkZDCstep);
      ZDCCconf->SetCalibrationMethod(1);
      //ZDCCconf->SetEqualizationMethod(2);
      ZDCCconf->SetTwistAndScalingMethod(2);
      ZDCCconf->SetChannelList(new TArrayS(10, ZDCchannels[1]));
      ZDCCconf->SetCalibrationDetectorName("ZDCC");
      ZDCCconf->SetEqualizationDetectorName("ZDCC");
      ZDCCconf->SetEventPlaneDetectorName("ZDCC");
      ZDCCconf->SetCorrelationDetectorNames("ZDCA","TPC");
      //ZDCCconf->SetEqualizationHistPath(gkZdcEqualizationPath);
      ZDCCconf->SetRecenteringHistPath(gkZdcRecenteringPath);
      ZDCCconf->SetCorrelationHistPath(gkZdcDiagonalizationPath);
      ZDCCconf->SetCalibrationBinning(ZDCbinning);
      EPmanager->AddEventPlaneConfiguration(ZDCCconf, AliEventPlaneManager::kZDC);

  }
  void AddFMD(AliEventPlaneManager* EPmanager){
const Char_t* gkFmdEqualizationPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20131010_0000-1/";
const Char_t* gkFmdRecenteringPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20131010_2111-1/";
const Char_t* gkFmdDiagonalizationPath   = "/hera/alice/jonderw/PbPb2010/EventPlane_20131010_3222-1/";

const Int_t gkFMDstep=0;
//    /////////////// Add FMD subdetectors ///////////////////
//
//    
//      AliEventPlaneConfiguration * FMD1conf = new AliEventPlaneConfiguration();
//      FMD1conf->SetDetectorType(AliEventPlaneManager::kFMD1);
//      FMD1conf->SetCalibrationStep(gkFMDstep);
//      FMD1conf->SetCalibrationMethod(1);
//      FMD1conf->SetEqualizationMethod(2);
//      FMD1conf->SetCalibrationDetectorName("FMD1");
//      FMD1conf->SetEqualizationDetectorName("FMD1");
//      FMD1conf->SetEventPlaneDetectorName("FMD1");
//      FMD1conf->SetCorrelationDetectorNames("FMD2I","TPC");
//      FMD1conf->SetEqualizationHistPath(gkFmdEqualizationPath);
//      FMD1conf->SetRecenteringHistPath(gkFmdRecenteringPath);
//      FMD1conf->SetCorrelationHistPath(gkFmdDiagonalizationPath);
//      EPmanager->AddEventPlaneConfiguration(FMD1conf);
//
//      AliEventPlaneConfiguration * FMD2Iconf = new AliEventPlaneConfiguration();
//      FMD2Iconf->SetDetectorType(AliEventPlaneManager::kFMD2I);
//      FMD2Iconf->SetCalibrationStep(gkFMDstep);
//      FMD2Iconf->SetCalibrationMethod(1);
//      FMD2Iconf->SetEqualizationMethod(2);
//      FMD2Iconf->SetCalibrationDetectorName("FMD2I");
//      FMD2Iconf->SetEqualizationDetectorName("FMD2I");
//      FMD2Iconf->SetEventPlaneDetectorName("FMD2I");
//      FMD2Iconf->SetCorrelationDetectorNames("FMD1","TPC");
//      FMD2Iconf->SetEqualizationHistPath(gkFmdEqualizationPath);
//      FMD2Iconf->SetRecenteringHistPath(gkFmdRecenteringPath);
//      FMD2Iconf->SetCorrelationHistPath(gkFmdDiagonalizationPath);
//      EPmanager->AddEventPlaneConfiguration(FMD2Iconf);
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

//TChain* CreateESDChain(const char* aDataDir, Int_t aRuns, Int_t offset)
//{
//  // creates chain of files in a given directory or file containing a list.
//  // In case of directory the structure is expected as:
//  // <aDataDir>/<dir0>/AliESDs.root
//  // <aDataDir>/<dir1>/AliESDs.root
//  // ...
//  
//  if (!aDataDir)
//    return 0;
//  
//  Long_t id, size, flags, modtime;
//  if (gSystem->GetPathInfo(aDataDir, &id, &size, &flags, &modtime))
//    {
//      printf("%s not found.\n", aDataDir);
//      return 0;
//    }
//  
//  TChain* chain = new TChain("esdTree");
//  TChain* chaingAlice = 0;
//  
//  if (flags & 2)
//    {
//      TString execDir(gSystem->pwd());
//      TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
//      TList* dirList            = baseDir->GetListOfFiles();
//      Int_t nDirs               = dirList->GetEntries();
//      gSystem->cd(execDir);
//      
//      Int_t count = 0;
//      
//      for (Int_t iDir=0; iDir<nDirs; ++iDir)
//	{
//	  TSystemFile* presentDir = (TSystemFile*) dirList->At(iDir);
//	  if (!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0)
//	    continue;
//	  
//	  if (offset > 0)
//	    {
//	      --offset;
//	      continue;
//	    }
//	  
//	  if (count++ == aRuns)
//	    break;
//	  
//	  TString presentDirName(aDataDir);
//	  presentDirName += "/";
//	  presentDirName += presentDir->GetName();	  
//	  chain->Add(presentDirName + "/AliESDs.root/esdTree");
//	  //  cerr<<presentDirName<<endl;
//	}
//      
//    }
//  else
//    {
//      // Open the input stream
//      ifstream in;
//      in.open(aDataDir);
//      
//      Int_t count = 0;
//      
//      // Read the input list of files and add them to the chain
//      TString esdfile;
//      while(in.good()) {
//	in >> esdfile;
//	if (!esdfile.Contains("root")) continue; // protection
//	
//	if (offset > 0)
//	  {
//	    --offset;
//	    continue;
//	  }
//	
//	if (count++ == aRuns)
//	  break;
//	
//	// add esd file
//	chain->Add(esdfile);
//      }
//      
//      in.close();
//    }
//  
//  return chain;
//}
//

