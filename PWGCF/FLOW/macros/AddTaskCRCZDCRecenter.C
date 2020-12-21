AliAnalysisTask * AddTaskCRCZDCRecenter(Double_t ptMin=0.2,                          //0.2
                             Double_t ptMax=5.,                          //50
                             TString analysisTypeUser="AOD",              //"AOD"
                             Int_t AODfilterBit=768,                      //"768"
                             TString sDataSet="2015",                     //"2015"
                             TString sIntRuns="high",                      //"high"
                             TString sIntRate="any",                      //"any"
                             TString EvTrigger="MB",                      //"MB"
                             Bool_t bSetStoreZDCQVecVtxPos=kTRUE,         //"kTRUE:
                             Bool_t bUseZDC=kTRUE,                       //"kTRUE"
                             TString ZDCCalibFileName="alien:///alice/cern.ch/user/j/jmargutt/15oHI_ZDCcalibVar_CenVtxCen_VtxRbR_Ecom_EcomVtx_Cen_VtxFit7.root",                    //"alien:///alice/cern.ch/user/j/jmargutt/10h_ZDCCalib_VtxCen.root"
                             TString sCorrWeight="TPCmVZuZDCm",           //"TPCmVZuZDCm"
                             Double_t etaMin=-0.8,                        //-0.8
                             Double_t etaMax=0.8,                         //0.8
                             TString Label="15_768",                            //"10_768"
                             TString sCentrEstimator="V0",                //"V0"
                             Double_t dVertexRange=10.,                   //10
                             Double_t dMinClusTPC=70,                     //70
                             Double_t dDCAxy=2.4,                       //1000
                             Double_t dDCAz=3.2,                        //1000
                             Double_t MaxFracSharedTPCCl=0.4,             //0.4
                             Double_t MaxFracSharedITSCl=1,            //1
                             Double_t MaxChi2PerClTPC=4.,                 //4
                             Double_t MaxChi2PerClITS=36.,               //36
                             TString sSelecCharge="",                     //""
                             Bool_t bPtDepDCAxyCut=kFALSE,                //kFALSE
                             Bool_t bRequireITSRefit=kFALSE,              //kFALSE
                             Bool_t bStoreExtraHistoForSubSampling=kFALSE,//kFALSE
                             Double_t DeltaEta=1,                       //1
                             Bool_t bRecZDCVtxRbR=kTRUE,                 //kFALSE swith on for filling phi eta dis vtx
                             TString PtWeightsFileName="alien:///alice/cern.ch/user/j/jmargutt/15o_FB768_hijing_eff_CorSec.root",                //"alien:///alice/cern.ch/user/j/jmargutt/10h_FB768_hijing_eff_CorSec.root"
                             TString sPhiEtaWeight="EtaPhiVtxRbR",                 //"EtaPhiVtxRbR"
                             Bool_t bRemoveSplitMergedTracks=kFALSE,      //kFALSE
                             Bool_t bUseTightPileUp=kFALSE,               //kFALSE
                             Int_t MinMulZN=10,                            //10
                             TString ZDCESEFileName="",                   //""
                             TString CenWeightsFileName="alien:///alice/cern.ch/user/j/jmargutt/15oCentrality2.root",  //"alien:///alice/cern.ch/user/j/jmargutt/10hCentrality.root"
                             TString ZDCRecenterFileName = "", //"alien:///alice/cern.ch/user/s/sqiu/15o_ZDCcalibVar_Step1.root",
                             Int_t bStepZDCRecenter = 3,
                             const char* suffix="") {
  // load libraries
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libOADB.so");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");

  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/OCDB -I$ALICE_ROOT/STEER/macros -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/TRD -I$ALICE_ROOT/ZDC -I$ALICE_ROOT/macros -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/OADB $ALICE_PHYSICS/OADB/macros -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGCF -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/TENDER -I$ALICE_PHYSICS/TENDER/Tender -I$ALICE_PHYSICS/TENDER/TenderSupplies -I$ALICE_PHYSICS/PARfiles -I$ALICE_PHYSICS/PWGCF/FLOW/macros I$ALICE_PHYSICS/PWGPP/ZDC -g ");

  // the manager is static, so get the existing manager via the static method
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    printf("No analysis manager to connect to!\n");
    return NULL;
  }

  // just to see if all went well, check if the input event handler has been connected
  if (!mgr->GetInputEventHandler()) {
    printf("This task requires an input event handler!\n");
    return NULL;
  }

  Int_t nCenBin = 10;
  Double_t centrMin=0.;
  Double_t centrMax=100.;
  Double_t CenBinWidth=10.;
  Bool_t bCalculateCRCInt=kFALSE;
  Bool_t bCalculateCRC2=kFALSE;
  Float_t MaxDevZN=10.;
  Int_t NumCenBins=100;
  Bool_t bCalculateCRC=kTRUE;
  if(analysisTypeUser == "TrackQA") bCalculateCRC=kFALSE;
  Bool_t bCalculateCRCVZ=kTRUE; // Control VZ QVector and Recenter. Will cause error if it is switched off when calculating e.g. CMESPPP()
  TString PhiEtaWeightsFileName="";
  Bool_t bCutsQA=kTRUE; 
  Bool_t bCalculateEbEFlow=kTRUE; 
  Bool_t bDivSigma=kFALSE;
  Bool_t bCalculateCRCZDC=kFALSE;
  Bool_t bCalculateCME=kTRUE;  
  Bool_t bUseVZERO=kFALSE;
  Int_t nHarmonic=2;
  Bool_t bMimicGlobalCuts=kFALSE;
  Bool_t bZDCMCCen=kTRUE;          // Default kTRUE use GetZNCentroidInPbPb, with correction from MC
  Bool_t bCorrSpecZDC=kFALSE;      // Whether using correction file for ZDC energy
  Bool_t bUsePhiEtaCuts=kFALSE;
  Bool_t bSetQAZDC=kTRUE;
  Int_t bCutTPCbound=0;
  Bool_t bCalculateFlow=kTRUE;
  Bool_t bCorrectForBadChannel=kFALSE;
  Bool_t bUsePileUp=kTRUE;
  Bool_t bSpecialVZERORingSelection=kFALSE;
  Bool_t bResetNegativeZDC=kFALSE;
  Bool_t bPhiExclZone=kFALSE;
  Bool_t bTestSin=kFALSE;
  Bool_t bZDCCut=kFALSE;
  Bool_t bUsePtWeights = (PtWeightsFileName.EqualTo("")?kFALSE:kTRUE);
  if(MinMulZN>=13) bZDCCut=kTRUE;
  Bool_t bUseCRCRecenter=kFALSE;
  Float_t ZDCGainAlpha=0.395;
  Int_t CRC2nEtaBins=5;
  Bool_t bCorrectPhiTracklets=kFALSE;
  Bool_t bStoreCalibZDCRecenter = kTRUE;
  Bool_t bStoreQAforDiffEventPlanes=kFALSE;
  Bool_t bFillZNCenDisRbR=kFALSE;
  Bool_t bRequireTOFSignal=kFALSE;             //kFALSE
  
  // define CRC suffix
  TString CRCsuffix = ":CRC";

  TString CentrName = "_";
  CentrName += (Int_t)centrMin;
  CentrName += "-";
  CentrName += (Int_t)centrMax;
  CRCsuffix += CentrName;

  TString pTName = "_";
  Int_t rt = (Int_t)(ptMin*10.);
  Int_t r = (Int_t)(ptMin);
  pTName += ( ptMin < 1. ? Form("0.%i",rt) : Form("%i.%i",r,rt-r*10));
  pTName += "-";
  rt = (Int_t)(ptMax*10.);
  r = (Int_t)(ptMax);
  pTName += ( ptMax < 1. ? Form("0.%i",rt) : Form("%i.%i",r,rt-r*10));
  CRCsuffix += pTName;

  if(!Label.EqualTo("")) {
    TString Appendix = "_";
    Appendix += Label;
    CRCsuffix += Appendix;
  }

  TGrid::Connect("alien://");                                //@Shi

  // create instance of the class: because possible qa plots are added in a second output slot,
  // the flow analysis task must know if you want to save qa plots at the time of class construction
  TString taskFEname = "FlowEventTask";
  taskFEname += CRCsuffix;
  taskFEname += suffix;
  
  // create instance of the class
  UInt_t seed=666;
  Bool_t bCandidates=kFALSE;
  AliAnalysisTaskCRCZDC* taskFE = new AliAnalysisTaskCRCZDC(taskFEname, "", bCutsQA, seed, bCandidates, bStepZDCRecenter);
  taskFE->SetCentralityRange(centrMin,centrMax);
  if(sCentrEstimator=="V0")  taskFE->SetCentralityEstimator(AliAnalysisTaskCRCZDC::kV0M);
  if(sCentrEstimator=="TRK") taskFE->SetCentralityEstimator(AliAnalysisTaskCRCZDC::kTRK);
  if(sCentrEstimator=="CL1") taskFE->SetCentralityEstimator(AliAnalysisTaskCRCZDC::kCL1);
  if(sCentrEstimator=="CL0") taskFE->SetCentralityEstimator(AliAnalysisTaskCRCZDC::kCL0);
  taskFE->SetRejectPileUp(bUsePileUp);
  taskFE->SetRejectPileUpTight(bUseTightPileUp);
  taskFE->SetUseMCCen(bZDCMCCen);
  taskFE->SetZDCGainAlpha(ZDCGainAlpha);
  taskFE->SetResetNegativeZDC(bResetNegativeZDC);
  taskFE->SetFillZNCenDisRbR(bFillZNCenDisRbR);  //@Shi add flag for run by run ZN centroid distribution. Do not turn on for large dataset when running on grid. It takes too much memory
  if (sDataSet == "2010") taskFE->SetDataSet(AliAnalysisTaskCRCZDC::k2010);
  if (sDataSet == "2011") taskFE->SetDataSet(AliAnalysisTaskCRCZDC::k2011);
  if (sDataSet == "2015") taskFE->SetDataSet(AliAnalysisTaskCRCZDC::k2015);
  if (sDataSet == "2015v6") taskFE->SetDataSet(AliAnalysisTaskCRCZDC::k2015v6);
  if (sDataSet == "2015pidfix") taskFE->SetDataSet(AliAnalysisTaskCRCZDC::k2015pidfix);
  if (sDataSet == "2018r") taskFE->SetDataSet(AliAnalysisTaskCRCZDC::k2018r);
  taskFE->SetQAOn(bCutsQA);
  // set the analysis type
  if (analysisTypeUser == "AOD" || analysisTypeUser == "AUTOMATIC") taskFE->SetAnalysisType(AliAnalysisTaskCRCZDC::kAUTOMATIC);
  if (analysisTypeUser == "MCAOD") taskFE->SetAnalysisType(AliAnalysisTaskCRCZDC::kMCAOD);
  if (analysisTypeUser == "MCkine") taskFE->SetAnalysisType(AliAnalysisTaskCRCZDC::kMCkine);
  if (analysisTypeUser == "MCESD") taskFE->SetAnalysisType(AliAnalysisTaskCRCZDC::kMCESD);
  if (analysisTypeUser == "TrackQA") taskFE->SetAnalysisType(AliAnalysisTaskCRCZDC::kTrackQA);
  if (analysisTypeUser == "Tracklets") taskFE->SetAnalysisType(AliAnalysisTaskCRCZDC::kTracklets);
  taskFE->SetCorrectPhiTracklets(bCorrectPhiTracklets);
  // set the trigger selection
  if (EvTrigger == "Cen")
    taskFE->SelectCollisionCandidates(AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB);
  if (EvTrigger == "SemiCen")
    taskFE->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral);
  if (EvTrigger == "MB")
    //taskFE->SelectCollisionCandidates(AliVEvent::kMB);
    taskFE->SelectCollisionCandidates(0);
  if (EvTrigger == "MB" && sDataSet.Contains("2015"))
    //taskFE->SelectCollisionCandidates(AliVEvent::kINT7);
    taskFE->SelectCollisionCandidates(0);
  if (EvTrigger == "Any")
    taskFE->SelectCollisionCandidates(AliVEvent::kAny);

  if(sDataSet=="2010" && !bZDCMCCen) {
    TString ZDCTowerEqFileName = "alien:///alice/cern.ch/user/j/jmargutt/Calib10hZDCEqTowerVtx.root";
    TFile* ZDCTowerEqFile = TFile::Open(ZDCTowerEqFileName,"READ");
    if(!ZDCTowerEqFile) {
      cout << "ERROR: ZDCTowerEqFile not found!" << endl;
      exit(1);
    }
  } // end of if(sDataSet=="2010" && !bZDCMCCen)
  if(sDataSet=="2015") {
    TString ZDCTowerEqFileName = "alien:///alice/cern.ch/user/j/jmargutt/15oHI_EZDCcalib.root";
    TFile* ZDCTowerEqFile = TFile::Open(ZDCTowerEqFileName,"READ");
    gROOT->cd();
    TList* ZDCTowerEqList = (TList*)(ZDCTowerEqFile->FindObjectAny("EZNcalib"));
    if(ZDCTowerEqList) {
      taskFE->SetTowerEqList(ZDCTowerEqList);
      cout << "ZDCTowerEq set (from " <<  ZDCTowerEqFileName.Data() << ")" << endl;
    } else {
      cout << "ERROR: ZDCTowerEqList not found!" << endl;
      exit(1);
    }
    delete ZDCTowerEqFile;
  }
  if(sDataSet=="2015pidfix") {
    TString ZDCTowerEqFileName = "alien:///alice/cern.ch/user/j/jmargutt/15oHIpidfix_EZDCcalib.root";
    TFile* ZDCTowerEqFile = TFile::Open(ZDCTowerEqFileName,"READ");
    gROOT->cd();
    TList* ZDCTowerEqList = (TList*)(ZDCTowerEqFile->FindObjectAny("EZNcalib"));
    if(ZDCTowerEqList) {
      taskFE->SetTowerEqList(ZDCTowerEqList);
      cout << "ZDCTowerEq set (from " <<  ZDCTowerEqFileName.Data() << ")" << endl;
    } else {
      cout << "ERROR: ZDCTowerEqList not found!" << endl;
      exit(1);
    }
    delete ZDCTowerEqFile;
  }
  if(sDataSet=="2018r") { //@Shi add Tower Eq file for 2018r
	TString ZDCTowerEqFileName = "alien:///alice/cern.ch/user/s/sqiu/18r_ZDCgainEqualization.root";
    TFile* ZDCTowerEqFile = TFile::Open(ZDCTowerEqFileName,"READ");
    gROOT->cd();
    TList* ZDCTowerEqList = (TList*)(ZDCTowerEqFile->FindObjectAny("EZNcalib"));
    if(ZDCTowerEqList) {
      taskFE->SetTowerEqList(ZDCTowerEqList);
      cout << "ZDCTowerEq set (from " <<  ZDCTowerEqFileName.Data() << ")" << endl;
    } else {
      cout << "ERROR: ZDCTowerEqList not found!" << endl;
      exit(1);
    }
    delete ZDCTowerEqFile;
  }
  if(bCorrectForBadChannel) {
    TString ZDCBadTowerFileName = "alien:///alice/cern.ch/user/j/jmargutt/ZDCCalibBadChannel.root";
    TFile* ZDCBadTowerFile = TFile::Open(ZDCBadTowerFileName,"READ");
    gROOT->cd();
    TList* ZDCBadTowerList = (TList*)(ZDCBadTowerFile->FindObjectAny("resp"));
    if(ZDCBadTowerList) {
      taskFE->SetBadTowerCalibList(ZDCBadTowerList);
      cout << "BadTowerCalibList set (from " <<  ZDCBadTowerFileName.Data() << ")" << endl;
    } else {
      cout << "ERROR: BadTowerCalibList not found!" << endl;
      exit(1);
    }
    delete ZDCBadTowerFile;
  }
  if(bCorrSpecZDC) {
    TString ZDCRecFileName = "alien:///alice/cern.ch/user/j/jmargutt/";
    if(bCorrectForBadChannel) ZDCRecFileName += "15o_ZDCSpectraCorr_BadCh_3.root";
    TFile* ZDCRecFile = TFile::Open(ZDCRecFileName,"READ");
    if(!ZDCRecFile) {
      cout << "ERROR: ZDC Spectra Calibration not found!" << endl;
      exit(1);
    }
    gROOT->cd();
    TList* ZDCRecList = (TList*)(ZDCRecFile->FindObjectAny("ZDCSpectraCorr"));
    if(ZDCRecList) {
      taskFE->SetZDCSpectraCorrList(ZDCRecList);
      cout << "ZDC Spectra Calibration set (from " <<  ZDCRecFileName.Data() << ")" << endl;
    }
    else {
      cout << "ERROR: ZDCSpectraCorrList not found!" << endl;
      exit(1);
    }
    delete ZDCRecFile;
  }
  if(sDataSet=="2015") {
    // load VZERO gain equalization
    TString VZEROGainEqFileName = "alien:///alice/cern.ch/user/j/jmargutt/15oHI_VZEROEqGain.root";
    TFile* VZEROGainEqFile = TFile::Open(VZEROGainEqFileName,"READ");
    if(!VZEROGainEqFile) {
      cout << "ERROR: VZERO gain equalization file not found!" << endl;
      exit(1);
    }
    gROOT->cd();
    TList* VZEROGainEqList = (TList*)(VZEROGainEqFile->FindObjectAny("VZEROEqGain"));
    if(VZEROGainEqList) {
      taskFE->SetVZEROGainEqList(VZEROGainEqList);
      cout << "VZERO gain equalization set (from " <<  VZEROGainEqFileName.Data() << ")" << endl;
    }
    else {
      cout << "ERROR: VZERO gain equalization list not found!" << endl;
      exit(1);
    }
    delete VZEROGainEqFile;
    // load VZERO Q-vector re-centering
    //TString VZEROQVecRecFileName = "alien:///alice/cern.ch/user/m/mhaque/jacopo/15oHI_VZEROQVecRec.root";
    TString VZEROQVecRecFileName = "alien:///alice/cern.ch/user/j/jmargutt/15oHI_VZEROQVecRec.root"; //@Shi
    TFile* VZEROQVecRecFile = TFile::Open(VZEROQVecRecFileName,"READ");
    if(!VZEROQVecRecFile) {
      cout << "ERROR: VZERO Q-vector re-centering file not found!" << endl;
      exit(1);
    }
    gROOT->cd();
    TList* VZEROQVecRecList = (TList*)(VZEROQVecRecFile->FindObjectAny("VZEROEqGain"));
    if(VZEROQVecRecList) {
      taskFE->SetVZEROQVecRecList(VZEROQVecRecList);
      cout << "VZERO Q-vector re-centering set (from " <<  VZEROQVecRecFileName.Data() << ")" << endl;
    }
    else {
      cout << "ERROR: Q-vector re-centering list not found!" << endl;
      exit(1);
    }
    delete VZEROQVecRecFile;
  }
  // set which rings of VZEROs to use
  if(bSpecialVZERORingSelection) taskFE->SetWhichVZERORings(1,2,7,8);

  //@Shi set ZDC recentering (begin)
  if (bStepZDCRecenter > 0 && bStepZDCRecenter <=2) {
	TFile* ZDCRecenterFile = TFile::Open(ZDCRecenterFileName, "READ");
	if(ZDCRecenterFile) {
	  TList* ZDCRecenterList = (TList*)(ZDCRecenterFile->FindObjectAny("Q Vectors")); // hardcoded TList AQ Vectors
	  if(ZDCRecenterList) {
		taskFE->SetZDCCalibList(ZDCRecenterList);
	  } else {
		cout << "ERROR: ZDCRecenterList do not exist!" << endl;
		exit(1);
	  }
	} else {
	  cout << "ERROR: if bStepZDCRecenter larger than 0, ZDCRecenterFile should exist!" << endl;
	  exit(1);
	}
	delete ZDCRecenterFile;
  }
  
  if (bStepZDCRecenter >= 3) {
	TFile* ZDCRecenterFileStep3CommonPart = TFile::Open("alien:///alice/cern.ch/user/s/sqiu/15o_ZDCcalibVar_Step3_commonPart.root", "READ");
	if(ZDCRecenterFileStep3CommonPart) {
	  TList* ZDCRecenterListStep3CommonPart = (TList*)(ZDCRecenterFileStep3CommonPart->FindObjectAny("Q Vectors"));
      if(ZDCRecenterListStep3CommonPart) {
	    taskFE->SetZDCCalibListStep3CommonPart(ZDCRecenterListStep3CommonPart);
	  } else {
	    cout << "ERROR: bStepZDCRecenter >= 3 ZDCRecenterListStep3CommonPart do not exist!" << endl;
	    exit(1);
	  }
	} else {
	  cout << "ERROR: if bStepZDCRecenter larger than 2, ZDCRecenterFileStep3CommonPart should exist!" << endl;
	  exit(1);
	}
	delete ZDCRecenterFileStep3CommonPart;
  }
  
  taskFE->SetStepZDCRecenter(bStepZDCRecenter);
  taskFE->SetStoreCalibZDCRecenter(bStoreCalibZDCRecenter);
  //@Shi set ZDC recentering (end)
  // add the task to the manager
  mgr->AddTask(taskFE);

  // define the event cuts object
  AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("EventCuts");
  cutsEvent->SetCheckPileup(kFALSE);
  // configure some event cuts, starting with centrality
  if(analysisTypeUser == "MCkine" || analysisTypeUser == "ESD") {
    // method used for centrality determination
    if(sCentrEstimator=="V0")  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
    if(sCentrEstimator=="TRK") cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kTPConly);
    if(sCentrEstimator=="CL1") cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kSPD1clusters);
    if (sDataSet == "2010" || sDataSet == "2011") {
      cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
    }
    if (sDataSet.Contains("2015")) {
      cutsEvent->SetCentralityPercentileRange(centrMin,centrMax,kTRUE);
    }
    // vertex-z cut
    cutsEvent->SetPrimaryVertexZrange(-dVertexRange,dVertexRange);
    cutsEvent->SetQA(bCutsQA);
  }
  else if (analysisTypeUser == "AOD" || analysisTypeUser == "TrackQA" || analysisTypeUser == "Tracklets" || analysisTypeUser == "MCAOD") {
    if (sDataSet == "2010" || sDataSet == "2011") {
      cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
    }
    if (sDataSet.Contains("2015")) {
      cutsEvent->SetCentralityPercentileRange(centrMin,centrMax,kTRUE);
    }
    // method used for centrality determination
    if(sCentrEstimator=="V0")  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
    if(sCentrEstimator=="TRK") cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kTPConly);
    if(sCentrEstimator=="CL1") cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kSPD1clusters);
    AliFlowTrackCuts* RefMultCuts = new AliFlowTrackCuts("RefMultCuts");
    RefMultCuts->SetParamType(AliFlowTrackCuts::kAODFilterBit);
    RefMultCuts->SetAODfilterBit(768);
    RefMultCuts->SetMinimalTPCdedx(-999999999);
    RefMultCuts->SetMaxDCAToVertexXY(1000.);
    RefMultCuts->SetMaxDCAToVertexZ(1000.);
    RefMultCuts->SetMinNClustersTPC(70.);
    RefMultCuts->SetMinChi2PerClusterTPC(0.1);
    RefMultCuts->SetMaxChi2PerClusterTPC(4.);
    RefMultCuts->SetPtRange(0.2,20.2);
    RefMultCuts->SetEtaRange(-0.8,0.8);
    RefMultCuts->SetAcceptKinkDaughters(kFALSE);
    cutsEvent->SetRefMultCuts(RefMultCuts);
    cutsEvent->SetRefMultMethod(AliFlowEventCuts::kTPConly);
    // vertex-z cut
    cutsEvent->SetPrimaryVertexZrange(-dVertexRange,dVertexRange);
    // enable the qa plots
    cutsEvent->SetQA(bCutsQA);
    // explicit multiplicity outlier cut
    cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
    if (sDataSet == "2011") cutsEvent->SetLHC11h(kTRUE);
    if (sDataSet == "2010") cutsEvent->SetLHC10h(kTRUE);
  }

  // pass these cuts to your flow event task
  taskFE->SetCutsEvent(cutsEvent);
  AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts("RP cuts");
  AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("POI cuts");

  if (analysisTypeUser == "MCkine") {
    // Track cuts for RPs
    cutsRP->SetParamType(AliFlowTrackCuts::kMC);
    cutsRP->SetCutMC(kTRUE);
    cutsRP->SetPtRange(ptMin,ptMax);
    cutsRP->SetEtaRange(etaMin,etaMax);
    cutsRP->SetQA(bCutsQA);
    if(bUseVZERO) {
      cutsRP->SetEtaRange(-10.,+10.);
      cutsRP->SetEtaGap(-1.,1.);
    }
    // Track cuts for POIs
    cutsPOI->SetParamType(AliFlowTrackCuts::kMC);
    cutsPOI->SetCutMC(kTRUE);
    cutsPOI->SetPtRange(ptMin,ptMax);
    cutsPOI->SetEtaRange(etaMin,etaMax);
    cutsPOI->SetQA(bCutsQA);
  }
  if (analysisTypeUser == "AOD" || analysisTypeUser == "MCAOD" || analysisTypeUser == "TrackQA" || analysisTypeUser == "Tracklets") {
    // Track cuts for RPs
    if(bUseVZERO) {
      if (sDataSet == "2011")
        cutsRP->SetParamType(AliFlowTrackCuts::kDeltaVZERO);
      if (sDataSet == "2010")
        cutsRP->SetParamType(AliFlowTrackCuts::kBetaVZERO);
      if (sDataSet.Contains("2015"))
        cutsRP->SetParamType(AliFlowTrackCuts::kKappaVZERO);
      cutsRP->SetEtaRange(-10.,+10.);
      cutsRP->SetEtaGap(-1.,1.);
      cutsRP->SetPhiMin(0.);
      cutsRP->SetPhiMax(TMath::TwoPi());
      // options for the reweighting
      cutsRP->SetVZEROgainEqualizationPerRing(bUseVZERO);
      cutsRP->SetApplyRecentering(bUseVZERO);
      cutsRP->SetDivSigma(bDivSigma);
      if (analysisTypeUser == "MCAOD")
        cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();
    } else {
      cutsRP->SetParamType(AliFlowTrackCuts::kAODFilterBit);
      cutsRP->SetAODfilterBit(AODfilterBit);
      cutsRP->SetMinimalTPCdedx(-999999999);
      cutsRP->SetMaxDCAToVertexXY(dDCAxy);
      cutsRP->SetMaxDCAToVertexZ(dDCAz);
      cutsRP->SetMinNClustersTPC(dMinClusTPC);
      cutsRP->SetMinChi2PerClusterTPC(0.1);
      cutsRP->SetMaxChi2PerClusterTPC(MaxChi2PerClTPC);
      cutsRP->SetCutChi2PerClusterITS(MaxChi2PerClITS);
      cutsRP->SetPtRange(ptMin,ptMax);
      cutsRP->SetEtaRange(etaMin,etaMax);
      cutsRP->SetAcceptKinkDaughters(kFALSE);
      cutsRP->SetMaxFracSharedTPCCluster(MaxFracSharedTPCCl);
      cutsRP->SetMaxFracSharedITSCluster(MaxFracSharedITSCl);
      if(bCutTPCbound==1) cutsRP->SetCutTPCSecbound(kTRUE,ptMin); // new cut for LHC15o
      if(bCutTPCbound==2) cutsRP->SetCutTPCSecboundVar(kTRUE); // new cut for LHC15o
      cutsRP->SetQA(bCutsQA);
    }
    // Track cuts for POIs
    cutsPOI->SetParamType(AliFlowTrackCuts::kAODFilterBit);
    cutsPOI->SetAODfilterBit(AODfilterBit);
    cutsPOI->SetMinimalTPCdedx(-999999999);
    cutsPOI->SetMaxDCAToVertexXYAOD(dDCAxy);
    cutsPOI->SetMaxDCAToVertexZAOD(dDCAz);
    cutsPOI->SetMinNClustersTPC(dMinClusTPC);
    cutsPOI->SetMinChi2PerClusterTPC(0.1);
    cutsPOI->SetMaxChi2PerClusterTPC(MaxChi2PerClTPC);
    cutsPOI->SetCutChi2PerClusterITS(MaxChi2PerClITS);
    if(bMimicGlobalCuts) {
      cutsPOI->SetMinNClustersTPC(50);
      cutsPOI->SetCutCrossedTPCRows(70,0.8);
      cutsPOI->SetRequireITSRefit(kTRUE);
      cutsPOI->SetMaxDCAToVertexXYPtDepAOD(kTRUE);
      cutsPOI->SetCutGoldenChi2(kTRUE);
      cutsPOI->SetCutChi2PerClusterITS(36.);
      cutsPOI->SetCutITSClusterGlobal(kTRUE);
    }
    if(bRequireITSRefit) {
      cutsPOI->SetRequireITSRefit(kTRUE);
      cutsPOI->SetCutChi2PerClusterITS(36.);
    }
    if(bPtDepDCAxyCut) cutsPOI->SetMaxDCAToVertexXYPtDepAOD(kTRUE);
    cutsPOI->SetPtRange(ptMin,ptMax);
    cutsPOI->SetEtaRange(etaMin,etaMax);
    cutsPOI->SetAcceptKinkDaughters(kFALSE);
    cutsPOI->SetMaxFracSharedTPCCluster(MaxFracSharedTPCCl);
    cutsPOI->SetMaxFracSharedITSCluster(MaxFracSharedITSCl);
    cutsPOI->SetRequireTOFSignal(bRequireTOFSignal);
    if(bCutTPCbound==1) cutsPOI->SetCutTPCSecbound(kTRUE,ptMin); // new cut for LHC15o
    if(bCutTPCbound==2) cutsPOI->SetCutTPCSecboundVar(kTRUE); // new cut for LHC15o
    cutsPOI->SetQA(bCutsQA);
  }
  if (analysisTypeUser == "ESD") {
    cutsRP = AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
    cutsPOI = AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
  }

  taskFE->SetCutsRP(cutsRP);
  taskFE->SetCutsPOI(cutsPOI);
  taskFE->SetSubeventEtaRange(-10.,-1.,1.,10.);
  if (analysisTypeUser == "MCkine")
    taskFE->SetSubeventEtaRange(-3.7,-1.7,2.8,5.1);

  // get the default name of the output file ("AnalysisResults.root")
  TString file = "AnalysisResults.root";

  // get the common input container from the analysis manager
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();

  // create a data container for the output of the flow event task
  TString taskFECname = "FlowEventContainer";
  taskFECname += CRCsuffix;
  taskFECname += suffix;
  AliAnalysisDataContainer *coutputFE = mgr->CreateContainer(taskFECname,
                                                             AliFlowEventSimple::Class(),
                                                             AliAnalysisManager::kExchangeContainer);
  // connect the input data to the flow event task
  mgr->ConnectInput(taskFE,0,cinput);
  // and connect the output to the flow event task
  mgr->ConnectOutput(taskFE,1,coutputFE);

  // QA OUTPUT CONTAINER
  TString taskFEQAname = file;
  taskFEQAname += ":CutsQA";
  taskFEQAname += CRCsuffix;
  taskFEQAname += suffix;
  AliAnalysisDataContainer* coutputFEQA = mgr->CreateContainer(taskFEQAname.Data(),
                                                               TList::Class(),
                                                               AliAnalysisManager::kOutputContainer,
                                                               taskFEQAname);
  // and connect the qa output container to the flow event.
  // this container will be written to the output file
  mgr->ConnectOutput(taskFE,2,coutputFEQA);
  
  AliAnalysisDataContainer* coutputFERecenter1;
  AliAnalysisDataContainer* coutputFERecenter2;
  if (bStepZDCRecenter >= 0) {
	  // OUTPUT CONTAINER TO SAVE ZDC RECENTERING
	  TString taskFERecenter1name = file;
	  taskFERecenter1name += ":RecenterList1";
	  taskFERecenter1name += CRCsuffix;
	  taskFERecenter1name += suffix;
	  coutputFERecenter1 = mgr->CreateContainer(taskFERecenter1name.Data(),
																   TList::Class(),
																   AliAnalysisManager::kOutputContainer,
																   taskFERecenter1name);
	  mgr->ConnectOutput(taskFE,3,coutputFERecenter1);
	  
	  TString taskFERecenter2name = file;
	  taskFERecenter2name += ":RecenterList2";
	  taskFERecenter2name += CRCsuffix;
	  taskFERecenter2name += suffix;
	  coutputFERecenter2 = mgr->CreateContainer(taskFERecenter2name.Data(),
																   TList::Class(),
																   AliAnalysisManager::kOutputContainer,
																   taskFERecenter2name);
	  mgr->ConnectOutput(taskFE,4,coutputFERecenter2);
	  
	  TString taskFERecenter3name = file;
	  taskFERecenter3name += ":RecenterList3";
	  taskFERecenter3name += CRCsuffix;
	  taskFERecenter3name += suffix;
	  coutputFERecenter3 = mgr->CreateContainer(taskFERecenter3name.Data(),
																   TList::Class(),
																   AliAnalysisManager::kOutputContainer,
																   taskFERecenter3name);
	  mgr->ConnectOutput(taskFE,5,coutputFERecenter3);
  }
  
  return taskFE;
}
