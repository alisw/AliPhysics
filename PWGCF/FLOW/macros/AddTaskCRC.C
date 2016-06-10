AliAnalysisTask * AddTaskCRC(Double_t ptMin=0.2,
                             Double_t ptMax=5.0,
                             TString analysisTypeUser="AOD",
                             Int_t AODfilterBit=768,
                             TString sDataSet="2010",
                             TString EvTrigger="MB",
                             Bool_t bCalculateEbEFlow=kFALSE,
                             Bool_t bUseCRCRecenter,
                             TString QVecWeightsFileName,
                             Bool_t bCalculateCME=kFALSE,
                             Bool_t bUseVZERO=kFALSE,
                             Bool_t bCalculateCRCVZ=kFALSE,
                             Bool_t bUseZDC=kFALSE,
                             TString ZDCCalibFileName,
                             Bool_t bDivSigma=kFALSE,
                             Bool_t bCalculateCRCZDC=kFALSE,
                             TString sCorrWeight="TPCmVZuZDCu",
                             Bool_t bZDCMCCen=kTRUE,
                             Float_t ZDCGainAlpha=0.395,
                             Bool_t bCutsQA=kFALSE,
                             TString Label="",
                             TString sCentrEstimator="V0",
                             Double_t dVertexRange=10.,
                             Double_t dDCAxy=2.4,
                             Double_t dDCAz=3.2,
                             Double_t dMinClusTPC=70,
                             Bool_t bCalculateFlow=kFALSE,
                             Int_t NumCenBins=100,
                             Bool_t bUsePtWeights=kFALSE,
                             TString PtWeightsFileName="",
                             Bool_t bUseEtaWeights=kFALSE,
                             TString EtaWeightsFileName="",
                             Bool_t bSetQAZDC=kFALSE,
                             Int_t MinMulZN=1,
                             TString ZDCESEFileName="",
                             Bool_t bCenFlattening=kTRUE,
                             TString CenWeightsFileName="",
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
 Int_t nHarmonic=1;
 Int_t CRC2nEtaBins=6;
 Bool_t bCalculateCRC2=kFALSE;
 Float_t MaxDevZN=10.;
 Bool_t bUsePhiEtaWeights=kFALSE;
 TString PhiEtaWeightsFileName="";
 
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
 
 Double_t etaMin=-0.8;
 Double_t etaMax=0.8;
 
 // create instance of the class: because possible qa plots are added in a second output slot,
 // the flow analysis task must know if you want to save qa plots at the time of class construction
 TString taskFEname = "FlowEventTask";
 taskFEname += CRCsuffix;
 taskFEname += suffix;
 // create instance of the class
 AliAnalysisTaskCRCZDC* taskFE = new AliAnalysisTaskCRCZDC(taskFEname, "", bCutsQA);
 taskFE->SetCentralityRange(centrMin,centrMax);
 taskFE->SetCentralityEstimator("V0M");
 taskFE->SetRejectPileUp(kTRUE);
 taskFE->SetUseMCCen(bZDCMCCen);
 taskFE->SetZDCGainAlpha(ZDCGainAlpha);
 taskFE->SetDataSet(sDataSet);
 taskFE->SetQAOn(bCutsQA);
 // set the analysis type
 TString analysisType = "AUTOMATIC";
 if (analysisTypeUser != "") analysisType = analysisTypeUser;
 if (analysisTypeUser == "AOD") analysisType = "AUTOMATIC";
 taskFE->SetAnalysisType(analysisType);
 // set the trigger selection
 if (EvTrigger == "Cen")
  taskFE->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
 else if (EvTrigger == "SemiCen")
  taskFE->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kSemiCentral);
 else if (EvTrigger == "MB")
  taskFE->SelectCollisionCandidates(AliVEvent::kMB);
 if (EvTrigger == "MB" && sDataSet == "2015")
   taskFE->SelectCollisionCandidates(AliVEvent::kINT7);
 else if (EvTrigger == "Any")
  taskFE->SelectCollisionCandidates(AliVEvent::kAny);
  // add the task to the manager
  mgr->AddTask(taskFE);
 
 // define the event cuts object
 AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("EventCuts");
  // configure some event cuts, starting with centrality
  if(analysisTypeUser == "MCkine" || analysisTypeUser == "MCAOD" || analysisTypeUser == "ESD") {
    // method used for centrality determination
    if(sCentrEstimator=="V0")  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
    if(sCentrEstimator=="TPC") cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kTPConly);
    if(sCentrEstimator=="CL1") cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kSPD1clusters);
    if (sDataSet == "2010" || sDataSet == "2011") {
      cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
    }
    if (sDataSet == "2015") {
      cutsEvent->SetCentralityPercentileRange(centrMin,centrMax,kTRUE);
    }
      cutsEvent->SetPrimaryVertexZrange(-dVertexRange,dVertexRange);
      cutsEvent->SetQA(bCutsQA);
  }
 else if (analysisTypeUser == "AOD") {
   if (sDataSet == "2010" || sDataSet == "2011") {
     cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
   }
   if (sDataSet == "2015") {
     cutsEvent->SetCentralityPercentileRange(centrMin,centrMax,kTRUE);
   }
  // method used for centrality determination
  if(sCentrEstimator=="V0")  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  if(sCentrEstimator=="TPC") cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kTPConly);
  if(sCentrEstimator=="CL1") cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kSPD1clusters);
  AliFlowTrackCuts* RefMultCuts = new AliFlowTrackCuts("RefMultCuts");
  RefMultCuts->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  RefMultCuts->SetAODfilterBit(AODfilterBit);
  RefMultCuts->SetMinimalTPCdedx(-999999999);
  RefMultCuts->SetMaxDCAToVertexXY(dDCAxy);
  RefMultCuts->SetMaxDCAToVertexZ(dDCAz);
  RefMultCuts->SetMinNClustersTPC(dMinClusTPC);
  RefMultCuts->SetPtRange(ptMin,ptMax);
  RefMultCuts->SetEtaRange(etaMin,etaMax);
  cutsEvent->SetRefMultCuts(RefMultCuts);
  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kTPConly);
  // vertex-z cut
  cutsEvent->SetPrimaryVertexZrange(-dVertexRange,dVertexRange);
  // explicit multiplicity outlier cut
   if (sDataSet == "2010" || sDataSet == "2011") {
     cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
     if (sDataSet == "2011")
       cutsEvent->SetLHC11h(kTRUE);
     else if (sDataSet == "2010")
       cutsEvent->SetLHC10h(kTRUE);
   }
   if (sDataSet == "2015") {
     cutsEvent->SetCutTPCmultiplicityOutliersAOD(kFALSE);
   }
   // enable the qa plots
   cutsEvent->SetQA(bCutsQA);
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
 if (analysisTypeUser == "AOD" || analysisTypeUser == "MCAOD") {
  // Track cuts for RPs
  if(bUseVZERO) {
   if (sDataSet == "2011")
    cutsRP->SetParamType(AliFlowTrackCuts::kDeltaVZERO);
   if (sDataSet == "2010")
    cutsRP->SetParamType(AliFlowTrackCuts::kBetaVZERO);
   if (sDataSet == "2015")
    cutsRP->SetParamType(AliFlowTrackCuts::kVZERO);
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
   cutsRP->SetPtRange(ptMin,ptMax);
   cutsRP->SetEtaRange(etaMin,etaMax);
   cutsRP->SetAcceptKinkDaughters(kFALSE);
   cutsRP->SetQA(bCutsQA);
  }
  // Track cuts for POIs
  cutsPOI->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  cutsPOI->SetAODfilterBit(AODfilterBit);
  cutsPOI->SetMinimalTPCdedx(-999999999);
  cutsPOI->SetMaxDCAToVertexXY(dDCAxy);
  cutsPOI->SetMaxDCAToVertexZ(dDCAz);
  cutsPOI->SetMinNClustersTPC(dMinClusTPC);
  cutsPOI->SetPtRange(ptMin,ptMax);
  cutsPOI->SetEtaRange(etaMin,etaMax);
  cutsPOI->SetAcceptKinkDaughters(kFALSE);
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

 //TString ParticleWeightsFileName = "ParticleWeights2D_FullLHC10h_2030.root";
 
 // create the flow analysis tasks
 TString taskCRCname = "AnalysisTask";
 taskCRCname += CRCsuffix;
 taskCRCname += suffix;
 AliAnalysisTaskCRC *taskQC = new AliAnalysisTaskCRC(taskCRCname, bUsePhiEtaWeights);
 // set number of centrality bins
 taskQC->SetnCenBin(nCenBin);
 taskQC->SetCenBinWidth(CenBinWidth);
 taskQC->SetDataSet(sDataSet);
 // set thei triggers
 if (EvTrigger == "Cen")
  taskQC->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
 else if (EvTrigger == "MB")
  taskQC->SelectCollisionCandidates(AliVEvent::kMB);
 if (EvTrigger == "MB" && sDataSet == "2015")
  taskQC->SelectCollisionCandidates(AliVEvent::kINT7);
 else if (EvTrigger == "Any")
  taskQC->SelectCollisionCandidates(AliVEvent::kAny);
 // and set the correct harmonic n
 taskQC->SetHarmonic(nHarmonic);
 // set standard flow settings
 taskQC->SetCalculateDiffFlow(kFALSE);
 taskQC->SetCalculateDiffFlowVsEta(kFALSE);
 taskQC->SetCalculateMixedHarmonics(kFALSE); // calculate all multi-partice mixed-harmonics correlators
 taskQC->SetStoreDistributions(kFALSE);
 taskQC->SetFillMultipleControlHistograms(kFALSE);
 taskQC->SetBookOnlyBasicCCH(kTRUE); // book only basic common control histograms
 //  CRC settings
 taskQC->SetStoreVarious(kTRUE);
 taskQC->SetCalculateCRC(kTRUE);
 taskQC->SetCalculateCRC2(bCalculateCRC2);
 taskQC->SetCalculateCRCVZ(bCalculateCRCVZ);
 taskQC->SetCalculateCRCZDC(bCalculateCRCZDC);
 taskQC->SetCalculateEbEFlow(bCalculateEbEFlow);
 taskQC->SetCRC2nEtaBins(CRC2nEtaBins);
 taskQC->SetCalculateCME(bCalculateCME);
 taskQC->SetCalculateFlowQC(bCalculateFlow);
 taskQC->SetCalculateFlowZDC(bCalculateFlow);
 taskQC->SetCalculateFlowVZ(bCalculateFlow);
 taskQC->SetFlowQCCenBin(NumCenBins);
 taskQC->SetUseVZERO(bUseVZERO);
 taskQC->SetUseZDC(kTRUE);
  if (ZDCCalibFileName != "" && bUseZDC) {
    taskQC->SetRecenterZDC(kTRUE);
  }
 taskQC->SetNUAforCRC(kTRUE);
 taskQC->SetCRCEtaRange(-0.8,0.8);
 taskQC->SetUseCRCRecenter(bUseCRCRecenter);
 taskQC->SetDivSigma(bDivSigma);
 taskQC->SetInvertZDC(bUseZDC);
 taskQC->SetCorrWeight(sCorrWeight);
 taskQC->SetQAZDCCuts(bSetQAZDC);
 taskQC->SetMinMulZN(MinMulZN);
 taskQC->SetMaxDevZN(MaxDevZN);
  if(bSetQAZDC) {
    TFile* ZDCESEFile = TFile::Open(ZDCESEFileName,"READ");
    if(!ZDCESEFile) {
      cout << "ERROR: ZDCESEFile not found!" << endl;
      exit(1);
    }
    TList* ZDCESEList = dynamic_cast<TList*>(ZDCESEFile->FindObjectAny("ZDCESE"));
    if(ZDCESEList) {
      taskQC->SetZDCESEList(ZDCESEList);
      cout << "ZDCESE set (from " <<  ZDCESEFileName.Data() << ")" << endl;
    }
    else {
      cout << "ERROR: ZDCESEList not found!" << endl;
      exit(1);
    }
  } // end of if(bSetQAZDC)
 
 if(bCenFlattening && sDataSet=="2011") {
  TFile* CenWeightsFile = TFile::Open(CenWeightsFileName,"READ");
  if(!CenWeightsFile) {
   cout << "ERROR: CenWeightsFile not found!" << endl;
   exit(1);
  }
  TCanvas* cav = dynamic_cast<TCanvas*>(CenWeightsFile->Get("Canvas_1"));
  TH1D* CenHist = dynamic_cast<TH1D*>(cav->GetPrimitive("Centrality"));
  if(CenHist) {
   taskQC->SetCenWeightsHist(CenHist);
   cout << "Centrality weights set (from " <<  CenWeightsFileName.Data() << ")" << endl;
  }
  else {
   cout << "ERROR: CenHist not found!" << endl;
   exit(1);
  }
 } // end of if(bCenFlattening)
 
 if(bUsePtWeights) {
  taskQC->SetUsePtWeights(bUsePtWeights);
  TFile* PtWeightsFile = TFile::Open(PtWeightsFileName,"READ");
  if(!PtWeightsFile) {
   cout << "ERROR: PtWeightsFile not found!" << endl;
   exit(1);
  }
  for(Int_t c=0; c<10; c++) {
   TH1D* PtWeightsHist = dynamic_cast<TH1D*>(PtWeightsFile->Get(Form("eff_unbiased_%d",c)));
   if(PtWeightsHist) {
    taskQC->SetPtWeightsHist(PtWeightsHist,c);
    cout << "Pt weights centr. "<<c<<" set (from " <<  PtWeightsFileName.Data() << ")" << endl;
   }
   else {
    cout << "ERROR: PtWeightsHist not found!" << endl;
    exit(1);
   }
  }
 } // end of if(bUsePtWeights)
  
  if(bUseEtaWeights) {
    taskQC->SetUseEtaWeights(bUseEtaWeights);
    TFile* EtaWeightsFile = TFile::Open(EtaWeightsFileName,"READ");
    if(!EtaWeightsFile) {
      cout << "ERROR: EtaWeightsFile not found!" << endl;
      exit(1);
    }
    for(Int_t c=0; c<10; c++) {
      for(Int_t b=0; b<21; b++) {
        for(Int_t k=0; k<2; k++) {
          TH1D* EtaWeightsHist = dynamic_cast<TH1D*>(EtaWeightsFile->Get(Form("Eta_Weights_cen%d_ptbin%d_ch%d",c,b+1,k)));
          if(EtaWeightsHist) {
            taskQC->SetEtaWeightsHist(EtaWeightsHist,c,b,k);
            printf("Eta_Weights_cen%d_ptbin%d_ch%d set (from %s)\n",c,b+1,k,EtaWeightsFileName.Data());
          }
          else {
            printf("ERROR: Eta_Weights_cen%d_ptbin%d_ch%d not found! \n",c,b+1,k);
            exit(1);
          }
        }
      }
    }
  } // end of if(bUseEtaWeights)
 
  if(bUseCRCRecenter) {
    TFile* QVecWeightsFile = TFile::Open(QVecWeightsFileName,"READ");
    if(!QVecWeightsFile) {
      cout << "ERROR: QVecWeightsFile not found!" << endl;
      exit(1);
    }
    TList* QVecWeightsList = dynamic_cast<TList*>(QVecWeightsFile->FindObjectAny("Q Vectors"));
    if(QVecWeightsList) {
      taskQC->SetQVecList(QVecWeightsList);
      cout << "Q Vector weights set (from " <<  QVecWeightsFileName.Data() << ")" << endl;
    }
    else {
      cout << "ERROR: QVecWeightsList not found!" << endl;
      exit(1);
    }
  } // end of if(bUseCRCRecenter)
 if(ZDCCalibFileName != "" && bUseZDC) {
  TFile* ZDCCalibFile = TFile::Open(ZDCCalibFileName,"READ");
  if(!ZDCCalibFile) {
   cout << "ERROR: ZDC calibration not found!" << endl;
   exit(1);
  }
  TList* ZDCCalibList = dynamic_cast<TList*>(ZDCCalibFile->FindObjectAny("Q Vectors"));
  if(ZDCCalibList) {
   taskQC->SetCRCZDCCalibList(ZDCCalibList);
   cout << "ZDC calibration set (from " <<  ZDCCalibFileName.Data() << ")" << endl;
  }
  else {
   cout << "ERROR: ZDCCalibList not found!" << endl;
   exit(1);
  }
 } // end of if(bUseZDC)
 taskQC->SetUsePhiEtaWeights(bUsePhiEtaWeights);
 if(bUsePhiEtaWeights) {
  TFile* PhiEtaWeightsFile = TFile::Open(PhiEtaWeightsFileName,"READ");
  if(!PhiEtaWeightsFile) {
   cout << "ERROR: PhiEtaWeightsFile not found!" << endl;
   exit(1);
  }
  TList* PhiEtaWeightsList = dynamic_cast<TList*>(PhiEtaWeightsFile->FindObjectAny("PhiEta Weights"));
  if(PhiEtaWeightsList) {
   taskQC->SetWeightsList(PhiEtaWeightsList);
   cout << "PhiEta weights set (from " <<  PhiEtaWeightsFileName.Data() << ")" << endl;
  }
  else {
   cout << "ERROR: PhiEtaWeightsList not found!" << endl;
   exit(1);
  }
 } // end of if(bUsePhiEtaWeights)

 // connect the task to the analysis manager
 mgr->AddTask(taskQC);
 
 // initialize output name
 TString outputQC = file;
 outputQC += CRCsuffix;
 outputQC += suffix;
 // create and connect the output containers
 AliAnalysisDataContainer *coutputQC = mgr->CreateContainer(outputQC.Data(),
                                                            TList::Class(),
                                                            AliAnalysisManager::kOutputContainer,
                                                            outputQC);
 // connect the output of the flow event task to the flow analysis task
 mgr->ConnectInput(taskQC, 0, coutputFE);
 // and connect the output of the flow analysis task to the output container
 // which will be written to the output file
 mgr->ConnectOutput(taskQC, 1, coutputQC);
 
 return taskQC;
}

