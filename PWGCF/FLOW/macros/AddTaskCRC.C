AliAnalysisTask * AddTaskCRC(Int_t nHarmonic,
                             Double_t centrMin,
                             Double_t centrMax,
                             Int_t nCenBin,
                             Double_t CenBinWidth,
                             Double_t ptMin=0.2,
                             Double_t ptMax=5.0,
                             Double_t etaMin=-0.8,
                             Double_t etaMax=0.8,
                             TString analysisTypeUser="AOD",
                             Int_t AODfilterBit=768,
                             TString sDataSet="2010",
                             TString EvTrigger="MB",
                             Bool_t bCalculateCME=kFALSE,
                             Bool_t bCalculateCRCPt=kFALSE,
                             Bool_t bUseCRCRecentering=kFALSE,
                             TString QVecWeightsFileName,
                             Bool_t bUsePhiEtaWeights,
                             TString PhiEtaWeightsFileName,
                             Bool_t bUseVZERO=kFALSE,
                             Bool_t bUseVZEROCalib=kFALSE,
                             Bool_t bUseZDC=kFALSE,
                             Bool_t bRecenterZDC=kFALSE,
                             TString sCorrWeight="TPCmVZuZDCu",
                             Bool_t bDivSigma=kTRUE,
                             Bool_t bZDCMCCen=kTRUE,
                             Bool_t bInvertZDC=kFALSE,
                             Bool_t bEventCutsQA=kFALSE,
                             Bool_t bTrackCutsQA=kFALSE,
                             TString Label="",
                             TString sCentrEstimator="V0M",
                             Double_t dVertexRange=10.,
                             Double_t dDCAxy=2.4,
                             Double_t dDCAz=3.2,
                             Double_t dMinClusTPC=70,
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
 
 // create instance of the class: because possible qa plots are added in a second output slot,
 // the flow analysis task must know if you want to save qa plots at the time of class construction
 TString taskFEname = "FlowEventTask";
 taskFEname += CRCsuffix;
 taskFEname += suffix;
 // create instance of the class
 Bool_t bCutsQA = (Bool_t)(bEventCutsQA || bTrackCutsQA);
 AliAnalysisTaskCRCZDC* taskFE = new AliAnalysisTaskCRCZDC(taskFEname, "", bCutsQA);
 taskFE->SetCentralityRange(centrMin,centrMax);
 taskFE->SetCentralityEstimator(sCentrEstimator);
 taskFE->SetUseMCCen(bZDCMCCen);
 taskFE->SetDataSet(sDataSet);
 taskFE->SetQAOn(bCutsQA);
 // set the analysis type
 TString analysisType = "AUTOMATIC";
 if (analysisTypeUser != "") analysisType = analysisTypeUser;
 if (analysisTypeUser == "AOD" || analysisTypeUser == "ESD") analysisType = "AUTOMATIC";
 taskFE->SetAnalysisType(analysisType);
 // add the task to the manager
 mgr->AddTask(taskFE);
 // set the trigger selection
 if (EvTrigger == "Cen")
  taskFE->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
 else if (EvTrigger == "MB")
  taskFE->SelectCollisionCandidates(AliVEvent::kMB);
 else if (EvTrigger == "Any")
  taskFE->SelectCollisionCandidates(AliVEvent::kAny);
 
 // define the event cuts object
 AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("EventCuts");
 // configure some event cuts, starting with centrality
 if(analysisTypeUser == "MCkine") {
  cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
  cutsEvent->SetQA(kFALSE);
 }
 else if (analysisTypeUser == "AOD") {
  cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
  // method used for centrality determination
  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kV0);
  // vertex-z cut
  cutsEvent->SetPrimaryVertexZrange(-dVertexRange,dVertexRange);
  // enable the qa plots
  cutsEvent->SetQA(bEventCutsQA);
  // explicit multiplicity outlier cut
  cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
  if (sDataSet == "2011")
   cutsEvent->SetLHC11h(kTRUE);
  else if (sDataSet == "2010")
   cutsEvent->SetLHC10h(kTRUE);
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
  cutsRP->SetQA(bTrackCutsQA);
  if(bUseVZERO) {
   cutsRP->SetEtaRange(-10.,+10.);
   cutsRP->SetEtaGap(-1.,1.);
  }
  // Track cuts for POIs
  cutsPOI->SetParamType(AliFlowTrackCuts::kMC);
  cutsPOI->SetCutMC(kTRUE);
  cutsPOI->SetPtRange(ptMin,ptMax);
  cutsPOI->SetEtaRange(etaMin,etaMax);
  cutsPOI->SetQA(bTrackCutsQA);
 }
 if (analysisTypeUser == "AOD") {
  // Track cuts for RPs
  if(bUseVZERO) {
   if (sDataSet == "2011")
    cutsRP->SetParamType(AliFlowTrackCuts::kDeltaVZERO);
   else if (sDataSet == "2010")
    cutsRP->SetParamType(AliFlowTrackCuts::kBetaVZERO);
   cutsRP->SetEtaRange(-10.,+10.);
   cutsRP->SetEtaGap(-1.,1.);
   cutsRP->SetPhiMin(0.);
   cutsRP->SetPhiMax(TMath::TwoPi());
   // options for the reweighting
   cutsRP->SetVZEROgainEqualizationPerRing(bUseVZEROCalib);
   cutsRP->SetApplyRecentering(bUseVZEROCalib);
   cutsRP->SetDivSigma(bDivSigma);
  } else {
   cutsRP->SetParamType(AliFlowTrackCuts::kAODFilterBit);
   cutsRP->SetAODfilterBit(AODfilterBit);
   cutsRP->SetMinimalTPCdedx(-999999999);
   cutsRP->SetPtRange(ptMin,ptMax);
   cutsRP->SetEtaRange(etaMin,etaMax);
   cutsRP->SetQA(bTrackCutsQA);
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
  cutsPOI->SetQA(bTrackCutsQA);
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
 if(bCutsQA) {
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
 }

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
 taskQC->SetCalculateCRCPt(bCalculateCRCPt);
 taskQC->SetCalculateCME(bCalculateCME);
 taskQC->SetUseVZERO(bUseVZERO);
 taskQC->SetUseZDC(bUseZDC);
 taskQC->SetRecenterZDC(bRecenterZDC);
 taskQC->SetNUAforCRC(kTRUE);
 taskQC->SetCRCEtaRange(-0.8,0.8);
 taskQC->SetUseCRCRecenter(bUseCRCRecentering);
 taskQC->SetDivSigma(bDivSigma);
 taskQC->SetInvertZDC(bInvertZDC);
 taskQC->SetCorrWeight(sCorrWeight);
 if(bUseCRCRecentering || bRecenterZDC) {
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
 } // end of if(bUseCRCRecentering)
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

