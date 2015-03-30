void AddTaskCRC(Double_t centrMin,
                Double_t centrMax,
                Double_t ptMin=0.2,
                Double_t ptMax=5.0,
                Double_t etaMin=-0.8,
                Double_t etaMax=0.8,
                TString analysisTypeUser="AOD",
                Int_t AODfilterBit=768,
                TString TPCMultOut="2010",
                TString EvTrigger="MB",
                Bool_t bCalculateCRCPt,
                Bool_t bCalculateCRCBck,
                Bool_t bEventCutsQA=kTRUE,
                Bool_t bTrackCutsQA=kTRUE,
                const char* suffix="") {
 // load libraries
 gSystem->Load("libGeom");
 gSystem->Load("libVMC");
 gSystem->Load("libXMLIO");
 gSystem->Load("libPhysics");
 gSystem->Load("libCore.so");
 gSystem->Load("libTree.so");
 gSystem->Load("libSTEERBase");
 gSystem->Load("libESD");
 gSystem->Load("libAOD");
 gSystem->Load("libANALYSIS");
 gSystem->Load("libANALYSISalice");
 gSystem->Load("libOADB.so");
 gSystem->Load("libPWGflowBase.so");
 gSystem->Load("libPWGflowTasks.so");
 
 gROOT->ProcessLine(".include $ALICE_ROOT/include");
 gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
 gSystem->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/OCDB -I$ALICE_ROOT/STEER/macros -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/OADB $ALICE_PHYSICS/OADB/macros -I$ALICE_PHYSICS/PWGGA -I$ALICE_PHYSICS/PWGCF -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/TENDER -I$ALICE_PHYSICS/TENDER/Tender -I$ALICE_PHYSICS/TENDER/TenderSupplies -I$ALICE_PHYSICS/PARfiles -I$ALICE_PHYSICS/PWGCF/FLOW/macros -g ");
 
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
 
 // set the analysis type automatically
 TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
 if(dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())) analysisType = "MC";
 // or manually, if specified
 if(analysisTypeUser != ""){
  analysisType = analysisTypeUser;
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
 
 // create instance of the class: because possible qa plots are added in a second output slot,
 // the flow analysis task must know if you want to save qa plots at the time of class construction
 TString taskFEname = "FlowEventTask";
 taskFEname += CRCsuffix;
 taskFEname += suffix;
 // create instance of the class
 Bool_t bCutsQA = (Bool_t)(bEventCutsQA || bTrackCutsQA);
 AliAnalysisTaskFlowEvent* taskFE = new AliAnalysisTaskFlowEvent(taskFEname, "", bCutsQA);
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
 if(analysisType == "MC") {
  cutsEvent->SetImpactParameterRange(centrMin,centrMax);
  cutsEvent->SetQA(kFALSE);
 }
 else if (analysisType == "AOD") {
  cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
  // method used for centrality determination
  cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
  cutsEvent->SetRefMultMethod(AliFlowEventCuts::kV0);
  // vertex-z cut
  cutsEvent->SetPrimaryVertexZrange(-10.,10.);
  // enable the qa plots
  cutsEvent->SetQA(bEventCutsQA);
  // explicit multiplicity outlier cut
  cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
  if (TPCMultOut == "2011")
   cutsEvent->SetLHC11h(kTRUE);
  else if (TPCMultOut == "2010")
   cutsEvent->SetLHC10h(kTRUE);
 }
 
 // pass these cuts to your flow event task
 taskFE->SetCutsEvent(cutsEvent);
 AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts("RP cuts");
 AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("POI cuts");
 
 if (analysisType == "MC") {
  // Track cuts for RPs
  cutsRP->SetParamType(AliFlowTrackCuts::kMC);
  cutsRP->SetCutMC(kTRUE);
  cutsRP->SetPtRange(ptMin,ptMax);
  cutsRP->SetEtaRange(etaMin,etaMax);
  cutsRP->SetQA(bTrackCutsQA);
  // Track cuts for POIs
  cutsPOI->SetParamType(AliFlowTrackCuts::kMC);
  cutsPOI->SetCutMC(kTRUE);
  cutsPOI->SetPtRange(ptMin,ptMax);
  cutsPOI->SetEtaRange(etaMin,etaMax);
  cutsPOI->SetQA(bTrackCutsQA);
 }
 else if (analysisType == "AOD") {
  // Track cuts for RPs
  cutsRP->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  cutsRP->SetAODfilterBit(768);
  cutsRP->SetMinimalTPCdedx(-999999999);
  cutsRP->SetPtRange(ptMin,ptMax);
  cutsRP->SetEtaRange(etaMin,etaMax);
  cutsRP->SetQA(bTrackCutsQA);
  // Track cuts for POIs
  cutsPOI->SetParamType(AliFlowTrackCuts::kAODFilterBit);
  cutsPOI->SetAODfilterBit(768);
  cutsPOI->SetMinimalTPCdedx(-999999999);
  cutsPOI->SetPtRange(ptMin,ptMax);
  cutsPOI->SetEtaRange(etaMin,etaMax);
  cutsPOI->SetQA(bTrackCutsQA);
 }
 
 taskFE->SetCutsRP(cutsRP);
 taskFE->SetCutsPOI(cutsPOI);
 
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
 
 // create an additional container for the QA output of the flow event task
 // the QA histograms will be stored in a sub-folder of the output file called 'QA'
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
 
 Bool_t UseParticleWeights = kFALSE;
 //TString ParticleWeightsFileName = "ParticleWeights2D_FullLHC10h_2030.root";
 
 // create the flow analysis tasks
 TString taskCRCname = "AnalysisTask";
 taskCRCname += CRCsuffix;
 taskCRCname += suffix;
 AliAnalysisTaskCRC *taskQC = new AliAnalysisTaskCRC(taskCRCname, UseParticleWeights);
 // set thei triggers
 if (EvTrigger == "Cen")
  taskQC->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
 else if (EvTrigger == "MB")
  taskQC->SelectCollisionCandidates(AliVEvent::kMB);
 else if (EvTrigger == "Any")
  taskQC->SelectCollisionCandidates(AliVEvent::kAny);
 // and set the correct harmonic n
 taskQC->SetHarmonic(1);
 // set standard flow settings
 taskQC->SetCalculateDiffFlow(kTRUE);
 taskQC->SetCalculateDiffFlowVsEta(kTRUE);
 taskQC->SetStoreDistributions(kTRUE);
 taskQC->SetFillMultipleControlHistograms(kFALSE);
 //  CRC settings
 taskQC->SetStoreVarious(kTRUE);
 taskQC->SetCalculateCRC(kTRUE);
 taskQC->SetCalculateCRCPt(bCalculateCRCPt);
 taskQC->SetCalculateCRCBck(bCalculateCRCBck);
 taskQC->SetNUAforCRC(kTRUE);
 taskQC->SetCRCEtaRange(-0.8,0.8);
 
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
 
}
