void AddTaskCRC(Int_t iHarmonic, Double_t centrMin, Double_t centrMax, Double_t pTMin, Double_t pTMax, Bool_t UseEtaGap = kTRUE, Double_t EtaGap = 0., Bool_t doQA = kTRUE, Bool_t bUse11h = kFALSE, const char* suffix = "")
{
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
 
 // specify charge/CEA in all names
 TString CRCsuffix = ":CRC";
 
 TString HarmName = "_n";
 HarmName += (Int_t)iHarmonic;
 CRCsuffix += HarmName;
 
 TString CentrName = "_";
 CentrName += (Int_t)centrMin;
 CentrName += "-";
 CentrName += (Int_t)centrMax;
 CRCsuffix += CentrName;
 
 TString pTName = "_";
 Int_t rt = (Int_t)(pTMin*10.);
 Int_t r = (Int_t)(pTMin);
 pTName += ( pTMin < 1. ? Form("0.%i",rt) : Form("%i.%i",r,rt-r*10));
 pTName += "-";
 rt = (Int_t)(pTMax*10.);
 r = (Int_t)(pTMax);
 pTName += ( pTMax < 1. ? Form("0.%i",rt) : Form("%i.%i",r,rt-r*10));
 CRCsuffix += pTName;
 
 TString EtaGapName = "_EG";
 if (!UseEtaGap) { EtaGapName += "NULL"; }
 else            { EtaGapName += Form("0.%i",(Int_t)(EtaGap*10.)); }
 CRCsuffix += EtaGapName;
 
 // specify some variables
 // Float_t centrMin = 20.;
 // Float_t centrMax = 30.;
 
 // create instance of the class: because possible qa plots are added in a second output slot,
 // the flow analysis task must know if you want to save qa plots at the time of class construction
 TString taskFEname = "FlowEventTask";
 taskFEname += CRCsuffix;
 taskFEname += suffix;
 // create instance of the class
 AliAnalysisTaskFlowEvent* taskFE = new AliAnalysisTaskFlowEvent(taskFEname, "", doQA);
 // add the task to the manager
 mgr->AddTask(taskFE);
 // set the trigger selection
 if (bUse11h) { taskFE->SelectCollisionCandidates(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral); }
 else         { taskFE->SelectCollisionCandidates(AliVEvent::kMB); }
 
 // define the event cuts object
 AliFlowEventCuts* cutsEvent = new AliFlowEventCuts("EventCuts");
 // configure some event cuts, starting with centrality
 cutsEvent->SetCentralityPercentileRange(centrMin,centrMax);
 // method used for centrality determination
 cutsEvent->SetCentralityPercentileMethod(AliFlowEventCuts::kV0);
 cutsEvent->SetRefMultMethod(AliFlowEventCuts::kV0);
 // vertex-z cut
 cutsEvent->SetPrimaryVertexZrange(-10.,10.);
 // enable the qa plots
 cutsEvent->SetQA(doQA);
 // explicit multiplicity outlier cut
 cutsEvent->SetCutTPCmultiplicityOutliersAOD(kTRUE);
 if (bUse11h) { cutsEvent->SetLHC11h(kTRUE); }
 else         { cutsEvent->SetLHC10h(kTRUE); }
 
 // pass these cuts to your flow event task
 taskFE->SetCutsEvent(cutsEvent);
 
 // create the track cuts for RPs
 //AliFlowTrackCuts* cutsRP = AliFlowTrackCuts::GetAODTrackCutsForFilterBit(768, "RP cuts");
 AliFlowTrackCuts* cutsRP = new AliFlowTrackCuts("RP cuts");
 cutsRP->SetMinimalTPCdedx(-999999999);
 cutsRP->SetAODfilterBit(768);
 cutsRP->SetParamType(AliFlowTrackCuts::kAODFilterBit);
 cutsRP->SetPtRange(pTMin,pTMax);
 cutsRP->SetEtaRange(-0.8,0.8);
 //   cutsRP->SetMinNClustersTPC(70);
 //   cutsRP->SetMinChi2PerClusterTPC(0.2);
 //   cutsRP->SetMaxChi2PerClusterTPC(4.0);
 //   cutsRP->SetRequireTPCRefit(kFALSE);
 //   cutsRP->SetMaxDCAToVertexXY(2.4);
 //   cutsRP->SetMaxDCAToVertexZ(3.2);
 //   cutsRP->SetDCAToVertex2D(kTRUE);
 //   cutsRP->SetAcceptKinkDaughters(kFALSE);
 //   cutsRP->SetMinimalTPCdedx(10.);
 cutsRP->SetQA(doQA);
 
 taskFE->SetCutsRP(cutsRP);
 
 // create the track cuts for POIs
 AliFlowTrackCuts* cutsPOI = new AliFlowTrackCuts("POI cuts");
 cutsPOI->SetMinimalTPCdedx(-999999999);
 cutsPOI->SetAODfilterBit(768);
 cutsPOI->SetParamType(AliFlowTrackCuts::kAODFilterBit);
 cutsPOI->SetPtRange(pTMin,pTMax);
 cutsPOI->SetEtaRange(-0.8,0.8);
 //   cutsPOI->SetMinNClustersTPC(70);
 //   cutsPOI->SetMinChi2PerClusterTPC(0.2);
 //   cutsPOI->SetMaxChi2PerClusterTPC(4.0);
 //   cutsPOI->SetRequireTPCRefit(kTRUE);
 //   cutsPOI->SetMaxDCAToVertexXY(2.4);
 //   cutsPOI->SetMaxDCAToVertexZ(3.2);
 //   cutsPOI->SetDCAToVertex2D(kTRUE);
 //   cutsPOI->SetAcceptKinkDaughters(kFALSE);
 //   cutsPOI->SetMinimalTPCdedx(10.);
 cutsPOI->SetQA(doQA);
 
 //cutsPOI->SetRequireCharge(kTRUE);
 //cutsPOI->SetMinNsigmaToVertex(4);
 //cutsPOI->SetRequireSigmaToVertex(kTRUE);
 //cutsPOI->SetMinNClustersITS(2);
 //cutsPOI->SetRequireITSRefit(kTRUE);
 //cutsPOI->SetMaxChi2PerClusterITS(36.);
 
 /*   cutsPOI->SetPID(AliPID::kPion, AliFlowTrackCuts::kTPCTOFNsigma, 0.9);
  cutsPOI->SetPriors((centrMin+centrMax)*0.5); // set priors and PID as a function of the centrality
  //cutsPOI->SetAllowTOFmismatchFlag(kTRUE);
  cutsPOI->SetRequireStrictTOFTPCagreement(kTRUE);
  cutsPOI->GetBayesianResponse()->ForceOldDedx();*/ // for 2010 data to use old TPC PID Response instead of the official one
 //example: francesco's tunig TPC Bethe Bloch for data:
 //cutsPOI->GetESDpid().GetTPCResponse().SetBetheBlochParameters(4.36414e-02,1.75977e+01,1.14385e-08,2.27907e+00,3.36699e+00);
 //cutsPOI->GetESDpid().GetTPCResponse().SetMip(49);
 
 //  cutsPOI->SetCharge(charge);
 
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
 if(doQA) {
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
 taskQC->SelectCollisionCandidates(AliVEvent::kMB);
 // and set the correct harmonic n
 taskQC->SetHarmonic(iHarmonic);
 // set standard flow settings
 taskQC->SetCalculateDiffFlow(kTRUE);
 taskQC->SetCalculateDiffFlowVsEta(kTRUE);
 taskQC->SetStoreDistributions(kTRUE);
 taskQC->SetFillMultipleControlHistograms(kFALSE);
 //  CRC settings
 taskQC->SetStoreVarious(kTRUE);
 taskQC->SetCalculateCRC(kTRUE);
 taskQC->SetCalculateCRCPt(kTRUE);
 taskQC->SetUseEtaGap(UseEtaGap);
 taskQC->SetEtaGap(EtaGap);
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
