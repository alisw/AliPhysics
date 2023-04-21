// ROOT6 modifications
#ifdef __CLING__
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliDielectronVarCuts.h>

#include <AliDielectronCutGroup.h>
#include <AliDielectronEventCuts.h>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#endif



AliAnalysisTask *AddTaskMLTreeMaker2018(TString taskname = "ESDExample",
	       				TString outputname = "AnalysisResults.root",
              				Bool_t isAOD=kTRUE,
              				Bool_t getFromAlien = kFALSE,
              				TString cFileName="Config_jjung_lowmass.C",	
					Int_t cutsetting = 0,
              				Int_t triggerMask = AliVEvent::kINT7,
              				Double_t centMin = 0.,
              				Double_t centMax = 100.,
              				Float_t PtMin =  0.4,
              				Float_t PtMax = 10.0,
              				Float_t EtaMin = -0.8,
              				Float_t EtaMax = +0.8,
					Bool_t DoTracks = kTRUE,
					Bool_t DoTrackQuality = kTRUE,

					Bool_t DoPairing =kTRUE,
              				Bool_t DoULS = kTRUE,
					Bool_t DoLS = kTRUE,	
					TString generatorNameForMCSignal  = "Pythia CC_0;Pythia B_0;Pythia BB_0;",
              				TString generatorNameForULSSignal = "Hijing_1;Pythia CC_0;Pythia B_0;Pythia BB_0;"

					) {				


AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
if (!mgr) {
  ::Error("AddTaskBalancePsiCentralityTrain",  "No analysis manager to connect to.");
  return NULL;
}

//Base Directory for GRID / LEGO Train
TString configBasePath= "/home/jerome/analysis/localLegotrain/005_tests_PbPbAOD/";
if (!gSystem->AccessPathName(cFileName))
{
  printf("file already present\n");
  configBasePath=Form("%s/",gSystem->pwd());
}
else if(getFromAlien && (!gSystem->Exec(Form("alien_cp alien:///alice/cern.ch/user/j/jjung/%s file:./",cFileName.Data()))) ){
  std::cout << "ALIEN???!  " << std::endl;
  configBasePath=Form("%s/",gSystem->pwd());
}

TString configFilePath(configBasePath+cFileName);

std::cout << "Configpath:  " << configFilePath << std::endl;

//Do we have an MC handler?
Bool_t hasMC = (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()!=0x0);

if (!gROOT->GetListOfGlobalFunctions()->FindObject(cFileName.Data()))
  gROOT->LoadMacro(configFilePath.Data());  //old root5


// Check the analysis type using the event handlers connected to the analysis manager.
//===========================================================================
if (!mgr->GetInputEventHandler()) {
  ::Error("AddTaskMLTreeMaker",  "This task requires an input event handler");
  return NULL;
}
TString analysisType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

if (analysisType!="ESD"){
  ::Warning("AddTaskMLTreeMaker",  "analysis type NOT ESD --> some variables might not be filled");
}
      

AliAnalysisTaskMLTreeMaker2018 *task = new AliAnalysisTaskMLTreeMaker2018(taskname);

   // ==========================================================================
  // user customization part

  task->SetEtaRange(EtaMin, EtaMax);
  task->SetPtRange(PtMin, PtMax);
  task->SetCentralityPercentileRange(centMin, centMax);
  task->SelectCollisionCandidates(triggerMask);
  task->isMC(hasMC);
  task->runningOnAOD(isAOD);


  task->SetupPairing(DoPairing);
  task->doULSpairs(DoULS);
  task->doLSpairs(DoLS);

  task->SetGeneratorName(generatorNameForMCSignal);
  task->SetGeneratorMCSignalName(generatorNameForMCSignal);
  task->SetGeneratorULSSignalName(generatorNameForULSSignal);

  task->SetupTracks(DoTracks);
  task->filterTrackQuality(DoTrackQuality);
  
  task->SetupTrackCuts((reinterpret_cast<AliDielectronCutGroup*>(gROOT->ProcessLine(Form("SetupAODtrackCuts(%i)",cutsetting)))));
  task->SetupEventCuts((reinterpret_cast<AliDielectronEventCuts*>(gROOT->ProcessLine(Form("GetEventCuts(%i)",cutsetting)))));


  gROOT->GetListOfSpecials()->Add(task);
  // #########################################################
  // #########################################################
  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  gROOT->ProcessLine(Form("AddSingleLegMCSignal(%s)",task->GetName()));//not task itself, task name
  gROOT->ProcessLine(Form("AddPairMCSignal(%s)"     ,task->GetName()));//not task itself, task name
  //Done in config
  //std::vector<bool> DielectronsPairNotFromSameMother = AddSingleLegMCSignal(task);
  //task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);



  task->SetFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);// for AOD analyses: TPC cuts + any SPD hit
  // ==========================================================================

  mgr->AddTask(task);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  //TString outputFileName = AliAnalysisManager::GetCommonFileName();

  
// AliAnalysisDataContainer *coutESD = mgr->CreateContainer("",  TTree::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
// AliAnalysisDataContainer *coutESD1 = mgr->CreateContainer("QAHist",  TH1::Class(), AliAnalysisManager::kOutputContainer, outputFileName.Data());
// mgr->ConnectInput(taskESD,  0,  mgr->GetCommonInputContainer());
// mgr->ConnectOutput(taskESD,  1,  coutESD);
// mgr->ConnectOutput(taskESD,  2,  coutESD1);
 
  AliAnalysisDataContainer *cout = mgr->CreateContainer("output", TList::Class(),AliAnalysisManager::kOutputContainer,outputname.Data());
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, cout);
 
 return task;
}
