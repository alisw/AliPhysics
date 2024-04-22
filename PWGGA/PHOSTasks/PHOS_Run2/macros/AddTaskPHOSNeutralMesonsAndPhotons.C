AliAnalysisPHOSNeutralMesonsAndPhotons* AddTaskPHOSNeutralMesonsAndPhotons(
  const char* name = "PHOSNeutralMesonsAndPhotons",
  const TString CollisionSystem = "PbPb",
  const Bool_t isMC = kFALSE,
  const UInt_t trigger = AliVEvent::kINT7,
  const AliCaloTriggerMimicHelper::phosTriggerType PHOSTrigType = AliCaloTriggerMimicHelper::kPHOSAny,
  const Double_t DispCut = 2.5,
  const Double_t CPVCut = 2.5,
  const Double_t TOFCut = 30.,
  const Double_t Emin = 0.2,
  const Bool_t useCorrE = kTRUE,
  const Bool_t useMinBiasLowEClustCut = kTRUE,
  const Bool_t doNonlinCorr = kTRUE,
  const TString period = "LHC18q",
  const Int_t NMix = 10,
  const TString QAOptionSet = "ClustersQA_CellsQA_TracksQA",
  const Bool_t doCentralityStudy = kTRUE,
  const TString CentralityEstimator = "V0M",
  const Float_t CentralityMin = 0.,
  const Float_t CentralityMax = 10.,
  const char* subname = "")
{

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSNeutralMesonsAndPhotons", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSNeutralMesonsAndPhotons", "This task requires an input event handler");
    return NULL;
  }

  TString taskname = TString(name);
  taskname += "_" + CollisionSystem;

  if (isMC) {
    taskname += "_MC";
  }
  TString TrigName = "";
  Bool_t MBTrigFlag = kFALSE;
  Bool_t PHOSTrigFlag = kFALSE;
  if (trigger == (UInt_t)AliVEvent::kAny)
    TrigName = "kAny";
  else if (trigger == (UInt_t)AliVEvent::kINT7)
    TrigName = "kINT7";
  else if (trigger == (UInt_t)AliVEvent::kCentral)
    TrigName = "kCentral";
  else if (trigger == (UInt_t)AliVEvent::kSemiCentral)
    TrigName = "kSemiCentral";
  else if (trigger == (UInt_t)AliVEvent::kMB)
    TrigName = "kMB";
  else if (trigger == (UInt_t)AliVEvent::kPHI7) {
    switch (PHOSTrigType) {
      case AliCaloTriggerMimicHelper::kPHOSAny:
        TrigName = "kPHI7Any";
        break;
      case AliCaloTriggerMimicHelper::kPHOSL0:
        TrigName = "kPHI7L0";
        break;
      case AliCaloTriggerMimicHelper::kPHOSL1low:
        TrigName = "kPHI7L1L";
        break;
      case AliCaloTriggerMimicHelper::kPHOSL1med:
        TrigName = "kPHI7L1M";
        break;
      case AliCaloTriggerMimicHelper::kPHOSL1high:
        TrigName = "kPHI7L1H";
        break;
      default:
        break;
    }
  } else
    TrigName = "UnknownTrigger";

  taskname += "_" + TrigName;

  UInt_t MBTriggerMask = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB;
  if (trigger & MBTriggerMask)
    MBTrigFlag = kTRUE;
  if (trigger & AliVEvent::kPHI7)
    PHOSTrigFlag = kTRUE;

  Int_t systemID = -1;
  if (CollisionSystem == "pp")
    systemID = 0;
  else if (CollisionSystem == "PbPb")
    systemID = 1;
  else if (CollisionSystem == "pPb" || CollisionSystem == "Pbp")
    systemID = 2;

  Bool_t doClusterQA = kFALSE;
  Bool_t doCellQA = kFALSE;

  if (QAOptionSet.Contains("ClustersQA")) {
    doClusterQA = kTRUE;
    taskname += "_ClustersQA";
  }
  if (QAOptionSet.Contains("CellsQA")) {
    doCellQA = kTRUE;
    taskname += "_CellsQA";
  }

  if (doCentralityStudy) {
    taskname += "_" + CentralityEstimator;
    taskname += Form("_%d_%d", (Int_t)CentralityMin, (Int_t)CentralityMax);
  }

  AliAnalysisPHOSNeutralMesonsAndPhotons* task = new AliAnalysisPHOSNeutralMesonsAndPhotons(taskname.Data());

  task->SelectCollisionCandidates(trigger);
  task->SetMBTrigger(MBTrigFlag);
  task->SetPHOSTrigger(PHOSTrigFlag, PHOSTrigType);
  task->SetMC(isMC);
  task->SetClusterCuts(useCorrE, Emin, DispCut, CPVCut, TOFCut, useMinBiasLowEClustCut);
  task->SetQAStudy(doClusterQA, doCellQA);
  task->SetCentralityStudy(doCentralityStudy, CentralityEstimator, CentralityMin, CentralityMax);
  task->SetNEventsForMixing(NMix);

  if (doNonlinCorr) {
    if (isMC){
      TF1* f1NonlinFunc = new TF1("f1NonlinFunc", "[2]*(1.+[0]/(1.+TMath::Power(x/[1],2)))", 0, 200);

      if (period.Contains("LHC18q") || period.Contains("LHC18r")){
        f1NonlinFunc->FixParameter(0, -0.06); // for core E at ZS 20 MeV with only MIP cut
        f1NonlinFunc->FixParameter(1, 0.7);   // for core E at ZS 20 MeV with only MIP cut
        f1NonlinFunc->FixParameter(2, 1.013); // for core E at ZS 20 MeV with only MIP cut
      }
      else{
        f1NonlinFunc->FixParameter(0, 0.);
        f1NonlinFunc->FixParameter(1, 1.);
        f1NonlinFunc->FixParameter(2, 1.);
      }

      task->SetUserNonlinearityFunction(f1NonlinFunc);
    }
    else{
      if (period.Contains("LHC18q") || period.Contains("LHC18r")){
        Double_t PHOSModWeights[4] = {1.025, 1., 1., 1.}; //Pi0 rec. mass in M1 is low
        task->SetPHOSModEnergyCorr(PHOSModWeights);
      }
    }
  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  TString outputFile = AliAnalysisManager::GetCommonFileName();
  TString outputList = taskname;
  if (TString(subname).Contains("_"))
    outputList += Form("%s", subname);
  else
    outputList += Form("_%s", subname);

  AliAnalysisDataContainer* coutput = mgr->CreateContainer(outputList.Data(),
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputFile.Data());
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}