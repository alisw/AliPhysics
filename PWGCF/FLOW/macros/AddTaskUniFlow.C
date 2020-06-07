///////////////////////////////////////////////////////////////////
//
// AddTaskUniFlow.C macro
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
//
//  See AliAnalysisTaskUniFlow(.cxx) for details & documentation
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class AliAnalysisTaskUniFlow;

AliAnalysisTaskUniFlow* AddTaskUniFlow(
    AliAnalysisTaskUniFlow::ColSystem colSys,
    TString sWeigthsFile = "",
    Bool_t bIsMC = kFALSE,
    const char* suffix = ""
)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { return NULL; }
  if (!mgr->GetInputEventHandler())	{ return NULL; }

  // checking weights
  Bool_t bUseWeights = kFALSE;
  if(!sWeigthsFile.IsNull()) { bUseWeights = kTRUE; }

  // crate a combined name for the task (required for LEGO trains)
  TString taskName = Form("UniFlow%s",suffix);

  AliAnalysisTaskUniFlow* task = new AliAnalysisTaskUniFlow(taskName.Data(), colSys, bUseWeights, bIsMC); // now we create an instance of your task
  if(!task) { return NULL; }

  // task default settings (ColSystem independent)
  task->SetAnalysisType(AliAnalysisTaskUniFlow::kAOD);
  task->SetRunMode(AliAnalysisTaskUniFlow::kFull);
  task->SetNumEventsAnalyse(50);
  task->SetTrigger(AliVEvent::kINT7);
  task->SetFlowFillWeights(kTRUE);
  task->SetUseWeigthsRunByRun(kTRUE);
  task->SetUseWeights3D(kFALSE);
  task->SetFillQAhistos(kTRUE);
  task->SetSampling(kFALSE);
  task->SetProcessPID(kTRUE);
  task->SetProcessPhi(kTRUE);
  task->SetProcessV0s(kTRUE);
  task->SetFlowRFPsPt(0.2,5.0);
  task->SetFlowPOIsPt(0.0,10.0);
  task->SetFlowEta(0.8);

  // task default settings dependent on ColSystem (colSys)
  if(colSys == AliAnalysisTaskUniFlow::ColSystem::kPbPb) {
    // Pb-Pb
    task->SetCentrality(AliAnalysisTaskUniFlow::kV0M);
    task->SetPVtxZMax(10.0);
    task->SetRejectAddPileUp(kFALSE);
    task->SetChargedNumTPCclsMin(70);
    task->SetChargedDCAzMax(0.0);
    task->SetChargedDCAxyMax(0.0);
    task->SetChargedTrackFilterBit(96);
    task->SetPIDUseAntiProtonOnly(kFALSE);
    task->SetPIDNumSigmasCombinedNoTOFrejection(kTRUE);
    task->SetPIDNumSigmasPionMax(3.0);
    task->SetPIDNumSigmasKaonMax(3.0);
    task->SetPIDNumSigmasProtonMax(3.0);
    task->SetPIDNumSigmasTPCRejectElectron(0.0);
    task->SetUseBayesPID(kTRUE);
    task->SetPIDBayesProbPionMin(0.95);
    task->SetPIDBayesProbKaonMin(0.85);
    task->SetPIDBayesProbProtonMin(0.85);
    task->SetV0sOnFly(kFALSE);
    task->SetV0sTPCRefit(kTRUE);
    task->SetV0sRejectKinks(kTRUE);
    task->SetV0sDaughterNumTPCClsMin(70);
    task->SetV0sDaughterNumTPCrossMin(70);
    task->SetV0sDaughterNumTPCFindMin(1);
    task->SetV0sDaughterNumTPCClsPIDMin(70);
    task->SetV0sDaughterRatioCrossFindMin(0.8);
    task->SetV0sUseCrossMassRejection(kTRUE);
    task->SetV0sCrossMassCutK0s(0.005);
    task->SetV0sCrossMassCutLambda(0.010);
    task->SetV0sDCAPVMin(0.1);
    task->SetV0sDCAPVMax(0.0);
    task->SetV0sDCAPVzMax(0.0);
    task->SetV0sDaughtersFilterBit(0);
    task->SetV0sDCADaughtersMin(0.0);
    task->SetV0sDCADaughtersMax(0.5);
    task->SetV0sDecayRadiusMin(5.0);
    task->SetV0sDecayRadiusMax(100.0);
    task->SetV0sDaughterEtaMax(0.8);
    task->SetV0sDaughterPtMin(0.0);
    task->SetV0sDaughterPtMax(0.0);
    task->SetV0sMotherRapMax(0.0);
    task->SetV0sK0sInvMassMin(0.4);
    task->SetV0sK0sInvMassMax(0.6);
    task->SetV0sLambdaInvMassMin(1.08);
    task->SetV0sLambdaInvMassMax(1.16);
    task->SetV0sK0sCPAMin(0.998);
    task->SetV0sLambdaCPAMin(0.998);
    task->SetV0sK0sNumTauMax(0.0);
    task->SetV0sLambdaNumTauMax(0.0);
    task->SetV0sK0sArmenterosAlphaMin(0.2);
    task->SetV0sLambdaArmenterosAlphaMax(0.0);
    task->SetV0sK0sPionNumTPCSigmaMax(3.0);
    task->SetV0sLambdaPionNumTPCSigmaMax(3.0);
    task->SetV0sLambdaProtonNumTPCSigmaMax(3.0);
    task->SetPhiInvMassMin(0.99);
    task->SetPhiInvMassMax(1.07);
  } else {
    // p-Pb & pp
    // NB: so far based on "previous" default ctor values
    task->SetCentrality(AliAnalysisTaskUniFlow::kV0A);
    task->SetPVtxZMax(10.0);
    task->SetRejectAddPileUp(kFALSE);
    task->SetChargedNumTPCclsMin(70);
    task->SetChargedDCAzMax(0.0);
    task->SetChargedDCAxyMax(0.0);
    task->SetChargedTrackFilterBit(96);
    task->SetPIDUseAntiProtonOnly(kFALSE);
    task->SetPIDNumSigmasCombinedNoTOFrejection(kTRUE);
    task->SetPIDNumSigmasPionMax(3.0);
    task->SetPIDNumSigmasKaonMax(3.0);
    task->SetPIDNumSigmasProtonMax(3.0);
    task->SetPIDNumSigmasTPCRejectElectron(0.0);
    task->SetUseBayesPID(kTRUE);
    task->SetPIDBayesProbPionMin(0.8);
    task->SetPIDBayesProbKaonMin(0.8);
    task->SetPIDBayesProbProtonMin(0.8);
    task->SetV0sOnFly(kFALSE);
    task->SetV0sTPCRefit(kTRUE);
    task->SetV0sRejectKinks(kTRUE);
    task->SetV0sDaughterNumTPCClsMin(70);
    task->SetV0sDaughterNumTPCrossMin(70);
    task->SetV0sDaughterNumTPCFindMin(1);
    task->SetV0sDaughterNumTPCClsPIDMin(70);
    task->SetV0sDaughterRatioCrossFindMin(0.8);
    task->SetV0sUseCrossMassRejection(kTRUE);
    task->SetV0sCrossMassCutK0s(0.005);
    task->SetV0sCrossMassCutLambda(0.010);
    task->SetV0sDCAPVMin(0.06);
    task->SetV0sDCAPVMax(0.0);
    task->SetV0sDCAPVzMax(0.0);
    task->SetV0sDaughtersFilterBit(0);
    task->SetV0sDCADaughtersMin(0.0);
    task->SetV0sDCADaughtersMax(1.0);
    task->SetV0sDecayRadiusMin(0.5);
    task->SetV0sDecayRadiusMax(200.0);
    task->SetV0sDaughterEtaMax(0.8);
    task->SetV0sDaughterPtMin(0.0);
    task->SetV0sDaughterPtMax(0.0);
    task->SetV0sMotherRapMax(0.0);
    task->SetV0sK0sInvMassMin(0.4);
    task->SetV0sK0sInvMassMax(0.6);
    task->SetV0sLambdaInvMassMin(1.08);
    task->SetV0sLambdaInvMassMax(1.16);
    task->SetV0sK0sCPAMin(0.97);
    task->SetV0sLambdaCPAMin(0.995);
    task->SetV0sK0sNumTauMax(0);
    task->SetV0sLambdaNumTauMax(0);
    task->SetV0sK0sArmenterosAlphaMin(0.2);
    task->SetV0sLambdaArmenterosAlphaMax(0.0);
    task->SetV0sK0sPionNumTPCSigmaMax(5.0);
    task->SetV0sLambdaPionNumTPCSigmaMax(5.0);
    task->SetV0sLambdaProtonNumTPCSigmaMax(5.0);
    task->SetPhiInvMassMin(0.99);
    task->SetPhiInvMassMax(1.07);
  }

  mgr->AddTask(task); // add your task to the manager

  // Create a common outpit file name with subfolder based on task name
  TString fileName = Form("%s:%s",AliAnalysisManager::GetCommonFileName(),taskName.Data());

  // Creating containers
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* cOutput1 = mgr->CreateContainer(TString("Flow_Refs_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput2 = mgr->CreateContainer(TString("Flow_Charged_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput3 = mgr->CreateContainer(TString("Flow_Pion_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput4 = mgr->CreateContainer(TString("Flow_Kaon_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput5 = mgr->CreateContainer(TString("Flow_Proton_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput6 = mgr->CreateContainer(TString("Flow_ChargedUnidentified_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput7 = mgr->CreateContainer(TString("Flow_K0s_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput8 = mgr->CreateContainer(TString("Flow_Lambda_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput9 = mgr->CreateContainer(TString("Flow_Phi_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput10 = mgr->CreateContainer(TString("QA_Events_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput11 = mgr->CreateContainer(TString("QA_Charged_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput12 = mgr->CreateContainer(TString("QA_PID_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput13 = mgr->CreateContainer(TString("QA_V0s_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput14 = mgr->CreateContainer(TString("QA_Phi_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
  AliAnalysisDataContainer* cOutput15 = mgr->CreateContainer(TString("Weights_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);

  // Connecting containers to task
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,cOutput1);
  mgr->ConnectOutput(task,2,cOutput2);
  mgr->ConnectOutput(task,3,cOutput3);
  mgr->ConnectOutput(task,4,cOutput4);
  mgr->ConnectOutput(task,5,cOutput5);
  mgr->ConnectOutput(task,6,cOutput6);
  mgr->ConnectOutput(task,7,cOutput7);
  mgr->ConnectOutput(task,8,cOutput8);
  mgr->ConnectOutput(task,9,cOutput9);
  mgr->ConnectOutput(task,10,cOutput10);
  mgr->ConnectOutput(task,11,cOutput11);
  mgr->ConnectOutput(task,12,cOutput12);
  mgr->ConnectOutput(task,13,cOutput13);
  mgr->ConnectOutput(task,14,cOutput14);
  mgr->ConnectOutput(task,15,cOutput15);

  if(bUseWeights) {
    TObjArray* taskContainers = mgr->GetContainers();
    if(!taskContainers) { printf("E-AddTaskUniFlow: Task containers does not exists!\n"); return NULL; }

    // check if the input weights are already loaded (e.g. in different subwagon)
    AliAnalysisDataContainer* weights = (AliAnalysisDataContainer*) taskContainers->FindObject("inputWeights");
    if(!weights) {
      // if it does not exists create it

      // in case of non-local run, establish connection to ALiEn for loading the weights
      if(sWeigthsFile.Contains("alien://")) { gGrid->Connect("alien://"); }

      TFile* weights_file = TFile::Open(sWeigthsFile.Data(),"READ");
      if(!weights_file) { printf("E-AddTaskUniFlow: Input file with weights not found!\n"); return NULL; }

      TList* weights_list = (TList*) weights_file->Get("weights");
      if(!weights_list) { printf("E-AddTaskUniFlow: Input list with weights not found!\n"); weights_file->ls(); return NULL; }

      AliAnalysisDataContainer* cInputWeights = mgr->CreateContainer("inputWeights",TList::Class(), AliAnalysisManager::kInputContainer);
      cInputWeights->SetData(weights_list);
      mgr->ConnectInput(task,1,cInputWeights);
    }
    else {
      // connect existing container
      mgr->ConnectInput(task,1,weights);
    }
  }

  // Monte Carlo
  if(bIsMC) {
    AliAnalysisDataContainer* cOutMC = mgr->CreateContainer(TString("MC_")+taskName, TList::Class(), AliAnalysisManager::kOutputContainer, fileName);
    mgr->ConnectOutput(task,16,cOutMC);

  }

  return task;
}
