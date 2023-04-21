AliAnalysisTaskHFEECorrelators* AddTaskHFEECorrelators(const char* ntracksData, const char* ntracksDet, const char* ntracksTrue, const Double_t R, const Double_t CandidateMinPt, const Double_t CandidateMaxRap, AliAnalysisTaskHFEECorrelators::ECandidateType_t ECandidateType = AliAnalysisTaskHFEECorrelators::kD0toKpi, AliAnalysisTaskHFEECorrelators::JetShapeType jetShapeType = AliAnalysisTaskHFEECorrelators::kData, Bool_t IncludeInclusive = kFALSE) {
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskHFEECorrelators", "No analysis manager found.");
    return 0;
  }
  Bool_t ismc = kFALSE;
  ismc = (mgr->GetMCtruthEventHandler()) ? kTRUE : kFALSE;

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AliAnalysisTaskHFEECorrelators", "This task requires an input event handler");
    return NULL;
  }
  TString wagonName1, wagonName2, wagonName3;
  TString tag = "";
  if (ECandidateType == AliAnalysisTaskHFEECorrelators::kD0toKpi) tag = "kD0toKpi";
  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kData || jetShapeType == AliAnalysisTaskHFEECorrelators::kDataInclusive) {
    wagonName1 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%s", ntracksData, tag.Data());
    wagonName2 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%sTree", ntracksData, tag.Data());
    wagonName3 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%sTreeSplittings", ntracksData, tag.Data());
  }
  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kDetSignal || jetShapeType == AliAnalysisTaskHFEECorrelators::kDetBackground || jetShapeType == AliAnalysisTaskHFEECorrelators::kDetReflection || jetShapeType == AliAnalysisTaskHFEECorrelators::kTrueDet || jetShapeType == AliAnalysisTaskHFEECorrelators::kDet) {
    wagonName1 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%s", ntracksDet, tag.Data());
    wagonName2 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%sTree", ntracksDet, tag.Data());
    wagonName3 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%sTreeSplittings", ntracksDet, tag.Data());
  }
  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kTrue) {
    wagonName1 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%s", ntracksTrue, tag.Data());
    wagonName2 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%sTree", ntracksTrue, tag.Data());
    wagonName3 = Form("AliAnalysisTaskHFEECorrelators_%s_TC%sTreeSplittings", ntracksTrue, tag.Data());
  }
  // Configure jet tagger task
  AliAnalysisTaskHFEECorrelators* task = new AliAnalysisTaskHFEECorrelators(wagonName1);

  task->SetECandidateType_t(ECandidateType);
  task->SetJetShapeType(jetShapeType);
  task->SetJetRadius(R);
  task->SetIncludeInclusive(IncludeInclusive);
  task->SetCandidateMinPt(CandidateMinPt);
  task->SetCandidateMaxY(CandidateMaxRap);

  // AliParticleContainer *trackContData=0x0;  //why not track containers?
  // AliParticleContainer *trackContDet=0x0;
  // AliParticleContainer *trackContTrue=0x0;

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kData || jetShapeType == AliAnalysisTaskHFEECorrelators::kDataInclusive) {
    AliHFTrackContainer* trackContData = new AliHFTrackContainer(ntracksData);
    task->AdoptParticleContainer(trackContData);
  }
  else if (jetShapeType == AliAnalysisTaskHFEECorrelators::kDetSignal || jetShapeType == AliAnalysisTaskHFEECorrelators::kDetBackground || jetShapeType == AliAnalysisTaskHFEECorrelators::kDetReflection || jetShapeType == AliAnalysisTaskHFEECorrelators::kTrueDet || jetShapeType == AliAnalysisTaskHFEECorrelators::kDet) {
    AliHFTrackContainer* trackContDet = new AliHFTrackContainer(ntracksDet);
    task->AdoptParticleContainer(trackContDet);
    AliMCParticleContainer* trackContTrue = new AliHFAODMCParticleContainer(ntracksTrue);
    trackContTrue->SetEtaLimits(-1.5, 1.5);
    trackContTrue->SetPtLimits(0, 1000);
    task->AdoptParticleContainer(trackContTrue);
  }
  else if (jetShapeType == AliAnalysisTaskHFEECorrelators::kTrue) {
    AliMCParticleContainer* trackContTrue = new AliHFAODMCParticleContainer(ntracksTrue);
    trackContTrue->SetEtaLimits(-1.5, 1.5);
    trackContTrue->SetPtLimits(0, 1000);
    task->AdoptParticleContainer(trackContTrue);
  }

  task->SetUseAliAnaUtils(kFALSE);

  mgr->AddTask(task);

  // Connnect input
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  // Connect output
  TString contName1(wagonName1);
  TString contName2(wagonName2);
  TString contName3(wagonName3);

  contName2 += "_Splittings";

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kData) {
    contName1 += "_Data";
    contName2 += "_Data";
    contName3 += "_Data";
  }

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kDetSignal) {
    contName1 += "_DetSignal";
    contName2 += "_DetSignal";
    contName3 += "_DetSignal";
  }

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kDetBackground) {
    contName1 += "_DetBackground";
    contName2 += "_DetBackground";
    contName3 += "_DetBackground";
  }

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kDetReflection) {
    contName1 += "_DetReflection";
    contName2 += "_DetReflection";
    contName3 += "_DetReflection";
  }

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kTrueDet) {
    contName1 += "_TrueDet";
    contName2 += "_TrueDet";
    contName3 += "_TrueDet";
  }

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kTrue) {
    contName1 += "_True";
    contName2 += "_True";
    contName3 += "_True";
  }
  if (IncludeInclusive) {
    contName1 += "_Inclusive";
    contName2 += "_Inclusive";
    contName3 += "_Inclusive";
  }

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kDataInclusive) {
    contName1 += "_DataInclusive";
    contName2 += "_DataInclusive";
    contName3 += "_DataInclusive";
  }

  if (jetShapeType == AliAnalysisTaskHFEECorrelators::kDet) {
    contName1 += "_Det";
    contName2 += "_Det";
    contName3 += "_Det";
  }

  TString outputfile = Form("%s", AliAnalysisManager::GetCommonFileName());
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(contName1.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 1, coutput1);
  AliAnalysisDataContainer* coutput2 = mgr->CreateContainer(contName2.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 2, coutput2);
  AliAnalysisDataContainer* coutput3 = mgr->CreateContainer(contName3.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer, outputfile);
  mgr->ConnectOutput(task, 3, coutput3);

  return task;
}