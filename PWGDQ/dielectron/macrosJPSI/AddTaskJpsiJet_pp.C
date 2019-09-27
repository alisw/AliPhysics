// Trigger Index : 0/ALL, 1/INT7, 2/EG1, 3/EG2, 4/DG1, 5/DG2, 6/MC
enum TriggerIndex{
  kALL, kINT7, kEG1, kEG2, kDG1, kDG2, kMC, kNTrigIndex
};
const char* TRIGGER_CLASS[kNTrigIndex] = {
  "INT7;EG1;EG2;DG1;DG2",
  "INT7", "EG1", "EG2", "DG1", "DG2", "MB"};
const char* TRIGGER_TAG[kNTrigIndex] = {
  "ALL",
  "MB", "EG1", "EG2", "DG1", "DG2", "MC"};

AliAnalysisTaskJpsiJet* AddTaskJpsiJet_pp(
    int trigIndex = int(kALL),
    Bool_t enableJetFinder = kTRUE,
    TString period = "16k/pass1",
    const char* suffix = ""){
  // Analysis Manager
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  // Analysis Task
  TString tag = suffix;// support subwagon
  if(tag == "") tag = TRIGGER_TAG[trigIndex];
  AliAnalysisTaskJpsiJet *task = new AliAnalysisTaskJpsiJet(Form("JpsiJet_PP13TeV_%s", tag.Data()));
  // Trigger
  switch(trigIndex){
    case kALL:
      task->SetTrigger(AliVEvent::kINT7 | AliVEvent::kEMCEGA);
      task->SetTriggerQA(kTRUE); // Multi-triggers in single task
      break;
    case kINT7: 
      task->SetTrigger(AliVEvent::kINT7);break;
    case kMC:
      task->SetTrigger(AliVEvent::kAny);
      task->SetMC(kTRUE);
      break;
    default:
      task->SetTrigger(AliVEvent::kEMCEGA);break;
  }
  TString trigClass = TRIGGER_CLASS[trigIndex];
  task->SetTriggerClasses(trigClass.Data());
  // Jet finder
  task->EnableJetFinder(enableJetFinder);
  //Event filter
  AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts", "Vertex Track && |vtxZ|<10 && ncontrib>1");
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxAny);
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10., 10.); // Unit: cm
  task->SetEventFilter(eventCuts);

  task->SetRejectPileup(kTRUE);

  // ITS Improver - by F. Fionda
  if(trigIndex == kMC && !period.Contains("16i")){
    // Period check - 16k_pass1 -> LHC16k/pass1
    period = "LHC" + period.ReplaceAll("_","/");
    //add improver task before the tree maker 
    TGrid::Connect("alien://"); // if not connected input files are not loaded by the improver task
    gSystem->Load("libPWGHFvertexingHF.so"); // load the needed library
    AliAnalysisTaskSEImproveITS *itsImpr =
      reinterpret_cast<AliAnalysisTaskSEImproveITS*>(
      gInterpreter->ExecuteMacro(
      Form("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskImproveITS.C(kFALSE, \"%s\",\"central\")", period.Data())));
    itsImpr->SetMimicData(kTRUE);
    itsImpr->SetAOD(kTRUE);
  }
  
  if(task) mgr->AddTask(task);

  // Output container
  TString containerName = mgr->GetCommonFileName();
	containerName += ":JpsiJetAnalysis";

  AliAnalysisDataContainer* cHistos = mgr->CreateContainer(Form("QAhistos_%s", tag.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data());

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1, cHistos);

  if(trigIndex == kMC){
    cHistos = mgr->CreateContainer("MChistos", TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data());
    mgr->ConnectOutput(task, 2, cHistos);
  }

  return task;
}
