void runLocal() {
  // header location
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");

  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisBG");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  // compile the class (locally) with debug symbols
  gInterpreter->LoadMacro("AliAnalysisTaskLegendreCoef.cxx++g");

  // load the addtask macro and create the task
  AliAnalysisTaskLegendreCoef *task = reinterpret_cast<AliAnalysisTaskLegendreCoef*>(gInterpreter->ExecuteMacro("macros/AddTaskLegendreCoef.C"));
  task->SetPileUpRead(kTRUE);
  // task->SetBuildBackground(kTRUE);
  TFile *f1 = new TFile("/home/alidock/localtest/AnalysisResults9599.root");
  TList* histlist = (TList*)f1->Get("LongFluctuations/EtaBG");
  if(!histlist) printf("error!!!!!!!! no list\n");
  if(!(TH2D*)histlist->FindObject("PosBGHistOut")) printf("error!!!!!!!! no hist\n");

  task->GetPosBackground((TH2D*)histlist->FindObject("PosBGHistOut"));
  task->GetNegBackground((TH2D*)histlist->FindObject("NegBGHistOut"));
  task->GetChargedBackground((TH2D*)histlist->FindObject("ChargedBGHistOut"));

  task->SetBuildLegendre(kTRUE);


  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1, 25);
  
  // if you want to run locally, we need to define some input
  TChain* chain = new TChain("aodTree");
  chain->Add("/home/alidock/longfluc/AliAOD.root");
  mgr->StartAnalysis("local", chain);
}
