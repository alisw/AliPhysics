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
  if(!task) return 0x0;

  // task->SetPileUpRead(kTRUE);
  task->SetMCRead(kTRUE);
  task->SetGeneratorName("Hijing");
  task->SetEtaMinLimit(-0.8);
  task->SetEtaMaxLimit(0.8);
  task->SetNEtaBins(16);
  task->SetFilterBit(768);
  task->SetChi2DoF(4);
  task->SetEfficiencyTree(kTRUE);//generates ttree for efficiency
  task->SetBuildBackground(kTRUE);

  // TFile *f1 = new TFile("AnalysisResultsbg.root");
  // TList* histlist = (TList*)f1->Get("LegCoef");
  // if(!histlist) printf("error!!!!!!!! no list\n");
  // if(!(TH2D*)histlist->FindObject("PosBGHistOut")) printf("error!!!!!!!! no hist\n");

  // task->GetPosBackground((TH2D*)histlist->FindObject("PosBGHistOut"));
  // task->GetNegBackground((TH2D*)histlist->FindObject("NegBGHistOut"));
  // task->GetChargedBackground((TH2D*)histlist->FindObject("ChargedBGHistOut"));
  // task->GetNeventsCentHist((TH1D*)histlist->FindObject("NeventsCentHist"));

  // task->SetBuildLegendre(kTRUE);


  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1, 25);
  
  // if you want to run locally, we need to define some input
  TChain* chain = new TChain("aodTree");
  chain->Add("/Users/raquelquishpe/analysis/XeMCAOD/AliAOD1.root");
  chain->Add("/Users/raquelquishpe/analysis/XeMCAOD/AliAOD3.root");
  // chain->Add("/Users/raquelquishpe/analysis/XeMCAOD/AliAOD4.root");
  // chain->Add("/Users/raquelquishpe/analysis/XeMCAOD/AliAOD5.root");
  // chain->Add("/Users/raquelquishpe/analysis/XeMCAOD/AliAOD6.root");


  mgr->StartAnalysis("local", chain);
}
