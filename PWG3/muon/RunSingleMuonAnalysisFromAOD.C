//--------------------------------------------------------------------------
// Base macro for submitting single muon analysis.
// 
// In case it is not run with full aliroot, it needs to have in the working directory:
//  - STEERBase.par
//  - AOD.par
//  - ANALYSIS.par
// 
// The macro reads AODs and outputs file:
// - outputDir/singleMuAnalysis.root
//--------------------------------------------------------------------------

void runSingleMuAnalysis(Char_t *inputDir=".", Char_t *outputDir=".") {
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");  
  gSystem->Load("libPWG3base.so");

  TString outFileName("singleMuAnalysis.root");
  outFileName.Prepend(Form("%s/",outputDir));

  //____________________________________________//
  AliTagAnalysis *TagAna = new AliTagAnalysis("AOD"); 

  AliRunTagCuts *runCuts = new AliRunTagCuts();
  AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts *evCuts = new AliEventTagCuts();

  TagAna->ChainLocalTags(inputDir);
  

  // Temporary workaround to avoid problems with AOD tags.
  TChain* chain = new TChain("aodTree");
  TString inFileName("AliAOD.root");
  inFileName.Prepend(Form("%s/",inputDir));
  chain->Add(inFileName);
  // When problems will be solved and/or you manage in getting a
  // Run*.Merged.AOD.tag.root, substitute with the following lines:

  //TChain* chain = 0x0;
  //chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);


  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  //____________________________________________//
  // 1st Pt task
  AliAnalysisTaskSingleMu *task1 = new AliAnalysisTaskSingleMu("SingleMu");
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cobjArray1", TObjArray::Class(),AliAnalysisManager::kOutputContainer,outFileName.Data());
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}
