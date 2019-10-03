//TString dataDir = "/alice/data/2013/LHC13b";
//TString dataPattern = "*/ESDs/pass3/AOD154/*/AliAOD.root";
//TString workingDir = "lhc13b";
//Int_t nRuns = 12;
//Int_t runList[] = {
//      195483, 195482, 195481, 195480, 195479, 195478, 195391, 195390, 195389, 195351, 
//      195346, 195344
//};

TString dataDir = "/alice/data/2013/LHC13c";
TString dataPattern = "*/ESDs/pass2/AOD154/*/AliAOD.root";
TString workingDir = "lhc13c";
Int_t nRuns = 13;
Int_t runList[] = {
      195677, 195675, 195673, 195644, 195633, 195596, 195593, 195592, 195568, 195567, 
      195566, 195531, 195529
};

//TString dataDir = "/alice/data/2013/LHC13d";
//TString dataPattern = "*/ESDs/muon_pass2/AOD134/*/AliAOD.root";
//
//TString workingDir = "lhc13d";
//Int_t nRuns = 20;
//Int_t runList[] = {
//      195873, 195872, 195871, 195869, 195867, 195831, 195830, 195829, 195827, 195826, 
//      195787, 195783, 195767, 195765, 195760, 195727, 195726, 195725, 195724, 195682
//};

//TString dataDir = "/alice/data/2013/LHC13e";
//TString dataPattern = "*/ESDs/muon_pass2/AOD134/*/AliAOD.root";
//TString workingDir = "lhc13e";
//Int_t nRuns = 27;
//Int_t runList[] = {
//195949, 195950, 195954, 195955, 195958, 195989, 195994, 196000, 196006, 196085,
//196089, 196090, 196091, 196105, 196107, 196185, 196187, 196194, 196199, 196200,
//196201, 196203, 196214, 196308, 196309, 196310, 196311
//};


//TString dataDir = "/alice/data/2013/LHC13f";
//TString dataPattern = "/ESDs/muon_pass2/AOD/*/AliAOD.root";
//TString workingDir = "lhc13f";
//Int_t nRuns = 64;
//Int_t runList[] = {
//    197388, 197387, 197386, 197349, 197348, 197342, 197341, 197302, 197299, 197298, 
//    197258, 197256, 197255, 197254, 197247, 197189, 197184, 197153, 197152, 197150, 
//    197148, 197147, 197145, 197144, 197143, 197142, 197139, 197138, 197099, 197098, 
//    197092, 197091, 197089, 197011, 197003, 196974, 196973, 196972, 196965, 196876, 
//    196869, 196774, 196773, 196772, 196722, 196721, 196720, 196702, 196701, 196648, 
//    196646, 196608, 196605, 196601, 196568, 196566, 196564, 196563, 196535, 196528, 
//    196477, 196475, 196474, 196433
//};

//TString dataDir = "/alice/data/2013/LHC13d";
//TString dataPattern = "*/pass2/AOD154/*/AliAOD.root";
//
//TString workingDir = "lhc13d_pass2";
//Int_t nRuns = 20;
//Int_t runList[] = {
//      195873, 195872, 195871, 195869, 195867, 195831, 195830, 195829, 195827, 195826, 
//      195787, 195783, 195767, 195765, 195760, 195727, 195726, 195725, 195724, 195682
//};

//TString dataDir = "/alice/data/2013/LHC13e";
//TString dataPattern = "*/pass2/AOD154/*/AliAOD.root";
//TString workingDir = "lhc13e_pass2";
//Int_t nRuns = 27;
//Int_t runList[] = {
//195949, 195950, 195954, 195955, 195958, 195989, 195994, 196000, 196006, 196085,
//196089, 196090, 196091, 196105, 196107, 196185, 196187, 196194, 196199, 196200,
//196201, 196203, 196214, 196308, 196309, 196310, 196311
//};


void runGridAod(){
  if (!TGrid::Connect("alien://")) return;

  AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
  mgr->SetDebugLevel(0);
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  

//  gROOT->LoadMacro("$ALICE_PHYSICS/../src/ANALYSIS/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AddTaskCentrality(1,1);

//    gROOT->LoadMacro("$ALICE_ROOT/PWGCF/Correlations/macros/cftree/AddTaskCFTree.C");
//    AliAnalysisTaskCFTree *taskCFTree = AddTaskCFTree(); 

  gROOT->LoadMacro("AliAnalysisTaskCFTreeNew.cxx+g");
  AliAnalysisTaskCFTreeNew *ana = new AliAnalysisTaskCFTreeNew();
  ana->SetEventSelectionBit(AliVEvent::kINT7 | AliVEvent::kMUS7 | AliVEvent::kMUSH7);
  ana->SetZVertex(7);
  ana->SetDphiCut(0.005);
  ana->SetTrackletEtaCut(1);
  ana->SetTrackEtaCut(1.),
  ana->SetPtMin(0.5);
  ana->SetTrackFilterBit(1<<8|1<<9);
  ana->SetStoreTracks(kTRUE);
  ana->SetStoreTracklets(kTRUE);
  ana->SetStoreMuons(kTRUE);
  ana->SetStoreOnlyEventsWithMuons(kTRUE);
  ana->SetApplyPhysicsSelectionCut(kTRUE);

  mgr->AddTask(ana);
  
  const char* outputFileName = "AnalysisResults.root";
  const char* folderName = "CorrelationTree";
  if (!outputFileName) outputFileName = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("histos", TList::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("events", TTree::Class(),AliAnalysisManager::kOutputContainer,Form("%s:%s", outputFileName, folderName));
  mgr->ConnectInput  (ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput (ana, 1, coutput1);
  mgr->ConnectOutput (ana, 2, coutput2);
  
  if (!mgr->InitAnalysis()) return;

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode("full");
  plugin->SetNtestFiles(2);
  plugin->SetAPIVersion("V1.1x");
  plugin->SetAliPhysicsVersion("vAN-20150319");
  plugin->SetGridDataDir(dataDir.Data());
  plugin->SetDataPattern(dataPattern.Data());
  plugin->SetGridWorkingDir(workingDir.Data());
  plugin->SetRunPrefix("000");
  for (Int_t i=0;i<nRuns;i++)  plugin->AddRunNumber(runList[i]);
  plugin->SetGridOutputDir("output");
  plugin->SetAnalysisSource("AliAnalysisTaskCFTreeNew.cxx");
  plugin->AddIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  plugin->SetAdditionalLibs("libOADB.so libCORRFW.so libPWGTools.so libPWGmuon.so libPWGCFCorrelationsBase.so AliAnalysisTaskCFTreeNew.cxx AliAnalysisTaskCFTreeNew.h");
  plugin->SetNrunsPerMaster(1);
  plugin->SetSplitMaxInputFileNumber(100);

  mgr->SetGridHandler(plugin);
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");

//  TChain *chain = new TChain("aodTree");
//  chain->AddFile("/alice/data/2013/LHC13d/000195873/pass2/AOD/001/AliAOD.root");
//  mgr->StartAnalysis("local",chain,100);

}


