void runGrid() {

  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  //  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  //  gSystem->Load("libPWGPP");

  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler();  
  if (!alienHandler) return;

  cout<<" Got alien handler"<<endl;
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  mgr->SetGridHandler(alienHandler);

  /*
  // For centrality Selection: Requires input event handler.
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  taskCentrality->SetPass(2); // remember to set the pass you are processing!!! 
  */
  // My Task
  gROOT->LoadMacro("AliAnalysisFBMultFluct.cxx++g");
  AliAnalysisTask *task = new AliAnalysisFBMultFluct("TaskFBGrid");

  //  task->SelectCollisionCandidates(AliVEvent::kMB) // Added for minbias trigger on 3Jun14
  //  task->SelectCollisionCandidates(AliVEvent::kCentral) // Added for minbias trigger on 3Jun14

  mgr->AddTask(task);

  /*
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);
  TChain *chain = new TChain("aodTree");
  //  chain->Add("/Users/sadhanadash/QA/AliESDs.root");
  chain->Add("/Users/rashmi/ALICE/InitFluct/alice/data/2010/LHC11h/000170388/ESDs/pass2_muon/AOD/132/AliAOD.root");
  //  chain->Add("/Users/rashmi/ALICE/InitFluct/alice/data/2011/LHC11h/000170388/ESDs/pass2_muon/AOD/132/AliESDs.root");
  chain->Add("AliESDs.root");
  */

  
  //  AliAODInputHandler* aodH = new AliAODInputHandler;
    AliVEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  cout<<" SetInputEventHandler"<<endl;
  
  //   TChain *chain = new TChain("aodTree");
   //alice/data/2010/LHC10h/000139510/ESDs/pass2/AOD160/0001
  //chain->Add("AODfile/139510/AliAOD.root");
  //  chain->Add("pass2/137161/1/AliAOD.root");
  //  chain->Add("pass2/137161/2/AliAOD.root");
  // chain->GetListOfFiles()->Print();

  // Create containers for input/output

  AliAnalysisDataContainer *coutput = 
    //     mgr->CreateContainer("coutput", TObjArray::Class(),
    //     mgr->CreateContainer("coutput", TList::Class(),
    //     AliAnalysisManager::kOutputContainer, "MC.AOD.root" );
    mgr->CreateContainer("coutput", TList::Class(),
			 AliAnalysisManager::kOutputContainer, "AOD.DATA.root" );

  AliAnalysisDataContainer *cTree = 
    mgr->CreateContainer("cTree",TTree::Class(),AliAnalysisManager::kOutputContainer,"AOD.DATA.root");
  
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput);
  mgr->ConnectOutput(task,2,cTree);
  /*
  mgr->CreateContainer("cinput",TFile::Class(),1,"AOD.137161.input.root");
  mgr->ConnectInput(task,2,cinput); 
  */
  mgr->SetDebugLevel(1);
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  cout<<" Starting Analysis now"<<endl;
  //  mgr->StartAnalysis("local",chain);
  mgr->StartAnalysis("grid");
}

