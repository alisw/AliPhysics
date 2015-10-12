class  AliAnalysisManager;
class  AliAnalysisAlien;

void runGridDeuteronTree(TString mode="terminate", TString fname="DeuteronTreeTest")
{

  AliLog::SetGlobalDebugLevel(5);

  //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gSystem->Load("libPWGTools.so");
  gSystem->Load("libPWGLFspectra.so");

  // Create and configure the alien handler plugin
  AliAnalysisGrid *alienHandler = CreateAlienHandler(mode,fname);
  if (!alienHandler) return;
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  Bool_t mc = kFALSE;
  //if(sub == 2)
  //mc = kTRUE;
  
  // Add PIDResponse task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(mc,kTRUE,kTRUE);

  // Add Physics Selec
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
    
  // Add centrality
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
  //taskCentrality->SetPass(2); // Change if changing pass


  // Add AliPPVsMultUtils task
  // AliAnalysisUtils* utils = new AliAnalysisUtils();

  // add DPhiCorrelations task from $ALICE_PHYSICS/../src/PWGCF/Correlations/DPhi/AliAnalysisTaskPhiCorrelations[.cxx .h]
  //gROOT->LoadMacro("AliAnalysisTaskPhiCorrelations.cxx+g");
  gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGLF/SPECTRA/Nuclei/deuteronpA/macros/AddTaskDeuteronTree.C");

  AliAnalysisDeuteronTree *taskDeuteronTree = AddTaskDeuteronTree("DeuteronTreeTest");

  taskDeuteronTree->SelectCollisionCandidates(AliVEvent::kINT7);

    //mgr->SetDebugLevel(5);

  if (!mgr->InitAnalysis())
    return;
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid",100);

  // local analysis
  //TChain *chain = new TChain("aodTree");
  //chain->AddFile("./traintest/__alice__data__2010__LHC10d__000126158__ESDs__pass2__AOD147_root_archive_AliAOD_2/0002/root_archive/AliAOD.root");
  //mgr->StartAnalysis("local",chain);

}


AliAnalysisGrid* CreateAlienHandler(TString mode="test",TString fname="testName"){
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(200);
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  //plugin->SetROOTVersion("v5-34-26");
  plugin->SetAliROOTVersion("v5-06-42");
  plugin->SetAliPhysicsVersion("vAN-20151006");
  // Declare input data to be processed.
  plugin->SetNrunsPerMaster(100);
  plugin->SetSplitMaxInputFileNumber(15); // 3 in the LEGO trains

    plugin->SetGridDataDir("/alice/data/2013/LHC13c/");
    plugin->SetGridWorkingDir(Form("%s/",fname.Data()));
    plugin->SetRunPrefix("000"); 
    plugin->SetDataPattern("/ESDs/pass2/*/AliESDs.root");
    plugin->SetAnalysisMacro("TaskDeuteronTreeTest.C");
    plugin->SetExecutable("TaskDeuteronTreeTest.sh");
    plugin->SetJDLName("TaskDeuteronTreeTest.jdl");
    Int_t runnumbers[]={195596,195592};

    for(Int_t irun=0;irun<2;irun++){
      Printf("%d %d",irun,runnumbers[irun]);
      plugin->AddRunNumber(runnumbers[irun]);
    }

 

    TString extraLibs;
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  
  //plugin->SetAdditionalLibs("libCORRFW.so libPWGTools.so libPWGCFCorrelationsBase.so libPWGCFCorrelationsDPhi.so AliAnalysisTaskPhiCorrelations.cxx AliAnalysisTaskPhiCorrelations.h");
  //plugin->SetAnalysisSource("AliAnalysisTaskPhiCorrelations.cxx");

  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  //plugin->SetDefaultOutputs(0);
  //plugin->SetOutputFiles("AnalysisResults.root.root");
  plugin->SetMergeViaJDL(kTRUE);
  plugin->SetOneStageMerging(kFALSE);
  //plugin->SetMaxMergeFiles(40);
  //plugin->SetMaxMergeStages(4);
  
  plugin->SetTTL(50000);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  //plugin->SetSplitMaxInputFileNumber();
  plugin->SetSplitMode("se");
  return plugin;
}
