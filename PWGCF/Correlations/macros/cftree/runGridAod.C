// example macro for configuring and running AliAnalysisTaskCFTree task on the grid

// copy this file to your local running directory, and copy env_grid.h to a
// directory named with the name of the period
// then everything can be configured within env_grid.h for a given period, and
// no changes should be necessary in this runGridAod.C


void runGridAod(char * period="lhc10d",char * runMode="test"){
  gROOT->LoadMacro(Form("%s/env_grid.h",period));

  if (!TGrid::Connect("alien://")) return;
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
  mgr->SetDebugLevel(0);
  AliAODInputHandler* aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  if(!isMC)
    {
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
      AliMultSelectionTask* multTask = AddTaskMultSelection();
      if(!is13TeV)
	multTask->SetSelectedTriggerClass(AliVEvent::kMB); 
    }

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); 
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,kTRUE);
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/Correlations/macros/cftree/AddTaskCFTree.C");
  AliAnalysisTaskCFTree *taskCFTree = AddTaskCFTree(); 
  //gROOT->LoadMacro("AliAnalysisTaskCFTree.cxx++g");
  //AliAnalysisTaskCFTree *ana = new AliAnalysisTaskCFTree();

  //event selection
  taskCFTree->SetIs13TeV(is13TeV);
  taskCFTree->SetClassBit(classbit);
  taskCFTree->SetEventSelectionBit(eventselectionbit);
  taskCFTree->SetZVertex(zvertex);

  //track selection
  taskCFTree->SetPtMin(ptmin);
  taskCFTree->SetTrackFilterBit(trackfilterbit);
  taskCFTree->SetTrackEtaCut(tracketacut);

  //tracklet selection
  taskCFTree->SetDphiCut(dphicut);
  taskCFTree->SetTrackletEtaCut(trackletetacut);

  //tree setting
  taskCFTree->SetStoreTracks(storeTracks);
  taskCFTree->SetStoreTracklets(storeTracklets);
  taskCFTree->SetStoreMuons(storeMuons);
  taskCFTree->SetStoreTrackInfo(storeTrackInfo);
  taskCFTree->SetStorePidInfo(storePidInfo);
  taskCFTree->SetApplyPhysicsSelectionCut(applyPhysicsSelectionCut);
  taskCFTree->SetStoreOnlyEventsWithMuons(storeOnlyEventsWithMuons);
  taskCFTree->SetStoreCutBitsInTrackMask(storeCutBitsInTrackMask);
  if(isMC){
    taskCFTree->SetStoreMcTracks(storeMcTracks);
    taskCFTree->SetStoreMcTracklets(storeMcTracklets);
    taskCFTree->SetStoreMcMuons(storeMcMuons);
  }

  if (!mgr->InitAnalysis()) return;

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  plugin->SetRunMode(runMode);
  plugin->SetNtestFiles(ntestfiles);
  plugin->SetKeepLogs(keeplogs);
  plugin->SetTTL(86000);
  //plugin->SetAPIVersion(APIVersion);
  //plugin->SetAliROOTVersion(AliROOTVersion);
  plugin->SetAliPhysicsVersion(AliPhysicsVersion);
  plugin->SetGridDataDir(dataDir.Data());
  plugin->SetDataPattern(dataPattern.Data());
  plugin->SetGridWorkingDir(Form("trees_pp/%s",workingDir.Data()));
  if(!isMC) plugin->SetRunPrefix("000");
  for (Int_t i=0;i<nRuns;i++) plugin->AddRunNumber(runList[i]);
  plugin->SetGridOutputDir("output");
  //plugin->SetAnalysisSource("AliAnalysisTaskCFTree.cxx");
  //plugin->SetAdditionalLibs("libpythia6_4_25.so libCORRFW.so libPWGTools.so libPWGCFCorrelationsBase.so AliAnalysisTaskCFTree.cxx AliAnalysisTaskCFTree.h");
  //plugin->SetAdditionalLibs("libpythia6_4_25.so libCORRFW.so libPWGTools.so libPWGCFCorrelationsBase.so");
  plugin->SetNrunsPerMaster(NrunsPerMaster);
  plugin->SetSplitMaxInputFileNumber(maxinputfilenumber);
  plugin->SetMergeViaJDL();
  plugin->SetMaxMergeStages(maxmergestage);
  plugin->SetCheckCopy(checkcopy);

  mgr->SetGridHandler(plugin);
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");

}

