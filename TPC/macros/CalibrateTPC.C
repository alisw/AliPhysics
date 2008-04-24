/*
  Simple calibration analysis

  
  //1. Load needed libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  // Setup analysis manager
  //
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  .L $ALICE_ROOT/TPC/macros/CalibrateTPC.C
  AliAnalysisManager * mgr = SetupCalibTask();
  //
  // Process data - chain
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("cahin.txt","esdTree",0,10)
  mgr->StartAnalysis("local",chain);
  

*/


AliAnalysisManager * SetupCalibTask() {
  //
  //
  //
  TStopwatch stopwatch;
  stopwatch.Start();

  AliAnalysisManager *mgr=new AliAnalysisManager("TestManager");

  AliESDInputHandler* esdH=new AliESDInputHandler;
  esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);  

  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("foo bar");
  
  AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);

  // ---*---*-----*-*-----*----------*---
  // ADD CALIB JOBS HERE!!!!!!!!!!!!!!!!
  task1->AddJob(new AliTPCcalibAlign);//"align","The kewl alignment job"));
  //  task1->AddJob(new AliTPCcalibTracks("resolution","I would have been called AliTPCcalibResolution in a bit more perfect world.",0,cuts));
  task1->AddJob(new AliTPCcalibTracksGain("resolution","I would have been called AliTPCcalibGain in a bit more perfect world.",cuts));
  //  task1->AddJob(new AliTPCcalibBase);
  // task1->AddJob(new AliTPCcalibV0);
  // -*----*----*---*-*------*-------**--
  // -------*--*---------*-----*-------*-

  mgr->AddTask(task1);

  AliAnalysisDataContainer *cinput1
    =mgr->CreateContainer("cchain1",TChain::Class(),
			  AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1
    =mgr->CreateContainer("asdofhaw",TObjArray::Class(),
			  AliAnalysisManager::kOutputContainer,
			  "CalibObjects.root");

  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus(); 
  
  stopwatch.Stop();
  stopwatch.Print();
  return mgr;
}
