/*
  Simple calibration analysis

  
  //1. Load needed libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //increased memstat
  gSystem->Load("$ROOTSYS/lib/libGui.so");
  gSystem->Load("$ROOTSYS/lib/libTree.so");
  gSystem->Load("$MEMSTAT/libMemStat.so");
  TMemStat memstat(100000000,10000000,kTRUE);
  memstat->AddStamp("aaaa");
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
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("chain.txt","esdTree",0,30);
  chain->Lookup();
  // memory
  mgr->SetNSysInfo(100); 
  AliSysInfo::AddCallBack(TMemStatManager::GetInstance()->fStampCallBack);
  //
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
  //
  //
  AliCDBManager::Instance()->SetRun(1) ;
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliTPCClusterParam * clusterParam = AliTPCcalibDB::Instance()->GetClusterParam();

  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("foo bar");
  
  AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);

  AliTPCcalibTracks *calibTracks =  new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts); 
  AliTPCcalibTracksGain *calibTracksGain =  new AliTPCcalibTracksGain("TPCGainTracks","TPCGainTracks",cuts); 
  calibTracks->SetDebugLevel(5);
  calibTracks->SetStreamLevel(5);
  calibTracksGain->SetDebugLevel(1);
  calibTracksGain->SetStreamLevel(1);
 // ---*---*-----*-*-----*----------*---
  // ADD CALIB JOBS HERE!!!!!!!!!!!!!!!!
  task1->AddJob(new AliTPCcalibAlign);//"align","The kewl alignment job"));
  task1->AddJob(calibTracksGain);
  task1->AddJob(calibTracks);
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
