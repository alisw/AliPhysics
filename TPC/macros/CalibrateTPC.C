/*
  Simple calibration analysis
  //
  //0. Setup memory chcecker if you want 
  //
  gSystem->Load("$ROOTSYS/lib/libGui.so");
  gSystem->Load("$ROOTSYS/lib/libTree.so");
  gSystem->Load("$MEMSTAT/libMemStat.so");
  TMemStat *memstat = new TMemStat(100000000,10000000,kTRUE);
  AliSysInfo::AddCallBack(TMemStatManager::GetInstance()->fStampCallBack);

  AliSysInfo::AddStamp("Start");  
  //1. Load needed libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  // Setup analysis manager
  //
  .L $ALICE_ROOT/TPC/macros/CalibrateTPC.C
  AliAnalysisManager * mgr = SetupCalibTask();
  //
  // Process data - chain
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("esd.txt","esdTree",0,1200);
  chain->Lookup();
  // memory
  mgr->SetNSysInfo(100); 
  //
  mgr->SetDebugLevel(1);
  mgr->StartAnalysis("local",chain);
  // delete manager
  //
  delete mgr;
  AliSysInfo::AddStamp("End");
  //
  // analyze memstat report
  //
  delete memstat;
  TMemStat draw("memstat.root");
  draw.MakeReport(0,0,"order 0 sortstat 3 sortstamp 0 sortdeep 10 stackdeep 15 maxlength 50")   
*/


AliAnalysisManager * SetupCalibTask() {
  //
  //
  //
  TStopwatch stopwatch;
  stopwatch.Start();
  //
  // set magnetic field form the cosmos - it should be provided by framework
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 2);
  AliTracker::SetFieldMap(field,0);
  TGeoManager::Import("geometry.root");
  //
  AliAnalysisManager *mgr=new AliAnalysisManager("TestManager");

  AliESDInputHandler* esdH=new AliESDInputHandler;
  esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);  
  //
  //
  AliCDBManager::Instance()->SetRun(1) ;
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliTPCClusterParam * clusterParam = AliTPCcalibDB::Instance()->GetClusterParam();

  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("TPC calibration task");
  
  AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);

  //
  AliTPCcalibTracks *calibTracks =  new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts); 
  AliTPCcalibTracksGain *calibTracksGain =  new AliTPCcalibTracksGain("calibTracksGain","Gain calibration using tracks",cuts); 
  AliTPCcalibAlign *calibAlign = new AliTPCcalibAlign("alignTPC","Alignment of the TPC sectors");
  AliTPCcalibLaser *calibLaser = new AliTPCcalibLaser("laserTPC","laserTPC");
  AliTPCcalibCosmic *calibCosmic = new AliTPCcalibCosmic("cosmicTPC","cosmicTPC");
  calibTracks->SetDebugLevel(0);
  calibTracks->SetStreamLevel(0);
  calibTracksGain->SetDebugLevel(0);
  calibTracksGain->SetStreamLevel(0);
  calibAlign->SetDebugLevel(20);
  calibAlign->SetStreamLevel(2);
  calibLaser->SetDebugLevel(20);
  calibLaser->SetStreamLevel(2);
  calibCosmic->SetDebugLevel(20);
  calibCosmic->SetStreamLevel(2);

  //
 // ---*---*-----*-*-----*----------*---
  // ADD CALIB JOBS HERE!!!!!!!!!!!!!!!!
  //task1->AddJob(calibAlign);
  //  task1->AddJob(calibLaser);
  task1->AddJob(calibCosmic);
  //task1->AddJob(calibTracksGain);
  //task1->AddJob(calibTracks);
  //  task1->AddJob(new AliTPCcalibBase);
  // task1->AddJob(new AliTPCcalibV0);
  // -*----*----*---*-*------*-------**--
  // -------*--*---------*-----*-------*-
  task1->SetDebugOuputhPath("/lustre_alpha/alice/miranov/rec/cosmic_jun2008/");
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
