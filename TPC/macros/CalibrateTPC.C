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
  AliAnalysisManager * mgr = SetupCalibTask("/V6/");
  //
  // Process data - chain
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("esd.txt","esdTree",0,100);
  chain->Lookup();
  // memory
  mgr->SetNSysInfo(100); 
  //
  mgr->SetDebugLevel(1);
  mgr->StartAnalysis("proof",chain);
  //mgr->StartAnalysis("local",chain);
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


AliAnalysisManager * SetupCalibTask(char * prefix ="/V12/") {
  //
  //
  //
  TStopwatch stopwatch;
  stopwatch.Start();
  //
  // set magnetic field form the cosmos - it should be provided by framework
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 2);
  AliTracker::SetFieldMap(field,0);
  TGeoManager::Import("/u/miranov/proof/geometry.root");
  //
  TFile f("/u/miranov/GainMap.root");
  AliTPCCalPad *gainMap = f.Get("GainMap");
  //
  // OCDB setup
  //
  AliCDBManager::Instance()->SetRun(1);
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliTPCcalibDB::Instance()->SetExBField(0);
  AliTPCClusterParam * param = AliTPCcalibDB::Instance()->GetClusterParam();
  AliTPCClusterParam::SetInstance(param);
  AliTPCcalibDB::Instance()->SetExBField(0);
  //


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
  
  AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(40, 0.4, 0.5, 0.13, 0.018);

  //
  AliTPCcalibTracks *calibTracks =  new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts); 


  AliTPCcalibTracksGain *calibTracksGain =  new AliTPCcalibTracksGain("calibTracksGain","Gain calibration using tracks",cuts); 
  AliTPCcalibAlign *calibAlign = new AliTPCcalibAlign("alignTPC","Alignment of the TPC sectors");
  AliTPCcalibLaser *calibLaser = new AliTPCcalibLaser("laserTPC","laserTPC");
  AliTPCcalibCosmic *calibCosmic = new AliTPCcalibCosmic("cosmicTPC","cosmicTPC");
  AliTPCcalibCalib *calibCalib = new AliTPCcalibCalib("calibTPC","calibTPC");
  TTimeStamp startTime(2008,9,0,0,0,0);
  TTimeStamp stopTime(2008,11,0,0,0,0);
  AliTPCcalibTime *calibTime = new AliTPCcalibTime("cosmicTime","cosmicTime",0, startTime.GetSec(), stopTime.GetSec(), 5*60, 5*60);

  calibCosmic->SetGainMap(gainMap);
  calibTracksGain->SetGainMap(gainMap);
  //
  calibTracks->SetDebugLevel(20);
  calibTracks->SetStreamLevel(20);
  calibTracksGain->SetDebugLevel(2);
  calibTracksGain->SetStreamLevel(20);
  calibAlign->SetDebugLevel(20);
  calibAlign->SetStreamLevel(10);
  calibLaser->SetDebugLevel(20);
  calibLaser->SetStreamLevel(20);
  calibCosmic->SetDebugLevel(20);
  calibCosmic->SetStreamLevel(2);
  calibCalib->SetDebugLevel(20);
  calibCalib->SetStreamLevel(10);

  //
 // ---*---*-----*-*-----*----------*---
  // ADD CALIB JOBS HERE!!!!!!!!!!!!!!!!
  task1->AddJob(calibCalib);
  task1->AddJob(calibAlign);
  task1->AddJob(calibLaser);
  task1->AddJob(calibCosmic);
  task1->AddJob(calibTime);

  task1->AddJob(calibTracksGain);
  task1->AddJob(calibTracks);
  // task1->AddJob(new AliTPCcalibV0);
  // -*----*----*---*-*------*-------**--
  // -------*--*---------*-----*-------*-
  TString path=gSystem->pwd();
  path+=prefix;
  gSystem->mkdir(path);
  task1->SetDebugOuputhPath(path.Data());

  mgr->AddTask(task1);

  AliAnalysisDataContainer *cinput1
    =mgr->CreateContainer("cchain1",TChain::Class(),
			  AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1
    =mgr->CreateContainer("TPCCalib",TObjArray::Class(),
			  AliAnalysisManager::kOutputContainer,
			  "CalibObjects.root");
  
  coutput1->SetSpecialOutput(kTRUE);
  //coutput1->SetFileName("CalibObjectFile.root");
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  //
  //mgr->SetSpecialOutputLocation(path->Data());

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus(); 
  
  stopwatch.Stop();
  stopwatch.Print();
  return mgr;
}


void Merge(){
  fstream finput;
  TString currentFile("");
  finput.open("mergelist.txt", ios_base::in);
  TFileMerger merger;
  merger.OutputFile("result.root");  
  finput >> currentFile;
  merger->AddFile(currentFile->Data());
  merger.Merge();
  TFile::Cp("result.root","last.root");
  //
  while(finput.good()) {
    TFileMerger  merger2;
    merger2.OutputFile("result.root");
    finput >> currentFile;    
    merger2.AddFile("last.root");
    merger2.AddFile(currentFile->Data());    
    merger2.Merge(); 
    TFile::Cp("result.root","last.root");
  }
}
