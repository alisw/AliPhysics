/*
  Simple test of the V0 finder
  //
  //0. Setup memory chcecker if you want 
  //
  gSystem->Load("$ROOTSYS/lib/libGui.so");
  gSystem->Load("$ROOTSYS/lib/libTree.so");
  gSystem->Load("$MEMSTAT/libMemStat.so");
  TMemStat *memstat = new TMemStat(100000000,10000000,kTRUE);
  AliSysInfo::AddCallBack(TMemStatManager::GetInstance()->fStampCallBack);
  AliSysInfo::AddStamp("Start");  
  //

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
  TChain * chain = tool.MakeChain("esd.txt","esdTree",0,50000);
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


AliAnalysisManager * SetupV0Task() {
  //
  //
  //
  TStopwatch stopwatch;
  stopwatch.Start();
  //
  AliAnalysisManager *mgr=new AliAnalysisManager("TestManager");

  AliESDInputHandler* esdH=new AliESDInputHandler;
  esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);  
  //
  //
  AliCDBManager::Instance()->SetRun(1) ;
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");

  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("TPC calibration task");
  
  AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(20, 0.4, 0.5, 0.13, 0.018);

  //
  AliTPCcalibV0 *calibV0 = new AliTPCcalibV0;
  calibV0->SetName("calibV0");
  calibV0->SetTitle("calibV0");
  calibV0->SetDebugLevel(20);
  calibV0->SetStreamLevel(2);
  //
  
  task1->AddJob(calibV0);
 
  TString path=gSystem->pwd();
  path+="/V0/";
  gSystem->mkdir(path);
  task1->SetDebugOuputhPath(path.Data());
  mgr->AddTask(task1);

  mgr->AddTask(task1);

  AliAnalysisDataContainer *cinput1
    =mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1
    =mgr->CreateContainer("TPCCalib",TObjArray::Class(),
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
