/*
  Simple calibration analysis
  //
  //0. Setup memory chcecker if you want 
  //
  TMemStat *memstat = new TMemStat("new,gnubuildin");
  AliSysInfo::AddCallBack(TMemStatManager::GetInstance()->fStampCallBack);

  AliSysInfo::AddStamp("Start");  
  //1. Load needed libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  // Setup analysis manager
  //
  TString path=gSystem->pwd();
  gROOT->Macro(Form("%s/ConfigOCDB.C",path->Data()));

  .L CalibrateTPC.C
  AliAnalysisManager * mgr = ( AliAnalysisManager *)SetupCalibTask("/V6/");  

  //
  // Process data - chain
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("esd.txt","esdTree",0,10);
  //chain->Lookup();
  // memory
  mgr->SetNSysInfo(1000); 
  //
  mgr->SetDebugLevel(10);
  mgr->StartAnalysis("proof",chain);
  //mgr->StartAnalysis("local",chain);//
  // delete manager
  //
  delete mgr;
  AliSysInfo::AddStamp("End");
  //
  // analyze memstat report
  //
  delete memstat;
  TMemStat mem;
  mem.MakeReport(0,0,"order 0 sortstat 3 sortstamp 0 sortdeep 10 stackdeep 15 maxlength 50")   
*/

char * prefix = "/V6/";
void SetupCalibTask(AliTPCAnalysisTaskcalib* task1);

TObject  * SetupCalibTask(char * tprefix ="/V12/") {
  //
  //
  //
  prefix=tprefix;
  TStopwatch stopwatch;
  stopwatch.Start();
  AliAnalysisManager *mgr=new AliAnalysisManager("TestManager");
  AliESDInputHandler* esdH=new AliESDInputHandler;
  esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH);  
  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("TPC calibration task");
  //
  SetupCalibTask(task1);
  //
  mgr->AddTask(task1);
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutput1 =mgr->CreateContainer("TPCCalib",TObjArray::Class(), AliAnalysisManager::kOutputContainer, "CalibObjects.root");  
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  //
  if (!mgr->InitAnalysis()) return 0;
  mgr->PrintStatus();   
  stopwatch.Stop();
  stopwatch.Print();
  return mgr;
}



void SetupCalibTask(AliTPCAnalysisTaskcalib* task1){
  //
  // Configure calibration task
  //
  AliTPCClusterParam * clusterParam = AliTPCcalibDB::Instance()->GetClusterParam();  
  AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(40, 0.4, 0.5, 0.13, 0.018);  
  //
  AliTPCcalibTracks *calibTracks =  new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts); 


  AliTPCcalibTracksGain *calibTracksGain =  new AliTPCcalibTracksGain("calibTracksGain","Gain calibration using tracks",cuts); 
  AliTPCcalibAlign *calibAlign = new AliTPCcalibAlign("alignTPC","Alignment of the TPC sectors");
  AliTPCcalibAlign *calibAlignAll = new AliTPCcalibAlign("alignTPCAll","Alignment of the TPC sectors- All");
  AliTPCcalibLaser *calibLaser = new AliTPCcalibLaser("laserTPC","laserTPC");
  AliTPCcalibCosmic *calibCosmic = new AliTPCcalibCosmic("cosmicTPC","cosmicTPC");
  AliTPCcalibCalib *calibCalib = new AliTPCcalibCalib("calibTPC","calibTPC");
  TTimeStamp startTime(2008,9,0,0,0,0);
  TTimeStamp stopTime(2008,11,0,0,0,0);
  AliTPCcalibTime *calibTime = new AliTPCcalibTime("calibTime","calibTime", startTime.GetSec(), stopTime.GetSec(), 5*60, 5*60);

  AliTPCcalibUnlinearity *calibUnlinearity = new AliTPCcalibUnlinearity("calibUnlinearity","calibUnlinearity");
  AliTPCcalibUnlinearity *calibUnlinearityAll = new AliTPCcalibUnlinearity("calibUnlinearityAll","calibUnlinearityAll");
  //
  //
  //
  AliExternalComparison *compCosmic = new AliExternalComparison("compCosmic","compCosmic");
  compCosmic->SetDefaultRange(1,160,20);
  //
  compCosmic->SetResolRange(0,-1.5,1.5,100);
  compCosmic->SetResolRange(1,-1,1,100);
  compCosmic->SetResolRange(2,-0.02,0.02,100);
  compCosmic->SetResolRange(3,-0.02,0.02,100);
  compCosmic->SetResolRange(4,-0.15,0.15,100);
  //
  compCosmic->SetParameterRange(0,-30,30,10);
  compCosmic->SetParameterRange(1,-250,250,10);
  compCosmic->SetParameterRange(2,-1,1,10);
  compCosmic->SetParameterRange(3,-1.5,1.5,20);
  compCosmic->SetParameterRange(4,-2,2,20);
  compCosmic->SetParameterRange(5,0,90,10);
  compCosmic->SetParameterRange(6,-3.14,3.14,20);
  //
  compCosmic->SetDistCut(3,8,0.06,0.06,0.2);
  compCosmic->SetPullDistCut(10,20,20,20,10);
  calibCosmic->SetComparison(compCosmic);
  //
  //
  calibTracks->SetDebugLevel(2);
  calibTracks->SetStreamLevel(20);
  calibTracksGain->SetDebugLevel(2);
  calibTracksGain->SetStreamLevel(20);
  calibAlign->SetDebugLevel(20);
  calibAlign->SetStreamLevel(10);
  calibAlignAll->SetDebugLevel(20);
  calibAlignAll->SetStreamLevel(10);
  calibLaser->SetDebugLevel(20);
  calibLaser->SetStreamLevel(20);
  calibCosmic->SetDebugLevel(20);
  calibCosmic->SetStreamLevel(20);
  calibCalib->SetDebugLevel(20);
  calibCalib->SetStreamLevel(20);
  calibTime->SetDebugLevel(20);
  calibTime->SetStreamLevel(10);
  //
  calibUnlinearity->SetDebugLevel(20);
  calibUnlinearity->SetStreamLevel(10);
  calibUnlinearityAll->SetDebugLevel(20);
  calibUnlinearityAll->SetStreamLevel(10);
  
  //  calibCalib->SetTriggerMask(-1,-1,kFALSE);       //accept everything 
  calibCalib->SetTriggerMask(-1,-1,kTRUE);       //accept everything  - except laser
  calibTracks->SetTriggerMask(-1,16,kTRUE);      //reject laser trigger, accept everything else
  calibTracksGain->SetTriggerMask(-1,16,kTRUE);  //reject laser trigger, accept everything else
  calibAlign->SetTriggerMask(-1,-1,kTRUE);       //accept everything 
  calibAlignAll->SetTriggerMask(-1,-1,kFALSE);   //accept everything 
  calibLaser->SetTriggerMask(-1,-1,kFALSE);       //accept only laser trigger
  calibCosmic->SetTriggerMask(-1,-1,kTRUE);      //reject laser trigger, accept everything else
  calibTime->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  calibUnlinearity->SetTriggerMask(-1,-1,kTRUE);   //reject laser 
  calibUnlinearityAll->SetTriggerMask(-1,-1,kFALSE);   //non reject laser 


  //
 // ---*---*-----*-*-----*----------*---
  // ADD CALIB JOBS HERE!!!!!!!!!!!!!!!!
  task1->AddJob(calibCalib);
  task1->AddJob(calibAlign);
  //task1->AddJob(calibAlignAll);
  //task1->AddJob(calibLaser);
  task1->AddJob(calibCosmic);
  //task1->AddJob(calibTime);

  //task1->AddJob(calibTracksGain);
  //task1->AddJob(calibTracks);
  //task1->AddJob(calibUnlinearity);
  //task1->AddJob(calibUnlinearityAll);
  // task1->AddJob(new AliTPCcalibV0);
  // -*----*----*---*-*------*-------**--
  // -------*--*---------*-----*-------*-
  TString path=gSystem->pwd();
  path+=prefix;
  gSystem->mkdir(path);
  task1->SetDebugOuputhPath(path.Data());

}




void CalibrateTPC(Float_t magf, Int_t first, Int_t last){
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  gROOT->Macro(Form("$ALICE_ROOT/TPC/macros/ConfigOCDB.C(%f)",magf));
  //
  // Setup analysis manager
  //
  //.L $ALICE_ROOT/TPC/macros/CalibrateTPC.C
  AliAnalysisManager * mgr = (AliAnalysisManager*)SetupCalibTask("/V3/");
  //
  // Process data - chain
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("esd.txt","esdTree",0,last-first,last);
  chain->Lookup();
  // memory
  mgr->SetNSysInfo(5000); 
  //
  mgr->SetDebugLevel(1);
  mgr->StartAnalysis("local",chain);
  
}
