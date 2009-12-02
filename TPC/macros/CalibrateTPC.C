/*
  Simple calibration analysis
  //
  //0. Setup memory chcecker if you want 
  //
  TMemStat *memstat = new TMemStat("new,gnubuildin");
  AliSysInfo::AddCallBack(TMemStatManager::GetInstance()->fStampCallBack);
  AliSysInfo::AddStamp("Start"); 

  gSystem->Load("$ROOTSYS/lib/libXrdClient.so");
  gSystem->Load("libNetx.so");

  gSystem->Setenv("alien_CLOSE_SE","ALICE::GSI::SE") 
  TGrid * alien =     TGrid::Connect("alien://",0,0,"t"); 
  gSystem->Setenv("alien_CLOSE_SE","ALICE::GSI::SE") 
 
  //1. Load needed libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //
  // Setup analysis manager
  //
  TString path=gSystem->pwd();
  gROOT->Macro(Form("%s/ConfigOCDB.C\(%d\)",path->Data(),0));

  .L CalibrateTPC.C
  AliAnalysisManager * mgr = ( AliAnalysisManager *)SetupCalibTask("/V6/");  

  //
  // Process data - chain
  // 
  //gEnv->SetValue("TFile.Recover", 0);  // dont try to recover anything
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chain = tool.MakeChain("esd.txt","esdTree",0,200000);
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


void SetupCalibTaskTrain1(TObject * task1);
void SetupCalibTaskTrain2(TObject * task2);

char * prefix = "/V6/";
// Global parameters to set
TTimeStamp startTime(2009,8,7,0,0,0);
TTimeStamp stopTime(2009,12,31,0,0,0);
Int_t debugLevel  = 2;
Int_t streamLevel = 20;

//

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

  //
  // Train 1 - to be run always on full statistic (
  //
  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
  //
  SetupCalibTaskTrain1(task1);
  //
  mgr->AddTask(task1);
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutput1 =mgr->CreateContainer("TPCCalib",TObjArray::Class(), AliAnalysisManager::kOutputContainer, "CalibObjectsTrain1.root");  
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  //
  //
  //
  AliTPCAnalysisTaskcalib *task2=new AliTPCAnalysisTaskcalib("CalibObjectsTrain2");
  //
  SetupCalibTaskTrain2(task2);
  //
  mgr->AddTask(task2);
  AliAnalysisDataContainer *cinput2 = mgr->GetCommonInputContainer();
  
  if (!cinput2) cinput2 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutput2 =mgr->CreateContainer("TPCCalib2",TObjArray::Class(), AliAnalysisManager::kOutputContainer, "CalibObjectsTrain2.root");  
  mgr->ConnectInput(task2,0,cinput2);
  mgr->ConnectOutput(task2,0,coutput2);



  //
  if (!mgr->InitAnalysis()) return 0;
  mgr->PrintStatus();   
  stopwatch.Stop();
  stopwatch.Print();
  return mgr;
}

void AddCalibCalib(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // calibCalib is a prefilter 
  // The current OCDB entries transformation are applied on cluster, tracks are refitted
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibCalib *calibCalib = new AliTPCcalibCalib("calibTPC","calibTPC");
  calibCalib->SetDebugLevel(debugLevel);
  calibCalib->SetStreamLevel(0);
  calibCalib->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  myTask->AddJob(calibCalib);

}
void AddCalibTimeGain(TObject* task){
  //
  //  Responsible: Alexander Kalweit
  //  Description:
  //  Parameters to set
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTimeGain *calibTimeGain = new AliTPCcalibTimeGain("calibTimeGain","calibTimeGain", startTime.GetSec(), stopTime.GetSec(), 30*60);
  //calibTimeGain->SetLowMemoryConsumption(kTRUE);
  //calibTimeGain->SetMIP(25.);
  calibTimeGain->SetUseCookAnalytical(kTRUE);
  calibTimeGain->SetUseMax(kFALSE);
  calibTimeGain->SetDebugLevel(debugLevel);
  calibTimeGain->SetStreamLevel(streamLevel);
  calibTimeGain->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  myTask->AddJob(calibTimeGain);
}

void AddCalibTime(TObject* task){
  //
  // Responsible: Dag Larsen
  // Description:
  //
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTime *calibTime = new AliTPCcalibTime("calibTime","calibTime",  startTime.GetSec(), stopTime.GetSec(), 20*60);
  calibTime->SetDebugLevel(debugLevel);
  calibTime->SetStreamLevel(streamLevel);
  calibTime->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  myTask->AddJob(calibTime);
}

void AddCalibTrigger(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  //             Export trees with summary information for the trigger efficieciency and purity study
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTrigger *calibTrigger = new AliTPCcalibTrigger();  
  calibTrigger->SetStreamLevel(20);  
  calibTrigger->SetTriggerMask(-1,-1,kFALSE);       //accept everything 
  myTask->AddJob(calibTrigger);
}

void AddCalibCosmic(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // Histogram residuals and pulls of the track parameters in bins of track parameters           
  // 
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibCosmic *calibCosmic = new AliTPCcalibCosmic("cosmicTPC","cosmicTPC");
  calibCosmic->SetDebugLevel(debugLevel);
  calibCosmic->SetStreamLevel(streamLevel);
  calibCosmic->SetTriggerMask(-1,-1,kTRUE);        //accept everything 
  myTask->AddJob(calibCosmic);
}

void AddCalibLaser(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // 
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibLaser *calibLaser = new AliTPCcalibLaser("laserTPC","laserTPC");
  calibLaser->SetDebugLevel(debugLevel);
  calibLaser->SetStreamLevel(streamLevel);
  calibLaser->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  myTask->AddJob(calibLaser);
}


void AddCalibAlign(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // 
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibAlign *calibAlign = new AliTPCcalibAlign("alignTPC","Alignment of the TPC sectors");
  calibAlign->SetDebugLevel(debugLevel);
  calibAlign->SetStreamLevel(streamLevel);
  calibAlign->SetTriggerMask(-1,-1,kTRUE);        //accept everything 
  myTask->AddJob(calibAlign);
}






void AddCalibPID(TObject* task){
  //
  // Responsible: Marian Ivanov, Alexander Kalweit
  // Description:
  // 
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibPID *calibPID06 = new AliTPCcalibPID("calibPID06","calibPID06");
  AliTPCcalibPID *calibPID08 = new AliTPCcalibPID("calibPID08","calibPID08");
  AliTPCcalibPID *calibPID10 = new AliTPCcalibPID("calibPID10","calibPID10");
  calibPID06->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibPID08->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibPID10->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibPID06->SetUpperTrunc(0.6);
  calibPID08->SetUpperTrunc(0.8);
  calibPID10->SetUpperTrunc(0.99);
  //
  calibPID06->SetUsePosNorm(2);
  calibPID08->SetUsePosNorm(2);
  calibPID10->SetUsePosNorm(2);
  //
  calibPID06->SetUseShapeNorm(kFALSE);
  calibPID08->SetUseShapeNorm(kFALSE);
  calibPID10->SetUseShapeNorm(kFALSE);
  //
  calibPID06->SetPadNorm(0);
  calibPID08->SetPadNorm(0);
  calibPID10->SetPadNorm(0);
  calibPID06->SetMIPvalue(50);
  calibPID08->SetMIPvalue(50);
  calibPID08->SetMIPvalue(50);
  //
  //
  //
  calibPID06->SetDebugLevel(debugLevel);
  calibPID06->SetStreamLevel(streamLevel);
  calibPID06->SetTriggerMask(-1,-1,kTRUE);        //accept everything 
  calibPID08->SetDebugLevel(debugLevel);
  calibPID08->SetStreamLevel(streamLevel);
  calibPID08->SetTriggerMask(-1,-1,kTRUE);        //accept everything 
  calibPID10->SetDebugLevel(debugLevel);
  calibPID10->SetStreamLevel(streamLevel);
  calibPID10->SetTriggerMask(-1,-1,kTRUE);        //accept everything 
  myTask->AddJob(calibPID06);
  myTask->AddJob(calibPID08);
  myTask->AddJob(calibPID10);
}




//
//

void SetupCalibTaskTrain1(TObject* task){
  //
  //
  //
  AliTPCClusterParam * clusterParam = AliTPCcalibDB::Instance()->GetClusterParam(); 
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AddCalibCalib(task);
  AddCalibTimeGain(task);
  AddCalibTime(task);
  AddCalibCosmic(task);
  //AddCalibTrigger(task);
  //
  TString path=gSystem->pwd();
  path+=prefix;
  gSystem->mkdir(path);
  myTask->SetDebugOuputhPath(path.Data());

}

void SetupCalibTaskTrain2(TObject* task){
  //
  //
  //
  AliTPCClusterParam * clusterParam = AliTPCcalibDB::Instance()->GetClusterParam(); 
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  // AddCalibCalib(task);
  //AddCalibAlign(task);
  AddCalibLaser(task);
  // AddCalibTracks()
  AddCalibPID(task);
  //
  TString path=gSystem->pwd();
  path+=prefix;
  gSystem->mkdir(path);
  myTask->SetDebugOuputhPath(path.Data());
}





//
// backup of old setups - to be removed soon
//




void SetupCalibTask(TObject* task1){
  //
  // Configure calibration task
  //
  AliTPCAnalysisTaskcalib* myTask =   (AliTPCAnalysisTaskcalib*)task1;
  AliTPCClusterParam * clusterParam = AliTPCcalibDB::Instance()->GetClusterParam();  
  AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(30, 0.4, 5, 0.13, 0.018);  
  //
  AliTPCcalibTracks *calibTracks =  new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts); 


  AliTPCcalibTracksGain *calibTracksGain =  new AliTPCcalibTracksGain("calibTracksGain","Gain calibration using tracks",cuts); 
  AliTPCcalibAlign *calibAlign = new AliTPCcalibAlign("alignTPC","Alignment of the TPC sectors");
  AliTPCcalibAlign *calibAlignAll = new AliTPCcalibAlign("alignTPCAll","Alignment of the TPC sectors- All");
  AliTPCcalibLaser *calibLaser = new AliTPCcalibLaser("laserTPC","laserTPC");
  AliTPCcalibCosmic *calibCosmic = new AliTPCcalibCosmic("cosmicTPC","cosmicTPC");
  AliTPCcalibTrigger *calibTrigger = new AliTPCcalibTrigger();
  AliTPCcalibCalib *calibCalib = new AliTPCcalibCalib("calibTPC","calibTPC");
  AliTPCcalibTime *calibTime = new AliTPCcalibTime("calibTime","calibTime",  startTime.GetSec(), stopTime.GetSec(), 20*60);

  AliTPCcalibUnlinearity *calibUnlinearity = new AliTPCcalibUnlinearity("calibUnlinearity","calibUnlinearity");
  AliTPCcalibUnlinearity *calibUnlinearityAll = new AliTPCcalibUnlinearity("calibUnlinearityAll","calibUnlinearityAll");
  //
  //
  AliTPCcalibPID *calibPID06 = new AliTPCcalibPID("calibPID06","calibPID06");
  AliTPCcalibPID *calibPID08 = new AliTPCcalibPID("calibPID08","calibPID08");
  AliTPCcalibPID *calibPID10 = new AliTPCcalibPID("calibPID10","calibPID10");
  calibPID06->SetUpperTrunc(0.6);
  calibPID08->SetUpperTrunc(0.8);
  calibPID10->SetUpperTrunc(0.99);
  //
  calibPID06->SetUsePosNorm(2);
  calibPID08->SetUsePosNorm(2);
  calibPID10->SetUsePosNorm(2);
  //
  calibPID06->SetUseShapeNorm(kFALSE);
  calibPID08->SetUseShapeNorm(kFALSE);
  calibPID10->SetUseShapeNorm(kFALSE);
  //
  calibPID06->SetPadNorm(0);
  calibPID08->SetPadNorm(0);
  calibPID10->SetPadNorm(0);
  calibPID06->SetMIPvalue(50*3);
  calibPID08->SetMIPvalue(50*3);
  calibPID08->SetMIPvalue(50*3);
  //
  //
  //
  calibPID06->SetStreamLevel(2);
  calibPID06->SetDebugLevel(2);
  calibPID08->SetStreamLevel(2);
  calibPID08->SetDebugLevel(2);
  calibPID10->SetStreamLevel(2);
  calibPID10->SetDebugLevel(2);

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
  //calibTime->SetDebugLevel(20);
  calibTrigger->SetStreamLevel(20);
  //
  calibUnlinearity->SetDebugLevel(20);
  calibUnlinearity->SetStreamLevel(10);
  calibUnlinearityAll->SetDebugLevel(20);
  calibUnlinearityAll->SetStreamLevel(10);
  
  //  calibCalib->SetTriggerMask(-1,-1,kFALSE);       //accept everything 
  calibTrigger->SetTriggerMask(-1,-1,kFALSE);       //accept everything 
  calibCalib->SetTriggerMask(-1,-1,kTRUE);       //accept everything  - except laser
  calibTracks->SetTriggerMask(-1,16,kTRUE);      //reject laser trigger, accept everything else
  calibTracksGain->SetTriggerMask(-1,16,kTRUE);  //reject laser trigger, accept everything else
  calibAlign->SetTriggerMask(-1,-1,kTRUE);       //accept everything 
  calibAlignAll->SetTriggerMask(-1,-1,kFALSE);   //accept everything 
  calibLaser->SetTriggerMask(-1,-1,kFALSE);       //accept only laser trigger
  calibCosmic->SetTriggerMask(-1,-1,kTRUE);      //reject laser trigger, accept everything else
  //calibTime->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  calibUnlinearity->SetTriggerMask(-1,-1,kTRUE);   //reject laser 
  calibUnlinearityAll->SetTriggerMask(-1,-1,kFALSE);   //non reject laser 


  //
 // ---*---*-----*-*-----*----------*---
  // ADD CALIB JOBS HERE!!!!!!!!!!!!!!!!
  //myTask->AddJob(calibCalib);
  //myTask->AddJob(calibAlign);
  //myTask->AddJob(calibAlignAll);
  //myTask->AddJob(calibLaser);
  //myTask->AddJob(calibCosmic);
  myTask->AddJob(calibTrigger);
  //myTask->AddJob(calibTime);
  //myTask->AddJob(calibPID06);
  //myTask->AddJob(calibPID08);
  //myTask->AddJob(calibPID10);
  //  myTask->AddJob(calibTracksGain);
  //myTask->AddJob(calibTracks);
  //myTask->AddJob(calibUnlinearity);
  //myTask->AddJob(calibUnlinearityAll);
  // myTask->AddJob(new AliTPCcalibV0);
  // -*----*----*---*-*------*-------**--
  // -------*--*---------*-----*-------*-
  TString path=gSystem->pwd();
  path+=prefix;
  gSystem->mkdir(path);
  myTask->SetDebugOuputhPath(path.Data());


  //
  // backup important parameters
  //
  AliTPCClusterParam * paramCl = AliTPCcalibDB::Instance()->GetClusterParam(); 
  TFile f("ClusterParam.root","recreate");
  paramCl->Write();
  f.Close();
}




void CalibrateTPC(Int_t first, Int_t last, Int_t run, const char*closeSE="ALICE::GSI::SE"){
  gSystem->Load("$ROOTSYS/lib/libXrdClient.so");
  gSystem->Load("libNetx.so");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  //gEnv->SetValue("TFile.Recover", 0);  // dont try to recover anything
  gSystem->Setenv("alien_CLOSE_SE",closeSE);
  TGrid * alien =     TGrid::Connect("alien://",0,0,"t"); 
  gSystem->Setenv("alien_CLOSE_SE",closeSE);

  gSystem->Exec("touch nonOK");
  gSystem->Exec("rm -f isOK");
  gROOT->Macro(Form("ConfigOCDB.C\(%d\)",run));
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
  TChain * chain = tool.MakeChain("esd.txt","esdTree",0,last-first,first);
  chain->Lookup();
  // memory
  mgr->SetNSysInfo(100); 
  //
  mgr->SetDebugLevel(0);
  mgr->StartAnalysis("local",chain);
  gSystem->Exec("touch isOK");
  gSystem->Exec("rm -f nonOK");
  printf("Kill ourself\n");
  printf("Kill ourself\n");
  printf("Kill ourself\n");
  printf("Kill ourself\n");
  gSystem->Exec(Form("kill -9 %d",gSystem->GetPid()));
}
