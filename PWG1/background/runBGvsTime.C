void Load() ;

void runBGvsTime(const char * incollection,const char *filename,UInt_t start, UInt_t end, Int_t cuts = 0, Bool_t mc = kFALSE, Bool_t useBINT = kFALSE, Int_t nev=123456789) {

  //  gDebug = 3;
  Load();

  AliLog::SetGlobalLogLevel(AliLog::kInfo);

  // Connect to the grid and create chain
  //  TGrid::Connect("alien://");
  TChain* analysisChain = 0;
  analysisChain = new TChain("esdTree");
  if (TString(incollection).Contains(".root")){
    analysisChain->Add(incollection);
  }
  else if (TString(incollection).Contains("xml")){
    TGrid::Connect("alien://");
    TAlienCollection * coll = TAlienCollection::Open (incollection);
    while(coll->Next()){
      analysisChain->Add(TString("alien://")+coll->GetLFN());
    }
  } else {
    ifstream file_collect(incollection);
    TString line;
    while (line.ReadLine(file_collect) ) {
      analysisChain->Add(line.Data());
    }
  }
  analysisChain->GetListOfFiles()->Print();

  //
  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  //  mgr->SetDebugLevel(3);
  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler; 
  //  AliESDInputHandler* esdH = new AliESDInputHandlerRP; // for RecPoints
  
//   esdH->SetInactiveBranches("AliESDACORDE FMD ALIESDTZERO ALIESDVZERO ALIESDZDC AliRawDataErrorLogs CaloClusters Cascades EMCALCells EMCALTrigger ESDfriend Kinks Kinks Cascades AliESDTZERO ALIESDACORDE MuonTracks TrdTracks CaloClusters");
  mgr->SetInputEventHandler(esdH);
	
//   AliMCEventHandler *mc = new AliMCEventHandler();
//   mc->SetReadTR(kFALSE);
//   mgr->SetMCtruthEventHandler(mc);
	
  if(mc) {
    AliMCEventHandler *mch = new AliMCEventHandler();
    mch->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mch);
  }
  // assign simple task
  //____________________________________________//
  // assign simple task
  AliAnalysisTaskBGvsTime * task = new AliAnalysisTaskBGvsTime("TaskBGvsTime");
  //  task->SetMC();
//   const Int_t mult_bins[] = {0,10,50,1000};
//   task->SetMultBins(4,mult_bins);
  const Int_t mult_bins[] = {0,1000};
  task->SetMultBins(2,mult_bins);
  task->SetTimes(start,end);
  task->SetBinWidth(300); // binw in secs

  if(mc) task->SetMC();

  if (cuts == 0 )       task->SetNoCuts();
  else if (cuts == 2 )  task->SetUsePhysicsSelection();
  else if (cuts == 3 )  {task->SetUsePhysicsSelection();task->SetUseZeroBin();}
  else if (cuts == 4 )  {task->SetUsePhysicsSelection();task->SetSkipV0();task->SetUseZeroBin();}
  else if (cuts == 5 )  {task->SetUsePhysicsSelection();task->SetSkipV0();task->SetSkipZeroBin();}

  if(useBINT) task->SetUseBI();

  mgr->AddTask(task);

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();	
  mgr->ConnectInput(task,0,cinput1);

  // Attach output
  cOutput = mgr->CreateContainer("BGvsTime", 
				 AliHistoListWrapper::Class(), AliAnalysisManager::kOutputContainer,filename);
  mgr->ConnectOutput(task, 1, cOutput);      
  cOutput = mgr->CreateContainer("PhysSel", 
				 AliPhysicsSelection::Class(), AliAnalysisManager::kOutputContainer,filename);
  mgr->ConnectOutput(task, 2, cOutput);      

	
  //____________________________________________//

  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  mgr->StartAnalysis("local",analysisChain,nev);




}

void Load() {

  //load the required aliroot libraries
  gSystem->Load("libANALYSIS") ;
  gSystem->Load("libANALYSISalice") ;
  gSystem->Load("libCORRFW") ;
  gSystem->Load("libITSbase") ;
  gSystem->Load("libPWG0base") ;

  //compile online the task class
  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_ROOT/PWG0/ -I$ROOTSYS/include");
  //  gROOT->LoadMacro("./AliAnalysisTaskCombPIDSpectra.cxx+");
  //  gROOT->LoadMacro("AliBackgroundSelection.cxx++g");
  gROOT->LoadMacro("AliHistoListWrapper.cxx++g");   
  gROOT->LoadMacro("AliAnalysisTaskBGvsTime.cxx++g");
}
