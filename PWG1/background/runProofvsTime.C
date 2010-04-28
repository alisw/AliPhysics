
void runProofvsTime(const char * dataset = "LHC09b12_7TeV_0.5T", TString dataSetPath ="/PWG0/jgrosseo/",const char * filename = "LHC09b12_7TeV_0.5T_bg.root",Int_t start=0, Int_t stop =1000, Int_t cuts = 0, Bool_t mc = kFALSE, Bool_t useBINT = kFALSE, Int_t nev =123456789) {
//void runProofvsTime(const char * dataset = "LHC09b12_7TeV_0.5T", TString dataSetPath ="/PWG0/jgrosseo/",const char * filename = "LHC09b12_7TeV_0.5T_bg.root",Int_t start=0, Int_t stop =1000, Int_t cuts = 0, Bool_t mc = kFALSE, Bool_t useBINT = kFALSE, Int_t nev =10000) {
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  //  TProof::Mgr("alicecaf")->SetROOTVersion("v5-24-00a_dbg");
  //TProof::Open("alicecaf", "valgrind=workers#4");
  //  TProof::Open("alicecaf");
  TProof::Open("mfloris@alicecaf.cern.ch");
  
  //  gSystem->AddIncludePath("-I${ALICE_ROOT}/include/ -I${ALICE_ROOT}/PWG0/ -I${ALICE_ROOT}/PWG0/dNdEta/");
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include/");
  gProof->UploadPackage("$ALICE_ROOT/STEERBase");
  gProof->EnablePackage("$ALICE_ROOT/STEERBase");
  gProof->UploadPackage("$ALICE_ROOT/ESD");
  gProof->EnablePackage("$ALICE_ROOT/ESD");
  gProof->UploadPackage("$ALICE_ROOT/AOD");
  gProof->EnablePackage("$ALICE_ROOT/AOD");
  gProof->UploadPackage("$ALICE_ROOT/ANALYSIS");
  gProof->EnablePackage("$ALICE_ROOT/ANALYSIS");
  gProof->UploadPackage("$ALICE_ROOT/ANALYSISalice");
  gProof->EnablePackage("$ALICE_ROOT/ANALYSISalice");
  gProof->UploadPackage("$ALICE_ROOT/CORRFW");
  gProof->EnablePackage("$ALICE_ROOT/CORRFW");
 

    // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  //  mgr->SetDebugLevel(3);
  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler; 

  mgr->SetInputEventHandler(esdH);
  
  if(mc) {
    AliMCEventHandler *mch = new AliMCEventHandler();
    mch->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mch);
  }
  // assign simple task
  gProof->Load(gSystem->ExpandPathName("$(ALICE_ROOT)/PWG1/background/AliHistoListWrapper.cxx++g"));   
  gProof->Load(gSystem->ExpandPathName("$(ALICE_ROOT)/PWG1/background/AliAnalysisTaskBGvsTime.cxx++g"));   
  //____________________________________________//
  // assign simple task
  AliAnalysisTaskBGvsTime * task = new AliAnalysisTaskBGvsTime("TaskBG");
  //  task->SetMC();
  const Int_t mult_bins[] = {0,1000};
  task->SetMultBins(2,mult_bins);
  //  task->SetBinWidth(2); // binw in secs
  task->SetBinWidth(300); // binw in secs
  task->SetTimes(start,stop);

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


	
  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",dataSetPath+dataset+"#esdTree",nev);

  if( cuts == 0 ){ 
    cout << "WARNING: disabled cuts" << endl;
  } else if (cuts == 2) {
    cout << "INFO: Using Physics Selection" << endl;
  }

}
