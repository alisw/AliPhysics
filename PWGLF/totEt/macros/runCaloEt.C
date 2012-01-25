//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runCaloEt(bool submit = false, // true or false 
	       const char *dataType="real", // "sim" or "real" etc.
	       const char *pluginRunMode="test", // "test" or "full" or "terminate"
	       const char *det = "EMCalDetail") // "PHOS" or "EMCAL"
{
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");

  gSystem->Load("libMinuit");

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I. -I$ALICE_ROOT/EMCAL -I$ALICE_ROOT/ANALYSIS");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");



  if (!submit) { 
    cout << "local - no submitting" << endl;
  }
  else { 
    cout << "submitting to grid" << endl;
  }
   
  gROOT->ProcessLine(".L AliAnalysisEtCuts.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisHadEtCorrections.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtCommon.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEt.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarlo.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarloPhos.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtMonteCarloEmcal.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructed.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructedPhos.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtReconstructedEmcal.cxx+g");  
  gROOT->ProcessLine(".L AliAnalysisEtSelectionContainer.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisEtSelectionHandler.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisTaskTransverseEnergy.cxx+g");
gROOT->ProcessLine(".L AliAnalysisEmEtMonteCarlo.cxx+g");
gROOT->ProcessLine(".L AliAnalysisEmEtReconstructed.cxx+g");
  gROOT->ProcessLine(".L AliAnalysisTaskTotEt.cxx+g");

  gInterpreter->GenerateDictionary("std::map<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;map")  ;
  gInterpreter->GenerateDictionary("std::pair<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;utility");

  char *kTreeName = "esdTree" ;
  TChain * chain   = new TChain(kTreeName,"myESDTree") ;
  
  if(submit){      
    gSystem->Load("libNetx") ; 
    gSystem->Load("libgapiUI");
    gSystem->Load("libRAliEn"); 
    TGrid::Connect("alien://") ;
  }
  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TotEtManager");
  
  TString detStr(det);
  TString taskName = "TaskTotEt" + detStr;
  TString dataStr(dataType);
  TString dataStrName(dataType);
  dataStrName.ReplaceAll("/",".");
  TString outputName = "Et.ESD." + dataStrName + "." + detStr + ".root";
  TString outputDir = "totEt" + dataStr;

  cout << " taskName " << taskName
       << " outputName " << outputName 
       << " outputDir (alien) " << outputDir << endl;

  if (submit) {
    gROOT->LoadMacro("CreateAlienHandlerCaloEtSim.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandlerCaloEtSim(outputDir, outputName, pluginRunMode);  
    if (!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler* handler = new AliMCEventHandler;
  Bool_t isMc = kTRUE;
  Bool_t isPb = kFALSE;
  if ( dataStr.Contains("PbPb") ) { isPb = kTRUE;}
  if ( dataStr.Contains("sim") ) {
    cout << " MC " << endl;
    if ( dataStr.Contains("PbPb") ) { // a la: simPbPb/LHC10e18a
      cout << " PbPb " << endl;
      TString fileLocation = "/home/dsilverm/data/E_T/" + dataStr + "/dir/AliESDs.root";
      cout << "fileLocation " << fileLocation.Data() << endl; 
      chain->Add(fileLocation.Data()); // link to local test file
    }
    else { // pp
      cout<<"adding pp simulation file"<<endl;
      chain->Add("/data/LHC10d15/1821/AliESDs.root");
      //chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
      //chain->Add("/home/dsilverm/data/E_T/sim/LHC10d1/117222/100/AliESDs.root"); // link to local test file
    }
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  else { // real data
    isMc = kFALSE;
    chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
    //chain->Add("/home/dsilverm/data/E_T/data/2010/LHC10b/000117222/ESDs/pass2/10000117222021.30/AliESDs.root"); // link to local test file
    cout << " not MC " << endl;
  }


  if(!isMc && detStr.Contains("EMC")){
    cout<<"You are running over EMCal data and using the tender supply"<<endl;
    gSystem->Load("libTENDER.so");
    gSystem->Load("libTENDERSupplies.so"); 
    gROOT->ProcessLine(".include $ALICE_ROOT/Tender/"); 
    gSystem->AddIncludePath("-I$ALICE_ROOT/ANALYSIS "); 
    //Tender Supplies
    gROOT->LoadMacro("CreateEMCALTender.C");
    AliAnalysisTaskSE *tender = CreateEMCALTender(kTRUE);
    mgr->AddTask(tender); 
  }

  if(isMc) cout<<"I am a MC"<<endl;
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  
  AliPhysicsSelectionTask *physicsSelectionTask = AddTaskPhysicsSelection(isMc);//isMC is true when processing monte carlo
  if(isPb){	 
    cout<<"Adding centrality selection task"<<endl;
    gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    gROOT->ProcessLine(".L AliCentralitySelectionTask.cxx++g");
    AliCentralitySelectionTask *centTask = AddTaskCentrality();
  }



  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();





  AliAnalysisTaskTotEt *task1 = new AliAnalysisTaskTotEt(taskName);
  task1->SetMcData(isMc);//necessary to tell the task to basically accept all MC events.
  mgr->AddTask(task1);

  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("out1", TList::Class(), AliAnalysisManager::kOutputContainer, outputName);
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,1,coutput1);


  
  mgr->SetDebugLevel(0);
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  if(submit){
    mgr->StartAnalysis("grid");
  }
  else{
    mgr->StartAnalysis("local",chain);
  }
  
  timer.Stop();
  timer.Print();
}
