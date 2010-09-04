//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runCaloEt(bool submit = true, // true or false 
	       const char *dataType="sim", // "sim" or "real" or "simPbPb"
	       const char *det = "EMCAL") // "PHOS" or "EMCAL"
{
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  if (!submit) { // assume the stuff you need is compiled
    cout << "local - no submitting" << endl;
    gSystem->Load("libPWG4totEt");
  }
  else { 
    cout << "submitting to grid" << endl;
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gROOT->ProcessLine(".L AliAnalysisEt.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtMonteCarlo.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtMonteCarloPhos.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtMonteCarloEmcal.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtReconstructed.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtReconstructedPhos.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtReconstructedEmcal.cxx+g");
    
    gROOT->ProcessLine(".L AliAnalysisTaskTotEt.cxx+g");
  }

  char *kTreeName = "esdTree" ;
  TChain * chain   = new TChain(kTreeName,"myESDTree") ;
  
  if(submit){      
    gSystem->Load("libNetx.so") ; 
    gSystem->Load("libgapiUI.so");
    gSystem->Load("libRAliEn.so"); 
    TGrid::Connect("alien://") ;
  }
  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TotEtManager");
  
  TString detStr(det);
  TString taskName = "TaskTotEt" + detStr;
  TString dataStr(dataType);
  TString outputName = "Et.ESD." + dataStr + "." + detStr + ".root";
  TString outputDir = "totEt" + dataStr;

  cout << " taskName " << taskName
       << " outputName " << outputName << endl;

  if(submit){
    gROOT->LoadMacro("CreateAlienHandlerCaloEtSim.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandlerCaloEtSim(outputDir, outputName);  
    if (!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler* handler = new AliMCEventHandler;
  if ( dataStr.Contains("sim") ) {
    cout << " MC " << endl;
    chain->Add("/home/dsilverm/data/E_T/sim/LHC10d1/117222/100/AliESDs.root"); // link to local test file
    if ( dataStr.Contains("PbPb") ) {
      cout << " PbPb " << endl;
      chain->Add("/home/dsilverm/data/E_T/sim/LHC10e11/191001/001/AliESDs.root"); // link to local test file
    }
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }
  else {
  chain->Add("/home/dsilverm/data/E_T/data/2010/LHC10b/000117222/ESDs/pass2/10000117222021.30/AliESDs.root"); // link to local test file
    cout << " not MC " << endl;
  }


  AliAnalysisTaskTotEt *task1 = new AliAnalysisTaskTotEt(taskName);
  mgr->AddTask(task1);

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
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
