//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runCaloEt(bool submit = false) {
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

    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gROOT->ProcessLine(".L AliAnalysisEt.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtMonteCarlo.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtMonteCarloPhos.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtMonteCarloEmcal.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtReconstructed.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtReconstructedPhos.cxx+g");
    gROOT->ProcessLine(".L AliAnalysisEtReconstructedEmcal.cxx+g");

    gROOT->ProcessLine(".L AliAnalysisTaskTotEt.cxx+g");


  char *kTreeName = "esdTree" ;
  TChain * chain   = new TChain(kTreeName,"myESDTree") ;

  
  if(submit){      
    gSystem->Load("libNetx.so") ; 
    gSystem->Load("libgapiUI.so");
    gSystem->Load("libRAliEn.so"); 
    TGrid::Connect("alien://") ;
  }
  chain->Add("/data/LHC10d15/1821/AliESDs.root");//CN changed


  
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TotEtManager");
  
  if(submit){
    gROOT->LoadMacro("CreateAlienHandlerPhosEtSim.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandlerPhosEtSim();  
    if (!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(handler);
  
  AliAnalysisTaskTotEt *task1 = new AliAnalysisTaskTotEt("TaskTotEt");
  mgr->AddTask(task1);

    // Create containers for input/output
    //AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1", TChain::Class(),AliAnalysisManager::kInputContainer);
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

    //AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("out1", TList::Class(), AliAnalysisManager::kOutputContainer,"Et.ESD.sim.root");
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("out1", TList::Class(), AliAnalysisManager::kOutputContainer,"Et.ESD.sim.root");

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
