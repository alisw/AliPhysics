//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runHadEt(bool submit = false, bool data = false, bool PbPb = false) {
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
   gROOT->ProcessLine(".L AliAnalysisEtCuts.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisHadEtCorrections.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisEtCommon.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisHadEt.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisHadEtMonteCarlo.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisHadEtReconstructed.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisEtSelectionContainer.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisEtSelectionHandler.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisTaskTransverseEnergy.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisTaskHadEt.cxx+g");


  char *kTreeName = "esdTree" ;
  TChain * chain   = new TChain(kTreeName,"myESDTree") ;
  if(submit){      
    gSystem->Load("libNetx.so") ; 
    gSystem->Load("libgapiUI.so");
    gSystem->Load("libRAliEn.so"); 
    TGrid::Connect("alien://") ;
  }
  //chain->Add("/data/LHC10h12/999/AliESDs.root");//Hijing Pb+Pb
  chain->Add("/data/LHC10d15/1821/AliESDs.root");//simulation p+p
  //chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TotEtManager");
  if(submit){
    if(PbPb){
      gROOT->LoadMacro("CreateAlienHandlerHadEtSimPbPb.C");
      AliAnalysisGrid *alienHandler = CreateAlienHandlerHadEtSimPbPb();  
    }
    else{
      gROOT->LoadMacro("CreateAlienHandlerHadEtSim.C");
      AliAnalysisGrid *alienHandler = CreateAlienHandlerHadEtSim();  
    }
    if (!alienHandler) return;
    mgr->SetGridHandler(alienHandler);
  }

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  AliMCEventHandler* handler = new AliMCEventHandler;
  if(!data){
    handler->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(handler);
  }

  bool isMc = true;
  if(data) isMC = false;
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");

  
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(isMc);
  if(PbPb){
   AliCentralitySelectionTask *centTask = AddTaskCentrality();
   mgr->AddTask(centTask);
  }

  if(data){
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("out1", TList::Class(), AliAnalysisManager::kOutputContainer,"event_stat.root");
    mgr->ConnectInput(physSelTask,0,cinput1);
    mgr->ConnectOutput(physSelTask,1,coutput1);
  }


    AliAnalysisTaskHadEt *task2 = new AliAnalysisTaskHadEt("TaskHadEt");
    mgr->AddTask(task2);
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("out2", TList::Class(), AliAnalysisManager::kOutputContainer,"Et.ESD.new.sim.root");
  mgr->ConnectInput(task2,0,cinput1);
  mgr->ConnectOutput(task2,1,coutput2);
  
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
