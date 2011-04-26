//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runHadEt(bool submit = false, bool data = false, Int_t dataset = 20100, Bool_t test = kTRUE) {
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
  bool PbPb = false;
  if(dataset ==20100){
    bool PbPb = true;
    if(data){
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.860/AliESDs.root");//Data Pb+Pb
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.870/AliESDs.root");//Data Pb+Pb
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.880/AliESDs.root");//Data Pb+Pb
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.890/AliESDs.root");//Data Pb+Pb
      chain->Add("/data/LHC10h/pass2_rev15/10000137366041.900/AliESDs.root");//Data Pb+Pb
    }
    else{
      chain->Add("/data/LHC10h8/137161/999/AliESDs.root");//Hijing Pb+Pb
      chain->Add("/data/LHC10h8/137161/111/AliESDs.root");//Hijing Pb+Pb
      chain->Add("/data/LHC10h8/137161/222/AliESDs.root");//Hijing Pb+Pb
      //chain->Add("/data/LHC11a4_bis/137161/999/AliESDs.root");//Hijing Pb+Pb
      //chain->Add("/data/LHC10h12/999/AliESDs.root");//Hijing Pb+Pb
    }
  } 
  else{
    if(data){
      chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
    }
    else{
      chain->Add("/data/LHC10d15/1821/AliESDs.root");//simulation p+p
    }
  }

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TotEtManager");
  if(submit){
      gROOT->LoadMacro("CreateAlienHandlerHadEt.C");
      AliAnalysisGrid *alienHandler = CreateAlienHandlerHadEt(dataset,data,test);//integer dataset, boolean isData, bool submit-in-test-mode
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

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(!data);
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");

  AliCentralitySelectionTask *centTask;
  
  if(PbPb){
    AliCentralitySelectionTask *centTask = AddTaskCentrality();
   if(!data){
     cout<<"Setting up centrality for MC"<<endl;
     centTask->SetMCInput();
   }
   else{
     cout<<"Setting up centrality for data"<<endl;
   }
  }

   if(dataset==20100){//PbPb 2.76 TeV
     if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a4_bis.PbPb.ForData.root","corrections.root",kTRUE);}
     else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11a4_bis.PbPb.ForSimulations.root","corrections.root",kTRUE);}
   }
   else{
     if(dataset==2009){//pp 900 GeV
       cout<<"Warning! You are using 7 TeV corrections for 900 GeV data!"<<endl;}
     if(dataset==20111){//pp 2.76 TeV
       cout<<"Warning! You are using 7 TeV corrections for 2.76eV data!"<<endl;}
     if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC10d4.pp.ForData.root","corrections.root",kTRUE);}
     else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC10d4.pp.ForSimulations.root","corrections.root",kTRUE);}
   }
   TString recoFile;
   TString mcFile;
   if(dataset==20100){
     recoFile = "ConfigHadEtReconstructedPbPb.C";
     mcFile = "ConfigHadEtMonteCarloPbPb.C";
   }
   if(dataset==2009){
     recoFile = "ConfigHadEtReconstructedpp900GeV.C";
     mcFile = "ConfigHadEtMonteCarlopp900GeV.C";
   }
   if(dataset==20111){
     recoFile = "ConfigHadEtReconstructedpp276TeV.C";
     mcFile = "ConfigHadEtMonteCarlopp276TeV.C";
   }
   if(dataset==2010){
     recoFile = "ConfigHadEtReconstructedpp7TeV.C";
     mcFile = "ConfigHadEtMonteCarlopp7TeV.C";
   }
   AliAnalysisTaskHadEt *task2 = new AliAnalysisTaskHadEt("TaskHadEt",!data,recoFile,mcFile);
    if(!data) task2->SetMcData();
    mgr->AddTask(task2);
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
