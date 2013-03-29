//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville

//by default this runs locally
//With the argument true this submits jobs to the grid
//As written this requires an xml script tag.xml in the ~/et directory on the grid to submit jobs
void runHadEt(bool submit = false, bool data = false, Int_t dataset = 20100, Bool_t test = kFALSE, Int_t material = 0, Bool_t altV0Scale = kFALSE, bool runCompiledVersion = kFALSE, int simflag = 0) {
    TStopwatch timer;
    timer.Start();
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libXMLIO.so");

    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");

    gSystem->Load("libANALYSIS");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libPWGUDbase.so");

    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I$ALICE_ROOT/PWGUD/base");
    if(runCompiledVersion){
      gSystem->Load("libPWGLFtotEt.so");
    }
    else{
      gROOT->ProcessLine(".L AliAnalysisEtCuts.cxx+g");
      gROOT->ProcessLine(".L AliAnalysisHadEtCorrections.cxx+g");
      gROOT->ProcessLine(".L AliAnalysisEtCommon.cxx+g");
      gROOT->ProcessLine(".L AliAnalysisHadEt.cxx+g");
      gROOT->ProcessLine(".L AliAnalysisHadEtMonteCarlo.cxx+g");
      gROOT->ProcessLine(".L AliAnalysisHadEtReconstructed.cxx+g");
      gROOT->ProcessLine(".L AliAnalysisTaskTransverseEnergy.cxx+g");
      gROOT->ProcessLine(".L AliAnalysisTaskHadEt.cxx+g");
    }

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
      cout<<"Yes I am analyzing the correct file"<<endl;
      //chain->Add("/data/LHC10dpass2/10000126403050.70/AliESDs.root");//data
      //chain->Add("/data/LHC11a/11000146856042.90/AliESDs.root");//data pp 2.76 TeV w/SDD pass 2
      chain->Add("/data/LHC11a/11000146860043.90/AliESDs.root");//data pp 2.76 TeV w/SDD pass 3
    }
    else{
      //chain->Add("/data/LHC10d15/1821/AliESDs.root");//simulation p+p
      chain->Add("/data/LHC11b10a/001/AliESDs.root");
    }
  }

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TotEtManager");
  if(submit){

    gROOT->LoadMacro("CreateAlienHandlerHadEt.C");
    AliAnalysisGrid *alienHandler = CreateAlienHandlerHadEt(dataset,data,test,material,altV0Scale,runCompiledVersion,simflag);//integer dataset, boolean isData, bool submit-in-test-mode, bool use alternatve V0 scaling
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


  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(!data);
  if(!physSelTask) { Printf("no physSelTask"); return; }
  AliPhysicsSelection *physSel = physSelTask->GetPhysicsSelection();
  //physSel->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL");// #3119 #769");
  //physSelTask->AddCollisionTriggerClass("kMB");// #3119 #769");


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
     gSystem->CopyFile("ConfigHadEtMonteCarloPbPb.C","ConfigHadEtMonteCarlo.C",kTRUE);
     gSystem->CopyFile("ConfigHadEtReconstructedPbPb.C","ConfigHadEtReconstructed.C",kTRUE);
     //centTask->SetPass(1);
   }
   else{
     if(dataset==2009){//pp 900 GeV
       gSystem->CopyFile("ConfigHadEtMonteCarlopp900GeV.C","ConfigHadEtMonteCarlo.C",kTRUE);
       gSystem->CopyFile("ConfigHadEtReconstructedpp900GeV.C","ConfigHadEtReconstructed.C",kTRUE);
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11b1a.pp.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11b1a.pp.ForSimulations.root","corrections.root",kTRUE);}
     }
     if(dataset==20111){//pp 2.76 TeV
       if(altV0Scale)gSystem->CopyFile("ConfigHadEtMonteCarlopp276TeVAlt.C","ConfigHadEtMonteCarlo.C",kTRUE);
       else{gSystem->CopyFile("ConfigHadEtMonteCarlopp276TeV.C","ConfigHadEtMonteCarlo.C",kTRUE);}
       gSystem->CopyFile("ConfigHadEtReconstructedpp276TeV.C","ConfigHadEtReconstructed.C",kTRUE);
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC11b10a.pp.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC11b10a.pp.ForSimulations.root","corrections.root",kTRUE);}
     }
     if(dataset==2010){//pp 7 TeV
       gSystem->CopyFile("ConfigHadEtMonteCarlopp7TeV.C","ConfigHadEtMonteCarlo.C",kTRUE);
       gSystem->CopyFile("ConfigHadEtReconstructedpp7TeV.C","ConfigHadEtReconstructed.C",kTRUE);
       if(data){gSystem->CopyFile("rootFiles/corrections/corrections.LHC10e20.pp.ForData.root","corrections.root",kTRUE);}
       else{gSystem->CopyFile("rootFiles/corrections/corrections.LHC10e20.pp.ForSimulations.root","corrections.root",kTRUE);}
     }
   }
   AliAnalysisTaskHadEt *task2 = new AliAnalysisTaskHadEt("TaskHadEt",!data);//,recoFile,mcFile);
   if(!data) task2->SetMcData();
   //Add thing here to select collision type!!
   //if(dataset!=20100){task2->SelectCollisionCandidates(AliVEvent::kMB ) ;}
   //if(dataset!=20100){task2->SelectCollisionCandidates(AliVEvent::kMB ) ;}
   mgr->AddTask(task2);
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
   AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("out2", TList::Class(), AliAnalysisManager::kOutputContainer,"Et.ESD.new.sim.root");
//    mgr->ConnectInput(task2,0,cinput1);
//    mgr->ConnectOutput(task2,1,coutput2);
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
