// $Id$
/*
 * Example macro to run locally an analysis task for comparing the offline
 * with the HLT esd tree.
 *
 * Its output is a root file containing the histograms defined in the
 * analysis task.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q -l compare_HLT_offline_local.C 2>&1 | tee task.log
 * </pre>
 *
 * @ingroup alihlt_tpc
 * @author zbyin@mail.ccnu.edu.cn, Kalliopi.Kanaki@ift.uib.no
 */

void compare_HLT_offline_local(){
 
  TStopwatch timer;
  timer.Start();

  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
 
  //----------- Loading the required libraries ---------//

  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libHLTbase.so");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");


  //-------------- Compile the analysis task ---------- //
  
  gROOT->LoadMacro("AliAnalysisTaskHLTTPC.cxx+"); 

  
  AliTagAnalysis *TagAna = new AliTagAnalysis("ESD"); 
  //TagAna->ChainLocalTags("../Tags");

//   AliRunTagCuts *runCuts = new AliRunTagCuts();
//   AliLHCTagCuts *lhcCuts = new AliLHCTagCuts();
//   AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
//   AliEventTagCuts *evCuts = new AliEventTagCuts();
//   evCuts->SetMultiplicityRange(11,12);  
  

  TChain *chain = 0x0;
  //chain = TagAna->QueryTags(runCuts,lhcCuts,detCuts,evCuts);
  chain = new TChain("esdTree");
 
  //gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  //chain = CreateESDChain("esd_run84254.txt", 2);
  
  //chain->Add("/afs/.alihlt.cern.ch/public/rec/79876/AliESDs.root");
  //chain->Add("/opt/HLT-public/rec/83683/rec/09000083683000.10/AliESDs.root");
  
  //chain->Add("/opt/HLT-public/rec/82762/ESDs/pass1/09000082762002.10/AliESDs.root");
  //chain->Add("/opt/HLT-public/rec/82762/ESDs/pass1/09000082762004.10/AliESDs.root");
  //chain->Add("/opt/HLT-public/rec/82762/ESDs/pass1/09000082762005.10/AliESDs.root");
  //chain->Add("/opt/HLT-public/rec/82762/ESDs/pass1/09000082762006.10/AliESDs.root");
  
  //chain->Add("/opt/HLT-public/rec/84254/alien/09000084254009.200/AliESDs.root");
  //chain->Add("/opt/HLT-public/rec/84254/alien/09000084254009.40/AliESDs.root");
  //chain->Add("/opt/HLT-public/rec/84254/alien/09000084254009.10/AliESDs.root");

  chain->Add("~/7TeV/115322/10000115322040.110/AliESDs.root");
  //chain->SetBranchStatus("*Calo*",0);

 
  //-------- Make the analysis manager ---------------//
 
  AliAnalysisManager *mgr  = new AliAnalysisManager("TestManager");
  AliESDInputHandler *esdH = new AliESDInputHandler;
  esdH->SetReadHLT();
  mgr->SetInputEventHandler(esdH);  
  mgr->SetNSysInfo(1000);
 
  //-------------- define the tasks ------------//
  
  AliAnalysisTaskHLTTPC *task1 = new AliAnalysisTaskHLTTPC("offhlt_comparison");
  mgr->AddTask(task1);

  AliAnalysisDataContainer *coutput1 =  mgr->CreateContainer("histograms", TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-TPC-comparison.root");  
  
  mgr->ConnectInput(task1,0,mgr->GetCommonInputContainer());
  //mgr->ConnectOutput (task1, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput(task1,1,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}
