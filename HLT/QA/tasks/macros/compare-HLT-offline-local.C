// $Id$
/*
 * Example macro to run locally an analysis task for comparing the offline
 * with the HLT esd tree.
 *
 * The output is a root file containing the histograms defined in the
 * analysis task. There is one output file per detector.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q -l compare_HLT_offline_local.C'("phos")' 2>&1 | tee task.log
 *   aliroot -b -q -l compare_HLT_offline_local.C'("phos tpc")' 2>&1 | tee task.log
 *   aliroot -b -q -l compare_HLT_offline_local.C 2>&1 | tee task.log
 * </pre>
 *
 * If no argument is specified, ALL detector tasks are run.
 *
 * @ingroup alihlt_tpc
 * @author zbyin@mail.ccnu.edu.cn, Kalliopi.Kanaki@ift.uib.no
 */

void compare_HLT_offline_local(const char* detectorTask="all"){
 
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

  
  Bool_t bAll=kFALSE, bTPC=kFALSE, bPHOS=kFALSE, bITS=kFALSE, bGLOBAL=kFALSE;
 
  TString allArgs = detectorTask;
  TString argument;
 
  TObjArray *pTokens = allArgs.Tokenize(" ");
  if(pTokens){
     for(int i=0; i<pTokens->GetEntries(); i++){
         argument=((TObjString*)pTokens->At(i))->GetString();
         if(argument.IsNull()) continue;

         if(argument.CompareTo("tpc", TString::kIgnoreCase)==0){
	    bTPC = kTRUE;
	    continue;
         }        
         if(argument.CompareTo("phos", TString::kIgnoreCase)==0){
  	    bPHOS = kTRUE;
	    continue;
         }         
	 if(argument.CompareTo("its", TString::kIgnoreCase)==0){
  	    bITS = kTRUE;
	    continue;
         }	
	 if(argument.CompareTo("global", TString::kIgnoreCase)==0){
  	    bGLOBAL = kTRUE;
	    continue;
         }        
	 if(argument.CompareTo("all",TString::kIgnoreCase)==0){
	    bTPC    = kTRUE;
	    bPHOS   = kTRUE;
	    bITS    = kTRUE;
	    bGLOBAL = kTRUE;
	    bAll    = kTRUE;
	    continue;
         }
         else break;
    }
  }
    
  
  //-------------- Compile the analysis tasks ---------- //
  if(bTPC)    gROOT->LoadMacro("AliAnalysisTaskHLTTPC.cxx+"); 
  if(bPHOS)   gROOT->LoadMacro("AliAnalysisTaskHLTPHOS.cxx+"); 
  if(bITS)    gROOT->LoadMacro("AliAnalysisTaskHLTITS.cxx+");
  if(bGLOBAL) gROOT->LoadMacro("AliAnalysisTaskHLT.cxx+");

  
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
  
  if(bTPC){ 
     AliAnalysisTaskHLTTPC *taskTPC = new AliAnalysisTaskHLTTPC("offhlt_comparison_TPC");
     mgr->AddTask(taskTPC);
     AliAnalysisDataContainer *coutput1 =  mgr->CreateContainer("tpc_histograms", TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-TPC-comparison.root");  
     mgr->ConnectInput(taskTPC,0,mgr->GetCommonInputContainer());
     //mgr->ConnectOutput (taskTPC, 0, mgr->GetCommonOutputContainer());
     mgr->ConnectOutput(taskTPC,1,coutput1);
  }

  if(bPHOS){
     AliAnalysisTaskHLTPHOS *taskPHOS = new AliAnalysisTaskHLTPHOS("offhlt_comparison_PHOS");
     mgr->AddTask(taskPHOS);
     AliAnalysisDataContainer *coutput2 =  mgr->CreateContainer("phos_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-PHOS-comparison.root");  
     mgr->ConnectInput(taskPHOS,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskPHOS,1,coutput2);
  }
  
  if(bITS){
     AliAnalysisTaskHLTITS *taskITS = new AliAnalysisTaskHLTITS("offhlt_comparison_ITS");
     mgr->AddTask(taskITS);
     AliAnalysisDataContainer *coutput3 =  mgr->CreateContainer("its_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-ITS-comparison.root");  
     mgr->ConnectInput(taskITS,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskITS,1,coutput3);
  }
 
  if(bGLOBAL){
     AliAnalysisTaskHLT *taskGLOBAL = new AliAnalysisTaskHLT("offhlt_comparison_GLOBAL");
     mgr->AddTask(taskGLOBAL);
     AliAnalysisDataContainer *coutput4 =  mgr->CreateContainer("global_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-GLOBAL-comparison.root");  
     mgr->ConnectInput(taskGLOBAL,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskGLOBAL,1,coutput4);
  }
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);

  timer.Stop();
  timer.Print();
}
