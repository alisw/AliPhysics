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
 *   aliroot -b -q -l compare_HLT_offline_local.C'("/home/blabla/AliESDs.root","phos")' 2>&1 | tee task.log
 *   aliroot -b -q -l compare_HLT_offline_local.C'("/home/blabla/AliESDs.root","phos tpc global")' 2>&1 | tee task.log
 *   aliroot -q compare-HLT-offline-local.C'("alien:///alice/data/2010/LHC10b/000115322/ESDs/pass1/10000115322040.20/AliESDs.root","global")' 2>&1 | tee log
 * </pre>
 * 
 * If alien:// is contained in the name of the file, then the macro connects to the grid to access the file.
 * 
 * In case you want to run over many ESD files, then prepare a list of them in a .txt file and they will be chained for the analysis.
 * The .txt file takes the place of the first argument in that case.
 *
 * @ingroup alihlt_qa
 * @author zbyin@mail.ccnu.edu.cn, Kalliopi.Kanaki@ift.uib.no
 */

void compare_HLT_offline_local(TString file, const char* detectorTask="global"){

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

  
  Bool_t bTPC=kFALSE, bPHOS=kFALSE, bITS=kFALSE, bGLOBAL=kFALSE, bEMCAL = kFALSE;
 
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
         else if(argument.CompareTo("emcal", TString::kIgnoreCase)==0){
	   bEMCAL = kTRUE;
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
	    bEMCAL  = kTRUE;
	    bITS    = kTRUE;
	    bGLOBAL = kTRUE;    
	    continue;
         }
         else break;
    }
  }
    
  
  //-------------- Compile the analysis tasks ---------- //
  if(bTPC) gROOT->LoadMacro("AliAnalysisTaskHLTTPC.cxx+"); 
  
  if(bPHOS) {
    AliHLTSystem * pHLT = AliHLTPluginBase::GetInstance();
    pHLT->LoadComponentLibraries("libHLTbase");
    pHLT->LoadComponentLibraries("libAliHLTUtil");
    pHLT->LoadComponentLibraries("libAliHLTGlobal");
    gROOT->LoadMacro("AliAnalysisTaskHLTCalo.cxx+"); 
    gROOT->LoadMacro("AliAnalysisTaskHLTPHOS.cxx+"); 
  }
  
  if(bEMCAL) {
    AliHLTSystem * pHLT = AliHLTPluginBase::GetInstance();
    pHLT->LoadComponentLibraries("libHLTbase");
    pHLT->LoadComponentLibraries("libAliHLTUtil");
    pHLT->LoadComponentLibraries("libAliHLTGlobal");
    gROOT->LoadMacro("AliAnalysisTaskHLTCalo.cxx+"); 
    gROOT->LoadMacro("AliAnalysisTaskHLTEMCAL.cxx+"); 
  }  
  
  if(bITS)    gROOT->LoadMacro("AliAnalysisTaskHLTITS.cxx+");
  if(bGLOBAL) gROOT->LoadMacro("AliAnalysisTaskHLT.cxx+");
  
  //if(!AliAnalysisGrid::CreateToken()) return NULL;
  
  if(file.Contains("alien")) TGrid::Connect("alien://");
    
  if(file.Contains("AliESDs.root")){
    TChain *chain = new TChain("esdTree"); 
    chain->Add(file);
  }
  
  //Constructs chain from filenames in *.txt
  //on the form $DIR/AliESDs.root
  else if(file.Contains(".txt")){
    gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
    chain=CreateESDChain(file.Data());
  }

  if(!chain){
    Printf("Chain is empty");
    return;
  }


   
  //-------- Make the analysis manager ---------------//
 
  AliAnalysisManager *mgr  = new AliAnalysisManager("TestManager");
  AliESDInputHandler *esdH = new AliESDInputHandler;
  esdH->SetReadHLT();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);  
  mgr->SetNSysInfo(1000);
 
  //-------------- define the tasks ------------//
  
  if(bTPC){ 
     AliAnalysisTaskHLTTPC *taskTPC = new AliAnalysisTaskHLTTPC("offhlt_comparison_TPC");
     mgr->AddTask(taskTPC);
     AliAnalysisDataContainer *coutput1 =  mgr->CreateContainer("tpc_histograms", TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-TPC-comparison.root");  
     mgr->ConnectInput(taskTPC,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskTPC,1,coutput1);
  }

  if(bPHOS){
     AliAnalysisTaskHLTPHOS *taskPHOS = new AliAnalysisTaskHLTPHOS("offhlt_comparison_PHOS");
     mgr->AddTask(taskPHOS);
     AliAnalysisDataContainer *coutput2 =  mgr->CreateContainer("phos_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-PHOS-comparison.root");  
     mgr->ConnectInput(taskPHOS,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskPHOS,1,coutput2);
  }

  if(bEMCAL){
     AliAnalysisTaskHLTEMCAL *taskEMCAL = new AliAnalysisTaskHLTEMCAL("offhlt_comparison_EMCAL");
     mgr->AddTask(taskEMCAL);
     AliAnalysisDataContainer *coutput5 =  mgr->CreateContainer("emcal_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-EMCAL-comparison.root");  
     mgr->ConnectInput(taskEMCAL,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskEMCAL,1,coutput5);
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
  mgr->StartAnalysis("local",chain,500);

  timer.Stop();
  timer.Print();
}
