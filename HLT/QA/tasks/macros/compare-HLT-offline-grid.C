// $Id$
/*
 * Example macro to run an HLT analysis task for comparing the offline
 * with the HLT esd tree on the GRID.
 *
 * The output is a root file containing the histograms defined in the
 * analysis task. There is one output file per detector.
 *
 * Usage:
 * <pre>
 *   aliroot -q compare-HLT-offline-grid.C'("000115322","/alice/data/2010/LHC10b","ESDcomparison","output","full","global")'
 * </pre>
 * - run number
 * - GRID input directory, where you define in which LHC period the run number belongs to
 * - GRID working directory, where the .xml, .jdl and the task are uploaded (you have to create it yourself in advance)
 * - GRID output directory with respect to the working one, where the output files of the task are located (you have to create it yourself in advance)
 * - run in full mode, i.e. completely on the GRID with all the chunks of the run processed
 * - specify the analysis task you want to run
 *
 * @ingroup alihlt_qa
 * @author zbyin@mail.ccnu.edu.cn, Kalliopi.Kanaki@ift.uib.no
 */

void compare_HLT_offline_grid(TString runNumber, TString dataDir, TString gridWorkingDir, TString gridOutputDir, const char* mode = "full", const char* detectorTask="global"){
 
  TStopwatch timer;
  timer.Start();

  //gSystem->Load("libTree");
  gSystem->Load("libCore");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
 
  //----------- Loading the required libraries ---------//

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libHLTbase");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");

  
  Bool_t bAll=kFALSE, bTPC=kFALSE, bPHOS=kFALSE, bEMCAL=kFALSE, bITS=kFALSE, bGLOBAL=kFALSE, bD0=kFALSE;
 
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
         if(argument.CompareTo("emcal", TString::kIgnoreCase)==0){
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
	 if(argument.CompareTo("D0", TString::kIgnoreCase)==0){
	   bD0 = kTRUE;
	   continue;
	 }  
	 if(argument.CompareTo("all",TString::kIgnoreCase)==0){
	    bTPC    = kTRUE;
	    bPHOS   = kTRUE;
	    bEMCAL   = kTRUE;
	    bITS    = kTRUE;
	    bGLOBAL = kTRUE;
	    bAll    = kTRUE;
	    continue;
         }
         else break;
    }
  }
    
  //-------- Make the analysis manager ---------------//
 
  AliAnalysisManager *mgr  = new AliAnalysisManager("TestManager");
  AliESDInputHandler *esdH = new AliESDInputHandler;
  esdH->SetReadHLT();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);  
  mgr->SetNSysInfo(1000);

  //To use Physics Selection
  AliPhysicsSelectionTask* physSelTask =AddTaskPhysicsSelection(kFALSE,kTRUE);
  
  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(runNumber, dataDir, gridWorkingDir, gridOutputDir, mode, detectorTask);
  if (!alienHandler) return;
  
  // Connect plugin to the analysis manager
  mgr->SetGridHandler(alienHandler);
 
  //-------------- Compile the analysis tasks ---------- //
  if(bTPC)    gROOT->LoadMacro("AliAnalysisTaskHLTTPC.cxx+"); 
  if(bPHOS && bEMCAL) {
    AliHLTSystem * pHLT = AliHLTPluginBase::GetInstance();
    pHLT->LoadComponentLibraries("libHLTbase");
    pHLT->LoadComponentLibraries("libAliHLTUtil");
    pHLT->LoadComponentLibraries("libAliHLTGlobal");
    gROOT->LoadMacro("AliAnalysisTaskHLTCalo.cxx+"); 
    gROOT->LoadMacro("AliAnalysisTaskHLTPHOS.cxx+");  
    gROOT->LoadMacro("AliAnalysisTaskHLTEMCAL.cxx+");  
    
  }
  else if(bPHOS) {
    AliHLTSystem * pHLT = AliHLTPluginBase::GetInstance();
    pHLT->LoadComponentLibraries("libHLTbase");
    pHLT->LoadComponentLibraries("libAliHLTUtil");
    pHLT->LoadComponentLibraries("libAliHLTGlobal");
    gROOT->LoadMacro("AliAnalysisTaskHLTCalo.cxx+"); 
    gROOT->LoadMacro("AliAnalysisTaskHLTPHOS.cxx+");  
  }
  else if(bEMCAL) {
     AliHLTSystem * pHLT = AliHLTPluginBase::GetInstance();
    pHLT->LoadComponentLibraries("libHLTbase");
    pHLT->LoadComponentLibraries("libAliHLTUtil");
    pHLT->LoadComponentLibraries("libAliHLTGlobal");
    gROOT->LoadMacro("AliAnalysisTaskHLTCalo.cxx+"); 
    gROOT->LoadMacro("AliAnalysisTaskHLTEMCAL.cxx+");  
  }
  if(bITS)    gROOT->LoadMacro("AliAnalysisTaskHLTITS.cxx+");
  if(bGLOBAL) gROOT->LoadMacro("AliAnalysisTaskHLT.cxx+");
  if(bD0)     gROOT->LoadMacro("AliAnalysisTaskD0Trigger.cxx+"); 
   
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
     taskGLOBAL->SelectCollisionCandidates();
     mgr->AddTask(taskGLOBAL);
     AliAnalysisDataContainer *coutput4 =  mgr->CreateContainer("global_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-GLOBAL-comparison.root");  
     mgr->ConnectInput(taskGLOBAL,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskGLOBAL,1,coutput4);
  }
  if(bD0){
    float cuts[7]={0.5,0.04,0.7,0.8,0.05,-0.00025,0.7};
    AliAnalysisTaskD0Trigger *taskD0 = new AliAnalysisTaskD0Trigger("offhlt_comparison_D0_Trigger",cuts);
    mgr->AddTask(taskD0);
    AliAnalysisDataContainer *coutput6 =  mgr->CreateContainer("D0_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-D0-comparison.root");  
    mgr->ConnectInput(taskD0,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskD0,1,coutput6);
  }
  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid");

  timer.Stop();
  timer.Print();
}
