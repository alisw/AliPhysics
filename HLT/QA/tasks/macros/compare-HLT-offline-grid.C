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
 * - GRID working directory, where the .xml, .jdl and the task are uploaded 
 * - GRID output directory with respect to the working one, where the output files of the task are located 
 * - run in full mode, i.e. completely on the GRID with all the chunks of the run processed
 * - specify the analysis task you want to run
 * - specify the path where the task is located, by default it takes $ALICE_ROOT/HLT/QA/tasks
 * - specify whether you are interested only in HLT triggered events
 * - specify how many events you want to analyze
 *
 * @ingroup alihlt_qa
 * @author Hege.Erdal@student.uib.no, Kalliopi.Kanaki@ift.uib.no
 */

void compare_HLT_offline_grid(TString runNumber, 
                              TString dataDir, 
			      TString gridWorkingDir, 
			      TString gridOutputDir, 
			      const char* mode = "full", 
			      const char* detectorTask="global",
			      TString taskFolder="$ALICE_ROOT/HLT/QA/tasks/",
			      TString beamType="p-p",
			      bool fUseHLTTrigger=kFALSE,
			      Long64_t nEvents=1234567890
			     )
{
 
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
  gSystem->AddIncludePath("-I$ALICE_ROOT/HLT/BASE -I. -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS");
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");

  
  Bool_t bAll=kFALSE, bTPC=kFALSE, bPHOS=kFALSE, bEMCAL=kFALSE, bITS=kFALSE, bGLOBAL=kFALSE, bD0=kFALSE, bCB=kFALSE;
 
  TString allArgs = detectorTask;
  TString argument;
 
  TObjArray *pTokens = allArgs.Tokenize(" ");
  if(pTokens){
     for(int i=0; i<pTokens->GetEntries(); i++){
         argument=((TObjString*)pTokens->At(i))->GetString();
         if(argument.IsNull()) continue;

         if(argument.CompareTo("phos", TString::kIgnoreCase)==0){
  	    bPHOS = kTRUE;
	    continue;
         }
         if(argument.CompareTo("emcal", TString::kIgnoreCase)==0){
  	    bEMCAL = kTRUE;
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
	 if(argument.CompareTo("cb", TString::kIgnoreCase)==0){
	   bCB = kTRUE;
	   continue;
	 }  
	 if(argument.CompareTo("all",TString::kIgnoreCase)==0){
	    bPHOS   = kTRUE;
	    bEMCAL  = kTRUE;
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
  //AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kFALSE,kTRUE);

  // Create and configure the alien handler plugin
  gROOT->LoadMacro("$ALICE_ROOT/HLT/QA/tasks/macros/CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(runNumber, dataDir, gridWorkingDir, gridOutputDir, mode, detectorTask);
  if (!alienHandler) return;
  
  // Connect plugin to the analysis manager
  mgr->SetGridHandler(alienHandler);
 
  //-------------- Compile the analysis tasks ---------- //
  if(bPHOS && bEMCAL) {
    gSystem->Load("libHLTbase");
    gSystem->Load("libAliHLTUtil");
    gSystem->Load("libAliHLTGlobal");
    TString strTask1("AliAnalysisTaskHLTCalo.cxx+");
    TString strTask2("AliAnalysisTaskHLTPHOS.cxx+");    
    TString strTask3("AliAnalysisTaskHLTEMCAL.cxx+");    
    gROOT->LoadMacro(taskFolder+strTask1); 
    gROOT->LoadMacro(taskFolder+strTask2);  
    gROOT->LoadMacro(taskFolder+strTask3);      
  }
  else if(bPHOS) {
    gSystem->Load("libHLTbase");
    gSystem->Load("libAliHLTUtil");
    gSystem->Load("libAliHLTGlobal");
    TString strTask1("AliAnalysisTaskHLTCalo.cxx+");
    TString strTask2("AliAnalysisTaskHLTPHOS.cxx+");    
    gROOT->LoadMacro(taskFolder+strTask1); 
    gROOT->LoadMacro(taskFolder+strTask2);  
  }
  else if(bEMCAL) {
    gSystem->Load("libHLTbase");
    gSystem->Load("libAliHLTUtil");
    gSystem->Load("libAliHLTGlobal");
    TString strTask1("AliAnalysisTaskHLTCalo.cxx+");
    TString strTask2("AliAnalysisTaskHLTEMCAL.cxx+");    
    gROOT->LoadMacro(taskFolder+strTask1); 
    gROOT->LoadMacro(taskFolder+strTask2);  
  }
  if(bGLOBAL){
     TString strTask("AliAnalysisTaskHLT.cxx+");   
     gROOT->LoadMacro(taskFolder+strTask);
  }
  if(bD0){
    TString strTask("AliAnalysisTaskD0Trigger.cxx+");
    gROOT->LoadMacro(taskFolder+strTask);     
  }
  if(bCB){
    TString strTask("AliAnalysisTaskHLTCentralBarrel.cxx+");
    gROOT->LoadMacro(taskFolder+strTask);
  } 
 
  //-------------- define the tasks ------------//
  
  if(bPHOS){
     AliAnalysisTaskHLTPHOS *taskPHOS = new AliAnalysisTaskHLTPHOS("offhlt_comparison_PHOS");
     mgr->AddTask(taskPHOS);
     AliAnalysisDataContainer *coutputPHOS =  mgr->CreateContainer("phos_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-PHOS-comparison.root");  
     mgr->ConnectInput(taskPHOS,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskPHOS,1,coutputPHOS);
  }
  if(bEMCAL){
     AliAnalysisTaskHLTEMCAL *taskEMCAL = new AliAnalysisTaskHLTEMCAL("offhlt_comparison_EMCAL");
     mgr->AddTask(taskEMCAL);
     AliAnalysisDataContainer *coutputEMCAL =  mgr->CreateContainer("emcal_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-EMCAL-comparison.root");  
     mgr->ConnectInput(taskEMCAL,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskEMCAL,1,coutputEMCAL);
  }
  
  if(bGLOBAL){
     AliAnalysisTaskHLT *taskGLOBAL = new AliAnalysisTaskHLT("offhlt_comparison_GLOBAL");
     taskGLOBAL->SetUseHLTTriggerDecision(fUseHLTTrigger);
     if(fUseHLTTrigger==kTRUE) printf("\n\nOnly HLT triggered events will be used to fill the distributions for task %s.\n\n", taskGLOBAL->GetName());
     //taskGLOBAL->SelectCollisionCandidates();
     mgr->AddTask(taskGLOBAL);
     if(fUseHLTTrigger==kFALSE)
       AliAnalysisDataContainer *coutputGLOBAL =  mgr->CreateContainer("global_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-GLOBAL-comparison.root");  
     else
       AliAnalysisDataContainer *coutputGLOBAL =  mgr->CreateContainer("global_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-GLOBAL-comparison_triggered.root");  
     mgr->ConnectInput(taskGLOBAL,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskGLOBAL,1,coutputGLOBAL);
  }
  if(bD0){
    float cuts[7]={0.5,0.04,0.7,0.8,0.05,-0.00025,0.7};
    AliAnalysisTaskD0Trigger *taskD0 = new AliAnalysisTaskD0Trigger("offhlt_comparison_D0_Trigger",cuts);
    mgr->AddTask(taskD0);
    AliAnalysisDataContainer *coutputD0 =  mgr->CreateContainer("D0_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-D0-comparison.root");  
    mgr->ConnectInput(taskD0,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskD0,1,coutputD0);
  }
  if(bCB){
     AliAnalysisTaskHLTCentralBarrel *taskCB = new AliAnalysisTaskHLTCentralBarrel("offhlt_comparison_CB");
     mgr->AddTask(taskCB);     
     taskCB->SetBeamType(beamType);
     if(beamType.Contains("Pb-Pb")){
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
        AliCentralitySelectionTask *taskCentrality = AddTaskCentrality(); 
     }   
     AliAnalysisDataContainer *coutputCB =  mgr->CreateContainer("esd_thnsparse",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-CentralBarrel-comparison.root");       
     mgr->ConnectInput(taskCB,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskCB,1,coutputCB);
  }
  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("grid", nEvents);

  timer.Stop();
  timer.Print();
}

void compare_HLT_offline_grid(){
  cout << " " << endl;
  cout << " Usage examples:" << endl;
  cout << "    compare-HLT-offline-grid.C'(runNumber, dataDir, gridWorkingDir, gridOutputDir, mode, taskOption, taskFolder, fUseHLTTrigger)' 2>&1 | tee log" << endl;
  cout << "    compare-HLT-offline-grid.C'(\"000115322\",\"/alice/data/2010/LHC10b\",\"ESDcomparison\",\"output\",\"full\",\"global\",\"./\",kTRUE)' 2>&1 | tee log" << endl;
  cout << " " << endl;
}
