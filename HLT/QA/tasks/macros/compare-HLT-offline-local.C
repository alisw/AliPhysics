// $Id$
/*
 * Example macro to run locally an analysis task for comparing the offline
 * with the HLT esd tree.
 *
 * The output is a root file containing the histograms defined in the
 * analysis task. There is one output file per detector.
 *
 * Run without arguments to get a few examples how to use the macro.
 * 
 * The arguments are:
 * - the input file or txt file containing a list of ESDs to be processed (CreateESDChain takes 20 files as a default argument)
 * - the task you want to use
 * - the path of the task location
 * - the beam type, "p-p" or "Pb-Pb", this is relevant ONLY for the central barrel task at the moment and is 
 *   used to select proper binning and axes ranges for the THnSparse objects that it fills
 * - options to make the central barrel task more flexible and lightweight; you can select if you want to 
 *   fill the THnSparse object with only event or track properties or only HLT data or only OFF
 *   possible options are: event-off event-hlt track-off track-hlt, all are turned on by default
 * - boolean variable for selecting events which contain an HLT trigger
 * - number of events to be analyzed
 *
 * If alien:// is placed before the input filename, then the macro connects to the grid to access the file.
 * 
 * In case you want to run over many ESD files, then prepare a list of them in a .txt file and they will be chained for the analysis.
 * The .txt file takes the place of the first argument in that case.
 *
 * @ingroup alihlt_qa
 * @author Kalliopi.Kanaki@ift.uib.no, Hege.Erdal@student.uib.no
 */

void compare_HLT_offline_local( TString file
                               ,const char* detectorTask="global"
			       ,TString taskFolder="$ALICE_ROOT/HLT/QA/tasks/"
			       ,TString beamType="p-p"
			       ,TString options="event-off event-hlt track-off track-hlt"
			       ,bool fUseHLTTrigger=kFALSE
			       ,Long64_t nEvents=1234567890
			      )
{

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
 
  gSystem->AddIncludePath("-I$ALICE_ROOT/HLT/BASE -I$ALICE_ROOT/PWG1/TPC -I. -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS");
  
  gSystem->Load("libTPCcalib.so");
  gSystem->Load("libTRDbase.so");
  gSystem->Load("libTRDrec.so");
  gSystem->Load("libITSbase.so");
  gSystem->Load("libITSrec.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libPWG1.so");
 
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    
  Bool_t bPHOS = kFALSE, bGLOBAL = kFALSE, bEMCAL = kFALSE, bPWG1 = kFALSE, bD0 = kFALSE, bCB = kFALSE;
 
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
      else if(argument.CompareTo("emcal", TString::kIgnoreCase)==0){
	bEMCAL = kTRUE;
	continue;
      }         
      if(argument.CompareTo("global", TString::kIgnoreCase)==0){
	bGLOBAL = kTRUE;
	continue;
      }      
      if(argument.CompareTo("pwg1", TString::kIgnoreCase)==0){
	bPWG1 = kTRUE;
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
      else break;
    }
  }
      
  //-------------- Compile the analysis tasks ---------- //
  
  if(bPHOS){
    //gSystem->Load("libHLTbase");
    gSystem->Load("libAliHLTUtil");
    gSystem->Load("libAliHLTGlobal");
    TString strTask1("AliAnalysisTaskHLTCalo.cxx+");
    TString strTask2("AliAnalysisTaskHLTPHOS.cxx+");    
    gROOT->LoadMacro(taskFolder+strTask1); 
    gROOT->LoadMacro(taskFolder+strTask2); 
    cout << "\n========= You are loading the following tasks --> "<< (taskFolder+strTask1).Chop()  << " and " <<  (taskFolder+strTask2).Chop() << endl;
  }
  
  if(bEMCAL){
    //gSystem->Load("libHLTbase");
    gSystem->Load("libAliHLTUtil");
    gSystem->Load("libAliHLTGlobal");
    TString strTask1("AliAnalysisTaskHLTCalo.cxx+");
    TString strTask2("AliAnalysisTaskHLTEMCAL.cxx+");
    gROOT->LoadMacro(taskFolder+strTask1); 
    gROOT->LoadMacro(taskFolder+strTask2); 
    cout << "\n========= You are loading the following tasks --> "<< (taskFolder+strTask1).Chop()  << " and " <<  (taskFolder+strTask2).Chop() << endl;
  }  
  
  if(bGLOBAL){
     TString strTask("AliAnalysisTaskHLT.cxx+");
     gROOT->LoadMacro(taskFolder+strTask);
     cout << "\n========= You are loading the following task --> "<< (taskFolder+strTask).Chop()  << endl;
  }
  if(bD0){
     TString strTask("AliAnalysisTaskD0Trigger.cxx+");
     gROOT->LoadMacro(taskFolder+strTask); 
     cout << "\n========= You are loading the following task --> "<< (taskFolder+strTask).Chop()  << endl;
  }
 
  if(bCB){
    TString strTask("AliAnalysisTaskHLTCentralBarrel.cxx+");
    gROOT->LoadMacro(taskFolder+strTask);
    cout << "\n========= You are loading the following task --> "<< (taskFolder+strTask).Chop()  << endl;
  }
  
  if(bPWG1) gROOT->LoadMacro("$ALICE_ROOT/HLT/QA/tasks/macros/AddTaskPerformance.C");
   
  if(file.Contains("alien")) TGrid::Connect("alien://");
    
  if(file.Contains("AliESDs.root")){
    TChain *chain = new TChain("esdTree"); 
    chain->Add(file);
  }
  
  //Constructs chain from filenames in *.txt
  //on the form $DIR/AliESDs.root
  else if(file.Contains(".txt")){
    gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
    chain=CreateESDChain(file.Data(),200);
  }
  
  else if(!file){
    printf("File %s does not exist or is corrupted.\n",file.Data());
    return;  
  }

  if(!chain){
    Printf("Chain is empty.\n");
    return;
  }
   
  //-------- Make the analysis manager ---------------//
 
  AliAnalysisManager *mgr  = new AliAnalysisManager("TestManager");
  AliESDInputHandler *esdH = new AliESDInputHandler;

  //For the PWG1 task, setting HLT is handled inside AliPerformanceTask.C
  if(!bPWG1)  esdH->SetReadHLT();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);  
  mgr->SetNSysInfo(1000);

  //AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kFALSE,kTRUE);
 
  //-------------- define the tasks ------------//
  
  if(bPHOS){
    AliAnalysisTaskHLTPHOS *taskPHOS = new AliAnalysisTaskHLTPHOS("offhlt_comparison_PHOS");
    taskPHOS->SetUseHLTTriggerDecision(fUseHLTTrigger);
    if(fUseHLTTrigger==kTRUE) printf("\n\nOnly HLT triggered events will be used to fill the distributions for task %s.\n\n", taskPHOS->GetName());
    mgr->AddTask(taskPHOS);
    if(fUseHLTTrigger==kFALSE)
       AliAnalysisDataContainer *coutputPHOS =  mgr->CreateContainer("phos_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-PHOS-comparison.root");  
    else 
       AliAnalysisDataContainer *coutputPHOS =  mgr->CreateContainer("phos_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-PHOS-comparison_triggered.root");      
    mgr->ConnectInput(taskPHOS,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskPHOS,1,coutputPHOS);
  }

  if(bEMCAL){
    AliAnalysisTaskHLTEMCAL *taskEMCAL = new AliAnalysisTaskHLTEMCAL("offhlt_comparison_EMCAL");
    taskEMCAL->SetUseHLTTriggerDecision(fUseHLTTrigger);
    if(fUseHLTTrigger==kTRUE) printf("\n\nOnly HLT triggered events will be used to fill the distributions for task %s.\n\n", taskEMCAL->GetName());
    mgr->AddTask(taskEMCAL);
    if(fUseHLTTrigger==kFALSE)
       AliAnalysisDataContainer *coutputEMCAL =  mgr->CreateContainer("emcal_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-EMCAL-comparison.root");  
    else
      AliAnalysisDataContainer *coutputEMCAL =  mgr->CreateContainer("emcal_histograms",TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-EMCAL-comparison_triggered.root");      
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
      AliAnalysisDataContainer *coutputGLOBAL =  mgr->CreateContainer("global_histograms",TList::Class(), AliAnalysisManager::kOutputContainer,"HLT-OFFLINE-GLOBAL-comparison_triggered.root");  
    mgr->ConnectInput(taskGLOBAL,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskGLOBAL,1,coutputGLOBAL);
  }

  if(bPWG1){
    Bool_t hasMC=kFALSE;  
    // -- Add Task for HLT and Offline
    AliPerformanceTask *HLTtpcQA = AddTaskPerformance(hasMC,kFALSE,kTRUE);
    AliPerformanceTask *tpcQA = AddTaskPerformance(hasMC,kFALSE); 
    if(!HLTtpcQA || !tpcQA) {
      Error("RunPerformanceTrain","AliPerformanceTask not created!");
      return;
    }
  }
  if(bD0){
    float cuts[7]={0.5,0.04,0.7,0.8,0.05,-0.00025,0.7};
    AliAnalysisTaskD0Trigger *taskD0 = new AliAnalysisTaskD0Trigger("offhlt_comparison_D0",cuts);
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
	taskCentrality->SetPass(1);
     }        
     taskCB->SetOptions(options);
     AliAnalysisDataContainer *coutputCB =  mgr->CreateContainer("esd_thnsparse", TList::Class(), AliAnalysisManager::kOutputContainer, "HLT-OFFLINE-CentralBarrel-comparison.root");  
     mgr->ConnectInput(taskCB,0,mgr->GetCommonInputContainer());
     mgr->ConnectOutput(taskCB,1,coutputCB);
  }
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain, nEvents);

  timer.Stop();
  timer.Print();
}

void compare_HLT_offline_local(){
  cout << " " << endl;
  cout << " Usage examples:" << endl;
  cout << "    compare-HLT-offline-local.C'(file, taskOption, taskFolder, beamType, options, fUseHLTTrigger, nEvents)' 2>&1 | tee log" << endl;
  cout << "    compare-HLT-offline-local.C'(\"AliESDs.root\",\"global\")' 2>&1 | tee log" << endl;
  cout << "    compare-HLT-offline-local.C'(\"AliESDs.root\",\"global\",\"./\", \"p-p\", \"event-off event-hlt track-off track-hlt\", kFALSE, nEvents)' 2>&1 | tee log" << endl;
  cout << "    compare-HLT-offline-local.C'(\"AliESDs.root\",\"global phos cb D0\", \"./\", \"Pb-Pb\", \"event-hlt\", kTRUE, nEvents)' 2>&1 | tee log" << endl;
  cout << "    compare-HLT-offline-local.C'(\"alien:///alice/data/2010/LHC10b/000115322/ESDs/pass1/10000115322040.20/AliESDs.root\",\"global\")' 2>&1 | tee log" << endl;
  cout << " " << endl;
}
