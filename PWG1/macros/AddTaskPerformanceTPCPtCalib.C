///////////////////////////////////////////////////////////////////////////////
// Macro to setup AliPerformanceTask for 
// TPC performance QA to run on PWG1 QA train. 
//
// Input: ESDs, ESDfriends (optional), Kinematics (optional), TrackRefs (optional)
// ESD and MC input handlers must be attached to AliAnalysisManager
// to run default configuration. 
//
// By default 1 performance component is added to 
// the task: 
// 1. AliPerformancePtCalib
// or AliPerformancePtCalibMC if bUseMCinfo = kTRUE (use MC info)

// Usage on the analysis train (default configuration):
// gSystem->Load("libANALYSIS");
// gSystem->Load("libANALYSISalice");
// gSystem->Load("libTPCcalib.so");
// gSystem->Load("libPWG1.so");
//
// gROOT->LoadMacro("$ALICE_ROOT/PWG1/macros/AddTaskPerformanceTPC.C");
// AliPerformanceTask *tpcQA = AddTaskPerformanceTPC("kTRUE","kTRUE"); 
// 
// Output:
// TPC.Performance.root file with TPC performance components is created.
//
// Each of the components contains THnSparse generic histograms which 
// have to be analysed (post-analysis) by using Analyse() function. 
// Each component contains such function.
//
//13.10.2009 -  J.Otwinowski@gsi.de
///////////////////////////////////////////////////////////////////////////////

//____________________________________________
AliPerformanceTask* AddTaskPerformanceTPCPtCalib(Bool_t bUseMCInfo=kTRUE, Bool_t bUseESDfriend=kTRUE)
{
  //
  // Add AliPerformanceTask with TPC performance components
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) { 
     Error("AddTaskPerformanceTPCPtCalib","AliAnalysisManager not set!");
    return NULL;
  }
  
  TString type = mgr->GetInputEventHandler()->GetDataType();
  if (!type.Contains("ESD")) {
     Error("AddTaskPerformanceTPCPtCalib", "ESD input handler needed!");
    return NULL;
  }
  
  AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
  if (!mcH && bUseMCInfo) {
    Error("AddTaskPerformanceTPCPtCalib", "MC input handler needed!");
    return NULL;
  }

  //
  // Create task
  //
  AliPerformanceTask *task = new AliPerformanceTask("Performance","TPC Performance PtCalib");
  if (!task) {
    Error("AddTaskPerformanceTPCPtCalib", "TPC performance task cannot be created!");
    return NULL;
  }
  task->SetUseMCInfo(bUseMCInfo);
  task->SetUseESDfriend(bUseESDfriend);

  //
  // Add task to analysis manager
  //
  mgr->AddTask(task);

  //
  // Create TPC-ESD track reconstruction cuts
  //
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetMaxDCAToVertexXY(3.0);
    pRecInfoCuts->SetMaxDCAToVertexZ(3.0);
    pRecInfoCuts->SetMinNClustersTPC(50);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetHistogramsOn(kFALSE); 
  } 
  else {
    Error("AddTaskPerformanceTPCPtCalib", "AliRecInfoCuts cannot be created!");
    return NULL;
  }
  //
  // Create TPC-MC track reconstruction cuts
  //
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinTrackLength(70);
  } 
  else {
    Error("AddTaskPerformanceTPCPtCalib", "AliMCInfoCuts cannot be created!");
    return NULL;
  }

  //
  // Create performance objects for TPC and set cuts 
  //
  enum { kTPC = 0, kTPCITS, kConstrained, kTPCInner, kTPCOuter, kTPCSec };


  AliPerformancePtCalib *ptCalib =  NULL;
  AliPerformancePtCalibMC *ptCalibMC = NULL;

  if(bUseMCInfo){
     ptCalibMC = new AliPerformancePtCalibMC("AliPerformancePtCalibMC","AliPerformancePtCalibMC");//,kTPC,kTRUE);
     if(!ptCalibMC) {
	Error("AddTaskPerformanceTPCPtCalib", "Cannot create AliPerformancePtCalibMC");
     }
    
     ptCalibMC->SetAliRecInfoCuts(pRecInfoCuts);
     ptCalibMC->SetReadTPCTracks(kTRUE);  
     ptCalibMC->SetTPCRefit(kFALSE) ;         
     ptCalibMC->SetITSRefit(kFALSE);          
     ptCalibMC->SetESDCuts(kTRUE);            
     ptCalibMC->SetDCACuts(kTRUE);             
     ptCalibMC->SetAcceptKinkDaughters(kFALSE);   
     ptCalibMC->SetRequireSigmaToVertex(kFALSE);
     ptCalibMC->SetfDCAToVertex2D(kFALSE)   ;

     // const Double_t esdCutvalues[6] ={};//set esd track cut values
     // ptCalibMC->SetESDcutValues(esdCutvalues);
     // ptCalibMC->SetEtaRange(0.9);
     // ptCalibMC->SetAliMCInfoCuts(pMCInfoCut);
     
  }
  else{

     ptCalib =  new AliPerformancePtCalib("AliPerformancePtCalib","AliPerformancePtCalib");//,kTPC,kFALSE);
     if(!ptCalib) {
	Error("AddTaskPerformanceTPCPtCalib", "Cannot create AliPerformancePtCalib");
     }
     ptCalib->SetAliRecInfoCuts(pRecInfoCuts);
     ptCalib->SetReadTPCTracks(kTRUE);  
     ptCalib->SetTPCRefit(kFALSE) ;         
     ptCalib->SetITSRefit(kFALSE);          
     ptCalib->SetESDCuts(kTRUE);            
     ptCalib->SetDCACuts(kTRUE);             
     ptCalib->SetAcceptKinkDaughters(kFALSE);   
     ptCalib->SetRequireSigmaToVertex(kFALSE);
     ptCalib->SetfDCAToVertex2D(kFALSE)   ;
     
     // const Double_t esdCutvalues[6] ={};
     // ptCalib->SetESDcutValues(esdCutvalues);
     // ptCalib->SetEtaRange(0.9);
     // ptCalib->SetAliMCInfoCuts(pMCInfoCut);
  }
 
     
  // add components to the performance task
  
  if(bUseMCInfo) task->AddPerformanceObject( ptCalibMC);
  else task->AddPerformanceObject( ptCalib );
  
  // Create containers for input
  //
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  
  //
  // Create containers for output
  //
  AliAnalysisDataContainer *coutput_tpcptcalib = mgr->CreateContainer("TPCPtCalib", TList::Class(), AliAnalysisManager::kOutputContainer, Form("TPCPtCalib.%s.root", task->GetName()));
  mgr->ConnectOutput(task, 0, coutput_tpcptcalib);

return task;  
}
