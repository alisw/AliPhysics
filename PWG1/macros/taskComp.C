/*

  Sequence hot to se the PWG1 analysis tasks:
  

  //1. Load libraries if needed:
  //
  gSystem->Load("/usr/local/grid/XRootd/GSI/lib/libXrdClient.so");  
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG1.so");

  AliLog::SetGlobalLogLevel(AliLog::kError);

  //2. Make a chain e.g.:
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  AliXRDPROOFtoolkit tool; 
  TChain * chainEsd = tool.MakeChain("esd.txt","esdTree",0,5000);
  chainEsd->Lookup();
  //

 
  //3. Make a analysis manager with attached task:
  .L $ALICE_ROOT/PWG1/Macros/taskComp.C
  Init();
  AliAnalysisManager *mgr = MakeManager();
  
  //4. Process task localy
  //
  mgr->SetNSysInfo(100);
  mgr->SetDebugLevel(1);
  mgr->StartAnalysis("local",chainEsd);

  //
  //4. Process task on proof
  //
  TProof::Open("");
  .L /u/miranov/macros/ProofEnableAliRoot.C
  ProofEnableAliRoot("/usr/local/grid/AliRoot/HEAD0108");
  gProof->Exec("gSystem->Load(\"libANALYSIS.so\")",kTRUE);
  gProof->Exec("gSystem->Load(\"libAOD.so\")",kTRUE);
  gProof->Exec("gSystem->Load(\"libANALYSISalice.so\")",kTRUE);
  gProof->Exec("gSystem->Load(\"libPWG0base.so\")",kTRUE);
  gProof->Exec("gSystem->Load(\"libPWG0dep.so\")",kTRUE);
  gProof->Exec("gSystem->Load(\"libPWG1.so\")",kTRUE);
  mgr->StartAnalysis("proof",chainEsd);
  //5. Get debug stream - if speciefied  
  TFile f("mcTaskDebug.root");
  TTree *treeCMP = (TTree*)f.Get("RC");

  //6. Read the analysis object
  TFile f("Output.root");
  TObjArray * array = (TObjArray*)f.Get("AliComparisonRes");
  AliComparisonRes * compObj = (AliComparisonRes*)array->FindObject("AliComparisonRes");
  //
  //7. Get debug streamer on PROOF
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chainTr = tool.MakeChain("cmp.txt","RC",0,1000000);
  chainTr->Lookup();
  chainTr->SetProof(kTRUE);
  TChain * chainTPC = tool.MakeChain("tpc.txt","Crefit",0,50);
  chainTPC->Lookup();
  chainTr->SetProof(kTRUE);



  TFile f("mcTaskDebug.root");
  TTree *treeCMP = (TTree*)f.Get("RC");

*/



void AddComparison( AliGenInfoTask * task);

void Init(){
  //
  // Init mag field and the geo manager
  // 
  TGeoManager::Import("/u/miranov/proof/geometry.root");
  AliGeomManager::LoadGeometry("/u/miranov/proof/geometry.root");
  
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));


}

AliAnalysisManager *  MakeManager(){
  //
  //
  //
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisComponentManager");
  mgr->SetDebugLevel(1);  
  cout << "Creating ESD event handler" << endl; 
  AliESDInputHandler* esdH = new AliESDInputHandler();
  // set the ESDfriend branch active (my modification of AliESDInputHandler)
  esdH->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdH); 
  
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  //


  AliGenInfoTask *genTask = new AliGenInfoTask("genTask");
  genTask->SetStreamLevel(10);
  genTask->SetDebugLevel(10); 
  genTask->SetDebugOuputhPath(gSystem->pwd());

 //  //AddComparison(genTask);
//   mgr->AddTask(genTask);
//   //
//   //
//   AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
//   mgr->ConnectInput(genTask,0,cinput1);
//   //
//   AliAnalysisDataContainer *coutput1
//     =mgr->CreateContainer("AliComparisonRes",TObjArray::Class(),
// 			   AliAnalysisManager::kOutputContainer,
//  			  "Output.root");
//   mgr->ConnectOutput(genTask,0,coutput1);

  
  //
  // TPC PID task
  //
  AliTPCtaskPID *pidTask = new AliTPCtaskPID("pidTask");
  mgr->AddTask(pidTask);
  AliAnalysisDataContainer *cinput2 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(pidTask,0,cinput2);
  //
  AliAnalysisDataContainer *coutput2
    =mgr->CreateContainer("tpcTaskPID", TObjArray::Class(),
			  AliAnalysisManager::kOutputContainer,
			  "OutputPID.root");
  mgr->ConnectOutput(pidTask,0,coutput2);

  //
  // TPC QA task
  //
  AliTPCtaskQA *qaTask = new AliTPCtaskQA("qaTask");
  mgr->AddTask(qaTask);
  AliAnalysisDataContainer *cinput3 = mgr->GetCommonInputContainer();
  mgr->ConnectInput(qaTask,0,cinput3);
  //
  AliAnalysisDataContainer *coutput3
    =mgr->CreateContainer("tpcTaskQA", TObjArray::Class(),
			  AliAnalysisManager::kOutputContainer,
			  "OutputQA.root");
  mgr->ConnectOutput(qaTask,0,coutput3);
  //
  if (!mgr->InitAnalysis()) return 0;
  return mgr;
}

void AddComparison( AliGenInfoTask * task){
  
  // Create ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetPtRange(0.15,200.0);
    pRecInfoCuts->SetMaxAbsTanTheta(1.0);
    pRecInfoCuts->SetMinNClustersTPC(10);
    pRecInfoCuts->SetMinNClustersITS(2);
    pRecInfoCuts->SetMinTPCsignalN(50);

	pRecInfoCuts->SetHistogramsOn(kFALSE); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliRecInfoCuts object");
  }

  // Create MC track reconstruction cuts
  AliMCInfoCuts  *pMCInfoCuts = new AliMCInfoCuts();
  if(pMCInfoCuts) {
    pMCInfoCuts->SetMinRowsWithDigits(50);
    pMCInfoCuts->SetMaxR(0.001);  
    pMCInfoCuts->SetMaxVz(0.001); 
    pMCInfoCuts->SetRangeTPCSignal(0.5,1.4); 
  } else {
    AliDebug(AliLog::kError, "ERROR: Cannot AliMCInfoCuts object");
  }

  //
  // Create comparison objects and set cuts 
  //

  // Resolution
  AliComparisonRes *pCompRes = new AliComparisonRes("AliComparisonRes","AliComparisonRes"); 
  if(!pCompRes) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonRes object");
  }
  pCompRes->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes->SetAliMCInfoCuts(pMCInfoCuts);

  // Efficiency
  AliComparisonEff *pCompEff =  new AliComparisonEff("AliComparisonEff","AliComparisonEff");
  if(!pCompEff) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonEff object");
  }
  pCompEff->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff->SetAliMCInfoCuts(pMCInfoCuts);

  // dE/dx
  AliComparisonDEdx *pCompDEdx = new AliComparisonDEdx("AliComparisonDEdx","AliComparisonDEdx");
  if(!pCompDEdx) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDEdx object");
  }
  pCompDEdx->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDEdx->SetAliMCInfoCuts(pMCInfoCuts);
  pCompDEdx->SetMCPtMin(0.5);
  pCompDEdx->SetMCAbsTanThetaMax(0.5);
  pCompDEdx->SetMCPdgCode(pMCInfoCuts->GetPiP()); // only pi+ particles

  // DCA
  AliComparisonDCA *pCompDCA = new AliComparisonDCA("AliComparisonDCA","AliComparisonDCA");
  if(!pCompDCA) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDCA object");
  }
  pCompDCA->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA->SetAliMCInfoCuts(pMCInfoCuts);
  //
  //
  //
  task->AddComparisonObject( pCompRes );
  task->AddComparisonObject( pCompEff );
  task->AddComparisonObject( pCompDEdx );
  task->AddComparisonObject( pCompDCA );  
}



