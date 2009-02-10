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

  //2. Make a chain e.g.:
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool; 
  TChain * chainEsd = tool.MakeChain("esd.txt","esdTree",0,5);
  chainEsd->Lookup();
  //

  //3. Make a analysis manager with attached task:
  .L $ALICE_ROOT/PWG1/Macros/taskComp.C
  Init();
  AliAnalysisManager *mgr = MakeManager();
  
  //4. Process task localy
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
  gProof->Exec("gROOT->Macro(Form(\"$ALICE_ROOT/TPC/macros/ConfigOCDB.C(%f)\",0))")
  mgr->StartAnalysis("proof",chainEsd);
  //5. Get debug stream - if speciefied  
  TFile f("mcTaskDebug.root");
  TTree *treeCMP = (TTree*)f.Get("RC");

  //6. Read the analysis object
  TFile f("Output.root");
  TObjArray * array = (TObjArray*)f.Get("AliComparisonRes");
  AliComparisonRes * compObj = (AliComparisonRes*)array->FindObject("AliComparisonRes");

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

  AliGenInfoTask *genTask = new AliGenInfoTask("genTask");
  genTask->SetStreamLevel(10);
  genTask->SetDebugLevel(10);
  AddComparison(genTask);
  mgr->AddTask(genTask);
  //
  //
  AliAnalysisDataContainer *cinput1 =
    mgr->CreateContainer("cchain1",TChain::Class(),						  AliAnalysisManager::kInputContainer);  
  mgr->ConnectInput(genTask,0,cinput1);
  //
  AliAnalysisDataContainer *coutput1
    =mgr->CreateContainer("AliComparisonRes",TObjArray::Class(),
			  AliAnalysisManager::kOutputContainer,
			  "Output.root");
  mgr->ConnectOutput(genTask,0,coutput1);
  //
  if (!mgr->InitAnalysis()) return 0;
  return mgr;
}

void AddComparison( AliGenInfoTask * task){
  
  Int_t magField = 5;  // 0 - 0.2 T, 1 = 0.4 T, 2  - 0.5 T

  // Create ESD track reconstruction cuts
  AliRecInfoCuts *pRecInfoCuts = new AliRecInfoCuts(); 
  if(pRecInfoCuts) {
    pRecInfoCuts->SetPtRange(0.15,200.0);
    pRecInfoCuts->SetMaxAbsTanTheta(1.0);
    pRecInfoCuts->SetMinNClustersTPC(10);
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
  AliComparisonRes *pCompRes = new AliComparisonRes(); 
  if(!pCompRes) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonRes object");
  }
  pCompRes->SetAliRecInfoCuts(pRecInfoCuts);
  pCompRes->SetAliMCInfoCuts(pMCInfoCuts);

  // Efficiency
  AliComparisonEff *pCompEff =  new AliComparisonEff();
  if(!pCompEff) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonEff object");
  }
  pCompEff->SetAliRecInfoCuts(pRecInfoCuts);
  pCompEff->SetAliMCInfoCuts(pMCInfoCuts);

  // dE/dx
  AliComparisonDEdx *pCompDEdx = new AliComparisonDEdx();
  if(!pCompDEdx) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDEdx object");
  }
  pCompDEdx->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDEdx->SetAliMCInfoCuts(pMCInfoCuts);
  pCompDEdx->SetMCPtMin(0.5);
  pCompDEdx->SetMCAbsTanThetaMax(0.5);
  pCompDEdx->SetMCPdgCode(pMCInfoCuts->GetPiP()); // only pi+ particles

  // DCA
  AliComparisonDCA *pCompDCA = new AliComparisonDCA();
  if(!pCompDCA) {
    AliDebug(AliLog::kError, "ERROR: Cannot create AliComparisonDCA object");
  }
  pCompDCA->SetAliRecInfoCuts(pRecInfoCuts);
  pCompDCA->SetAliMCInfoCuts(pMCInfoCuts);
  //
  //
  //
  //task->SetMagField(magField);
  task->AddComparisonObject( pCompRes );
  task->AddComparisonObject( pCompEff );
  task->AddComparisonObject( pCompDEdx );
  task->AddComparisonObject( pCompDCA );  
}

/*
  chaiEsd->SetAlias("dp2K","(Kinks[].fParamMother.fP[2]-Kinks[].fParamDaughter.fP[2])");
  chaiEsd->SetAlias("sp2K","sqrt(Kinks[].fParamMother.fC[5]+Kinks[].fParamDaughter.fC[5])");
  chaiEsd->SetAlias("dp3K","(Kinks[].fParamMother.fP[3]-Kinks[].fParamDaughter.fP[3])");
  chaiEsd->SetAlias("sp3K","sqrt(Kinks[].fParamMother.fC[9]+Kinks[].fParamDaughter.fC[9])");

*/






