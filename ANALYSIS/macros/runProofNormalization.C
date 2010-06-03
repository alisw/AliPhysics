
void runProofNormalization(const char * dataset = "LHC09b12_7TeV_0.5T", TString dataSetPath ="/PWG0/jgrosseo/",const char * filename = "LHC09b12_7TeV_0.5T_norm.root", Bool_t isMC = 1,Int_t nev =123456789) {

  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  TProof::Open("alicecaf");
  
  //  gSystem->AddIncludePath("-I${ALICE_ROOT}/include/ -I${ALICE_ROOT}/PWG0/ -I${ALICE_ROOT}/PWG0/dNdEta/");
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include/");
  gProof->UploadPackage("$ALICE_ROOT/STEERBase");
  gProof->EnablePackage("$ALICE_ROOT/STEERBase");
  gProof->UploadPackage("$ALICE_ROOT/ESD");
  gProof->EnablePackage("$ALICE_ROOT/ESD");
  gProof->UploadPackage("$ALICE_ROOT/AOD");
  gProof->EnablePackage("$ALICE_ROOT/AOD");
  gProof->UploadPackage("$ALICE_ROOT/ANALYSIS");
  gProof->EnablePackage("$ALICE_ROOT/ANALYSIS");
  gProof->UploadPackage("$ALICE_ROOT/ANALYSISalice");
  gProof->EnablePackage("$ALICE_ROOT/ANALYSISalice");
  gProof->UploadPackage("$ALICE_ROOT/CORRFW");
  gProof->EnablePackage("$ALICE_ROOT/CORRFW");
 

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  //  mgr->SetDebugLevel(3);
  // Add ESD handler
  AliESDInputHandler* esdH = new AliESDInputHandler; 

  mgr->SetInputEventHandler(esdH);
	
  if(isMC) {
    AliMCEventHandler *mc = new AliMCEventHandler();
    mc->SetReadTR(kFALSE);
    mgr->SetMCtruthEventHandler(mc);
  }
  // assign simple task
//   gProof->Load("AliCollisionNormalization.cxx++g");   
//   gProof->Load("AliCollisionNormalizationTask.cxx++g");   
  //____________________________________________//
  // assign simple task
  AliCollisionNormalizationTask * task = new AliCollisionNormalizationTask("TaskNormalization");
  //  task->SetMC();
  task->SetMC(isMC);
  mgr->AddTask(task);


  gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC,1,!isMC); // Use Physics Selection. Enable computation of BG if is not MC
  //  task->SelectCollisionCandidates(); /// This should be disabled, at least for MC: we need all the events
  physSelTask->GetPhysicsSelection()->SetBin0Callback("TaskNormalization");

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();	
  mgr->ConnectInput(task,0,cinput1);


  
  // Attach output
  cOutput = mgr->CreateContainer("Norm", TList::Class(), AliAnalysisManager::kOutputContainer,filename);
  mgr->ConnectOutput(task, 1, cOutput);      
	
  if (!mgr->InitAnalysis()) return;
	
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",dataSetPath+dataset+"#esdTree",nev);

}
