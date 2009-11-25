void runLocal() {
/*
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
*/
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWG1");

  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  
  //gROOT->LoadMacro("AliAnalysisTaskGlobalQA.cxx++g");
  AliAnalysisTask *task = new AliAnalysisTaskGlobalQA();
  mgr->AddTask(task);
  
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  chain=CreateESDChain("list.txt",1);

  // Create containers for input/output
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput = 
     mgr->CreateContainer("coutput", TObjArray::Class(),
     AliAnalysisManager::kOutputContainer, "GlobalQA.root" );
  mgr->ConnectOutput(task,1,coutput);

  mgr->SetDebugLevel(0);
  
  if (!mgr->InitAnalysis()) return;
  //mgr->PrintStatus();

  mgr->StartAnalysis("local",chain);
}

void runProof() {
    TProof::Open("belikov@localhost"); 
  /*
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
  */
    gProof->UploadPackage("STEERBase");
    gProof->EnablePackage("STEERBase");
    gProof->UploadPackage("ESD");
    gProof->EnablePackage("ESD");
    gProof->UploadPackage("AOD");
    gProof->EnablePackage("AOD");
    gProof->UploadPackage("ANALYSIS");
    gProof->EnablePackage("ANALYSIS");
    gProof->UploadPackage("ANALYSISalice");
    gProof->EnablePackage("ANALYSISalice");

    gROOT->ProcessLine(".include $ALICE_ROOT/include");

  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  gProof->Load("$ALICE_ROOT/PWG1/global/AliAnalysisTaskGlobalQA.cxx++g");
  AliAnalysisTask *task = new AliAnalysisTaskGlobalQA();
  mgr->AddTask(task);
  
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  AliAnalysisDataContainer *coutput = 
     mgr->CreateContainer("coutput", TObjArray::Class(),
     AliAnalysisManager::kOutputContainer, "GlobalQA.root" );
  mgr->ConnectOutput(task,1,coutput);

  if (!mgr->InitAnalysis()) return;

  mgr->StartAnalysis("proof","/COMMON/COMMON/LHC09d5_0.9TeV_0T",3000);

}

