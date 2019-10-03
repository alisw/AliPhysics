void runLocal() {
/*
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
*/
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTender");
  gSystem->Load("libPWGPP");

  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  
  //gROOT->LoadMacro("AliAnalysisTaskGlobalQA.cxx++g");
  AliAnalysisTask *task = new AliAnalysisTaskGlobalQA();
  mgr->AddTask(task);
  
  AliESDInputHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);
  
  AliMCEventHandler *mc = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mc);

  gROOT->LoadMacro("$ALICE_PHYSICS/PWG0/CreateESDChain.C");
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
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libPhysics");
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

    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");

  gProof->Load("$ALICE_PHYSICS/PWGPP/global/AliAnalysisTaskGlobalQA.cxx++g");
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

