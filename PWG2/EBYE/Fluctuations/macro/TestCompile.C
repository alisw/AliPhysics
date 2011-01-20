void TestCompile()
{

  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG2ebye.so");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  printf("Library is Loaded \n");


  gROOT->LoadMacro("../AliEbyEEventBase.cxx++g");  
  gROOT->LoadMacro("../AliEbyEMultiplicityFluctuationAnalysis.cxx++g"); 
  gROOT->LoadMacro("../AliEbyEMFAnalysisTask.cxx++g"); 
  gROOT->LoadMacro("../AliEbyEFluctuationAnalysis.cxx++g"); 
  gROOT->LoadMacro("../AliEbyEFluctuationTask.cxx++g"); 
  gROOT->LoadMacro("../AliEbyEMFAnalysisTaskT.cxx++g"); 
  gROOT->LoadMacro("../AliEbyEChargeFluctuationAnalysis.cxx++g");  
  gROOT->LoadMacro("../AliEbyECFAnalysisTask.cxx++g");  

  gROOT->LoadMacro("AddFluctuationTask.C");

  AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  AddTaskMF();
  AddTaskMFT();
  AddTaskCF();
  AddTaskF();

  mgr->InitAnalysis();
  mgr->PrintStatus();



}
