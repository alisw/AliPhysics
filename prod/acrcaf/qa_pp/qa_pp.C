
void AddTaskQAsym  (Int_t runNumber);
void AddTaskVZEROQA(Int_t runNumber);

void qa_pp(Int_t runNumber) {
  TStopwatch timer;
  timer.Start();

  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  /// Select ROOT version
  TProof::Mgr("proof02@alicecaf:31093")->SetROOTVersion("v5-24-00a");
  // Login to CAF
  TProof::Open("proof02@alicecaf:31093");

  // Enable AliRoot
  gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-18-12-AN/AF-v4-18-12-AN.par");
  gProof->EnablePackage("AF-v4-18-12-AN.par");

  
  // Enable analysis libs
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPHOSUtils.so");
  gSystem->Load("libEMCALUtils.so");
  //  gSystem->Load("libPWG4PartCorrBase.so");
  //  gSystem->Load("libPWG4PartCorrDep.so");

  gProof->Exec("gSystem->Load(\"libANALYSIS.so\");",        kTRUE);
  gProof->Exec("gSystem->Load(\"libANALYSISalice.so\");",   kTRUE);
  gProof->Exec("gSystem->Load(\"libPHOSUtils.so\");",       kTRUE);
  gProof->Exec("gSystem->Load(\"libEMCALUtils.so\");",      kTRUE);
  // gProof->Exec("gSystem->Load(\"libPWG4PartCorrBase.so\");",kTRUE);
  // gProof->Exec("gSystem->Load(\"libPWG4PartCorrDep.so\");", kTRUE);

  gProof->Load("AliAnalysisTaskQASym.cxx++g");
  gProof->Load("AliAnaVZEROQA.cxx++g");

  gProof->UploadPackage("PWG4PartCorrBase.par");
  gProof->EnablePackage("PWG4PartCorrBase");

  gProof->UploadPackage("PWG4PartCorrDep.par");
  gProof->EnablePackage("PWG4PartCorrDep");

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisQAManager");
  AliESDInputHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);  
  mgr->SetDebugLevel(10);

  // Wagons
  gROOT->LoadMacro("AddTaskQAsym.C");
  AddTaskQAsym(runNumber);

  gROOT->LoadMacro("AddTaskVZEROQA.C");
  AddTaskVZEROQA(runNumber);

  gROOT->LoadMacro("AddTaskCalorimeterQA.C");
  AliAnalysisTaskParticleCorrelation *taskQAcalo = AddTaskCalorimeterQA("ESD", kFALSE, kTRUE);

  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("proof",
		     Form("/ALIREC/aliprod/run%d",runNumber));

  timer.Stop();
  timer.Print();
}
