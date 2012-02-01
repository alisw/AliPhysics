void runAODProof(const char * proofMode = "full")
{

   gEnv->SetValue("XSec.GSI.DelegProxy", "2");
   //  TProof::Open("rbertens@alice-caf.cern.ch");

   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libSTEERBase.so");
   gSystem->Load("libESD.so");
   gSystem->Load("libAOD.so");
   gSystem->Load("libANALYSIS.so");
   gSystem->Load("libOADB.so");
   gSystem->Load("libANALYSISalice.so");
   gSystem->AddIncludePath("-I$ALICE_ROOT/include");


   AliAnalysisAlien * handler = new AliAnalysisAlien("test");
   handler->SetOverwriteMode();
   handler->SetRunMode(proofMode);
   handler->SetProofReset(0);
   handler->SetAliROOTVersion("v4-21-29-AN");
   handler->SetProofCluster(Form("%s@alice-caf.cern.ch", gSystem->Getenv("CAFUSER")));
   //handler->SetProofCluster("rbertens@skaf.saske.sk");
   //   handler->SetProofDataSet("/alice/data/LHC10h_000138653_p2_AOD049#aodTree");
   handler->SetProofDataSet("/alice/sim/LHC11a10a_000138653_AOD048#aodTree");
   handler->SetNproofWorkersPerSlave(1);
   handler->SetAliRootMode("default");
   handler->SetAdditionalLibs("AliSpectraAODHistoManager.cxx AliSpectraAODHistoManager.h AliSpectraAODEventCuts.cxx AliSpectraAODEventCuts.h AliSpectraAODTrackCuts.cxx AliSpectraAODTrackCuts.h AliAnalysisTaskSpectraAOD.cxx AliAnalysisTaskSpectraAOD.h");
   handler->SetAnalysisSource("AliSpectraAODHistoManager.cxx+ AliSpectraAODEventCuts.cxx+ AliSpectraAODTrackCuts.cxx+ AliAnalysisTaskSpectraAOD.cxx+");
   handler->SetFileForTestMode("filelist.txt"); // list of local files for testing

   //  handler->SetAliRootMode("");
   handler->SetClearPackages();


   AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
   mgr->SetGridHandler(handler);
   AliAODInputHandler* aodH = new AliAODInputHandler();
   mgr->SetInputEventHandler(aodH);

   gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
   gROOT->LoadMacro("AliSpectraAODEventCuts.cxx+g");
   gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
   gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");

   // Add PID task
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   AliAnalysisTask * taskPID = AddTaskPIDResponse();
   mgr->AddTask(taskPID);

   AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD("TaskAODExercise");
   task->SetIsMC(1);
   mgr->AddTask(task);

   // Set the cuts
   AliSpectraAODEventCuts * vcuts = new AliSpectraAODEventCuts("Event Cuts");
   AliSpectraAODTrackCuts  * tcuts = new AliSpectraAODTrackCuts("Track Cuts");
   tcuts->SetTrackType(6);
   tcuts->SetEta(1.);
   tcuts->SetP(.7);
   // vcuts->SetCentralityCutMin(0.0) // default
   // vcuts->SetCentralityCutMax(0.0) // default
   task->SetEventCuts(vcuts);
   task->SetTrackCuts(tcuts);
   task->SetNSigmaForIdentification(3.);

   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("chistpt", AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, "Pt.AOD.1.root");
   AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("cvcutpt", AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, "Pt.AOD.1.root");
   AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer("ctcutpt", AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, "Pt.AOD.1.root");

   mgr->ConnectInput(task, 0, cinput);
   mgr->ConnectOutput(task, 1, coutputpt1);
   mgr->ConnectOutput(task, 2, coutputpt2);
   mgr->ConnectOutput(task, 3, coutputpt3);
   mgr->SetDebugLevel(2);

   if (!mgr->InitAnalysis()) return;
   mgr->PrintStatus();
   mgr->StartAnalysis("proof");
}
