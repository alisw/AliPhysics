void runAODProof(Int_t c, const char * proofMode = "full")
{

   gEnv->SetValue("XSec.GSI.DelegProxy", "2");

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
//   handler->SetProofCluster(Form("%s@skaf.saske.sk",gSystem->Getenv("CAFUSER")));

   // Set handler for Real DATA:
   if (c == 1)
   {
     handler->SetProofDataSet("/alice/data/LHC10h_000138653_p2_AOD049#aodTree");
 //     handler->SetProofDataSet("/alice/sim/LHC11a10a_000139107_AOD048#aodTree|alice/sim/LHC11a10a_000138653_AOD048#aodTree");
   }
   if (c == 2)
   {
      handler->SetProofDataSet("/alice/sim/LHC11a10a_000139107_AOD048#aodTree|/alice/sim/LHC11a10a_000138653_AOD048#aodTree|/alice/sim/LHC11a10a_000139110_AOD048#aodTree|/alice/sim/LHC11a10a_000138662_AOD048#aodTree|/alice/sim/LHC11a10a_000138666_AOD048#aodTree|/alice/sim/LHC11a10a_000138795_AOD048#aodTree");      
   }
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
   Bool_t isMC = kFALSE;
   if (c == 2 || c == 3) isMC = kTRUE;   
   AliAnalysisTask * taskPID = AddTaskPIDResponse(isMC);
   mgr->AddTask(taskPID);


   AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD("TaskAODExercise");
   mgr->AddTask(task);

   // Set the cuts
   AliSpectraAODEventCuts * vcuts = new AliSpectraAODEventCuts("Event Cuts");
   AliSpectraAODTrackCuts  * tcuts = new AliSpectraAODTrackCuts("Track Cuts");
   tcuts->SetTrackType(6);
   tcuts->SetEta(.8);
   //   tcuts->SetDCA(.1);
   tcuts->SetPt(1.2);   
   vcuts->SetCentralityCutMax(5.);  // example min max cuts
   vcuts->SetCentralityCutMin(0.);
   task->SetEventCuts(vcuts);
   task->SetTrackCuts(tcuts);
   task->SetNSigmaForIdentification(3.); // FIXME
   task->SetYCut(.5);
   // check for MC or real data
   if (c == 2 || c == 3)
   {
      task->SetIsMC(kTRUE);


      AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
      AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("chistpt", AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, "Pt.AOD.1._MC.root");
      AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("cvcutpt", AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, "Pt.AOD.1._MC.root");
      AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer("ctcutpt", AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, "Pt.AOD.1._MC.root");
   }
   if (c == 1)
   {
      AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
      AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("chistpt", AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, "Pt.AOD.1._data_ptcut.root");
      AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("cvcutpt", AliSpectraAODEventCuts::Class(),    AliAnalysisManager::kOutputContainer, "Pt.AOD.1._data_ptcut.root");
      AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer("ctcutpt", AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, "Pt.AOD.1._data_ptcut.root");

   }
   mgr->ConnectInput(task, 0, cinput);
   mgr->ConnectOutput(task, 1, coutputpt1);
   mgr->ConnectOutput(task, 2, coutputpt2);
   mgr->ConnectOutput(task, 3, coutputpt3);
   mgr->SetDebugLevel(2);

   if (!mgr->InitAnalysis()) return;
   mgr->PrintStatus();
   if (c == 3)
   {
      mgr->StartAnalysis("local");
   }
   if (c != 3)
   {
      mgr->StartAnalysis("proof");
   }
}
