void runAOD(Bool_t isMC = 0)
{

//gEnv->SetValue("XSec.GSI.DelegProxy", "2");
//TProof::Open("pverstee@alice-caf.cern.ch");

   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libSTEERBase.so");
   gSystem->Load("libESD.so");
   gSystem->Load("libAOD.so");
   gSystem->Load("libANALYSIS.so");
   gSystem->Load("libANALYSISalice.so");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");
   gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");

   ifstream flist;
   flist.open("filelist.txt");
   TChain * chain = new TChain("aodTree");
   TString line;
   while (line.ReadLine(flist))
   {
      gSystem->ExpandPathName(line);
      chain->AddFile(line.Data());
      cout << "Adding file " << line.Data() << endl;
   }
//  gProof->UploadPackage("AF-v4-19-04-AN");
//  gProof->EnablePackage("AF-v4-19-04-AN");
   gSystem->AddIncludePath("-I$ALICE_ROOT/include");
   gStyle->SetPalette(1);



   AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
   AliAODInputHandler* aodH = new AliAODInputHandler();
   mgr->SetInputEventHandler(aodH);

   gROOT->LoadMacro("AliSpectraAODTrackCuts.cxx+g");
   gROOT->LoadMacro("AliSpectraAODVertexCuts.cxx+g");
   gROOT->LoadMacro("AliSpectraAODHistoManager.cxx+g");
   gROOT->LoadMacro("AliAnalysisTaskSpectraAOD.cxx+g");

   gROOT->LoadMacro("$(ALICE_ROOT)/ANALYSIS/macros/AddTaskPIDResponse.C");
   AliAnalysisTask * taskPID = AddTaskPIDResponse();
   mgr->AddTask(taskPID);


   AliAnalysisTaskSpectraAOD *task = new AliAnalysisTaskSpectraAOD("TaskAODExercise");
   task->SetIsMC(isMC);
   mgr->AddTask(task);

   // Set the cuts
   AliSpectraAODVertexCuts * vcuts = new AliSpectraAODVertexCuts("Vertex Cuts");
   AliSpectraAODTrackCuts  * tcuts = new AliSpectraAODTrackCuts("Tracks Cuts");
   tcuts->SetTrackType(6);
   tcuts->SetEta(1.);
   task->SetVertexCuts(vcuts);
   task->SetTrackCuts(tcuts);


   AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
   AliAnalysisDataContainer *coutputpt1 = mgr->CreateContainer("chistpt", AliSpectraAODHistoManager::Class(),  AliAnalysisManager::kOutputContainer, "Pt.AOD.1.root");
   AliAnalysisDataContainer *coutputpt2 = mgr->CreateContainer("cvcutpt", AliSpectraAODVertexCuts::Class(),    AliAnalysisManager::kOutputContainer, "Pt.AOD.1.root");
   AliAnalysisDataContainer *coutputpt3 = mgr->CreateContainer("ctcutpt", AliSpectraAODTrackCuts::Class(),     AliAnalysisManager::kOutputContainer, "Pt.AOD.1.root");

   mgr->ConnectInput(task, 0, cinput);
   mgr->ConnectOutput(task, 1, coutputpt1);
   mgr->ConnectOutput(task, 2, coutputpt2);
   mgr->ConnectOutput(task, 3, coutputpt3);
   mgr->SetDebugLevel(2);

   if (!mgr->InitAnalysis()) return;
   mgr->PrintStatus();
   mgr->StartAnalysis("local", chain, 100);
}
