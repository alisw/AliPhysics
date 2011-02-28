void AliAnalysisTaskSigma1385()
{
   TStopwatch timer;
   timer.Start();

   printf("*** Connect to AliEn ***\n");
   TGrid::Connect("alien://");

   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   //____________________________________________________//
   //_____________Setting up required packages___________//
   //____________________________________________________//
   gSystem->Load("libSTEERBase.so");
   gSystem->Load("libSTEER.so");
   gSystem->Load("libESD.so");
   gSystem->Load("libAOD.so");
   gSystem->Load("libANALYSIS.so");
   gSystem->Load("libANALYSISalice.so");
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/TOF");
   //ANALYSIS PART
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
   gROOT->LoadMacro("AliAnalysisTaskSigma1385.cxx+g");

   const char *collectionfile = "sigmaLHC10d1.xml";


   TAlienCollection * myCollection  = TAlienCollection::Open(collectionfile);

   if (!myCollection) {
      Error("AliAnalysisTaskSigma1385", Form("Cannot create an AliEn collection from %s", collectionfile)) ;
      return kFALSE ;
   }

   Info("AliAnalysisTaskSigma1385", Form("Creating the analysis chain %s", "esdTree")) ;
   TChain* chain = new TChain("esdTree");

   Info("AliAnalysisTaskSigma1385", "Preparing the file list") ;
   myCollection->Reset() ;
   while (myCollection->Next()) {
      char esdFile[255] ;
      sprintf(esdFile, "%s", myCollection->GetTURL("")) ;
      Info("AliAnalysisTaskSigma1385", Form("Adding %s", esdFile)) ;
      chain->Add(esdFile) ;
   }


   //____________________________________________//
   // Make the analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
   AliVEventHandler* esdH = new AliESDInputHandler;
   mgr->SetInputEventHandler(esdH);


   //____________________________________________//
   // 1st Pt task
   //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
   AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kTRUE);
   AliAnalysisTaskSigma1385 *task1 = new AliAnalysisTaskSigma1385("TaskSigma1385");


   mgr->AddTask(task1);

   task1->SetDataType("SIM");
   task1->SelectCollisionCandidates();

   if (!(strcmp(task1->GetDataType(), "SIM"))) {

      AliMCEventHandler *mc = new AliMCEventHandler();
      mc->SetReadTR(kFALSE);
      mgr->SetMCtruthEventHandler(mc);
   }
   mgr->SetDebugLevel(10);



   //------ input data ------
   AliAnalysisDataContainer *cinput0 = mgr->GetCommonInputContainer();

   // ----- output data -----

   AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("chist1", TList::Class(), AliAnalysisManager::kOutputContainer, "Sigma.ESD.root");

   //____________________________________________//

   mgr->ConnectInput(task1, 0, cinput0);
   mgr->ConnectOutput(task1, 1, coutput1);
   if (!mgr->InitAnalysis()) return;
   mgr->PrintStatus();
   mgr->StartAnalysis("local", chain);

   timer.Stop();
   timer.Print();
}
