Bool_t SetupPar(const char *parfile);

void runGridAODpp()
{
   // Load common libraries
  gSystem->Load("libCore");
   gSystem->Load("libTree");
   gSystem->Load("libGeom");
   gSystem->Load("libVMC");
   gSystem->Load("libPhysics");
   gSystem->Load("libMinuit");
   gSystem->Load("libGui");
   gSystem->Load("libXMLParser");
   //gSystem->Load("libSTEER");
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libCDB");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");

   //lib necessary for dielectron
   gSystem->Load("libCORRFW");

   gSystem->Load("libTender"); 
   gSystem->Load("libProof");
   gSystem->Load("libRAWDatabase");
   gSystem->Load("libSTEER");
   gSystem->Load("libTOFbase");

   gSystem->Load("libTRDbase");
   gSystem->Load("libVZERObase");
   gSystem->Load("libPWGHFhfe");
   //gSystem->Load("libTenderSupplies");

   // Load common libraries


   // Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWGHF/");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWGHF/hfe");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWG/FLOW");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWG/FLOW/Base");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWG/FLOW/Tasks");

   // Create and configure the alien handler plugin
   gROOT->LoadMacro("CreateAlienHandlerAODpp.C");
   AliAnalysisGrid *alienHandler = CreateAlienHandlerAODpp();  
   if (!alienHandler) return;

   // Create the analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("tpctofpp");

   // Connect plug-in to the analysis manager
   mgr->SetGridHandler(alienHandler);
/*
   AliESDInputHandler* esdH = new AliESDInputHandler();
   esdH->SetReadFriends(kFALSE);
   mgr->SetInputEventHandler(esdH);
*/
   // Read AODs
   aodH = new AliAODInputHandler();
   mgr->SetInputEventHandler(aodH);

   // // Read MC info  // not for real data
   // mcHandler = new AliMCEventHandler();
   // mcHandler->SetReadTR(kFALSE);
   // mgr->SetMCtruthEventHandler(mcHandler);

   //AOD Output Hanlder for Filter:
   /*
	  AliAODHandler* aodHandler   = new AliAODHandler();
	  mgr->SetOutputEventHandler(aodHandler);
	  aodHandler->SetOutputFileName("AliAODs.root");
	*/

   // gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
   // AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();



   //==== Physics Selection ====
   /*  gROOT->LoadMacro("AddTask_tender_PhysicsSelection.C");
	   AddTask_tender_PhysicsSelection();

   //==== Add tender ====
   gROOT->LoadMacro("AddTaskTender.C");
   AddTaskTender();
	*/

   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   //AliAnalysisTaskPIDResponse* aodPIDresponse = AddTaskPIDResponse();
   AddTaskPIDResponse();

   gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/AddTaskHFEtpctof.C");
   AliAnalysisTask* anaTask = AddTaskHFEtpctof();



   // Enable debug printouts
   mgr->SetDebugLevel(10);

   if (!mgr->InitAnalysis())
	  return;

   mgr->PrintStatus();
   // Start analysis in grid.
   mgr->StartAnalysis("grid");
};
