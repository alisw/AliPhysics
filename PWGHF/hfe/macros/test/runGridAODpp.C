Bool_t SetupPar(const char *parfile);

void runGridAODpp()
{
   // Load common libraries
  gSystem->Load("libCore.so");  
   gSystem->Load("libTree.so");
   gSystem->Load("libGeom.so");
   gSystem->Load("libVMC.so");
   gSystem->Load("libPhysics.so");
   gSystem->Load("libMinuit.so"); 
   gSystem->Load("libGui.so");
   gSystem->Load("libXMLParser.so");
   //gSystem->Load("libSTEER.so");
   gSystem->Load("libSTEERBase.so");
   gSystem->Load("libESD.so");
   gSystem->Load("libCDB.so");
   gSystem->Load("libAOD.so");
   gSystem->Load("libANALYSIS.so");
   gSystem->Load("libANALYSISalice.so");

   //lib necessary for dielectron
   gSystem->Load("libCORRFW.so");

   gSystem->Load("libTender"); 
   gSystem->Load("libProof.so");
   gSystem->Load("libRAWDatabase.so");
   gSystem->Load("libSTEER.so");
   gSystem->Load("libTOFbase.so");

   gSystem->Load("libTRDbase.so");
   gSystem->Load("libVZERObase.so");
   gSystem->Load("libPWGHFhfe.so"); 
   //gSystem->Load("libTenderSupplies.so"); 

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

   // gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
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
