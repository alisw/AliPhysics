Bool_t SetupPar(const char *parfile);

void runGridAODPbPb()
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
   gROOT->LoadMacro("CreateAlienHandlerAODPbPb.C");
   AliAnalysisGrid *alienHandler = CreateAlienHandlerAODPbPb();  
   if (!alienHandler) return;

   // Create the analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("tpctofv2");

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


   //===== ADD PID RESPONSE: ===
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   //AliAnalysisTaskPIDResponse* aodPIDresponse = AddTaskPIDResponse();
   AddTaskPIDResponse();

   //===== ADD VZERO event plane: ===
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
   AddTaskVZEROEPSelection();

   //===== ADD TPC event plane: ===
   gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/AddTaskEventPlaneTPC.C");
   AddTaskEventPlaneTPC(kTRUE,0.,kFALSE);

   //===== ADD TASK::
   gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/AddTaskHFEFlowTPCTOFEPSP.C");
   AddTaskHFEFlowTPCTOFEPSP(131073,16,kFALSE,kFALSE,110,60,80,4,2,100,200,30.,50., 1,2,80,kTRUE,kTRUE,0,-2.0,5.0);

   //===== ADD TASK::
   //gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/AddTaskHFEreducedEvent.C");
   //AddTaskHFEreducedEvent();

   //===== ADD TASK::
   //gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/AddTaskHFEnpePbPb.C");
   //AddTaskHFEnpePbPb();


   // Enable debug printouts
   mgr->SetDebugLevel(10);

   if (!mgr->InitAnalysis())
	  return;

   mgr->PrintStatus();
   // Start analysis in grid.
   mgr->StartAnalysis("grid");
};
