//Replace User Task with your Add Task and appropriate parameters

Bool_t SetupPar(const char *parfile);

void runGridPbPb()
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
   gSystem->Load("libSTEERBase");
   gSystem->Load("libESD");
   gSystem->Load("libCDB");
   gSystem->Load("libAOD");
   gSystem->Load("libANALYSIS");
   gSystem->Load("libANALYSISalice");

   //lib necessary for dielectron
   gSystem->Load("libCORRFW");

   gSystem->Load("libTender"); 
   gSystem->Load("libTenderSupplies"); 
   gSystem->Load("libProof");
   gSystem->Load("libRAWDatabase");
   gSystem->Load("libSTEER");
   gSystem->Load("libTOFbase");

   gSystem->Load("libTRDbase");
   gSystem->Load("libVZERObase");
   gSystem->Load("libPWGDQdielectron"); 
   gSystem->Load("libTenderSupplies"); 

   // Load common libraries

   // Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWGDQ/");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWGDQ/dielectron");

   // Create and configure the alien handler plugin
   gROOT->LoadMacro("CreateAlienHandlerPbPb.C");
   AliAnalysisGrid *alienHandler = CreateAlienHandlerPbPb();  
   if (!alienHandler) return;

   // Create the analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("lowmassAnalysis");

   // Connect plug-in to the analysis manager
   mgr->SetGridHandler(alienHandler);

   AliESDInputHandler* esdH = new AliESDInputHandler();
   esdH->SetReadFriends(kFALSE);
   mgr->SetInputEventHandler(esdH);



   //==== Physics Selection ====
    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

   //==== Add tender ====

//   gROOT->LoadMacro("AddTaskTender.C");
//   AddTaskTender();

   //===== ADD PID RESPONSE: ===

   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   AddTaskPIDResponse();

   //===== ADD CENTRALITY: ===
   gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskCentrality.C");
   AddTaskCentrality();
   

   //===== ADD TASK::

   gROOT->LoadMacro("$ALICE_ROOT/PWGDQ/dielectron/macrosLMEE/AddTask_cbaumann_LMEEPbPb2011.C");
   AddTask_cbaumann_LMEEPbPb2011(kFALSE,kFALSE,kTRUE);



   // Enable debug printouts
   mgr->SetDebugLevel(10);

   if (!mgr->InitAnalysis())
	  return;

   mgr->PrintStatus();
   // Start analysis in grid.
   mgr->StartAnalysis("grid");
};
