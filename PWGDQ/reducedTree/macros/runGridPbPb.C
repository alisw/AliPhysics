//Replace User Task with your Add Task and appropriate parameters

void runGridPbPb()
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
   gSystem->Load("libSTEERBase.so");
   gSystem->Load("libESD.so");
   gSystem->Load("libCDB.so");
   gSystem->Load("libAOD.so");
   gSystem->Load("libANALYSIS.so");
   gSystem->Load("libANALYSISalice.so");

   //lib necessary for dielectron
   gSystem->Load("libCORRFW.so");

   gSystem->Load("libTENDER"); 
   gSystem->Load("libTENDERSupplies"); 
   gSystem->Load("libProof.so");
   gSystem->Load("libRAWDatabase.so");
   gSystem->Load("libSTEER.so");
   gSystem->Load("libTOFbase.so");

   gSystem->Load("libTRDbase.so");
   gSystem->Load("libVZERObase.so");
   gSystem->Load("libPWGDQdielectron.so"); 
   //gSystem->Load("/hera/alice/iarsene/trainTrunk/util/dielectron/libPWGDQdielectron.so");
   gSystem->Load("libTENDERSupplies.so"); 

   // Load common libraries

   // Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/PWGDQ/");

   gROOT->ProcessLine(".include $ALICE_ROOT/PWGDQ/dielectron");
   //gROOT->ProcessLine(".include /hera/alice/iarsene/trainTrunk/util/dielectron/dielectron");

   // Create and configure the alien handler plugin
   gROOT->LoadMacro("CreateAlienHandlerPbPb.C");
   AliAnalysisGrid *alienHandler = CreateAlienHandlerPbPb();  
   if (!alienHandler) return;

   // Create the analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("jpsiAnalysis");

   // Connect plug-in to the analysis manager
   mgr->SetGridHandler(alienHandler);

   AliESDInputHandler* esdH = new AliESDInputHandler();
   esdH->SetReadFriends(kFALSE);
   mgr->SetInputEventHandler(esdH);

   //===== ADD PID RESPONSE: ===
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
   AddTaskPIDResponse();

   //==== Add tender ====
   //   gROOT->LoadMacro("AddTaskTender.C");
   //   AddTaskTender();

   //==== Physics Selection ====
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

   //===== ADD CENTRALITY: ===
   gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
   AddTaskCentrality();
   
   //===== ADD TASK::
   //  gROOT->LoadMacro("$ALICE_ROOT/PWGDQ/dielectron/macrosJPSI/AddTask_ReducedTree.C");
   //gROOT->LoadMacro("AddTask_ReducedTree.C");
   gROOT->LoadMacro("AddTask_iarsene_dst.C");
   //AddTask_ReducedTree();   
   AddTask_iarsene_dst();


   // Enable debug printouts
   mgr->SetDebugLevel(10);

   if (!mgr->InitAnalysis())
	  return;

   mgr->PrintStatus();
   // Start analysis in grid.
   mgr->StartAnalysis("grid");
};
