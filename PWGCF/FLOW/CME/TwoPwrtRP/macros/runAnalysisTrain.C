//Replace User Task with your Add Task and appropriate parameters

#include<TSystem.h>
#include<TROOT.h>

void runAnalysisTrain(const Char_t* infile="")
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

   gSystem->Load("libTender");
   gSystem->Load("libTenderSupplies");
   gSystem->Load("libProof.so");
   gSystem->Load("libRAWDatabase.so");
   gSystem->Load("libSTEER.so");
   gSystem->Load("libTOFbase.so");

   gSystem->Load("libTRDbase.so");
   gSystem->Load("libVZERObase.so");
   gSystem->Load("libPWGDQdielectron.so");
   gSystem->Load("libPWGDQreducedTree.so");
   gSystem->Load("libPWGPPevcharQn.so");
   gSystem->Load("libPWGCFflowCMETwoPwrtRP.so");
   //gSystem->Load("/hera/alice/iarsene/trainTrunk/util/dielectron/libPWGDQdielectron.so");

  // Load common libraries

   // Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGDQ/");
   gROOT->ProcessLine(".include $ROOTSYS/include");

   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGDQ/dielectron");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGDQ/ReducedTreeMaker");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGPP/EVCHAR/FlowEventCorrections/QnCorrections/");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGCF/FLOW");
   gROOT->ProcessLine(".include $ALICE_PHYSICS/PWGCF/FLOW/CME/TwoPwrtRP");
   //gROOT->ProcessLine(".include /hera/alice/iarsene/trainTrunk/util/dielectron/dielectron");

   // Create and configure the alien handler plugin
   gROOT->LoadMacro("CreateAlienHandler.C");
   AliAnalysisGrid *alienHandler = CreateAlienHandler();
   if (!alienHandler) return;

   // Create the analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("ReducedTreeMaker");

   // Connect plug-in to the analysis manager
   //mgr->SetGridHandler(alienHandler);

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
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    //AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE, kTRUE, AliVEvent::kMB)
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

   //===== ADD CENTRALITY: ===
   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
   AliMultSelectionTask * multtask = AddTaskMultSelection(kFALSE);
   multtask->SetSelectedTriggerClass(AliVEvent::kMB);
   //task->SetUseDefaultCalib(kTRUE); // data
   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
   AddTaskCentrality();

   //===== ADD TASK::
   //  gROOT->LoadMacro("$ALICE_ROOT/PWGDQ/dielectron/macrosJPSI/AddTask_ReducedTree.C");
   //gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/dielectron/ReducedTreeMaker/macros/AddTask_iarsene_dst.C");
   gROOT->LoadMacro("AddTask_jonderw_TwoPwrtRP.C");
   //AddTask_ReducedTree();
   //AddTask_iarsene_dst();
   AddTask_jonderw_TwoPwrtRP(kTRUE);


   // Enable debug printouts
   mgr->SetDebugLevel(10);

   if (!mgr->InitAnalysis())
	  return;

   TChain* chain = makeChain(infile);

   mgr->PrintStatus();
   // Start analysis in grid.
   mgr->StartAnalysis("local",chain);
   //mgr->StartAnalysis("grid");
};

//_______________________________________________________________________________
TChain* makeChain(const Char_t* filename) {
  //
  // make a chain using the trees from the list of filenames in "filename"
  //
  TChain *chain=new TChain("esdTree");
  ifstream in;
  in.open(filename);

  TString line;
  //loop over file
  while(in.good()) {
    in >> line;
    if (!line.IsNull()) {
      cout << "Adding file: " << line <<endl;
      chain->AddFile(line);
    }
  }
  return chain;
}
