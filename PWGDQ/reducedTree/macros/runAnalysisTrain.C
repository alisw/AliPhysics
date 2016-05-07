//Replace User Task with your Add Task and appropriate parameters

#include<TSystem.h>
#include<TROOT.h>

void runAnalysisTrain(const Char_t* infile, const Char_t* runmode = "local", const Char_t* inputType="ESD", Bool_t hasMC = kFALSE,
                                     Int_t reducedEventType = -1, Bool_t writeTree = kFALSE,
                                     Int_t nEntries=1234567890, Int_t firstEntry=0)
{
   //
   // infile: list of input files if mode is local, or list of runs for job submission if mode is grid
   //
   // Load common libraries
/*   gSystem->Load("libCore.so");  
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
   gSystem->Load("libVZERObase.so");*/
   gSystem->Load("libPWGDQdielectron.so"); 
   gSystem->Load("libPWGDQreducedTree.so"); 
   //gSystem->Load("/hera/alice/iarsene/trainTrunk/util/dielectron/libPWGDQdielectron.so");

  // Load common libraries

   // Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ROOTSYS/include");

   // Create the analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("ReducedTreeAnalysis");
   
   // Create and configure the alien handler plugin
   TString runmodestr(runmode); runmodestr.ToLower();
   AliAnalysisGrid *alienHandler = NULL;
   if(runmodestr.Contains("grid")) {
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/CreateAlienHandler.C");
     alienHandler = CreateAlienHandlerPbPb(infile);  
     if (!alienHandler) {
        cout << "runAnalysisTrain.C ::      Could not create the alien handler. Check it out!" << endl;
        return;
     }
     mgr->SetGridHandler(alienHandler);
   }

   // Create the input handler
   TString inputTypeStr(inputType); inputTypeStr.ToLower();
   AliInputEventHandler* inputHandler = NULL; 
   if(inputTypeStr.Contains("esd")) {                               // ESDs
     inputHandler = new AliESDInputHandler();
     ((AliESDInputHandler*)inputHandler)->SetReadFriends(kFALSE);
   }
   if(inputTypeStr.Contains("reducedevent")) {              // AliReducedEventInfo
      inputHandler = new AliReducedEventInputHandler();
      ((AliReducedEventInputHandler*)inputHandler)->SetInputEventType(AliReducedEventInputHandler::kReducedEventInfo);
   }
   if(inputTypeStr.Contains("baseevent")) {                    // AliReducedBaseEvent
      inputHandler = new AliReducedEventInputHandler();
      ((AliReducedEventInputHandler*)inputHandler)->SetInputEventType(AliReducedEventInputHandler::kReducedBaseEvent);
   }
   mgr->SetInputEventHandler(inputHandler);
   
   // Add the MC handler if needed
   AliMCEventHandler *mc = NULL;
   if(hasMC) {
      mc = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mc);
   }
   
   //==== Add tender ====
   //   gROOT->LoadMacro("AddTaskTender.C");
   //   AddTaskTender();

   if(inputTypeStr.Contains("esd")) {         // no need if we run over reduced events
      //==== Physics Selection ====
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();

      //===== ADD CENTRALITY: ===
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
      AddTaskCentrality();
      /* // For Run-2 Pb-Pb there is a new way to get centrality. Below are the lines of code presented by David Chinellato at the PF in 11 december 2015
       gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
       AddTaskMultSelection();
       // we also need to add something like the lines below, in the analysis task:
       // AliMultSelection* MultSelection = (AliMultSelection*)esdEvent->FindListObject("MultSelection");
       // Float_t lPerc = 300;  // nonsense
       // if(MultSelection) {
       //   lPerc = MultSelection->GetMultiplicityPercentile("V0M");
       //   // Quality check
       //   Int_t lEvSelCode = MultSelection->GetEvSelCode();
       //   if(lEvSelCode>0) lPerc = lEvSelCode; // disregard!
       // }
       // else {
       //   // If this happens, re-check if AliMultSelectionTask ran before your task!
       //   AliInfo("Didn't find MultSelection!");
       // }
       */

      //===== ADD PID RESPONSE: ===
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
      AddTaskPIDResponse();
   }
   
   //===== Add the reduced event producer task, if needed
   if(inputTypeStr.Contains("esd")) {         // if we run over reduced events this not needed anymore
      gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/AddTask_iarsene_dst.C");
      AddTask_iarsene_dst(reducedEventType, writeTree);
   }
   //===== Add consumer tasks for the reduced events
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/AddTask_iarsene_testTask.C");
   if(inputTypeStr.Contains("esd"))
      AddTask_iarsene_testTask(kTRUE, AliAnalysisTaskReducedEventProcessor::kUseOnTheFlyReducedEvents);
   if(inputTypeStr.Contains("reducedevent") || inputTypeStr.Contains("baseevent"))
      AddTask_iarsene_testTask(kTRUE, AliAnalysisTaskReducedEventProcessor::kUseEventsFromTree);

   // Enable debug printouts
   //mgr->SetDebugLevel(10);

   if (!mgr->InitAnalysis()) return;

   TChain* chain = NULL;
   if(!runmodestr.Contains("grid"))
      chain = makeChain(infile, inputType);
   
   mgr->PrintStatus();
   // Start analysis
   if(runmodestr.Contains("local"))
      mgr->StartAnalysis("local", chain, nEntries, firstEntry);
   if(runmodestr.Contains("grid"))
      mgr->StartAnalysis("grid", nEntries, firstEntry);
};

//_______________________________________________________________________________
TChain* makeChain(const Char_t* filename, const Char_t* inputType) {
  //
  // make a chain using the trees from the list of filenames in "filename"
  //
  TString itStr(inputType); itStr.ToLower();
  TChain *chain = NULL;
  if(itStr.Contains("reducedevent") || itStr.Contains("baseevent"))
    chain=new TChain("DstTree");
  if(itStr.Contains("esd"))
    chain=new TChain("esdTree");
  
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
