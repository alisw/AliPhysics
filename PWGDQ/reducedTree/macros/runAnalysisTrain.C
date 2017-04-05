// Grid running parameters
TString gGridRunMode = "full";
TString gRootVersion = "v5-34-30-alice7-1";
TString gAlirootVersion = "v5-08-18c-1";
TString gAliphysicsVersion = "vAN-20161208-1";
//TString gGridDataDir = "/alice/data/2015/LHC15o";
TString gGridDataDir = "/alice/data/2016/LHC16q";
//TString gGridDataDir = "/alice/cern.ch/user/i/iarsene/work/outputDst";
//TString gGridDataPattern = "*/pass1/*/AliESDs.root";
//TString gGridDataPattern = "*/pass1/PWGDQ/DQ_PbPb/231_20161009-2048/*/dstTree.root";
TString gGridDataPattern = "*/pass1/*/AliESDs.root";
TString gGridWorkingDir = "work";
TString gGridOutputDir = "testTrees_runBLALA";
Int_t gGridMaxInputFileNumber = 40;

//______________________________________________________________________________________________________________________________________
void runAnalysisTrain(const Char_t* infile, const Char_t* runmode = "local", const Char_t* inputType="ESD", Bool_t hasMC = kFALSE,
                                     Int_t reducedEventType = -1, Bool_t writeTree = kFALSE, TString tasks="dst", TString prod = "LHC10h",
                      Int_t nEntries=1234567890, Int_t firstEntry=0, TString pathForMacros="$ALICE_PHYSICS/PWGDQ/reducedTree/macros")
{
   //
   // infile: list of input files if mode is local, or list of runs for job submission if mode is grid
   //
   // Load common libraries
   gSystem->Load("libPWGDQdielectron.so"); 
   gSystem->Load("libPWGDQreducedTree.so"); 

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
     // Detect which tasks are being run and add the needed grid output files
     TObjArray* arr = tasks.Tokenize(";");
     TString outputFiles = "";
     if(writeTree) outputFiles = "dstTree.root";
     if((!writeTree && arr->GetEntries()>0) ||
        (writeTree && arr->GetEntries()>1)) outputFiles += " dstAnalysisHistograms.root";
        
     alienHandler = CreateAlienHandlerPbPb(infile, gGridRunMode, gGridDataDir, gGridDataPattern, gGridMaxInputFileNumber, 
                                           gGridWorkingDir, gGridOutputDir, outputFiles,
                                           gRootVersion, gAlirootVersion, gAliphysicsVersion);  
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
   if(inputTypeStr.Contains("aod")) {                               // AODs
      inputHandler = new AliAODInputHandler();
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
   if(hasMC && (inputTypeStr.Contains("esd") || inputTypeStr.Contains("aod"))) {
      mc = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mc);
   }
   
   //==== Add tender ====
   //   gROOT->LoadMacro("AddTaskTender.C");
   //   AddTaskTender();

   if(inputTypeStr.Contains("esd") || inputTypeStr.Contains("aod")) {         // no need if we run over reduced events
      //==== Physics Selection ====
      gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
      Bool_t applyPileupCuts = kTRUE;
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(hasMC, applyPileupCuts);

      //===== ADD CENTRALITY: ===
      if(!prod.CompareTo("LHC10h") || !prod.CompareTo("LHC11h")) {         // Run-1 Pb-Pb
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
        AddTaskCentrality();
      }
      else  { // Run-2
      //if(!prod.CompareTo("LHC15o") || !prod.CompareTo("LHC16l") || !prod.CompareTo("LHC16q") || !prod.CompareTo("LHC16t")) {         // Run-2 Pb-Pb
         gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
         AliMultSelectionTask* multTask = AddTaskMultSelection();
         if(hasMC && !prod.CompareTo("LHC15o"))
           multTask->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
      }      

      //===== ADD PID RESPONSE: ===
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
      Bool_t tuneOnData = kTRUE;
      Int_t recoPass = 1;
      if(hasMC) AddTaskPIDResponse(hasMC, kTRUE, tuneOnData, recoPass);
      else AddTaskPIDResponse();
   }
   
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/AddTask_TrainTreeAnalysis.C");
   AddTask_TrainTreeAnalysis((!runmodestr.CompareTo("grid") ? kTRUE : kFALSE), prod, reducedEventType, writeTree, tasks, pathForMacros);
   
   // Enable debug printouts
   //mgr->SetDebugLevel(10);

   if (!mgr->InitAnalysis()) return;

   TChain* chain = NULL;
   if(!runmodestr.Contains("grid"))
      chain = makeChain(infile, inputType);
   
   TProof* proof=0x0;
   if(runmodestr.Contains("proof")) {
      proof = TProof::Open("");
      chain->SetProof();
   }
   
   mgr->PrintStatus();
   // Start analysis
   if(runmodestr.Contains("local"))
      mgr->StartAnalysis("local", chain, nEntries, firstEntry);
   if(runmodestr.Contains("proof"))
      mgr->StartAnalysis("proof", chain, nEntries, firstEntry);
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
  if(itStr.Contains("aod"))
     chain=new TChain("aodTree");
  
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
  cout << "Number of events in chain: " << chain->GetEntries() << endl;
  return chain;
}
