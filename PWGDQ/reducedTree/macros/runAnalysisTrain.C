// Grid running parameters
TString gGridRunMode = "full";
TString gRootVersion = "v5-34-30-alice7-6";
TString gAlirootVersion = "v5-09-03a-1";
TString gAliphysicsVersion = "vAN-20170804-1";
//TString gGridDataDir = "/alice/data/2015/LHC15o";
//TString gGridDataDir = "/alice/data/2016/LHC16q";
TString gGridDataDir = "/alice/data/2016/LHC16s";
//TString gGridDataDir = "/alice/sim/2017/LHC17d13_cent/";
//TString gGridDataDir = "/alice/cern.ch/user/i/iarsene/work/outputDst";
//TString gGridDataPattern = "*/pass1/*/AliESDs.root";
//TString gGridDataPattern = "*/pass1/PWGDQ/DQ_PbPb/231_20161009-2048/*/dstTree.root";
TString gGridDataPattern = "*/pass1_FAST/*/AliESDs.root";
//TString gGridDataPattern = "*/AOD/*/AliAOD.root";
//TString gGridDataPattern = "*/pass1_CENT_wSDD/AOD/*/AliAOD.root";
TString gGridWorkingDir = "work";
TString gGridOutputDir = "20170805_dstTrees_LHC16rs";
Int_t gGridMaxInputFileNumber = 50;

TChain* makeChain(const Char_t* filename, const Char_t* inputType);

//______________________________________________________________________________________________________________________________________
void runAnalysisTrain(const Char_t* infile, const Char_t* runmode = "local", const Char_t* inputType="ESD", Bool_t hasMC = kFALSE,
                                     Int_t reducedEventType = -1, Bool_t writeTree = kFALSE, TString tasks="dst", TString prod = "LHC10h",
                      Int_t nEntries=-1, Int_t firstEntry=0, TString pathForMacros="$ALICE_PHYSICS/PWGDQ/reducedTree/macros")
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
     // Detect which tasks are being run and add the needed grid output files
     TObjArray* arr = tasks.Tokenize(";");
     TString outputFiles = "";
     if(writeTree) outputFiles = "dstTree.root";
     if((!writeTree && arr->GetEntries()>0) ||
        (writeTree && arr->GetEntries()>1)) outputFiles += " dstAnalysisHistograms.root";
#ifdef __CLING__
     // ROOT6 version
     std::stringstream creatalienhandleradd;
     creatalienhandleradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWGDQ/reducedTree/macros/CreateAlienHandler.C(";
     creatalienhandleradd << "\"" << infile << "\", ";
     creatalienhandleradd << "\"" << gGridRunMode.Data() << "\", ";
     creatalienhandleradd << "\"" << gGridDataDir.Data() << "\", ";
     creatalienhandleradd << "\"" << gGridDataPattern.Data() << "\", ";
     creatalienhandleradd         << gGridMaxInputFileNumber << ", ";
     creatalienhandleradd << "\"" << gGridWorkingDir.Data() << "\", ";
     creatalienhandleradd << "\"" << gGridOutputDir.Data() << "\", ";
     creatalienhandleradd << "\"" << outputFiles.Data() << "\", ";
     creatalienhandleradd << "\"" << gRootVersion.Data() << "\", ";
     creatalienhandleradd << "\"" << gAlirootVersion.Data() << "\", ";
     creatalienhandleradd << "\"" << gAliphysicsVersion.Data() << "\")";
     std::string creatalienhandleraddstr = creatalienhandleradd.str();
     std::cout << "Calling Add macro using command string: " << creatalienhandleraddstr << std::endl;
     alienHandler = (AliAnalysisGrid*)gROOT->ProcessLine(creatalienhandleraddstr.c_str());
#else
     // ROOT5 version
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/CreateAlienHandler.C");
     alienHandler = CreateAlienHandlerPbPb(infile, gGridRunMode, gGridDataDir, gGridDataPattern, gGridMaxInputFileNumber,
                                           gGridWorkingDir, gGridOutputDir, outputFiles,
                                           gRootVersion, gAlirootVersion, gAliphysicsVersion);
#endif
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
     Bool_t applyPileupCuts = kTRUE;
#ifdef __CLING__
     // ROOT6 version
     std::stringstream physseladd;
     physseladd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/OADB/macros/AddTaskPhysicsSelection.C(";
     physseladd << hasMC << ", ";
     physseladd << applyPileupCuts << ")";
     std::string physseladdstr = physseladd.str();
     std::cout << "Calling Add macro using command string: " << physseladdstr << std::endl;
     AliPhysicsSelectionTask* physSelTask = (AliPhysicsSelectionTask*)gROOT->ProcessLine(physseladdstr.c_str());
#else
     // ROOT5 version
     gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
     AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(hasMC, applyPileupCuts);
#endif
      //===== ADD CENTRALITY: ===
      if(!prod.CompareTo("LHC10h") || !prod.CompareTo("LHC11h")) {         // Run-1 Pb-Pb
#ifdef __CLING__
        // ROOT6 version
        std::stringstream centradd;
        centradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/OADB/macros/AddTaskCentrality.C";
        std::string centraddstr = centradd.str();
        std::cout << "Calling Add macro using command string: " << centraddstr << std::endl;
        gROOT->ProcessLine(centraddstr.c_str());
#else
        // ROOT5 version
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
        AddTaskCentrality();
#endif
      }
      else  { // Run-2
      //if(!prod.CompareTo("LHC15o") || !prod.CompareTo("LHC16l") || !prod.CompareTo("LHC16q") || !prod.CompareTo("LHC16t")) {         // Run-2 Pb-Pb
#ifdef __CLING__
        // ROOT6 version
        std::stringstream multsel;
        multsel << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C";
        std::string multselstr = multsel.str();
        std::cout << "Calling Add macro using command string: " << multselstr << std::endl;
        AliMultSelectionTask* multTask = (AliMultSelectionTask*)gROOT->ProcessLine(multselstr.c_str());
#else
        // ROOT5 version
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
        AliMultSelectionTask* multTask = AddTaskMultSelection();
#endif
         if(hasMC && !prod.CompareTo("LHC15o"))
           multTask->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
         if(hasMC && (!prod.CompareTo("LHC16q") || !prod.CompareTo("LHC16t")))
            multTask->SetAlternateOADBforEstimators("LHC16q-DefaultMC-HIJING");
         if(hasMC && (!prod.CompareTo("LHC16r")))
            multTask->SetAlternateOADBforEstimators("LHC16r-DefaultMC-EPOSLHC");
         if(hasMC && (!prod.CompareTo("LHC16s")))
            multTask->SetAlternateOADBforEstimators("LHC16s-DefaultMC-EPOSLHC");
      }      

     //===== ADD PID RESPONSE: ===
     Bool_t tuneOnData = kTRUE;
     TString recoPass = "1";
#ifdef __CLING__
     // ROOT6 version
     std::stringstream pidresp;
     if (hasMC) {
       pidresp << ".x " << gSystem->Getenv("ALICE_ROOT") << "/ANALYSIS/macros/AddTaskPIDResponse.C(";
       pidresp << hasMC << ", ";
       pidresp << kTRUE << ", ";
       pidresp << tuneOnData << ", ";
       pidresp << recoPass << ")";
     } else {
       pidresp << ".x " << gSystem->Getenv("ALICE_ROOT") << "/ANALYSIS/macros/AddTaskPIDResponse.C";
     }
     std::string pidrespstr = pidresp.str();
     std::cout << "Calling Add macro using command string: " << pidrespstr << std::endl;
     gROOT->ProcessLine(pidrespstr.c_str());
#else
     // ROOT5 version
     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
     if(hasMC) AddTaskPIDResponse(hasMC, kTRUE, tuneOnData, recoPass);
     else AddTaskPIDResponse();
#endif
   }

#ifdef __CLING__
  // ROOT6 version
  std::stringstream trainadd;
  trainadd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWGDQ/reducedTree/macros/AddTask_TrainTreeAnalysis.C(";
  trainadd << (!runmodestr.CompareTo("grid") ? kTRUE : kFALSE) << ", ";
  trainadd << "\"" << prod.Data() << "\", ";
  trainadd << reducedEventType << ", ";
  trainadd << writeTree << ", ";
  trainadd << "\"" << tasks.Data() << "\", ";
  trainadd << "\"" << pathForMacros.Data() << "\")";
  std::string trainaddstr = trainadd.str();
  std::cout << "Calling Add macro using command string: " << trainaddstr << std::endl;
  gROOT->ProcessLine(trainaddstr.c_str());
#else
  // ROOT5 version
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/AddTask_TrainTreeAnalysis.C");
  AddTask_TrainTreeAnalysis((!runmodestr.CompareTo("grid") ? kTRUE : kFALSE), prod, reducedEventType, writeTree, tasks, pathForMacros);
#endif

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
   if(nEntries==-1) nEntries=1234567890;
   if(runmodestr.Contains("local"))
      mgr->StartAnalysis("local", chain, nEntries, firstEntry);
   if(runmodestr.Contains("proof"))
      mgr->StartAnalysis("proof", chain, nEntries, firstEntry);
   if(runmodestr.Contains("grid"))
      mgr->StartAnalysis("grid", nEntries, firstEntry);

  // piping additional histograms from tree file into analysis output
  TObjArray* taskArr = tasks.Tokenize(";");
  if ((!writeTree && taskArr->GetEntries()>0) ||
      (writeTree && taskArr->GetEntries()>1)) {
    ifstream in;
    in.open(infile);
    TObjArray* histArr = new TObjArray();

    // loop over keys in file
    TString line;
    Int_t nFile = 0;
    while(in.good()) {
      in >> line;
      if (!line.IsNull()) {
        nFile++;
        TFile* tmpFile = TFile::Open(Form("%s",line.Data()));
        cout << "Looking for histograms in " << tmpFile->GetName() << ": " << endl;
        TIter next(tmpFile->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
          TClass *cl = gROOT->GetClass(key->GetClassName());
          if (!cl->InheritsFrom("TH1")) continue;
          TH1 *h = (TH1*)key->ReadObj();
          cout << " " << h->GetName();
          if (nFile==1) {
            cout << " -> added" << endl;
            histArr->Add(h);
          } else {
            if (histArr->FindObject(h->GetName())) {
              cout << " -> merged" << endl;
              ((TH1*)histArr->FindObject(h->GetName()))->Add(h);
            } else {
              cout << " -> added" << endl;
              histArr->Add(h);
            }
          }
        }
      }
    }

    // write histograms to output file
    TObjArray* outputs = mgr->GetOutputs();
    for (Int_t i=0; i<outputs->GetEntries(); i++) {
      TFile* tmpOut = (TFile*)mgr->OpenFile((AliAnalysisDataContainer*)outputs->At(i), "UPDATE", kTRUE);
      if (!tmpOut) continue;
      TString tmpOutName = tmpOut->GetName();
      if (tmpOutName.Contains("dstTree")) continue;
      if (!histArr->GetEntries()) continue;
      cout << "Writing histograms to " << tmpOut->GetName() << endl;
      for (Int_t j=0; j<histArr->GetEntries(); j++) histArr->At(j)->Write();
      tmpOut->Close();
    }
  }
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
