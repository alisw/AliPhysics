//--------------------------------------------------------------------------
// Base macro for submitting muon QA analysis.
//
// It needs the following libraries:
//  - libSTEERBase
//  - libESD
//  - libAOD
//  - libANALYSIS
//  - libANALYSISalice
//  - libCORRFW
//  - libPWGmuon
//  - libPWGPPMUONlite
// The macro reads ESDs and store outputs in standard output file (AnalysisResults.root)
//
// Author: Philippe Pillot / Cynthia Hadjidakis 
//--------------------------------------------------------------------------

enum {kLocal, kInteractif_xml, kInteractif_ESDList};

void RunMuonQA(TString inputFileName = "AliESDs.root", Bool_t isMC = kFALSE, 
	       Bool_t selectEvent = kTRUE, Bool_t selectMatched = kTRUE, 
	       Bool_t applyAccCut = kTRUE, Short_t selectCharge = 0)
{
  TStopwatch timer;
  timer.Start();
  
  // Check runing mode
  Int_t mode = GetMode(inputFileName);
  if(mode < 0){
    Error("RunMuonQA","Please provide either an ESD root file or a collection of ESDs.");
    return;
  }
  
  // Load common libraries
  //gSystem->Load("libVMC");
  //gSystem->Load("libTree");
  //gSystem->Load("libPhysics");
  //gSystem->Load("libMinuit");
  //gSystem->Load("libXMLParser");
  //gSystem->Load("libGui");
  //gSystem->Load("libSTEERBase");
  //gSystem->Load("libESD");
  //gSystem->Load("libAOD");
  //gSystem->Load("libANALYSIS");
  //gSystem->Load("libANALYSISalice");
  //gSystem->Load("libCORRFW");
  //gSystem->Load("libPWGmuon");
  //gSystem->Load("libPWGPPMUONlite");
  
  // Create input chain
  TChain* chain = CreateChain(inputFileName);
  if (!chain) return;
  
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("MuonQAAnalysis");
  
  // ESD input handler
  AliESDInputHandler* esdH = new AliESDInputHandler();
  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // event selection
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelection = AddTaskPhysicsSelection(isMC);
  if(!physicsSelection) {
    Error("RunMuonQA","AliPhysicsSelectionTask not created!");
    return;
  }
  
  // Muon QA analysis
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskMuonQA.C");

  // Create and configure task
  AliAnalysisTaskMuonQA *muonQA = AddTaskMuonQA(selectEvent,selectMatched,applyAccCut,selectCharge);
  if (!muonQA) {
    Error("RunMuonQA", "Muon QA task not created!");
    return;
  }

  // Enable debug printouts
  // mgr->SetDebugLevel(2);
  
  // start local analysis
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }
  
  timer.Stop();
  timer.Print();
}

//______________________________________________________________________________
Int_t GetMode(TString inputFileName)
{
  if ( inputFileName.EndsWith(".xml") ) return kInteractif_xml;
  else if ( inputFileName.EndsWith(".txt") ) return kInteractif_ESDList;
  else if ( inputFileName.EndsWith(".root") ) return kLocal;
  return -1;
}

//______________________________________________________________________________
TChain* CreateChainFromCollection(const char *xmlfile)
{
  // Create a chain from the collection of tags.
  TGridCollection* coll = gGrid->OpenCollection(xmlfile);
  if (!coll) {
    ::Error("CreateChainFromTags", "Cannot create an AliEn collection from %s", xmlfile);
    return NULL;
  }
  
  TGridResult* tagResult = coll->GetGridResult("",kFALSE,kFALSE);
  AliTagAnalysis *tagAna = new AliTagAnalysis("ESD");
  tagAna->ChainGridTags(tagResult);
  
  AliRunTagCuts      *runCuts = new AliRunTagCuts();
  AliLHCTagCuts      *lhcCuts = new AliLHCTagCuts();
  AliDetectorTagCuts *detCuts = new AliDetectorTagCuts();
  AliEventTagCuts    *evCuts  = new AliEventTagCuts();
  
  // Check if the cuts configuration file was provided
  if (!gSystem->AccessPathName("ConfigureCuts.C")) {
    gROOT->LoadMacro("ConfigureCuts.C");
    ConfigureCuts(runCuts, lhcCuts, detCuts, evCuts);
  }
  
  TChain *chain = tagAna->QueryTags(runCuts, lhcCuts, detCuts, evCuts);
  if (!chain || !chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChainFromFile(const char *rootfile)
{
  // Create a chain using the root file.
  TChain* chain = new TChain("esdTree");
  chain->Add(rootfile);
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChainFromESDList(const char *esdList)
{
  // Create a chain using tags from the run list.
  TChain* chain = new TChain("esdTree");
  ifstream inFile(esdList);
  TString inFileName;
  if (inFile.is_open()) {
    while (! inFile.eof() ) {
      inFileName.ReadLine(inFile,kFALSE);
      if(!inFileName.EndsWith(".root")) continue;
      chain->Add(inFileName.Data());
    }
  }
  inFile.close();
  if (!chain->GetNtrees()) return NULL;
  chain->ls();
  return chain;
}

//______________________________________________________________________________
TChain* CreateChain(TString inputFileName)
{
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
  Int_t mode = GetMode(inputFileName);
  if(mode == kInteractif_xml) return CreateChainFromCollection(inputFileName.Data());
  else if (mode == kInteractif_ESDList) return CreateChainFromESDList(inputFileName.Data());
  else if (mode == kLocal) return CreateChainFromFile(inputFileName.Data());
  else return NULL;
}

