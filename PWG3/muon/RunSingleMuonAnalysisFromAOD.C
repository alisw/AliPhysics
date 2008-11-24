//--------------------------------------------------------------------------
// Base macro for submitting single muon analysis.
// 
// In case it is not run with full aliroot, it needs to have in the working directory:
//  - STEERBase.par
//  - ESD.par
//  - AOD.par
//  - ANALYSIS.par
//  - ANALYSISalice.par
//  - PWG3muon.par
//
// The inputPath is either:
//  - The directory containing the AOD file in local mode
//  - The xml file with the list AODs in the alien catalogue in grid mode
//  - The proof dataset in proof mode
// 
// The macro reads AODs and outputs file:
//  - outputDir/singleMuAnalysis.root
//--------------------------------------------------------------------------

enum analysisMode {kMlocal, kMgridInteractive, kMgridBatch, kMproof};
TString modeName[4] = {"local", "local", "grid", "proof"};

void RunSingleMuonAnalysisFromAOD(Int_t mode=kMlocal, Char_t *inputPath=".", Char_t *outputDir=".", Char_t *aodFilename = "AliAODs.root", Long64_t nRuns = -1, Long64_t offset = 0) {
  TStopwatch timer;
  timer.Start();

  // Check if user is running root or aliroot
  TString checkString = gSystem->Getenv("ALICE_ROOT");
  checkString.Append("/lib/tgt_linux/libSTEERBase.so");
  TString foundLib = gSystem->GetLibraries(checkString.Data());
  Bool_t isFullAliroot = (foundLib.Length()==0) ? kFALSE : kTRUE;

  // Load libraries
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  if(mode==kMproof)
    TProof::Open("alicecaf.cern.ch");

  if(isFullAliroot){
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libAOD");
    gSystem->Load("libESD");  
    gSystem->Load("${ALICE_ROOT}/lib/tgt_${ALICE_TARGET}/libPWG3muon.so");
  }
  else {
    const Int_t nNeededPar = 6;
    TString parList[nNeededPar] = {"STEERBase", "ESD", "AOD", "ANALYSIS", "ANALYSISalice", "PWG3muon"};
    if(mode==kMproof){
      gProof->UploadPackage("AF-v4-15");
      gProof->EnablePackage("AF-v4-15");
      if(!SetupPar("PWG3muon")) return;
    }
    else {
      for(Int_t ipar=0; ipar<nNeededPar; ipar++){
	if(!SetupPar(parList[ipar].Data())) return;
      }
    }
  }

  // Connect to alien
  //
  if(mode==kMgridInteractive || mode==kMgridBatch)
    TGrid::Connect("alien://"); 


  TString outFileName("singleMuAnalysis.root");
  outFileName.Prepend(Form("%s/",outputDir));  

  // Get the chain.
  TChain* chain = 0x0;
  if(mode!=kMproof) chain = CreateChain(mode, inputPath, aodFilename);

  //____________________________________________//
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TestManager");
  AliVEventHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  //____________________________________________//
  // Single muon task
  AliAnalysisTaskSingleMu *task1 = new AliAnalysisTaskSingleMu("SingleMu");
  mgr->AddTask(task1);
  // Create containers for input/output
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree1", TTree::Class(),AliAnalysisManager::kOutputContainer,outFileName.Data());
  
  //____________________________________________//
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  if(mode==kMproof)
    mgr->StartAnalysis(modeName[mode].Data(), inputPath, nRuns, offset);
  else 
    mgr->StartAnalysis(modeName[mode].Data(),chain);

  timer.Stop();
  timer.Print();
}


//______________________________________________________________________________
Bool_t SetupPar(char* pararchivename)
{
  if (pararchivename) {
    FileStat_t fs;
    char pararchivenameFull[1024];
    sprintf(pararchivenameFull, "%s.par", pararchivename);
    if(gSystem->GetPathInfo(pararchivenameFull, fs)){
      Error("SetupPar", "PAR Archive %s not found!\nPlease either copy it in the current directory\nor run full aliroot", pararchivenameFull);
      return kFALSE;
    }
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);
    printf("Current directory = %s\n",gSystem->pwd());

    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");

      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
	Error("runProcess","Cannot Build the PAR Archive! - Abort!");
	return kFALSE;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
	
    gSystem->ChangeDirectory("../");
    return kTRUE;
  }
  return kFALSE;
}


//______________________________________________________________________________
TChain* CreateChain(Int_t mode, Char_t* inputPath, Char_t* aodFilename = "AliAOD.root")
{
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
  TChain *chain = 0x0;
  if(mode==kMgridInteractive || mode==kMgridBatch){
    AliTagAnalysis *analysis = new AliTagAnalysis();
    chain = analysis->GetChainFromCollection(inputPath,"aodTree");
  }
  else{
    chain = new TChain("aodTree");
    TString inFileName(aodFilename);
    inFileName.Prepend(Form("%s/",inputPath));
    chain->Add(inFileName);
  }
  //if (chain) chain->ls();
  return chain;
}
