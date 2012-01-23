//--------------------------------------------------------------------------
// Base macro for submitting trigger chamber efficiency determination.
// 
// In case it is not run with full aliroot, it needs to have in the working directory:
//  - STEERBase.par
//  - ESD.par
//  - AOD.par
//  - ANALYSIS.par
//  - ANALYSISalice.par
//  - PWG3muon.par
// 
// The macro reads ESDs and outputs file:
// - MUON.TriggerEfficiencyMap.root
//
// To display the trigger chamber efficiency:
// > aliroot
// > AliMUONTriggerEfficiencyCells effCells("MUON.TriggerEfficiencyMap.root")
// > effCells.DisplayEfficiency()
//--------------------------------------------------------------------------

enum anaModes {kMlocal, kMgridLocal, kMgrid};
TString modeName[3] = {"local", "local", "grid"};
const TString kDefaultLocalInputDir="$ALICE_ROOT/MUON/test_out.100";
const TString kDefaultXmlFile="wn.xml";

void AnalysisTrigChEff(Int_t mode=kMlocal)
{
// Example of running analysis train
  TStopwatch timer;
  timer.Start();

  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");

  // if root setup Par
  TString checkString=gSystem->Getenv("ALICE_ROOT");
  checkString.Append("/lib/tgt_linux/libMUONbase.so");
  TString foundLib=gSystem->GetLibraries(checkString.Data());

  if(foundLib.Length()==0){
    // Common packages
    SetupPar("STEERBase");
    SetupPar("ESD");
    SetupPar("AOD");
    SetupPar("ANALYSIS");
    SetupPar("ANALYSISalice");
  }
  SetupPar("PWG3muon");

  // Analysis using standard AliRoot libraries
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG3muon.so");

  //
  // Connect to alien
  //
  if(mode==kMgridLocal || mode==kMgrid)
    TGrid::Connect("alien://"); 

  //
  // Create the chain
  //
  TChain* chain = CreateChain(mode);

  TString outFileName("MUON.TriggerEfficiencyMap.root");
  //
  //
  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("Analysis Train", "A test setup for the analysis train");
  // ESD input handler
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  //esdHandler->SetInactiveBranches("FMD CaloCluster");
  // Monte Carlo handler
  //AliMCEventHandler* mcHandler = new AliMCEventHandler();

  mgr->SetInputEventHandler(esdHandler);
  //mgr->SetMCtruthEventHandler(mcHandler);
  
  // Trigger chamber efficiency analysis
  AliAnalysisTaskTrigChEff* taskTrigChEff = new AliAnalysisTaskTrigChEff("TaskTrigChEff");
  mgr->AddTask(taskTrigChEff);
  
  // Create containers for input/output
  // Top container for ESD input 
  AliAnalysisDataContainer *cin_esd = mgr->GetCommonInputContainer();

  // Output histograms list for single muons analysis
  AliAnalysisDataContainer *cout_trigChEff = mgr->CreateContainer("triggerChamberEff", TList::Class(),
								 AliAnalysisManager::kOutputContainer, outFileName.Data());

  mgr->ConnectInput  (taskTrigChEff, 0, cin_esd);
  mgr->ConnectOutput (taskTrigChEff, 0, cout_trigChEff);

  //
  // Run the analysis
  //
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis(modeName[mode].Data(), chain);
  }
  timer.Stop();
  timer.Print();
}

//______________________________________________________________________________
TChain* CreateChain(Int_t mode)
{
  printf("*******************************\n");
  printf("*** Getting the Chain       ***\n");
  printf("*******************************\n");
  TChain *chain = 0x0;
  if(mode==kMgridLocal || mode==kMgrid){
    AliTagAnalysis *analysis = new AliTagAnalysis();
    chain = analysis->GetChainFromCollection(kDefaultXmlFile.Data(),"esdTree");
  }
  else{
    chain = new TChain("esdTree");
    TString inFileName("AliESDs.root");
    inFileName.Prepend(Form("%s/",kDefaultLocalInputDir.Data()));
    chain->Add(inFileName.Data());
  }
  if (chain) chain->ls();
  return chain;
}

//______________________________________________________________________________
void SetupPar(char* pararchivename)
{
  if (pararchivename) {
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
	return -1;
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
  } 
}
