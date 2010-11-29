/** 
 * Script to set-up a train 
 * 
 * @param esd           ESD file 
 * @param mc            Whether to do MC or not
 * @param nEvents       Number of events
 * @param nCutBins      Bins to cut away 
 * @param correctionCut 
 *
 * @ingroup pwg2_forward_analysis_scripts
 */
void RunManager(const char* esddir, Bool_t mc=kFALSE, Int_t nEvents=1000,
		Int_t nCutBins=1, Float_t correctionCut=0.1, 
		Bool_t proof=false)
{
  // --- Libraries to load -------------------------------------------
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  
  gSystem->Load("libSTEERBase");
  
  gSystem->Load("libESD") ;
  gSystem->Load("libAOD") ;
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  
  gSystem->Load("libPhysics");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWG2forward");
  gSystem->Load("libPWG2forward2");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof) { 
    TProof::Open("workers=2");
    const char* pkgs[] = { "STEERBase", "ESD", "AOD", "ANALYSIS", 
			   "ANALYSISalice", "PWG2forward", "PWG2forward2", 0};
    const char** pkg = pkgs;
    while (*pkg) { 
      gProof->UploadPackage(Form("${ALICE_ROOT}/%s.par",*pkg));
      gProof->EnablePackage(*pkg);    
      pkg++;
    }
  }
  
  // --- Our data chain ----------------------------------------------
  TChain* chain = new TChain("esdTree");

  // --- Get list of ESDs --------------------------------------------
  // Open source directory, and make sure we go back to were we were 
  TString oldDir(gSystem->WorkingDirectory());
  TSystemDirectory d(esddir, esddir);
  TList* files = d.GetListOfFiles();
  gSystem->ChangeDirectory(oldDir);

  // Sort list of files and check if we should add it 
  files->Sort();
  TIter next(files);
  TSystemFile* file = 0;
  while ((file = static_cast<TSystemFile*>(next()))) {
    if (file->IsDirectory()) continue;
    TString name(file->GetName());
    if (!name.EndsWith(".root")) continue;
    if (!name.Contains("AliESDs")) continue;
    TString esd(Form("%s/%s", file->GetTitle(), name.Data()));
    Info("RunManager", "Adding %s to chain", esd.Data());
    chain->Add(esd);
  }  

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Analysis Train", 
						    "FMD analysis train");

  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  esdHandler->SetInactiveBranches("AliESDACORDE "
				  "AliRawDataErrorLogs "
				  "CaloClusters "
				  "Cascades "
				  "EMCALCells "
				  "EMCALTrigger "
				  "Kinks "
				  "Cascades "
				  "MuonTracks "
				  "TrdTracks "
				  "CaloClusters "
				  "HLTGlobalTrigger");
  mgr->SetInputEventHandler(esdHandler);      
       
  // Monte Carlo handler
  // AliMCEventHandler* mcHandler = new AliMCEventHandler();
  // mgr->SetMCtruthEventHandler(mcHandler);
  // mcHandler->SetReadTR(readTR);    
  
  // AOD output handler
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAODs.root");

  // --- Add tasks ---------------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskFMD.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliAnalysisTask* task = AddTaskFMD(nCutBins, correctionCut);
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());

  task = AddTaskPhysicsSelection(mc, kTRUE, kTRUE);
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  
  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("RunManager", "Failed to initialize analysis train!");
    return;
  }
  // Some informative output 
  mgr->PrintStatus();
  // mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1 && !proof) 
    mgr->SetUseProgressBar(kTRUE);

  // Run the train 
  t.Start();
  mgr->StartAnalysis(proof ? "proof" : "local", chain, nEvents);
  t.Stop();
  t.Print();
}
//
// EOF
//
