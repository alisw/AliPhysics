//
// Run an analysis job 
// 
void RunManager(const char* esddir=".", 
		Int_t       nEvents=1000,
		Float_t     cmsGeV=900.,
		const char* col="p-p",
		Float_t     bkG=5., 
		Bool_t      proof=false)
{
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
  gSystem->Load("libPWGLFforward");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof) { 
    TProof::Open("workers=2");
    const char* pkgs[] = { "STEERBase", "ESD", "AOD", "ANALYSIS", 
			   "ANALYSISalice", "PWGLFforward", 0};
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
						    "The Analysis Train");
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

  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  //mcHandler->SetReadTR(readTR);    
  
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAODs.root");
    
  // --- Add our task -----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis/AddTaskFMD.C");
  AliFMDAnalysisTaskSE *fmdtask = AddTaskFMD(energy, col, bkG);
  
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
