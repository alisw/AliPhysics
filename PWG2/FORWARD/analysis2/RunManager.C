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
void RunManager(const char* esd, Bool_t mc=kFALSE, Int_t nEvents=1000,
		Int_t nCutBins=1, Float_t correctionCut=0.1, 
		Bool_t proof=false)
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
  gSystem->Load("libPWG2forward");
  gSystem->Load("libPWG2forward2");

  if (proof) { 
    gEnv->SetValue("Proof.GlobalPackageDirs", 
		   Form("%s:%s", 
			gEnv->GetValue("Proof.GlobalPackageDirs", "."), 
			gSystem->Getenv("ALICE_ROOT")));
    Info("RunManager", "PAR path=%s", 
	 gEnv->GetValue("Proof.GlobalPackageDirs", "."));
    TProof::Open("workers=1");
    const char* pkgs[] = { "STEERBase", "ESD", "AOD", "ANALYSIS", 
			   "ANALYSISalice", "PWG2forward", "PWG2forward2", 0};
    const char** pkg = pkgs;
    while (*pkg) { 
      gProof->UploadPackage(Form("${ALICE_ROOT}/%s.par",*pkg));
      gProof->EnablePackage(*pkg);    
      pkg++;
    }
    gProof->ShowPackages();
  }
  
  //You can expand this chain if you have more data :-)
  TChain* chain = new TChain("esdTree");
  chain->Add(esd);
  
  //Creating the manager and handlers
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
  aodHandler->SetOutputFileName("foo.root");
    
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/AddTaskFMD.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AliAnalysisTask* task = AddTaskFMD(nCutBins, correctionCut);
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());

  task = AddTaskPhysicsSelection(mc, kTRUE, kTRUE);
  mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
  
  // Run the analysis
    
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("RunManager", "Failed to initialize analysis train!");
    return;
  }
  // Some informative output 
  mgr->PrintStatus();
  mgr->SetUseProgressBar(kTRUE);

  // Write train to file - a test 
#if 0
  TDirectory* savDir = gDirectory;
  TFile* file = TFile::Open("analysisTrain.root", "RECREATE");
  mgr->Write();
  file->Close();
  savDir->cd();
#endif

  // Run the train 
  t.Start();
  mgr->StartAnalysis(proof ? "proof" : "local", chain, nEvents);
  t.Stop();
  t.Print();
}
//
// EOF
//
