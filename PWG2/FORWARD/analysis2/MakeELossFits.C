/**
 * @file 
 * 
 * @ingroup pwg2_forward_scripts
 */
/** 
 * Run a pass on ESD data to produce the energ loss fits 
 * 
 *
 * @ingroup pwg2_forward_scripts
 */
void MakeELossFits(const char* esddir, 
		   Int_t       nEvents=1000, 
		   Int_t       proof=0,
		   Bool_t      mc=false,
		   Bool_t      cent=false)
{
  // --- Libraries to load -------------------------------------------
  gROOT->Macro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadLibs.C");

  // --- Check for proof mode, and possibly upload pars --------------
  if (proof > 0) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/LoadPars.C");
    LoadPars(proof);
  }
  
  // --- Our data chain ----------------------------------------------
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/FORWARD/analysis2/scripts/MakeChain.C");
  TChain* chain = MakeChain("ESD",esddir, true);
  // If 0 or less events is select, choose all 
  if (nEvents <= 0) nEvents = chain->GetEntries();
  Info("MakeELossFits", "Will analyse %d events", nEvents);

  // --- Set the macro path ------------------------------------------
  gROOT->SetMacroPath(Form("%s:$(ALICE_ROOT)/PWG2/FORWARD/analysis2:"
			   "$ALICE_ROOT/ANALYSIS/macros",
			   gROOT->GetMacroPath()));

  // --- Creating the manager and handlers ---------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Forward ELoss Train", 
						    "Forward energy loss");
  AliAnalysisManager::SetCommonFileName("forward_eloss.root");

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
  
  // AOD output handler
  AliAODHandler* aodHandler   = new AliAODHandler();
  mgr->SetOutputEventHandler(aodHandler);
  aodHandler->SetOutputFileName("AliAODs.root");

  // --- Add tasks ---------------------------------------------------
  gROOT->LoadMacro("AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(mc, kTRUE, kFALSE);
  
  gROOT->LoadMacro("AddTaskCentrality.C");
  AddTaskCentrality();
  
  AliFMDEnergyFitterTask* task = new AliFMDEnergyFitterTask("fmdEnergyFitter");
  mgr->AddTask(task);
  
  // --- Make the output container and connect it --------------------
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  AliAnalysisDataContainer* histOut = 
    mgr->CreateContainer("Forward", TList::Class(), 
			 AliAnalysisManager::kOutputContainer,outputfile);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, histOut);

  // --- Set parameters on the algorithms ----------------------------
  // Set the number of SPD tracklets for which we consider the event a
  // low flux event
  task->GetEventInspector().SetLowFluxCut(1000); 
  // Set the maximum error on v_z [cm]
  task->GetEventInspector().SetMaxVzErr(0.2);
  // Set the eta axis to use - note, this overrides whatever is used
  // by the rest of the algorithms - but only for the energy fitter
  // algorithm. 
  task->GetEnergyFitter().SetEtaAxis(200, -4, 6);
  // Set maximum energy loss to consider 
  task->GetEnergyFitter().SetMaxE(15); 
  // Set number of energy loss bins 
  task->GetEnergyFitter().SetNEbins(100);
  // Set whether to use increasing bin sizes 
  task->GetEnergyFitter().SetUseIncreasingBins(true);
  // Set whether to do fit the energy distributions 
  task->GetEnergyFitter().SetDoFits(kTRUE);
  // Set whether to make the correction object 
  task->GetEnergyFitter().SetDoMakeObject(kTRUE);
  // Set the low cut used for energy
  task->GetEnergyFitter().SetLowCut(0.4);
  // Set the number of bins to subtract from maximum of distributions
  // to get the lower bound of the fit range
  task->GetEnergyFitter().SetFitRangeBinWidth(4);
  // Set the maximum number of landaus to try to fit (max 5)
  task->GetEnergyFitter().SetNParticles(5);
  // Set the minimum number of entries in the distribution before
  // trying to fit to the data
  task->GetEnergyFitter().SetMinEntries(1000);
  // --- Set limits on fits the energy -------------------------------
  // Maximum relative error on parameters 
  AliFMDCorrELossFit::ELossFit::fgMaxRelError = .12;
  // Least weight to use 
  AliFMDCorrELossFit::ELossFit::fgLeastWeight = 1e-5;
  // Maximum value of reduced chi^2 
  AliFMDCorrELossFit::ELossFit::fgMaxChi2nu   = 10;
  
  // --- Run the analysis --------------------------------------------
  TStopwatch t;
  if (!mgr->InitAnalysis()) {
    Error("RunManager", "Failed to initialize analysis train!");
    return;
  }
  // Some informative output 
  mgr->PrintStatus();
  // mgr->SetDebugLevel(3);
  if (mgr->GetDebugLevel() < 1 && proof <= 0) mgr->SetUseProgressBar(kTRUE);

  // Run the train 
  t.Start();
  Printf("=== RUNNING ANALYSIS ==================================");
  mgr->StartAnalysis(proof > 0 ? "proof" : "local", chain, nEvents);
  t.Stop();
  t.Print();
}
//
// EOF
//
