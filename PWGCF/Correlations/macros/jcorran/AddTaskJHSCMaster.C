//_____________________________________________________________________
AliAnalysisTask *AddTaskJHSCMaster (TString taskName = "JHSCMaster", UInt_t period = 0, double ptmin = 0.5, Bool_t removebadarea = kFALSE)
{
// Load Custom Configuration and parameters.
  enum
  { lhc15o = 0, lhc18q = 1, lhc18r = 2 };
  TString speriod[3] = { "15o", "18q", "18r" };	// Needed string to load correct map config based on string.
  cout << "AddTaskJHSCMaster:: period=" << period << "\t ptmin=" << ptmin <<
    endl;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager ();

  const int Nsets = 4;		// Number of configurations.
  TString configNames[Nsets] = {
    "hybrid",			// 0 = default.
    "global",			// 1 for filtering.
    "SPD",			// 2 for centrality.
    "pileup"			// 3 for outliers.
  };

// Loading correction map.
  TString MAPfilenames[Nsets];
  TString MAPdirname =
    "alien:///alice/cern.ch/user/a/aonnerst/legotrain/NUAError/";
  AliJCorrectionMapTask *cmaptask =
    new AliJCorrectionMapTask ("JCorrectionMapTask");

  if (period == lhc18q || period == lhc18r)
    {				// Correction map for 2018 PbPb datasets.
      cmaptask->EnableCentFlattening (Form ("alien:///alice/cern.ch/user/j/jparkkil/legotrain/Cent/CentWeights_LHC%s_pass13.root", speriod[period].Data ()));	// Centrality flattening.
      cmaptask->EnableEffCorrection (Form ("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC18l8-0-Lists.root", speriod[period].Data ()));	// Efficiency correction.
    }
  if (period == lhc15o)
    {				// Correction map for 2015 PbPb dataset.
      cmaptask->EnableEffCorrection (Form ("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data/Eff--LHC%s-LHC16g-0-Lists.root", speriod[period].Data ()));	// Efficiency correction.
    }
  for (int i = 0; i < Nsets; i++)
    {
      MAPfilenames[i] = Form ("%sPhiWeights_LHC%s_Error_pt%02d_s_%s.root", MAPdirname.Data (), speriod[period].Data (), Int_t (ptmin * 10), configNames[i].Data ());	// Azimuthal correction.
      cmaptask->EnablePhiCorrection (i, MAPfilenames[i]);	// i is index for set file correction ->SetPhiCorrectionIndex(i);
    }
  mgr->AddTask ((AliAnalysisTask *) cmaptask);	// Loading of the correction map added to the analysis manager.

  Int_t hybridCut = 768;
  Int_t globalCut = 96;
  UInt_t selEvt;		// Trigger.
  if (period == lhc15o)
    {
      selEvt = AliVEvent::kINT7;
    }				// Minimum bias.
  else if (period == lhc18q || period == lhc18r)
    {
      selEvt = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral;
    }				// Minimum bias + central + semicentral.

// Setting the JCatalystTasks for the event and track selection.
  AliJCatalystTask *fJCatalyst[Nsets];	// One tracking per configuration with FLUC_CUT_OUTLIERS on.

  for (int i = 0; i < Nsets; i++) {
      fJCatalyst[i] = new AliJCatalystTask (Form ("JCatalystTask_%s", configNames[i].Data ()));
      cout << fJCatalyst[i]->GetJCatalystTaskName() << endl;
      // Set the correct flags to use.
      if (i != 3) fJCatalyst[i]->AddFlags (AliJCatalystTask::FLUC_CUT_OUTLIERS);
      if (period == lhc18q || period == lhc18r)	fJCatalyst[i]->AddFlags (AliJCatalystTask::FLUC_CENT_FLATTENING);
      // Set the trigger and centrality selection.
      fJCatalyst[i]->SelectCollisionCandidates (selEvt);
      if (i == 2) {			// SPD clusters for the systematics.
	  	fJCatalyst[i]->SetCentDetName ("CL1");
	  } else {
     	fJCatalyst[i]->SetCentDetName ("V0M"); // V0M in the default analysis and other systematics.
      }
      // Set the filtering and kinematic cuts.
      if (i == 1) {			// Global tracks for the systematics.
	    fJCatalyst[i]->SetTestFilterBit (globalCut);
	  } else {			// Hybrid tracks for the default analysis.
	    fJCatalyst[i]->SetTestFilterBit (hybridCut);
	  }
      fJCatalyst[i]->SetEtaRange (-0.8, 0.8);
      fJCatalyst[i]->SetPtRange (ptmin, 5.0);
      fJCatalyst[i]->SetPhiCorrectionIndex (i);

  }				// End for.
// Configuration of the analysis task itself.
  //-------- JHSC Wagons -------
  AliJHSCTask *myTask[Nsets];
  for (int i = 0; i < Nsets; i++) {
  	myTask[i] = new AliJHSCTask (Form("%s_s_%s", taskName.Data (), configNames[i].Data ()));
  	myTask[i]->SetJCatalystTaskName (fJCatalyst[i]->GetJCatalystTaskName());
	myTask[i]->SetDebugLevel(0);
  }

  // Must add the tasks
  for (int i = 0; i < Nsets; i++) {
  	mgr->AddTask ((AliAnalysisTask *) fJCatalyst[i]);
  	mgr->AddTask ((AliAnalysisTask *) myTask[i]);
  }
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer ();

  // Connect input/output
  for (int i = 0; i < Nsets; i++) {
      mgr->ConnectInput (fJCatalyst[i], 0, cinput); // MUST
      mgr->ConnectInput (myTask[i], 0, cinput); // MUST
      AliAnalysisDataContainer *jHist = mgr->CreateContainer (Form ("%scontainer", myTask[i]->GetName ()), TList::Class (), AliAnalysisManager::kOutputContainer, Form ("%s:%s",
				    AliAnalysisManager::GetCommonFileName (),myTask[i]->GetName ()));
      mgr->ConnectOutput (myTask[i], 1, jHist);
      //mgr->ConnectOutput (fJCatalyst[i], 1, jHist);
    }

  return myTask[0];
}
