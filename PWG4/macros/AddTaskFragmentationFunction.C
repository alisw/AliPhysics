
/*************************************************************************************************
***  Add Fragmentation Function Task ***
**************************************************************************************************
The fragmentation function task expects an ESD filter and jet finder running before this task. 
Or it runs on delta-AODs filled with filtered tracks and jets before.

** Parameters **
(char) recJetsBranch: branch in AOD for (reconstructed) jets
(char) genJetsBranch: branch in AOD for (generated) jets
(char) jetType: "AOD"   jets from recJetsBranch
                "AODMC" jets from genJetsBranch
                "KINE"  jets from PYCELL
                 +"b" (e.g. "AODb") jets with acceptance cuts
(char) trackType: "AOD"     reconstructed tracks from AOD filled by ESD filter (choose filter mask!)
                  "AODMC"   MC tracks from AOD filled by kine filter
                  "KINE"    kine particles from MC event 
                  +"2" (e.g. "AOD2")  charged tracks only
                  +"b" (e.g. "AOD2b") with acceptance cuts
(UInt_t) filterMask: select filter bit of ESD filter task

***************************************************************************************************/




AliAnalysisTaskFragmentationFunction *AddTaskFragmentationFunction(UInt_t iFlag=1, UInt_t filterMask=32, Int_t eventClass=0){
        
        AliAnalysisTaskFragmentationFunction *ff=0;

        // UA1 Jets
        // only reconstructed (default)
	if(iFlag&(1<<0)) ff = AddTaskFragmentationFunction("jetsAOD_UA1", "", "", "", "", filterMask, 0.4,0,1000., eventClass);
        // charged MC tracks and jets
	if(iFlag&(1<<1)) ff = AddTaskFragmentationFunction("jetsAOD_UA1", "", "jetsAODMC2_UA1", "AODMC", "AODMC2", filterMask, 0.4,0,1000., eventClass);
        // charged MC tracks and jets with acceptance cuts
	if(iFlag&(1<<2)) ff = AddTaskFragmentationFunction("jetsAOD_UA1", "", "jetsAODMC2_UA1", "AODMCb", "AODMC2b", filterMask, 0.4,0,1000., eventClass);
        // kine tracks in acceptance, pythia jets in acceptance
	if(iFlag&(1<<3)) ff = AddTaskFragmentationFunction("jetsAOD_UA1", "", "", "KINEb", "KINEb", filterMask, 0.4,0,1000., eventClass);
        // reconstructed charged tracks after cuts, MC jets in acceptance 
	if(iFlag&(1<<4)) ff = AddTaskFragmentationFunction("jetsAOD_UA1", "", "jetsMC2b", "AODMCb", "AOD2b", filterMask, 0.4,0,1000., eventClass);
	// reconstruction efficiency: pointing with rec jet axis into gen tracks 
	if(iFlag&(1<<5)) ff = AddTaskFragmentationFunction("jetsAOD_UA1", "", "jetsAODMC2_UA1", "AODb", "AODMC2b", filterMask, 0.4,0,1000., eventClass);



     	// Jet background subtracted
	// anti-kt, pt cut 1 GeV
	if(iFlag&(1<<20)) ff = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "", "", "", filterMask, 0.4,2,1000.,eventClass, "_Skip00");
	// anti-kt, pt cut 2 GeV
	if(iFlag&(1<<21)) ff = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "", "", "", filterMask, 0.4,2,2000.,eventClass, "_Skip00");
	// anti-kt, pt cut 150 MeV
	if(iFlag&(1<<22)) ff = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "", "", "", filterMask, 0.2,2,150.,eventClass, "_Skip00");

	
	// Jet background subtracted
	if(iFlag&(1<<23)) ff = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "", "", "", filterMask, 0.4,2,150.,eventClass, "_Skip00");
        // charged MC tracks and jets
	if(iFlag&(1<<24)) ff = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "jetsAODMC2_FASTJET", "AODMC", "AODMC2", filterMask, 0.4,2,150.,eventClass, "_Skip00");
        // charged MC tracks and jets with acceptance cuts
	if(iFlag&(1<<25)) ff = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "jetsAODMC2_FASTJET", "AODMCb", "AODMC2b", filterMask, 0.4,2,150., eventClass, "_Skip00");

       if(iFlag&(1<<26)) ff = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "", "", "", filterMask, 0.4,1,150.,eventClass, "_Skip00");

       if(iFlag&(1<<27)) ff = AddTaskFragmentationFunction("clustersAOD_ANTIKT", "", "", "", "", filterMask, 0.4,3,150.,eventClass, "_Skip00");	

      // SISCONE 
      if(iFlag&(1<<28)) ff = AddTaskFragmentationFunction("jetsAOD_SISCONE", "", "", "", "", filterMask, 0.4,1,150.,eventClass);
      if(iFlag&(1<<29)) ff = AddTaskFragmentationFunction("jetsAOD_SISCONE", "", "", "", "", filterMask, 0.4,2,150.,eventClass);
      if(iFlag&(1<<30)) ff = AddTaskFragmentationFunction("jetsAOD_SISCONE", "", "", "", "", filterMask, 0.4,3,150.,eventClass);

	return ff;
}

// _______________________________________________________________________________________

AliAnalysisTaskFragmentationFunction *AddTaskFragmentationFunction(
        const char* recJetsBranch,
	const char* recJetsBackBranch,
	const char* genJetsBranch,
	const char* jetType,
	const char* trackType,
	UInt_t filterMask,
        Float_t radius,
        int kBackgroundMode,
        float PtTrackMin,
        Int_t eventClass=0,
        TString BrOpt="",
        TString BrOpt2="",
        Float_t radiusBckg=0.4)
{
   // Creates a fragmentation function task,
   // configures it and adds it to the analysis manager.

   //******************************************************************************
   //*** Configuration Parameter **************************************************
   //******************************************************************************

   // space for configuration parameter: histo bin, cuts, ...
   // so far only default parameter used

   Int_t debug = -1; // debug level, -1: not set here

   //******************************************************************************


   
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
	  ::Error("AddTaskFragmentationFunction", "No analysis manager to connect to.");
	  return NULL;
   }
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
	 ::Error("AddTaskFragmentationFunction", "This task requires an input event handler");
	  return NULL;
   }

   TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
   Printf("Data Type: %s", type.Data());

   TString branchRecBackJets(recJetsBackBranch);
   TString branchRecJets(recJetsBranch);
   TString branchGenJets(genJetsBranch);
   TString typeJets(jetType);
   TString typeTracks(trackType);

   if(branchRecBackJets.Length()==0) branchRecBackJets = "noRecBackJets";
   if(branchRecJets.Length()==0) branchRecJets = "noRecJets";
   if(branchGenJets.Length()==0) branchGenJets = "noGenJets";
   if(typeTracks.Length()==0) typeTracks = "trackTypeUndef";
   if(typeJets.Length()==0)   typeJets   = "jetTypeUndef";
   
   // Create the task and configure it.
   //===========================================================================

   AliAnalysisTaskFragmentationFunction *task = new AliAnalysisTaskFragmentationFunction(
        Form("Fragmentation Function %s %s %s %s", branchRecJets.Data(), branchGenJets.Data(), typeJets.Data(), typeTracks.Data()));
   
   if(debug>=0) task->SetDebugLevel(debug);
   
   Printf("Rec Jets %s", branchRecJets.Data());
   Printf("Back Rec Jets %s", branchRecBackJets.Data());
   Printf("Gen Jets %s", branchGenJets.Data());
   Printf("Jet Type %s", typeJets.Data());
   Printf("Track Type %s", typeTracks.Data());
   
   // attach the filter mask and options
   TString cAdd = "";
   cAdd += Form("%02d",(int)((radius+0.01)*10.));
   cAdd += Form("_B%d",(int)((kBackgroundMode)));
   cAdd += Form("_Filter%05d",filterMask);
   cAdd += Form("_Cut%05d",(int)((PtTrackMin)));
   cAdd += Form("%s",BrOpt.Data());
   cAdd += Form("%s",BrOpt2.Data());
   Printf("%s",cAdd.Data());

   TString cAddb = "";
   cAddb += Form("%02d",(int)((radiusBckg+0.01)*10.));
   cAddb += Form("_B%d",(int)((kBackgroundMode)));
   cAddb += Form("_Filter%05d",filterMask);
   cAddb += Form("_Cut%05d",(int)((PtTrackMin)));
   cAddb += Form("%s",BrOpt.Data());
   cAddb += Form("%s",BrOpt2.Data());
   Printf("%s",cAddb.Data());

   TString cAddmc = "";
   cAddmc += Form("%02d",(int)((radius+0.01)*10.));
   cAddmc += Form("_B%d",(int)((kBackgroundMode)));
   cAddmc += Form("_Filter%05d",filterMask);
   cAddmc += Form("_Cut%05d",(int)((PtTrackMin)));
   Printf("%s",cAddmc.Data());


   if(branchRecJets.Contains("AOD")&&branchRecJets.Contains("jets")&&!branchRecJets.Contains("MC"))branchRecJets = branchRecJets + cAdd;
   if(branchRecJets.Contains("AOD")&&branchRecJets.Contains("cluster")&&!branchRecJets.Contains("MC"))branchRecJets = branchRecJets + cAdd;

   if(branchRecBackJets.Contains("back")&&branchRecBackJets.Contains("cluster")&&!branchRecBackJets.Contains("MC"))branchRecBackJets = branchRecBackJets + cAddb; 

   if(branchGenJets.Contains("AOD")&&branchGenJets.Contains("MC"))branchGenJets = branchGenJets + cAddmc;

   Printf("Gen jets branch %s: ", branchGenJets.Data());
   Printf("Rec jets branch %s: ", branchRecJets.Data());
   Printf("Jet backg branch %s: ", branchRecBackJets.Data());

   if(!branchRecJets.Contains("noRecJets")) task->SetBranchRecJets(branchRecJets);
   if(!branchRecBackJets.Contains("noRecBackJets")) task->SetBranchRecBackJets(branchRecBackJets);
   if(!branchGenJets.Contains("noGenJets")) task->SetBranchGenJets(branchGenJets);


   if(typeTracks.Contains("AODMC2b"))      task->SetTrackTypeGen(AliAnalysisTaskFragmentationFunction::kTrackAODMCChargedAcceptance);
   else if(typeTracks.Contains("AODMC2"))  task->SetTrackTypeGen(AliAnalysisTaskFragmentationFunction::kTrackAODMCCharged);
   else if(typeTracks.Contains("AODMC"))   task->SetTrackTypeGen(AliAnalysisTaskFragmentationFunction::kTrackAODMCAll);
   else if(typeTracks.Contains("KINE2b"))  task->SetTrackTypeGen(AliAnalysisTaskFragmentationFunction::kTrackKineChargedAcceptance);
   else if(typeTracks.Contains("KINE2"))   task->SetTrackTypeGen(AliAnalysisTaskFragmentationFunction::kTrackKineCharged);
   else if(typeTracks.Contains("KINE"))    task->SetTrackTypeGen(AliAnalysisTaskFragmentationFunction::kTrackKineAll);
   else if(typeTracks.Contains("AODb"))    task->SetTrackTypeGen(AliAnalysisTaskFragmentationFunction::kTrackAODCuts);
   else if(typeTracks.Contains("AOD"))     task->SetTrackTypeGen(AliAnalysisTaskFragmentationFunction::kTrackAOD);
   else if(typeTracks.Contains("trackTypeUndef")) task->SetTrackTypeGen(0); // undefined
   else Printf("trackType %s not found", typeTracks.Data());

   if(typeJets.Contains("AODMCb"))         task->SetJetTypeGen(AliAnalysisTaskFragmentationFunction::kJetsGenAcceptance);
   else if(typeJets.Contains("AODMC"))     task->SetJetTypeGen(AliAnalysisTaskFragmentationFunction::kJetsGen);
   else if(typeJets.Contains("KINEb"))     task->SetJetTypeGen(AliAnalysisTaskFragmentationFunction::kJetsKineAcceptance);
   else if(typeJets.Contains("KINE"))      task->SetJetTypeGen(AliAnalysisTaskFragmentationFunction::kJetsKine);
   else if(typeJets.Contains("AODb"))      task->SetJetTypeGen(AliAnalysisTaskFragmentationFunction::kJetsRecAcceptance);
   else if(typeJets.Contains("AOD"))       task->SetJetTypeGen(AliAnalysisTaskFragmentationFunction::kJetsRec);
   else if(typeJets.Contains("jetTypeUndef")) task->SetJetTypeGen(0); // undefined
   else Printf("jetType %s not found", typeJets.Data());
   
   if(typeJets.Contains("AODMCb")) task->SetJetTypeRecEff(AliAnalysisTaskFragmentationFunction::kJetsGenAcceptance); // kJetsRecAcceptance
   else if(typeJets.Contains("AODb")) task->SetJetTypeRecEff(AliAnalysisTaskFragmentationFunction::kJetsRecAcceptance); 
   else task->SetJetTypeRecEff(0);

   task->SetFilterMask(filterMask);
   task->SetEventClass(eventClass);
  
   // Set default parameters 
   // Cut selection 
   
   if((int)PtTrackMin == 2000)      task->SetTrackCuts(2.0, -0.9, 0.9, 0., 2*TMath::Pi());
   else if((int)PtTrackMin == 1000) task->SetTrackCuts(1.0, -0.9, 0.9, 0., 2*TMath::Pi());
   else                             task->SetTrackCuts();  // default : pt > 0.150 GeV, |eta|<0.9, full phi acc


   task->SetJetCuts();         // default: jet pt > 5 GeV, |eta|<0.5, full phi acc
   task->SetDiJetCuts();       // default: type of cut = 1 (cut in deltaPhi), deltaPhi = 0., cdf = 0.5, fraction of pt = 0.6
   task->SetKindSlices();      // default: kindSlice = 1 (inv mass)
   if(radius <= 0.2) task->SetFFRadius(0.2); // R = 0.2   
   else              task->SetFFRadius();    // default: R = 0.4
   task->SetFFBckgRadius();    // default: R = 0.7
   task->SetBckgSubMethod();   // default: subMethod = O, 1 = leading jet removed for rho extraction, 2 = 2 leading jets removed
   task->SetIJMode(0);         // default: ijMode = 1
   task->SetQAMode();          // default: qaMode = 3
   task->SetFFMode();          // default: ffMode = 1
   task->SetDJMode(0);          // default: djMode = 1
   task->SetEffMode(0);         // default: effMode = 1
   task->SetPhiCorrMode(0);     // default: phiCorrMode = 1
   task->SetHighPtThreshold();  // default: pt > 5 Gev
   task->UseRecEffRecJetPtBins(); // efficiency in bins of rec/gen jet pt - default: kTRUE  

   task->SetBckgMode(1);        // default: bgMode = 1 
   task->SetBckgType();
   task->SetBranchRecBackClusters(Form("clustersAOD_KT04_B0_Filter%05d_Cut00150_Skip00",filterMask));

   // Define histo bins
   task->SetFFHistoBins(23, 5, 120, 480, 0., 120.,70,  0., 7.,22,  0.,  1.1);

   task->SetQAJetHistoBins();
   task->SetQATrackHistoBins();
   task->SetIJHistoBins();
   task->SetDiJetHistoBins();
   task->SetQADiJetHistoBins();

   mgr->AddTask(task);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   AliAnalysisDataContainer *coutput_FragFunc = mgr->CreateContainer(
      Form("fracfunc_%s_%s_%s_%s_cl%d", branchRecJets.Data(), branchGenJets.Data(), typeTracks.Data(), typeJets.Data(), eventClass),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:PWG4_FragmentationFunction_%s_%s_%s_%s_cl%d", 
         AliAnalysisManager::GetCommonFileName(), branchRecJets.Data(), branchGenJets. Data(), typeTracks.Data(), typeJets.Data(), eventClass));

   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (task, 1, coutput_FragFunc);
   
   return task;
}
