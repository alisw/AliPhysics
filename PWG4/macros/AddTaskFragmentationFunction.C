
/*************************************************************************************************
***  Add Fragmentation Function Task ***
**************************************************************************************************
The fragmenation function task expects an ESD filter and jet finder running before this task. 
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




AliAnalysisTaskFragmentationFunction *AddTaskFragmentationFunction(UInt_t iFlag=1, UInt_t filterMask=16){
        
        AliAnalysisTaskFragmentationFunction *ff=0;

        // only reconstructed (default)
	if(iFlag&(1<<0)) ff = AddTaskFragmentationFunction("jets", "", "", "", filterMask);
        // MC tracks in acceptance, MC jets in acceptance
	if(iFlag&(1<<1)) ff = AddTaskFragmentationFunction("jets", "AODMC2b", "AODMCb", "AODMCb", filterMask);
        // kine tracks in acceptance, pythia jets in acceptance
	if(iFlag&(1<<2)) ff = AddTaskFragmentationFunction("jets", "", "KINEb", "KINEb", filterMask);
        // reconstructed charged tracks after cuts, MC jets in acceptance 
	if(iFlag&(1<<3)) ff = AddTaskFragmentationFunction("jets", "AODMC2b", "AODMCb", "AOD2b", filterMask);
	
	return ff;
}

// _______________________________________________________________________________________

AliAnalysisTaskFragmentationFunction *AddTaskFragmentationFunction(
        const char* recJetsBranch,
	const char* genJetsBranch,
	const char* jetType,
	const char* trackType,
	UInt_t filterMask)
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

   TString branchRecJets(recJetsBranch);
   TString branchGenJets(genJetsBranch);
   TString typeJets(jetType);
   TString typeTracks(trackType);

   if(branchRecJets.Length()==0) branchRecJets = "noRecJets";
   if(branchGenJets.Length()==0) branchGenJets = "noGenJets";
   if(typeTracks.Length()==0) typeTracks = "trackTypeUndef";
   if(typeJets.Length()==0)   typeJets   = "jetTypeUndef";
   
   // Create the task and configure it.
   //===========================================================================

   AliAnalysisTaskFragmentationFunction *task = new AliAnalysisTaskFragmentationFunction(
        Form("Fragmenation Function %s %s %s %s", branchRecJets.Data(), branchGenJets.Data(), typeJets.Data(), typeTracks.Data()));
   
   if(debug>=0) task->SetDebugLevel(debug);
   
   Printf("Rec Jets %s", branchRecJets.Data());
   Printf("Gen Jets %s", branchGenJets.Data());
   Printf("Jet Type %s", typeJets.Data());
   Printf("Track Type %s", typeTracks.Data());
   
   if(!branchRecJets.Contains("noRecJets")) task->SetBranchRecJets(branchRecJets);
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
   
   task->SetFilterMask(filterMask);
  
   // set default parameter 
   task->SetTrackCuts();       // default : pt > 0.150 GeV, |eta|<0.9, full phi acc
   task->SetJetCuts();         // default: jet pt > 5 GeV, |eta|<0.5, full phi acc
   task->SetDijetCuts();       // default: to be defined
   task->SetFFRadius();        // default: R = 0.4
   
   task->SetHighPtThreshold(); // default: pt > 5 Gev
   task->SetFFHistoBins();
   task->SetQAJetHistoBins();
   task->SetQATrackHistoBins();

   mgr->AddTask(task);

   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   AliAnalysisDataContainer *coutput_FragFunc = mgr->CreateContainer(
      Form("fracfunc_%s_%s_%s_%s", branchRecJets.Data(), branchGenJets.Data(), typeTracks.Data(), typeJets.Data()),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s:PWG4_FragmentationFunction_%s_%s_%s_%s", 
         AliAnalysisManager::GetCommonFileName(), branchRecJets.Data(), branchGenJets. Data(), typeTracks.Data(), typeJets.Data()));

   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput (task, 1, coutput_FragFunc);
   
   return task;
}
