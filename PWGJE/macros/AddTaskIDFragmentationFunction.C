R__ADD_INCLUDE_PATH($ALICE_PHYSICS);
#include <PWG/EMCAL/macros/AddTaskMCTrackSelector.C>
#include <PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C>
#include <PWGJE/EMCALJetTasks/macros/AddTaskRhoNew.C>

/*************************************************************************************************
***  Add Fragmentation Function Task ***
**************************************************************************************************
The ID fragmentation function task expects an ESD filter and jet finder running before this task. 
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

Typical parameters to run on 11a1* (MC_pp@7TeV):
"clustersAOD_ANTIKT", "", "clustersAODMC2_ANTIKT", "AODMCb", "AODMC2b", AliAnalysisManager::GetGlobalInt("kHighPtFilterMask", gDebug), 
-0.4, 0, 1000*AliAnalysisManager::GetGlobalDbl("kTrackPtCut", gDebug), 0, "_Skip00", "", "_Skip00", 0.4, -1, 0, 0, kFALSE,
"PWGJE_taskPID_Jets", "", "PWGJE_taskPID_Jets_Inclusive", "" 

***************************************************************************************************/

void SetEfficiencyFunctionsFastSimulation(AliAnalysisTaskIDFragmentationFunction *task = 0x0, TString fastSimParamFile = "", Double_t effFactor = 1.0, Double_t resFactor = 1.0, Int_t ffChange = 0) {
  if (!task || !task->GetUseFastSimulations())
    return;

  TF1** effFunctions = new TF1*[2*AliPID::kSPECIES];
  
  //effFactor: Multiplied with the parametrized efficiency
	//resFactor: Multiplied with the observed resolution
  //ffChange: 0 for normal Fragmentation, 1 for low-Pt enhancement, 2 for low-Pt depletion
  
  const Double_t kObsResolution = 0.002;  //Only change if the observed resolution is changing (external study)
  
  TString parameterString = AliPieceWisePoly::ReadFSParameters(fastSimParamFile.Data(), effFunctions);                                                        
	
  task->SetEfficiencyFunctions(effFunctions);
  task->SetFastSimulationParameters(parameterString);    // Set string with joined parameters additionally, because the functions are transient and will not be copied. Functions are initialized again at the start of the run
  task->SetFastSimEffFactor(effFactor);
  task->SetFastSimRes(kObsResolution);   
  task->SetFastSimResFactor(resFactor);
  task->SetFFChange(static_cast<AliAnalysisTaskIDFragmentationFunction::FragmentationFunctionChange>(ffChange));
  return;
}

void postConfig(AliAnalysisTaskIDFragmentationFunction* task, TString namesOfInclusivePIDtasks, TString namesOfJetPIDtasks,
                TString namesOfJetUEPIDtasks, TString namesOfJetUEMethods) {
	
	TObjArray* nameArrayInclusive = namesOfInclusivePIDtasks.Tokenize(";");
	Int_t numInclusivePIDtasks = nameArrayInclusive->GetEntriesFast();
	TString* namesOfInclusiveTasks = new TString[numInclusivePIDtasks];
	for (Int_t i=0;i<numInclusivePIDtasks;i++) {
		namesOfInclusiveTasks[i] = (((TObjString*)(nameArrayInclusive->At(i)))->GetString());
	}
	task->SetNamesOfInclusivePIDtasks(numInclusivePIDtasks, namesOfInclusiveTasks);  

	
	TObjArray* nameArrayJets = namesOfJetPIDtasks.Tokenize(";");
	Int_t numJetPIDtasks = nameArrayJets->GetEntriesFast();
	TString* namesOfJetTasks = new TString[numJetPIDtasks];
	for (Int_t i=0;i<numJetPIDtasks;i++) {
		namesOfJetTasks[i] = (((TObjString*)(nameArrayJets->At(i)))->GetString());
	}
	task->SetNamesOfJetPIDtasks(numJetPIDtasks, namesOfJetTasks);  
	
	
	TObjArray* nameArrayJetsUE = namesOfJetUEPIDtasks.Tokenize(";");
	Int_t numJetUEPIDtasks = nameArrayJetsUE->GetEntriesFast();
    TObjArray* nameArrayJetsUEMethods = namesOfJetUEMethods.Tokenize(";");
    //TODO: Check if length is same for both, if not exit (or take shorter one with warning?)
	TString* namesOfJetUETasks = new TString[numJetUEPIDtasks];
    TString* namesOfJetUEMethodsArray = new TString[numJetUEPIDtasks];
	for (Int_t i=0;i<numJetUEPIDtasks;i++) {
		namesOfJetUETasks[i] = (((TObjString*)(nameArrayJetsUE->At(i)))->GetString());
        namesOfJetUEMethodsArray[i] = (((TObjString*)(nameArrayJetsUEMethods->At(i)))->GetString());
	}
	task->SetNamesOfJetUEPIDtasks(numJetUEPIDtasks, namesOfJetUETasks, namesOfJetUEMethodsArray);  
  
	
  printf("PID framework:\n");
  printf("Jet PID tasks: ");
  for (Int_t i = 0; i < numJetPIDtasks; i++)
    printf("%s ", task->GetNamesOfJetPIDtasks()[i].Data());
  printf("\n");
  printf("Jet Underlying Event PID tasks: ");
  for (Int_t i = 0; i < numJetUEPIDtasks; i++) {
    printf("%s ", task->GetNamesOfJetUEPIDtasks()[i].Data());
    cout << ";  Underlying event subtraction method: " << task->GetUEMethods()[i] << endl; 
  }
  printf("\n");  
  printf("Inclusive PID task: ");
  for (Int_t i = 0; i < numInclusivePIDtasks; i++) {
    printf("%s ", task->GetNamesOfInclusivePIDtasks()[i].Data());
  }
}

// _______________________________________________________________________________________

AliAnalysisTaskIDFragmentationFunction *AddTaskIDFragmentationFunction(
  const char* recJetsBranch = "clustersAOD_ANTIKT",
  const char* recJetsBackBranch = "",
  const char* genJetsBranch = "clustersAODMC2_ANTIKT",
  const char* jetType = "AODMCb",
  const char* trackType = "AODMC2b",
  Float_t radius = -0.4,
  Int_t PtTrackMin = 150.0,
  Int_t eventClass=-1,
  Int_t FFMaxTrackPt = -1,
  Int_t FFMinNTracks = 0,
  UInt_t filterMaskTracks = 0,
  Bool_t onlyConsiderLeadingJets = kFALSE,
  Float_t MC_pThard_cut = -1.,
  TString namesOfInclusivePIDtasks = "",
  TString namesOfJetPIDtasks = "",
  TString namesOfJetUEPIDtasks = "",
  TString namesOfJetUEMethods = "",
  AliAnalysisTaskEmcal::BeamType iBeamType = AliAnalysisTaskEmcal::kpp,
  Int_t nJetContainer = 1,
  AliEmcalJet::JetAcceptanceType* jetAcceptanceRegion = 0x0,
  AliJetContainer::EJetType_t* typeOfJet = 0x0,
  TString correctionTaskFileName = "",
  TString period = "lhc13c",
  TString fastSimParamFile = "")
{
   // Creates a fragmentation function task,
   // configures it and adds it to the analysis manager.

   //******************************************************************************
   //*** Configuration Parameter **************************************************
   //******************************************************************************

	// space for configuration parameter: histo bin, cuts, ...
	// so far only default parameter used
  
     // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if (!mgr) {
    ::Error("AddTaskIDFragmentationFunction", "No analysis manager to connect to.");
    return NULL;
   }

	Int_t debug = -1; // debug level, -1: not set here
	UInt_t kPhysSel = AliVEvent::kMB | AliVEvent::kINT8 | AliVEvent::kINT7;
	AliEmcalTrackSelection::ETrackFilterType_t trackFilterType = AliEmcalTrackSelection::kHybridTracks2011woNoRefit;
	
	Bool_t subtractBackgroundPt = kFALSE;
	UInt_t kNumberExcludeJetsInBackground = 2;
  
	Double_t trackPtCut = PtTrackMin/1000.0;
	Double_t ghostArea = 0.005;
	Double_t clustPtCut = 0.5;
	Double_t jetMinPt = 1.0;
	Float_t kEtaMin = -0.9;
	Float_t kEtaMax = 0.9;
	
	Int_t kMinMultiplicity = -1;
	Int_t kMaxMultiplicity = -1;
	
  //Calculated parameters
	Bool_t useJets = nJetContainer > 0;
	
	Bool_t useFastSimulation = fastSimParamFile != "";
	
	Double_t absRadiusCut = TMath::Abs(radius);
	
	AliTrackContainer::SetDefTrackCutsPeriod(period);  
	
  if (!jetAcceptanceRegion) {
    jetAcceptanceRegion = new AliEmcalJet::JetAcceptanceType[nJetContainer];
    for (Int_t i=0;i<nJetContainer;i++) {
      jetAcceptanceRegion[i] = AliEmcalJet::kTPC;
    }
  }
  
  if (!typeOfJet) {
    typeOfJet = new AliJetContainer::EJetType_t[nJetContainer];
    for (Int_t i=0;i<nJetContainer;i++) {
      typeOfJet[i] = AliJetContainer::kChargedJet;
    }
  }
  
  Bool_t bDoFullJets = kFALSE;
	for (Int_t i=0;i<nJetContainer;i++) {
		bDoFullJets = bDoFullJets || (typeOfJet[i] == AliJetContainer::kFullJet);
  }

   //******************************************************************************
	
	Bool_t isMC = mgr->GetMCtruthEventHandler() != 0x0;

	if (isMC) {
    std::cout << "MCtruthEventHandler found" << std::endl;
    AddTaskMCTrackSelector("mcparticles", "usedefault", kFALSE, kFALSE, 1, kFALSE);   
	}	
  
	//Create jet Finder
	if (useJets) {
		AliEmcalJetTask** jetFinderTask = new AliEmcalJetTask*[nJetContainer];
			
		if (bDoFullJets) {
			AliEmcalCorrectionTask* correctionTask = AliEmcalCorrectionTask::AddTaskEmcalCorrectionTask();
			correctionTask->SelectCollisionCandidates(kPhysSel);
			correctionTask->SetUserConfigurationFilename(correctionTaskFileName.Data());
			correctionTask->Initialize();
		}
		
		for (Int_t i=0;i<nJetContainer;++i) {
			Bool_t isChargedJet = (typeOfJet[i] == AliJetContainer::kChargedJet);
			const char* particleName = "usedefault";
			const char* clustName = isChargedJet ? "" : "usedefault";
			Double_t clustCut = isChargedJet ? 0.0 : clustPtCut;
			jetFinderTask[i] = AddTaskEmcalJet(particleName, clustName, AliJetContainer::antikt_algorithm, absRadiusCut, typeOfJet[i], trackPtCut, clustCut, ghostArea, AliJetContainer::pt_scheme, "Jet", jetMinPt, kFALSE, kFALSE);
			jetFinderTask[i]->SelectCollisionCandidates(kPhysSel);
			jetFinderTask[i]->SetForceBeamType(iBeamType);
		}

		AliEmcalJetTask* jetFinderTaskMC = 0x0;
		if (isMC) {
			AliEmcalJetTask* jetFinderTaskMC = AddTaskEmcalJet("mcparticles", "", AliJetContainer::antikt_algorithm, absRadiusCut, AliJetContainer::kChargedJet, trackPtCut, clustPtCut, ghostArea, AliJetContainer::pt_scheme, "Jet", jetMinPt, kFALSE, kFALSE );  
			jetFinderTaskMC->SelectCollisionCandidates(kPhysSel);
			jetFinderTaskMC->SetForceBeamType(iBeamType);    
		}
	}  
   
   // Check the analysis type using the event handlers connected to the analysis manager.
   //==============================================================================
   if (!mgr->GetInputEventHandler()) {
   ::Error("AddTaskIDFragmentationFunction", "This task requires an input event handler");
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

   AliAnalysisTaskIDFragmentationFunction *task = new AliAnalysisTaskIDFragmentationFunction(
        Form("Fragmentation_Function_%s_%s_%s_%s_filterMaskTracks%d", branchRecJets.Data(), branchGenJets.Data(), typeJets.Data(),
             typeTracks.Data(), filterMaskTracks));
   
   if (debug>=0) 
		 task->SetDebugLevel(debug);

   task->SetUseNewCentralityEstimation(kTRUE);
   task->SetEventSelectionMask(kPhysSel);
   task->SetEventClass(eventClass);
	 
   // Set default parameters 
   // Cut selection 
   task->SetFFRadius(radius); 
	 task->SelectCollisionCandidates(kPhysSel);
   
   task->SetOnlyLeadingJets(onlyConsiderLeadingJets); // default: kFALSE
   
   task->SetMCPtHardCut(MC_pThard_cut);
	//Set Centrality Estimator (V0M ist default, V0A should be used for pPb). NoCentrality doesn't store the centrality information
	if (iBeamType != AliAnalysisTaskEmcal::kpA)
		task->SetCentralityEstimator("V0M");
	else
		task->SetCentralityEstimator("V0A");
   
   //TODO For PbPb: Search for Centrality estimator  
  
   task->SetUseFastSimulations(useFastSimulation);	
   if (useFastSimulation) {
    SetEfficiencyFunctionsFastSimulation(task, fastSimParamFile);
   }
// 	 task->SetDoGroomedJets(bDoJetGrooming);
// 	 task->SetStudyTransversalJetStructure(bStudyTransversalJetStructure);
// 	 task->SetMinMaxMultiplicity(kMinMultiplicity,kMaxMultiplicity);

   mgr->AddTask(task);
	 
	 	 
//   -------------------------------------------------------
//   Init track and cluster EMCAL containers
//   -------------------------------------------------------
	 
	//Setting up EMCAL jet framework
	 
	TString trackName;
	TString clusName;

	if (type == "ESD") {
		trackName = "Tracks";
		clusName = "CaloClusters";
	}
	else if (type == "AOD") {
		trackName = "tracks";
		clusName = "caloClusters";
	}
	else {
		trackName = "";
	} 
		
	task->SetVzRange(-10,10);
	task->SetNeedEmcalGeom(bDoFullJets);
	
	//Add Track Container to jet Task

	AliTrackContainer* trackCont = task->AddTrackContainer(trackName);
	trackCont->SetEtaLimits(kEtaMin,kEtaMax);
	trackCont->SetName("ReconstructedTracks");
	task->SetNameTrackContainer("ReconstructedTracks");
	
	trackCont->SetParticlePtCut(trackPtCut);
	trackCont->SetTrackFilterType(trackFilterType);
	trackCont->SetAODFilterBits(filterMaskTracks);
	
	AliClusterContainer* clusCont = 0x0;
	if (bDoFullJets) {
		clusCont = task->AddClusterContainer(clusName);
		clusCont->SetClusECut(clustPtCut);
		clusCont->SetExoticCut(kTRUE);
		clusCont->SetClusECut(0.);
		clusCont->SetClusPtCut(0.);
		clusCont->SetClusNonLinCorrEnergyCut(0.);
		clusCont->SetClusHadCorrEnergyCut(clustPtCut);
		clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr); 
	}
	 
	 AliJetContainer** jetContainer = new AliJetContainer*[nJetContainer];
	
	if (useJets) {
		TString jetContainerNames[nJetContainer];
		for (Int_t i=0;i<nJetContainer;++i) {
			Bool_t isChargedJet = (typeOfJet[i] == AliJetContainer::kChargedJet);
			jetContainer[i] = task->AddJetContainer(typeOfJet[i], AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, absRadiusCut, jetAcceptanceRegion[i], trackCont, isChargedJet ? 0x0 : clusCont);
			jetContainer[i]->SetMinPt(jetMinPt);
//       jetContainer[i]->SetEtaLimits(kEtaMin, kEtaMax); //Normally not needed, handled by the *fid region
			jetContainer[i]->SetPhiLimits(-10.0, 10.0);   //Full acceptance  
			jetContainer[i]->PrintCuts();
			jetContainerNames[i] = jetContainer[i]->GetName();
		}
		
		task->SetNameJetContainer(jetContainerNames[0]);
	}
		
	AliJetContainer* jetMCContainer = 0x0;
		
	AliTrackContainer* trackContainerEfficiency = 0x0;
	
	if (isMC) {
		AliMCParticleContainer* mcPartCont = task->AddMCParticleContainer("mcparticles");
		
		mcPartCont->SetParticlePtCut(trackPtCut);
		mcPartCont->SetParticleEtaLimits(kEtaMin,kEtaMax);
		mcPartCont->SetParticlePhiLimits(-10.0,10.0);
		mcPartCont->SetCharge(AliParticleContainer::kCharged);
		task->SetNameMCParticleContainer(mcPartCont->GetName());  
		
		trackContainerEfficiency = task->AddTrackContainer(trackName);
		trackContainerEfficiency->SetAODFilterBits(filterMaskTracks);
		trackContainerEfficiency->SetFilterHybridTracks(kTRUE);
		trackContainerEfficiency->SetParticlePtCut(0.0);
		trackContainerEfficiency->SetParticleEtaLimits(-100.0,100.0);
		trackContainerEfficiency->SetParticlePhiLimits(-10.0,10.0);
		trackContainerEfficiency->SetName("ReconstructedTracksEfficiency");
		task->SetNameTrackContainerEfficiency(trackContainerEfficiency->GetName());
		
		if (useJets) {
			jetMCContainer = task->AddJetContainer(AliJetContainer::kChargedJet, AliJetContainer::antikt_algorithm, AliJetContainer::pt_scheme, absRadiusCut, AliJetContainer::kTPCfid, mcPartCont,0x0);
			task->SetNameMCParticleJetContainer(jetMCContainer->GetName());
			jetMCContainer->SetMinPt(jetMinPt);
			jetMCContainer->SetEtaLimits(kEtaMin,kEtaMax);
			jetMCContainer->SetPhiLimits(-10.0,10.0);
		}
	}
		
	// Background
	if (subtractBackgroundPt && useJets && jetContainer[0]) {
		if (iBeamType != AliAnalysisTaskEmcal::kpp) {
			TString sRhoChName = "Rho";
			TString sRhoFuName = "Rho_Scaled";

			AliEmcalJetTask *pKtChJetTask = AddTaskEmcalJet("usedefault", "", AliJetContainer::kt_algorithm, absRadiusCut, AliJetContainer::kChargedJet, trackPtCut, 0, ghostArea, AliJetContainer::pt_scheme, "Jet", 1.0, kFALSE, kFALSE);
			pKtChJetTask->SelectCollisionCandidates(kPhysSel);
			
			AliAnalysisTaskRho* pRhoTask = AddTaskRhoNew("usedefault", "usedefault", sRhoChName, 0.4);
			pRhoTask->SetExcludeLeadJets(kNumberExcludeJetsInBackground);
			pRhoTask->SelectCollisionCandidates(kPhysSel);

//       if (bDoFullJets) {
//         TString sFuncPath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
//         TString sFuncName = "LHC11h_HadCorr20_ClustersV2";
//         pRhoTask->LoadRhoFunction(sFuncPath, sFuncName);
//       }
			
			for (Int_t i=0;i<nJetContainer;++i) {
				jetContainer[i]->SetRhoName(sRhoChName);
				jetContainer[i]->SetPercAreaCut(0.6);
			}
			if (jetMCContainer) {
				jetMCContainer->SetRhoName(sRhoChName);
				jetMCContainer->SetPercAreaCut(0.6);
			}
			
		}
	}
		
   // Create ONLY the output containers for the data produced by the task.
   // Get and connect other common input/output containers via the manager as below
   //==============================================================================

   TString strList(Form("idfracfunc_%s_%s_%s_%s_cl%d", branchRecJets.Data(), branchGenJets.Data(), typeTracks.Data(), typeJets.Data(), eventClass));
   
   TString strDir(Form("%s:PWGJE_IDFragmentationFunction_%s_%s_%s_%s_cl%d", 
           AliAnalysisManager::GetCommonFileName(), branchRecJets.Data(), branchGenJets. Data(), 
           typeTracks.Data(), typeJets.Data(), eventClass));


   if(FFMaxTrackPt>0){
     strList += Form("_FFMaxPt%d", FFMaxTrackPt);
     strDir  += Form("_FFMaxPt%d", FFMaxTrackPt);
   }
   if(FFMinNTracks>0){
     strList += Form("_minNTr%d",FFMinNTracks);
     strDir  += Form("_minNTr%d",FFMinNTracks);
   }

   if(radius<0){
     strList += "_trackRefs";
     strDir  += "_trackRefs";
   }
   if(filterMaskTracks){
     strList += Form("_TrackFilter%05d",filterMaskTracks);
     strDir  += Form("_TrackFilter%05d",filterMaskTracks);
   }
   
   AliAnalysisDataContainer *coutput_FragFunc = mgr->CreateContainer(strList,TList::Class(),
                     AliAnalysisManager::kOutputContainer,
                     strDir);   
   
   mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  if (mgr->GetCommonOutputContainer()) {
    //Not present for local runs
    mgr->ConnectOutput (task,  0, mgr->GetCommonOutputContainer());
  }
   mgr->ConnectOutput (task, 1, coutput_FragFunc);   
   
   postConfig(task, namesOfInclusivePIDtasks, namesOfJetPIDtasks, namesOfJetUEPIDtasks, namesOfJetUEMethods);
   
   return task;
}

