/// \file AddTaskGammaJetCorrelation.C
/// \ingroup CaloTrackCorrMacros
/// \brief gamma-jet correlation configuration.
///
/// Configuration of the gamma-jet correlation analysis
/// based on AddTaskIsoPhoton by Gustavo Conesa & Marie Germain.
///
/// \author Adam Matyja <Adam.Matyja@cern.ch>, INP-PAN-Krakow.
///

/// Global name to be composed of the settings, used to set the AOD branch name
TString kGammaJetCorrelationName = "";

///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param isoCone : A float setting the isolation cone size
/// \param isoPth : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param maxLambda0Cut : A float setting the maximum value of the shower shape of the clusters for the correlation analysis
/// \param maxNLMcut : maximum value of shower shape parameter
/// \param timecut : activate time cut
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle
/// \param simulation : A bool identifying the data as simulation
/// \param eventsel : reject bad events (pile-up ...)
/// \param exotic : reject exotic clusters
/// \param nonlin : A bool to set the use of the non linearity correction
/// \param collision : A string with the colliding system
/// \param trigger : A string with the trigger class, abbreviated, defined in method belowSetTriggerMaskFromName()
/// \param firedTrigger : In case of events with 2 L1 triggers, specify which one
/// \param clustersArray : A string with the array of clusters not being the default (default is empty string)
/// \param mix : A bool to switch the correlation mixing analysis
/// \param tm : A bool to select neutral clusters as triggers
/// \param minCen : An int to select the minimum centrality, -1 means no selection
/// \param maxCen : An int to select the maximum centrality, -1 means no selection
/// \param jetBranchName : Name of branch with reconstructed jets
/// \param jetBkgBranchName : Name of branch with reconstructed background jets
/// \param jetMinPt : Minimum jet pT.
/// \param minDeltaPhi : Minimum cut on photon-jet azimuthal angle
/// \param maxDeltaPhi : Maximum cut on photon-jet azimuthal angle
/// \param minPtRatio : Minimum cut on jet/photon pT ratio
/// \param maxPtRatio : Maximum cut on jet/photon pT ratio
/// \param debug : An int to define the debug level of all the tasks
/// \param printSettings : A bool to enable the print of the settings per task
/// \param scaleFactor : Scale factor in case for pT-hard simulation bins. Not useful with train.
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskGammaJetCorrelation (
                                const Float_t  isoCone       = 0.4,
								const Float_t  isoPth        = 0.5,
								const Double_t maxLambda0Cut = 0.5,
								const Int_t    maxNLMcut     = 2,
								const Bool_t   timecut       = kFALSE,
								const TString  calorimeter   = "EMCAL",
								const Bool_t   simulation    = kFALSE,
								const Bool_t   eventsel      = kFALSE,
								const Bool_t   exotic        = kTRUE,
								const Bool_t   nonlin        = kFALSE,
								const TString  collision     = "pp",
								const TString  trigger       = "MB",
								const TString  firedTrigger  = "EG1",
								const TString  clustersArray = "V1",
								const Bool_t   mix           = kTRUE,
								const Bool_t   tm            = kTRUE,
								const Int_t    minCen        = -1,
								const Int_t    maxCen        = -1,
								const TString  jetBranchName = "clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00",
								const TString  jetBkgBranchName = "jeteventbackground_clustersAOD_KT04_B0_Filter00768_Cut00150_Skip00",
                                const Float_t  jetMinPt      = 0,
                                const Float_t  minDeltaPhi   = 1.5,
								const Float_t  maxDeltaPhi   = 4.5,
								const Float_t  minPtRatio    = 0,
								const Float_t  maxPtRatio    = 5,   
								const Int_t    debug         = -1,
								const Bool_t   printSettings = kFALSE,
								const Double_t scaleFactor   = -1      )
{
  // Get the pointer to the existing analysis manager via the static access method.
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) 
  {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  
  if (!mgr->GetInputEventHandler()) 
  {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }
  
  // Make sure the B field is enabled for track selection, some cuts need it
  ((AliInputEventHandler*)mgr->GetInputEventHandler())->SetNeedField(kTRUE);
  
  // input jets
  TString     deltaAODJetName   = "AliAOD.Jets.root"; //Jet input AOD name
  if(deltaAODJetName.Length()!=0)
  {
    // External file with Jets
    // aodHandler->AddFriend(deltaAODJetName.Data());
    mgr->RegisterExtraFile(deltaAODJetName.Data());
    cout<<"Jet file registered "<<endl;
    cout<<"Extra files: "<<mgr->GetExtraFiles()<<endl;
  }

  Bool_t useKinematics = kFALSE;
  useKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Name for containers
  
  //kGammaJetCorrelationName = Form("%s_Trig%s_Cl%s_TM%d_R%1.1f_Pt%1.1f",calorimeter.Data(), trigger.Data(),clustersArray.Data(),tm,cone,pth);
  //  kGammaJetCorrelationName = Form("%s_Trig%s_Cl%s_TM%d",calorimeter.Data(), trigger.Data(),clustersArray.Data(),tm);
  kGammaJetCorrelationName = Form("%s_Trig%s_Fired%s_Cl%s_TM%d_l02%1.2f",calorimeter.Data(), trigger.Data(),firedTrigger.Data(),clustersArray.Data(),tm,maxLambda0Cut);//<<<---changed here

  if(collision=="PbPb" && maxCen>=0) kGammaJetCorrelationName+=Form("Cen%d_%d",minCen,maxCen);
    
  printf("<<<< NAME: %s >>>>>\n",kGammaJetCorrelationName.Data());
  
  // #### Configure analysis ####
    
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  printf("SCALE FACTOR %e\n",scaleFactor);
  maker->SetScaleFactor(scaleFactor); // for MC, negative (not scaled) by default
  
  // General frame setting and configuration
  maker->SetReader   (ConfigureReader(mgr->GetInputEventHandler()->GetDataType(),calorimeter,useKinematics,simulation,eventsel,nonlin,timecut,collision,trigger,firedTrigger,clustersArray,jetBranchName,jetBkgBranchName,mix,minCen,maxCen,debug,printSettings)   ); 
  maker->SetCaloUtils(ConfigureCaloUtils(clustersArray,collision,nonlin,exotic,simulation,timecut,debug,printSettings)); 
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important

  // Isolation settings
  Int_t partInCone = AliIsolationCut::kNeutralAndCharged; // kOnlyCharged;
  Int_t thresType  = AliIsolationCut::kPtThresIC;//  AliIsolationCut::kSumPtFracIC ; 
  //Float_t isoCone = -1;
  //Float_t isoPth  = -1;
  
  maker->AddAnalysis(ConfigurePhotonAnalysis(calorimeter,tm,simulation,maxLambda0Cut,maxNLMcut,debug,printSettings), n++); // Photon cluster selection
  maker->AddAnalysis(ConfigureIsolationAnalysis(calorimeter,collision,"Photon", partInCone,thresType, isoCone, isoPth,tm,kFALSE,simulation,debug,printSettings), n++); // Photon isolation   
  maker->AddAnalysis(ConfigurePhotonJetAnalysis(calorimeter,isoCone,jetMinPt,minDeltaPhi,maxDeltaPhi,minPtRatio,maxPtRatio,simulation,debug,printSettings), n++);// photon-jet correlation analysis

  maker->SetAnaDebug(debug)  ;
  maker->SwitchOnHistogramsMaker()  ;
  //if(kData.Contains("delta")) maker->SwitchOffAODsMaker() ;
  //else                        maker->SwitchOnAODsMaker()  ;
  
  if(printSettings) maker->Print("");
  
  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, calorimeter.Data());
  // Create task
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("CaloTrackCorr%s",kGammaJetCorrelationName.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(debug);
  task->SetBranches("ESD:AliESDRun.,AliESDHeader"); 
  task->SetAnalysisMaker(maker);
  mgr->AddTask(task);
  
  //Create containers
  
  //if(outputfile.Length()==0) outputfile = AliAnalysisManager::GetCommonFileName(); 
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(kGammaJetCorrelationName, TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s",outputfile.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",kGammaJetCorrelationName.Data()), TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             "AnalysisParameters.root");
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  //if(!kData.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
    
  if(!mix)
  {    
    UInt_t mask =  SetTriggerMaskFromName(trigger);
    task->SelectCollisionCandidates(mask);
  } 
  
  return task;
}

///
/// Configure the class handling the events and cluster/tracks filtering.
///
AliCaloTrackReader * ConfigureReader(TString inputDataType = "AOD",TString calorimeter = "EMCAL",Bool_t useKinematics = kFALSE,
				     Bool_t simulation = kFALSE,Bool_t   eventsel      = kFALSE,Bool_t nonlin = kTRUE, Bool_t timecut = kFALSE,
				     TString  collision     = "pp",TString trigger="MB",TString firedTrigger="EG1",
				     TString  clustersArray = "V1", TString  jetBranchName = "jets", TString  jetBkgBranchName = "jets",
				     Bool_t  mix           = kFALSE,
				     Float_t minCen = -1, Float_t maxCen = -1,
				     Int_t debug = -1,Bool_t printSettings = kFALSE)
{
  Bool_t useTender=kTRUE;

  if(simulation)
  {
    if (!useKinematics && inputDataType=="AOD") useKinematics = kTRUE; //AOD primary should be available ...
  }
  
  cout<<"********* ACCESS KINE? "<<useKinematics<< endl;



  AliCaloTrackReader * reader = 0;
  if     (inputDataType == "AOD") reader = new AliCaloTrackAODReader();
  else if(inputDataType == "ESD") reader = new AliCaloTrackESDReader();
  else printf("AliCaloTrackReader::ConfigureReader() - Data not known InputData=%s\n",inputDataType.Data());
  
  reader->SetDebug(debug);//10 for lots of messages

  if(useTender){
    reader->SwitchOffTriggerPatchMatching();
    reader->SwitchOffBadTriggerEventsRemoval();
  }

  if(simulation)
  {
    // Event rejection cuts for jet-jet simulations, do not use in other
    reader->SetPtHardAndJetPtComparison(kTRUE);
    reader->SetPtHardAndJetPtFactor(4);
    
    reader->SetPtHardAndClusterPtComparison(kTRUE);
    reader->SetPtHardAndClusterPtFactor(1.5);
  }

  //Delta AOD?
  //reader->SetDeltaAODFileName("");
  //if(kOutputAOD) reader->SwitchOnWriteDeltaAOD()  ;
  
  // MC settings
  if(useKinematics){
    if(inputDataType == "ESD"){
      reader->SwitchOnStack();          
      reader->SwitchOffAODMCParticles(); 
    }
    else if(inputDataType == "AOD"){
      reader->SwitchOffStack();          
      reader->SwitchOnAODMCParticles(); 
    }
  }  
  
  //------------------------
  // Detector input filling
  //------------------------
  
  //Min cluster/track E
  reader->SetEMCALEMin(0.3); 
  //reader->SetEMCALEMin(0.);// <<<----changed here 
  reader->SetEMCALEMax(1000); 
  reader->SetPHOSEMin(0.3);
  reader->SetPHOSEMax(1000);
  //  reader->SetCTSPtMin(0.2);
  reader->SetCTSPtMin(0.15);// <<<----changed here 
  reader->SetCTSPtMax(1000);

  //-----------------------------------------------------------------
  // Jet part
  //-----------------------------------------------------------------
  reader->SwitchOnNonStandardJets();
  //reader->SetInputNonStandardJetBranchName("clustersAOD_ANTIKT04_B0_Filter00768_Cut00150_Skip00");
  //reader->SetInputNonStandardJetBranchName("clustersAOD_ANTIKT04_B1_Filter00768_Cut00150_Skip00");//in PbPb
  //reader->SetInputNonStandardJetBranchName("clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00");//in pp 2.76 LHC11a,7 LHC11c
  //reader->SetInputNonStandardJetBranchName("clustersAOD_ANTIKT04_B0_Filter00768_Cut00150_Skip02");//in pp 7 LHC13e4 MC
  //reader->SetInputNonStandardJetBranchName("clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00");//in pp 7 LHC12a15f MC
  reader->SetInputNonStandardJetBranchName(jetBranchName.Data());
  if(jetBkgBranchName.Length()!=0) {
    reader->SwitchOnBackgroundJets();
    reader->SetInputBackgroundJetBranchName(jetBkgBranchName.Data());
  } else {
    reader->SwitchOffBackgroundJets();
  }
  //reader->SetInputBackgroundJetBranchName("jeteventbackground_clustersAOD_KT04_B0_Filter00768_Cut00150_Skip00");//in pp 7 LHC13e4 MC


  // Time cuts
  if(simulation) 
  {
    reader->SwitchOffUseTrackTimeCut();
    reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
    reader->SwitchOffUseEMCALTimeCut() ;   
    reader->SwitchOffUseParametrizedTimeCut();
  }
  else
  {
    reader->SwitchOffUseParametrizedTimeCut();
    if(timecut)
     {
	printf("Set time cut \n");
	reader->SwitchOnUseEMCALTimeCut();
	//Absolute window
	reader->SetEMCALTimeCut(-30.,30.); // default is -25ns-20ns
      }
    else
      {
	printf("Off time cuts time cut \n");
	reader->SwitchOffUseEMCALTimeCut();
	//Absolute window
	reader->SetEMCALTimeCut(-1.e6,1.e6);
      }

    //if(kCalibT)
    //{ 
      //printf("Set time cut parameters for run %d\n",kRunNumber);
      //reader->SetEMCALTimeCut(-20,20); 
      //reader->SwitchOnUseParametrizedTimeCut();
      //if     (kRunNumber >= 151636 && kRunNumber <= 155384 )
      //{
      //  printf("Set time parameters for LHC11c\n");
      //  reader->SetEMCALParametrizedMinTimeCut(0,-5  ); reader->SetEMCALParametrizedMinTimeCut(1,-1 ); reader->SetEMCALParametrizedMinTimeCut(2, 1.87); reader->SetEMCALParametrizedMinTimeCut(3, 0.4);   
      //  reader->SetEMCALParametrizedMaxTimeCut(0, 3.5); reader->SetEMCALParametrizedMaxTimeCut(1, 50); reader->SetEMCALParametrizedMaxTimeCut(2, 0.15); reader->SetEMCALParametrizedMaxTimeCut(3, 1.6);   
      //}
      //else if(kRunNumber >= 156447 && kRunNumber <= 159635 )
      //{
      //  printf("Set time parameters for LHC11d\n");
      //  reader->SetEMCALParametrizedMinTimeCut(0,-5);  reader->SetEMCALParametrizedMinTimeCut(1,-1 );  reader->SetEMCALParametrizedMinTimeCut(2, 3.5 ); reader->SetEMCALParametrizedMinTimeCut(3, 1.  );   
      //  reader->SetEMCALParametrizedMaxTimeCut(0, 5);  reader->SetEMCALParametrizedMaxTimeCut(1, 50);  reader->SetEMCALParametrizedMaxTimeCut(2, 0.45); reader->SetEMCALParametrizedMaxTimeCut(3, 1.25);   
      //}
      //else 
      //{
      //  printf("*** Fixed time cut 20 ns *** \n");
      //  reader->SetEMCALTimeCut(-20,20);
      //}
//      
//    }
//    else
//    {
//      reader->SwitchOffUseEMCALTimeCut();
//      reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
//      reader->SwitchOffUseParametrizedTimeCut(); 
//    }

  }  
  
  reader->SwitchOnFiducialCut();
  //reader->SwitchOffFiducialCut();// <<<----changed here 
  //reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;
  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360) ;// <<<--- changed here

  // Tracks

  //  reader->SwitchOffCTS();//here changed 0n->off
  reader->SwitchOnCTS();//here changed 0n->off

  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();
  
  reader->SwitchOffUseTrackTimeCut();
  reader->SetTrackTimeCut(0,50);

  if(inputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
    //Hybrids 2011
    AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
    reader->SetTrackCuts(esdTrackCuts);
    AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
    reader->SetTrackComplementaryCuts(esdTrackCuts2);
    
    //Hybrids 2010
    //AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001006);
    //reader->SetTrackCuts(esdTrackCuts);
    //AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10041006);
    //reader->SetTrackComplementaryCuts(esdTrackCuts2);
  }
  else if(inputDataType=="AOD")
  {
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SetTrackStatus(AliVTrack::kITSrefit);    
    //reader->SetTrackFilterMask(128);           // Filter bit, not mask, use if off hybrid
  }
  
  // Calorimeter
  
  reader->SetEMCALClusterListName(clustersArray);
  if(clustersArray == "" && !useTender)
  {
    printf("**************** Standard EMCAL clusters branch analysis **************** \n");
    reader->SwitchOnClusterRecalculation();
    // Check in ConfigureCaloUtils that the recalibration and bad map are ON 
  }
  else 
  {
    printf("**************** Input for analysis is Clusterizer %s **************** \n", clustersArray.Data());
    reader->SwitchOffClusterRecalculation();
  }  
  
  if(!nonlin) reader->SwitchOffClusterELinearityCorrection();
  else        reader->SwitchOnClusterELinearityCorrection();

  if(calorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();  
    reader->SwitchOnEMCAL();
  }
  if(calorimeter == "PHOS") { 
    reader->SwitchOnPHOSCells();  
    reader->SwitchOnPHOS();
  }
  
  // for case data="deltaAOD", no need to fill the EMCAL/PHOS cluster lists
  //if(kData.Contains("delta"))
  //{
  //  reader->SwitchOffEMCAL();
  //  reader->SwitchOffPHOS();
  //  reader->SwitchOffEMCALCells(); 
  //  reader->SwitchOffPHOSCells(); 
  //}
  
  //-----------------
  // Event selection
  //-----------------
  
  //if(!useKinematics) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
  if(!useKinematics) {
    if(collision =="pPb" && trigger=="EMCEGA") {
      reader->SetFiredTriggerClassName(firedTrigger);
    }
  }

  // For mixing with AliAnaParticleHadronCorrelation switch it off
  if(mix)
  {
    reader->SwitchOffEventTriggerAtSE();
    UInt_t mask =  SetTriggerMaskFromName(trigger);
    reader->SetEventTriggerMaks(mask); // Only for mixing and SwitchOffEventTriggerAtSE();
    //reader->SetMixEventTriggerMaks(AliVEvent::kMB); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
    reader->SetMixEventTriggerMaks(AliVEvent::kAnyINT); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT

    printf("---Trigger selection done in AliCaloTrackReader!!!\n");
  }
  else 
    reader->SwitchOnEventTriggerAtSE();
  
  reader->SetZvertexCut(10.);                // Open cut
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
  reader->SwitchOnRejectNoTrackEvents();//<<<--- changed here, new added

  if(eventsel)
  {
    reader->SwitchOnPileUpEventRejection();   // remove pileup by default  
    reader->SwitchOnV0ANDSelection() ;        // and besides v0 AND
  }
  else 
  {
    reader->SwitchOffPileUpEventRejection();  // remove pileup by default   
    reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
  }
    
  if(collision=="PbPb") 
  {
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(10);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
    reader->SetCentralityBin(minCen,maxCen); // Accept all events, if not select range
    
    // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
    reader->SetEventPlaneMethod("V0");
  }
  
  if(printSettings) reader->Print("");
  
  return reader;
}

///
/// Configure the class handling the calorimeter clusters specific methods
///
AliCalorimeterUtils* ConfigureCaloUtils( TString  clustersArray = "V1",TString  collision     = "pp",Bool_t nonlin = kTRUE,Bool_t exotic = kTRUE ,Bool_t simulation = kFALSE,Bool_t timecut = kFALSE,Int_t debug = -1,Bool_t printSettings = kFALSE)
{
  Bool_t useTender=kTRUE;

  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(debug);
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  //cu->SetNumberOfCellsFromEMCALBorder(0);// <<<----changed here 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  // Search of local maxima in cluster
  if(collision=="pp" || collision=="pPb")
  {
    cu->SetLocalMaximaCutE(0.1);
    cu->SetLocalMaximaCutEDiff(0.03);
  }
  else 
  {
    cu->SetLocalMaximaCutE(0.2);
    cu->SetLocalMaximaCutEDiff(0.03);
  }
  
  cu->SwitchOffClusterPlot();

  //if(kRecalTM) cu->SwitchOnRecalculateClusterTrackMatching(); // Done in clusterization
  //else         
  cu->SwitchOffRecalculateClusterTrackMatching();
  
  if(useTender)
    cu->SwitchOffBadChannelsRemoval() ;
  else
    cu->SwitchOnBadChannelsRemoval() ;

  //EMCAL settings
  if(!simulation)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();

  if(!useTender){  
    AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
    
    if(!simulation)
      {
	cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
	//if(clustersArray == "" && !kTender) cu->SwitchOnRunDepCorrection(); 
	if(clustersArray == "") cu->SwitchOnRunDepCorrection(); 
      }
    
    cu->SwitchOnEMCALOADB();//FIX ME!!!

    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
    ConfigureEMCALRecoUtils(recou,
			    simulation,                             
			    exotic,
			    nonlin,
			    kFALSE, // e calib
			    kFALSE, // bad map
			    kFALSE); // time calib
    //kCalibE, 
    //kBadMap,
    //kCalibT);   
    recou->SetExoticCellDiffTimeCut(1e10);
    if(timecut) recou->SetExoticCellDiffTimeCut(50.);
  }//end tender


  if( nonlin ) 
  { 
//    printf("ConfigureCaloUtils() - Apply non linearity to EMCAL\n");
//    //CAREFUL only for the latest simulation
//    recou->SetNonLinearityFunction(AliEMCALRecoUtils::kBeamTestCorrected);
//    recou->SetNonLinearityParam(0,9.81039e-01);
//    recou->SetNonLinearityParam(1,1.13508e-01);
//    recou->SetNonLinearityParam(2,1.00173e+00); 
//    recou->SetNonLinearityParam(3,9.67998e-02);
//    recou->SetNonLinearityParam(4,2.19381e+02);
//    recou->SetNonLinearityParam(5,6.31604e+01);
//    recou->SetNonLinearityParam(6,1);
    printf("*** SET cluster non linearity correction ***\n");
    cu->SwitchOnCorrectClusterLinearity();
  }
    
  if(!useTender){
    printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
    printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
  }

  cu->SetNumberOfSuperModulesUsed(10);
  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  if(printSettings) cu->Print("");
  
  return cu;
}

///
/// Configure the task doing the first photon cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection ...
///
AliAnaPhoton* ConfigurePhotonAnalysis(TString calorimeter = "EMCAL",Bool_t tm = kFALSE,Bool_t simulation = kFALSE,Double_t maxLambda0Cut=0.5,Int_t maxNLMcut=2,Int_t debug = -1,Bool_t printSettings = kFALSE)
{
  AliAnaPhoton *ana = new AliAnaPhoton();
  ana->SetDebug(debug); //10 for lots of messages
  
  // cluster selection cuts
  
  ana->SwitchOffFiducialCut();

  ana->SetCalorimeter(calorimeter);
  
  if(calorimeter == "PHOS")
  {
    ana->SetNCellCut(2);// At least 3 cells
    ana->SetMinPt(0.3);
    ana->SetMinDistanceToBadChannel(2, 4, 5);
    ana->SetTimeCut(-1e10,1e10); // open cut
  }
  else 
  {//EMCAL
    ana->SetNCellCut(1);// At least 2 cells
    ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
    ana->SetMaxEnergy(1000); 
    //ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off //<<<---modified here
    // restrict to less than 100 ns when time calibration is on 
    ana->SetMinDistanceToBadChannel(2, 4, 6); 
    // Not useful if M02 cut is already strong
    ana->SetNLMCut(1, maxNLMcut) ;//[1,2]
    //ana->SetNLMCut(1, 10) ;//<<<----changed here 
    //ana->SetNLMCut(1, 1) ;//<<<----changed here 
  }
  
  if(tm)
  {
    ana->SwitchOnTrackMatchRejection() ;
    ana->SwitchOffTMHistoFill() ;
  }
  else
  {
    ana->SwitchOffTrackMatchRejection() ;
    ana->SwitchOnTMHistoFill() ;
  }
  
  
  //PID cuts (shower shape)
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
  //Not used in bayesian
  
  //EMCAL
  caloPID->SetEMCALLambda0CutMax(maxLambda0Cut);//0.27 was before//0.50//<<<----changed here
  caloPID->SetEMCALLambda0CutMin(0.10);
  
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
    
  //PHOS
  caloPID->SetPHOSDispersionCut(2.5);
  caloPID->SetPHOSRCut(2.);
  //if(kInputData=="AOD") caloPID->SetPHOSRCut(2000.); // Open cut since dX, dZ not stored
      
  ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection <<<--- changed here
  //ana->SwitchOffFillShowerShapeHistograms();  // Filled before photon shower shape selection
  ana->SwitchOffFillPileUpHistograms();

  // Input / output delta AOD settings
  
  //if(!kData.Contains("delta")) 
  //{
  ana->SetOutputAODName(Form("Photon%s",kGammaJetCorrelationName.Data()));
  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  //  //ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
  //}
  //else 
  ana->SetInputAODName(Form("Photon%s",kGammaJetCorrelationName.Data()));
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaPhoton_TM%d_",tm));
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms(20);
  ana->FillNPrimaryHistograms(20);
  
  ConfigureMC(ana,simulation);
  
  if(printSettings) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the trigger cluster isolation
///
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString calorimeter = "EMCAL",
						    TString  collision     = "pp",
						    TString particle="Photon", 
                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                    Float_t cone = 0.3,
                                                    Float_t pth  = 0.3,
						    Bool_t tm = kFALSE,
                                                    Bool_t multi      = kFALSE,Bool_t simulation = kFALSE,
						    Int_t debug = -1,
						    Bool_t printSettings = kFALSE)
{
  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
  ana->SetDebug(debug);
    
  ana->SwitchOnFiducialCut();
  //Avoid borders of EMCal
  if(calorimeter=="EMCAL")
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;

  // Same Eta as EMCal, cut in phi if EMCAL was triggering
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    //if(trigger.Contains("EMC"))
    //  ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
    //else
      ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;    
  }
  
  ana->SetMinPt(10);//<<---changed here
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s",particle.Data(),kGammaJetCorrelationName.Data()));
  ana->SetAODObjArrayName(Form("IC%s_%s",particle.Data(),kGammaJetCorrelationName.Data())); 
  
  ana->SetCalorimeter(calorimeter);
  
  if(!tm)  ana->SwitchOnTMHistoFill();
  else      ana->SwitchOffTMHistoFill();
  
  //if(particle=="Photon")ana->SwitchOnSSHistoFill();
  //else                  ana->SwitchOffSSHistoFill();
  
  ana->SwitchOffSSHistoFill();
  ana->SwitchOffFillPileUpHistograms();

  //Do settings for main isolation cut class
  AliIsolationCut * ic =  ana->GetIsolationCut();	
  ic->SetDebug(debug);
  
  if(cone >0 && pth > 0)
  {
    ic->SetPtThreshold(pth);
    ic->SetConeSize(cone);
  }
  else
  {
    if(collision=="pp") 
    {
      ic->SetPtThreshold(1.);//<<---changed here was 0.5,1
      ic->SetConeSize(0.3);//<<---changed here was 0.4
    }
    if(collision=="pPb")
      {
	ic->SetPtThreshold(1.0);
	ic->SetConeSize(0.3);
      }

    if(collision=="PbPb")
    {
      ic->SetPtThreshold(3.);
      //ic->SetPtThreshold(1.);
      ic->SetConeSize(0.3);
    }
  }
  
  ic->SetPtFraction(0.1);
  ic->SetSumPtThreshold(1.0) ;
  ic->SetParticleTypeInCone(partInCone);
  ic->SetICMethod(thresType);
  
  //Do or not do isolation with previously produced AODs.
  //No effect if use of SwitchOnSeveralIsolation()
  ana->SwitchOffReIsolation();
  
  //Multiple IC
  if(multi) 
  {
    ic->SetConeSize(1.);    // Take all for first iteration
    ic->SetPtThreshold(100);// Take all for first iteration
    ana->SwitchOnSeveralIsolation() ;
    ana->SetAODObjArrayName(Form("MultiIC%sTM%d",particle.Data(),tm));
     
    ana->SetNCones(4);
    ana->SetNPtThresFrac(4);
    ana->SetConeSizes(0,0.3);       ana->SetConeSizes(1,0.4);
    ana->SetConeSizes(2,0.5);       ana->SetConeSizes(3,0.6);
    ana->SetPtThresholds(0, 0.5);   ana->SetPtThresholds(1, 1);     ana->SetPtThresholds(2, 2);
    ana->SetPtFractions (0, 0.05) ; ana->SetPtFractions (1, 0.1);   ana->SetPtFractions (2, 0.2) ;  ana->SetPtFractions (3, 0.3) ;
    ana->SetSumPtThresholds(0, 1) ; ana->SetSumPtThresholds(1, 3) ; ana->SetSumPtThresholds(2, 5);  ana->SetSumPtThresholds(3, 7)  ;
    
    ana->SwitchOffTMHistoFill();
    ana->SwitchOffSSHistoFill();
  }
  else      
    ana->SwitchOffSeveralIsolation() ;
  
  AliCaloPID* caloPID = ana->GetCaloPID();
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
  
  //Set Histograms name tag, bins and ranges
  
  if(!multi)ana->AddToHistogramsName(Form("AnaIsol%s_TM%d_",particle.Data(),tm));
  else      ana->AddToHistogramsName(Form("AnaMultiIsol%s_TM%d_",particle.Data(),tm));

  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below
  
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  ConfigureMC(ana,simulation);
  
  if(printSettings) ic ->Print("");
  if(printSettings) ana->Print("");
  
  return ana;
  
}

///
/// Configure the selection of MC events
///
void ConfigureMC(AliAnaCaloTrackCorrBaseClass* ana,Bool_t simulation = kFALSE)
{
  if(simulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else            ana->SwitchOffDataMC() ;

  //Set here generator name, default pythia
  //ana->GetMCAnalysisUtils()->SetMCGenerator("");
}  

///
/// Set common histograms binning and ranges
///
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges,TString calorimeter = "EMCAL")
{
  histoRanges->SetHistoPtRangeAndNBins(-0.25, 99.75, 200) ; // Energy and pt histograms
  
  if(calorimeter=="EMCAL")
  {
    histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
    histoRanges->SetHistoXRangeAndNBins(-600,90,200); // QA
    histoRanges->SetHistoYRangeAndNBins(100,450,100); // QA

    histoRanges->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  }
  else
  {
    histoRanges->SetHistoPhiRangeAndNBins(260*TMath::DegToRad(), 320*TMath::DegToRad(), 60) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
  }
  
  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 4.9, 500);
  
  // Invariant mass histoRangeslysis
  histoRanges->SetHistoMassRangeAndNBins(0., 1., 200) ;
  histoRanges->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  
  // check if time calibration is on
  //histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  histoRanges->SetHistoTimeRangeAndNBins(-400.,400,400);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistodRRangeAndNBins(0.,0.15,150);//QA

  // QA, electron, charged
  histoRanges->SetHistoPOverERangeAndNBins(0,10.,100);
  histoRanges->SetHistodEdxRangeAndNBins(0.,200.,200);
  
  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,100);
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  histoRanges->SetHistoNClusterCellRangeAndNBins(0,500,500);
  histoRanges->SetHistoZRangeAndNBins(-400,400,200);
  histoRanges->SetHistoRRangeAndNBins(400,450,25);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,5000,500);
  
  // Isolation
  histoRanges->SetHistoPtInConeRangeAndNBins(0, 50 , 250);
  histoRanges->SetHistoPtSumRangeAndNBins   (0, 100, 250);
}

///
/// Set the trigger requested for the analysis,
/// depending on a string given
///
UInt_t SetTriggerMaskFromName(TString trigger)
{
  if(trigger=="EMC7")
  {
    printf("CaloTrackCorr trigger EMC7\n");
    return AliVEvent::kEMC7;
  }
  else if (trigger=="INT7")
  {
    printf("CaloTrackCorr trigger INT7\n");
    return AliVEvent::kINT7;
  }
  else if(trigger=="EMC1")
  {
    printf("CaloTrackCorr trigger EMC1\n");
    return AliVEvent::kEMC1;
  }
  else if(trigger=="MB")
  {
    printf("CaloTrackCorr trigger MB\n");
    return AliVEvent::kMB;
  }  
  else if(trigger=="PHOS")
  {
    printf("CaloTrackCorr trigger PHOS\n");
    return AliVEvent::kPHI7;
  }  
  else if(trigger=="PHOSPb")
  {
    printf("CaloTrackCorr trigger PHOSPb\n");
    return AliVEvent::kPHOSPb;
  }
  else if(trigger=="AnyINT")
  {
    printf("CaloTrackCorr trigger AnyINT\n");
    return AliVEvent::kAnyINT;
  }  
  else if(trigger=="INT")
  {
    printf("CaloTrackCorr trigger AnyINT\n");
    return AliVEvent::kAny;
  }
  else if(trigger=="EMCEGA")
  {
    printf("CaloTrackCorr trigger EMC Gamma\n");
    return AliVEvent::kEMCEGA;
  } 
  else if(trigger=="EMCEJE")
  {
    printf("CaloTrackCorr trigger EMC Jet\n");
    return AliVEvent::kEMCEJE;
  }
  else if(trigger=="Central")
  {
    printf("CaloTrackCorr trigger Central\n");
    return AliVEvent::kCentral;
  } 
  else if(trigger=="CentralEGA")
  {
    printf("CaloTrackCorr trigger Central+EMCEGA\n");
    return (AliVEvent::kCentral | AliVEvent::kEMCEGA);
  } 
  else if(trigger=="SemiCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral\n");
    return AliVEvent::kSemiCentral;
  }
  else if(trigger=="SemiOrCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral Or Central\n");
    return (AliVEvent::kSemiCentral | AliVEvent::kCentral);
  }  
  else if(trigger=="SemiOrCentralOrAnyINT")
    {
      printf("CaloTrackCorr trigger SemiCentral Or Central Or AnyINT\n");
      return (AliVEvent::kSemiCentral | AliVEvent::kCentral | AliVEvent::kAnyINT);
    }

}

///
/// Configure the task doing the trigger cluster-jet correlation
///
AliAnaParticleJetFinderCorrelation* ConfigurePhotonJetAnalysis(TString calorimeter = "EMCAL",Float_t gammaConeSize = 0.3, Float_t  jetMinPt  = 0, 
							       Float_t  minDeltaPhi   = 1.5,Float_t  maxDeltaPhi   = 4.5,
							       Float_t  minPtRatio    = 0,Float_t  maxPtRatio    = 5,
							       Bool_t simulation = kFALSE,Int_t debug = -1,Bool_t printSettings = kFALSE)
{
  AliAnaParticleJetFinderCorrelation *ana = new AliAnaParticleJetFinderCorrelation();
  ana->SetDebug(debug);
  TString particle="Photon";
  ana->SetInputAODName(Form("%s%s",particle.Data(),kGammaJetCorrelationName.Data()));

  ana->SwitchOffFiducialCut();

  ana->SetConeSize(0.4); //was 1 - cone to calculate FF
  ana->SelectIsolated(kTRUE); // do correlation with isolated photons <<---changed here
  ana->SetMakeCorrelationInHistoMaker(kFALSE);
  ana->SetPtThresholdInCone(0.150);//<<---- change here
  //ana->SetDeltaPhiCutRange(TMath::Pi()/2.,TMath::Pi()*3./2.);//Mostly Open Cuts 
  ana->SetDeltaPhiCutRange(minDeltaPhi,maxDeltaPhi);  // Delta phi cut for correlation
  ana->SetJetConeSize(0.4);//jet cone size / check the reco jet name
  ana->SetJetMinPt(jetMinPt);//min jet pt
  ana->SetJetAreaFraction(0.8);//min area fraction was 0.6
  ana->SetMinPt(0.3);//min cluster pt repeated from reader
  ana->SetGammaConeSize(gammaConeSize);//isolation cone repeated from isolation ana
  //ana->SetRatioCutRange(0.01,5.); //Mostly Open Cuts //0.01-5//<<---- change here
  ana->SetRatioCutRange(minPtRatio,maxPtRatio); // Delta pt cut for correlation

  ana->UseJetRefTracks(kTRUE); //Working now
  //Set Histograms bins and ranges
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below 0,100,200
  //ana->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;

  //  ana->SwitchOnNonStandardJetFromReader();
  ana->SwitchOnBackgroundJetFromReader();
  //background subtraction for photons
  //ana->SwitchOnBackgroundSubtractionGamma();
  ana->SwitchOffBackgroundSubtractionGamma();

  ana->SwitchOnSaveGJTree();
  ana->SwitchOnMostOpposite();
  //ana->SwitchOnMostEnergetic();


  //if(useKinematics) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  if(printSettings) 
    ana->Print("");

  return ana;

}
