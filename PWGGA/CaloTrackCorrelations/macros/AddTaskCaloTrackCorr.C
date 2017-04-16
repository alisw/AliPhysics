/// \file AddTaskCaloTrackCorr.C
/// \ingroup CaloTrackCorrMacros
/// \brief Example of configuration of CaloTrackCorrelation package.
///
/// Example of configuration of different analysis combinations
/// of the package CaloTrackCorrelations.
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

/// Global variables to be accessed by the different methods
Bool_t  kPrint         = kFALSE;    ///< Print setted parameters when configuring
Bool_t  kSimulation    = kFALSE;    ///< Declare the analysis simulation
Bool_t  kUseKinematics = kFALSE;    ///< Use the MC information
Bool_t  kOutputAOD     = kFALSE;    ///< Create output AOD with generated particle AOD objects
Bool_t  kEventSelection= kFALSE;    ///< Remove bad events
Bool_t  kExotic        = kTRUE;     ///< Remove exotic clusters
Bool_t  kNonLinearity  = kFALSE;    ///< Correct cluster non linearity
Int_t   kYears         = 2011;      ///< Declare the year of the data
TString kCollisions    = "pp";      ///< Declare the collision type of the data
TString kTrig          = "EMC7" ;   ///< Set the trigger type to analyze in data
TString kClusterArray  = "";        ///< Name of branch with clusters, default none, standard branch
TString kData          = "";        ///< Declare data MC or deltaAOD
TString kInputDataType = "ESD";     ///< Declare data ESD/AOD
TString kCalorimeter   = "EMCAL";   ///< Use main analysis detector EMCal or PHOS or CTS
Bool_t  kTM            = kTRUE;     ///< Remove matched clusters to tracks
Bool_t  kRecalTM       = kTRUE;     ///< Recalculate track-cluster matching
Int_t   kMinCen        = -1;        ///< Set the minimum centrality to be analyzed
Int_t   kMaxCen        = -1;        ///< Set the maximum centrality to be analyzed
TString kName          = "";        ///< Name of the analysis, used in created AOD branches and histo container
Int_t   kDebug         = -1;        ///< Do the analysis with this debug level
Bool_t  kQA            = kFALSE;    ///< Execute the calorimeter QA analaysis
Bool_t  kHadronAN      = kFALSE;    ///< Execute the hadron selection and correlation analysis
Bool_t  kCalibE        = kTRUE;     ///< Calibrate energy of clusters
Bool_t  kCalibT        = kTRUE;     ///< Calibrate time of clusters
Bool_t  kBadMap        = kTRUE;     ///< Reject bad cells/clusters
Bool_t  kTender        = kFALSE;    ///< Declare that tender was executed
Bool_t  kMix           = kFALSE;    ///< Do mixing analysis
Int_t   kRunNumber     = -1;        ///< Declare the run number

///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskCaloTrackCorr(const TString  data          = "",
                                                          const TString  calorimeter   = "EMCAL", 
                                                          const Bool_t   simulation    = kFALSE,
                                                          const Bool_t   eventsel      = kFALSE,
                                                          const Bool_t   exotic        = kTRUE,
                                                          const Bool_t   nonlin        = kFALSE,
                                                          TString        outputfile    = "",
                                                          const Int_t    year          = 2010,
                                                          const TString  col           = "pp", 
                                                          const TString  trigger       = "MB", 
                                                          const TString  clustersArray = "V1",
                                                          const Bool_t   mix           = kTRUE,
                                                          const Bool_t   recaltm       = kTRUE,
                                                          const Bool_t   tm            = kTRUE,
                                                          const Int_t    minCen        = -1,
                                                          const Int_t    maxCen        = -1,
                                                          const Bool_t   qaan          = kFALSE,
                                                          const Bool_t   hadronan      = kFALSE,
                                                          const Bool_t   calibE        = kTRUE,
                                                          const Bool_t   badmap        = kTRUE,
                                                          const Bool_t   calibT        = kTRUE,
                                                          const Bool_t   tender        = kFALSE,
                                                          const Bool_t   outputAOD     = kFALSE, 
                                                          const Bool_t   printSettings = kFALSE,
                                                          const Double_t scaleFactor   = -1, 
                                                          const Int_t    runNumber     = -1
                                                          )
{
  // Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
  
  kPrint         = printSettings;
  kSimulation    = simulation;
  kYears         = year;
  kCollisions    = col;
  kExotic        = exotic;
  kNonLinearity  = nonlin;
  kTrig          = trigger;
  kClusterArray  = clustersArray;
  kData          = data;
  kCalorimeter   = calorimeter;
  kOutputAOD     = outputAOD;
  kTM            = tm;
  kRecalTM       = recaltm;
  kMinCen        = minCen;
  kMaxCen        = maxCen;
  kEventSelection= eventsel;
  kQA            = qaan;
  kHadronAN      = hadronan;
  kCalibE        = calibE;
  kCalibT        = calibT;
  kBadMap        = badmap;
  kTender        = tender;
  kMix           = mix;
  kRunNumber     = runNumber;
  
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
  
  kInputDataType = "AOD";
  if(!kData.Contains("delta"))
    kInputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  if(kSimulation) 
  { 
    kUseKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE; 
    if (!kUseKinematics && data=="AOD" && kInputDataType != "ESD") kUseKinematics = kTRUE; //AOD primary should be available ... 
  } 
  
  cout<<"********* ACCESS KINE? "<<kUseKinematics<<endl;
  
  // Name for containers
  
  kName = Form("%s_Trig%s_Cl%s_TM%d",kCalorimeter.Data(), kTrig.Data(),kClusterArray.Data(),kTM);

  if(kCollisions=="PbPb" && kMaxCen>=0) kName+=Form("Cen%d_%d",kMinCen,kMaxCen);
    
  printf("<<<< NAME: %s >>>>>\n",kName.Data());
  
  // #### Configure analysis ####
    
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  
  maker->SetScaleFactor(scaleFactor); // for MC, negative (not scaled) by default
  
  // General frame setting and configuration
  maker->SetReader   (ConfigureReader()   ); 
  maker->SetCaloUtils(ConfigureCaloUtils()); 
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important

  // Isolation settings
  Int_t partInCone = AliIsolationCut::kNeutralAndCharged; // kOnlyCharged;
  Int_t thresType  = AliIsolationCut::kPtThresIC;//  AliIsolationCut::kSumPtFracIC ; 
  Float_t cone = -1;
  Float_t pth  = -1;
  
  // Photon analysis
  
  maker->AddAnalysis(ConfigurePhotonAnalysis(), n++); // Photon cluster selection

    maker->AddAnalysis(ConfigurePi0Analysis(), n++); // Invariant mass of photon clusters

  // Invariant mass analysis Put here to tag selected photons as decay
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kIMCalo), n++); // Pi0 event by event selection, invariant mass and photon tagging from decay    
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Eta", AliAnaPi0EbE::kIMCalo), n++); // Eta event by event selection, invariant mass and photon tagging from decay
  
  // Photon analysis
  maker->AddAnalysis(ConfigureIsolationAnalysis("Photon", partInCone,thresType,cone, pth), n++); // Photon isolation   
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kFALSE), n++); // Gamma hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Photon",kTRUE,partInCone,thresType, cone, pth) , n++); // Isolated gamma hadron correlation
  //maker->AddAnalysis(ConfigureIsolationAnalysis("Photon", partInCone,thresType,kTRUE), n++); // Photon multi isolation, leave it the last   

  
  // Split cluster analysis
  if(kCalorimeter == "EMCAL")
  {
    maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0", AliAnaPi0EbE::kSSCalo), n++); // Pi0 event by event selection, cluster splitting
    maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SS", partInCone,thresType,cone, pth), n++);       // Pi0 isolation, cluster splits
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS" ,kFALSE), n++); // Pi0 hadron correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SS" ,kTRUE,partInCone,thresType, cone, pth) , n++); // Isolated pi0 hadron correlation
    //maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SS",  partInCone,thresType,kTRUE), n++); // Pi0 multi isolation, split cluster  
    maker->AddAnalysis(ConfigureInClusterIMAnalysis(kTRUE , kTRUE ), n++);
  }
  
  // Invariant mass analysis
  maker->AddAnalysis(ConfigurePi0EbEAnalysis("Pi0SideBand", AliAnaPi0EbE::kIMCalo), n++); // Pi0 event by event selection, and photon tagging from decay    
  maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0", partInCone,thresType,cone, pth), n++);         // Pi0 isolation, invariant mass   
  maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0SideBand", partInCone,thresType,cone, pth), n++); // Pi0 isolation, side band   
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0"   ,kFALSE), n++); // Pi0 hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0"   ,kTRUE,partInCone,thresType, cone, pth) , n++); // Isolated pi0 hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SideBand" ,kFALSE), n++); // Pi0 hadron correlation
  maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Pi0SideBand" ,kTRUE,partInCone,thresType, cone, pth) , n++); // Isolated pi0 hadron correlation
  //maker->AddAnalysis(ConfigureIsolationAnalysis("Pi0",    partInCone,thresType,kTRUE), n++); // Pi0 multi isolation, invariant mass, leave it the last   

  if(kHadronAN)
  {
    maker->AddAnalysis(ConfigureChargedAnalysis(), n++);                                // track selection
    maker->AddAnalysis(ConfigureIsolationAnalysis("Hadron",AliIsolationCut::kOnlyCharged,thresType,cone, pth), n++); // track isolation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Hadron",kFALSE), n++);       // track-track correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("Hadron",kTRUE,partInCone,thresType, cone, pth) , n++);       // Isolated track-track correlation
    //maker->AddAnalysis(ConfigureIsolationAnalysis("Hadron",partInCone,thresType,kTRUE), n++);// Hadron multi isolation  
  }
  
  // Analysis with ghost triggers, only for Min Bias like events
  if( kTrig.Contains("INT") || kTrig.Contains("Central") || kTrig.Contains("MB")  )
  {
    maker->AddAnalysis(ConfigureRandomTriggerAnalysis(), n++); 
    maker->AddAnalysis(ConfigureIsolationAnalysis(Form("RandomTrigger%s",kCalorimeter.Data()), partInCone,thresType,cone, pth), n++); // Ghost trigger isolation   
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis(Form("RandomTrigger%s",kCalorimeter.Data()),kFALSE), n++); // Ghost trigger hadron correlation
    maker->AddAnalysis(ConfigureHadronCorrelationAnalysis(Form("RandomTrigger%s",kCalorimeter.Data()),kTRUE,partInCone,thresType, cone, pth) , n++); // Isolated ghost hadron correlation
    //maker->AddAnalysis(ConfigureIsolationAnalysis(Form("RandomTrigger%s",kCalorimeter.Data()), partInCone,thresType,kTRUE), n++); // Ghost multi isolation   
    
    if(kHadronAN)
    {
      maker->AddAnalysis(ConfigureRandomTriggerAnalysis("CTS"), n++);                                // track selection
      maker->AddAnalysis(ConfigureIsolationAnalysis("RandomTriggerCTS",AliIsolationCut::kOnlyCharged,thresType,cone, pth), n++); // track isolation
      maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("RandomTriggerCTS",kFALSE), n++);       // track-track correlation
      maker->AddAnalysis(ConfigureHadronCorrelationAnalysis("RandomTriggerCTS",kTRUE,partInCone,thresType, cone, pth) , n++);       // Isolated track-track correlation
      //maker->AddAnalysis(ConfigureIsolationAnalysis("RandomTriggerCTS",AliIsolationCut::kOnlyCharged,thresType,kTRUE), n++); // Ghost multi isolation   
    }
  }
  
  if(kQA)  maker->AddAnalysis(ConfigureQAAnalysis(),n++);
  
  maker->SetAnaDebug(kDebug)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(kData.Contains("delta")) maker->SwitchOffAODsMaker() ;
  else                        maker->SwitchOnAODsMaker()  ;
  
  if(kPrint) maker->Print("");
  
  if(kSimulation) maker->SwitchOffDataControlHistograms();
  
  if(simulation)
  {
    // Calculate the cross section weights, apply them to all histograms 
    // and fill xsec and trial histo. Sumw2 must be activated.
    //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionCalculation(); 
    //maker->SwitchOnSumw2Histograms();
    
    // For recent productions where the cross sections and trials are not stored in separate file
    //maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionFromEventHeader() ;
    
    // Just fill cross section and trials histograms.
    maker->GetReader()->GetWeightUtils()->SwitchOnMCCrossSectionHistoFill(); 
    
    // Add control histogram with pT hard to control aplication of weights 
    maker->SwitchOnPtHardHistogram();
  }

  
  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, kCalorimeter.Data());
 
  // Create task
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("CaloTrackCorr%s",kName.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(kDebug);
  task->SetBranches("ESD:AliESDRun.,AliESDHeader"); 
  task->SetAnalysisMaker(maker);
  mgr->AddTask(task);
  
  //Create containers
  
  if(outputfile.Length()==0) outputfile = AliAnalysisManager::GetCommonFileName(); 
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(kName, TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s",outputfile.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",kName.Data()), TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             "AnalysisParameters.root");
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  if(!kData.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
  
  if(!kMix)
  {    
    UInt_t mask =  SetTriggerMaskFromName();
    task->SelectCollisionCandidates(mask);
  } 
  
  return task;
}

///
/// Configure the class handling the events and cluster/tracks filtering.
///
AliCaloTrackReader * ConfigureReader()
{
  AliCaloTrackReader * reader = 0;
  if     (kInputDataType == "ESD"&& kData=="MC" ) 
    reader = new AliCaloTrackMCReader();
  else if(kInputDataType=="AOD" || kData.Contains("AOD"))   
    reader = new AliCaloTrackAODReader();
  else if(kInputDataType=="ESD")            
    reader = new AliCaloTrackESDReader();
  else 
    printf("AliCaloTrackReader::ConfigureReader() - Data combination not known kData=%s, kInputData=%s\n",kData.Data(),kInputDataType.Data());
  
  reader->SetDebug(kDebug);//10 for lots of messages
  
  /*
   if(kSimulation)
   {
     // Event rejection cuts for jet-jet simulations
     reader->SetPtHardAndJetPtComparison(kTRUE);
     reader->SetPtHardAndJetPtFactor(4);
   
     reader->SetPtHardAndClusterPtComparison(kTRUE);
     reader->SetPtHardAndClusterPtFactor(1.5);
   }
   */
  
  //Delta AOD?
  //reader->SetDeltaAODFileName("");
  if(kOutputAOD) reader->SwitchOnWriteDeltaAOD()  ;
  
  // MC settings
  if(kUseKinematics){
    if(kInputDataType == "ESD"){
      reader->SwitchOnStack();          
      reader->SwitchOffAODMCParticles(); 
    }
    else if(kInputDataType == "AOD"){
      reader->SwitchOffStack();          
      reader->SwitchOnAODMCParticles(); 
    }
  }  
  
  //------------------------
  // Detector input filling
  //------------------------
  
  //Min cluster/track E
  reader->SetEMCALEMin(0.3); 
  reader->SetEMCALEMax(1000); 
  reader->SetPHOSEMin(0.3);
  reader->SetPHOSEMax(1000);
  reader->SetCTSPtMin(0.2);
  reader->SetCTSPtMax(1000);

  // Time cuts
  if(kSimulation) 
  {
    reader->SwitchOffUseTrackTimeCut();
    reader->SwitchOffUseParametrizedTimeCut();
    reader->SwitchOffUseEMCALTimeCut();
    reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
  }
  else
  {
    if(kCalibT)
    { 
      printf("Set time cut parameters for run %d\n",kRunNumber);
      reader->SwitchOnUseEMCALTimeCut();
      reader->SwitchOnUseParametrizedTimeCut();
      
      //Absolute window
      reader->SetEMCALTimeCut(-25,20);
      
      //Parametrization
      if     (kRunNumber >= 151636 && kRunNumber <= 155384 )
      {
        printf("Set time parameters for LHC11c");
        reader->SetEMCALParametrizedMinTimeCut(0,-5  ); reader->SetEMCALParametrizedMinTimeCut(1,-1 ); reader->SetEMCALParametrizedMinTimeCut(2, 1.87); reader->SetEMCALParametrizedMinTimeCut(3, 0.4);   
        reader->SetEMCALParametrizedMaxTimeCut(0, 3.5); reader->SetEMCALParametrizedMaxTimeCut(1, 50); reader->SetEMCALParametrizedMaxTimeCut(2, 0.15); reader->SetEMCALParametrizedMaxTimeCut(3, 1.6);   
      }
      else if(kRunNumber >= 156447 && kRunNumber <= 159635 )
      {
        printf("Set time parameters for LHC11d");
        reader->SetEMCALParametrizedMinTimeCut(0,-5);  reader->SetEMCALParametrizedMinTimeCut(1,-1 );  reader->SetEMCALParametrizedMinTimeCut(2, 3.5 ); reader->SetEMCALParametrizedMinTimeCut(3, 1.  );   
        reader->SetEMCALParametrizedMaxTimeCut(0, 5);  reader->SetEMCALParametrizedMaxTimeCut(1, 50);  reader->SetEMCALParametrizedMaxTimeCut(2, 0.45); reader->SetEMCALParametrizedMaxTimeCut(3, 1.25);   
      }
      else
      {
        reader->SwitchOffUseParametrizedTimeCut();
      }
    }
    else
    {
      reader->SwitchOffUseParametrizedTimeCut();
      reader->SwitchOffUseEMCALTimeCut();
      reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
    }
  }
  
  reader->SwitchOnFiducialCut();
  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;

  // Tracks
  reader->SwitchOnCTS();
  reader->SwitchOffRejectNoTrackEvents();

  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();
  
  reader->SwitchOffUseTrackTimeCut();
  reader->SetTrackTimeCut(0,50);
  
  reader->SwitchOffUseTrackDCACut();
  //reader->SetTrackDCACut(0,0.0105);
  //reader->SetTrackDCACut(1,0.035);
  //reader->SetTrackDCACut(2,1.1);
  
  if(kInputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
    //AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
    //reader->SetTrackCuts(esdTrackCuts);
    //reader->SwitchOnConstrainTrackToVertex();
    
    if(kYears>2010)
    {
      //Hybrids 2011
      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
      reader->SetTrackCuts(esdTrackCuts);
      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
      reader->SetTrackComplementaryCuts(esdTrackCuts2);
    }
    else
    {
      //Hybrids 2010
      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001006);
      reader->SetTrackCuts(esdTrackCuts);
      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10041006);
      reader->SetTrackComplementaryCuts(esdTrackCuts2);
    }
  }
  else if(kInputDataType=="AOD")
  {
    //reader->SetTrackFilterMask(128);           // Filter bit, not mask, use if off hybrid
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SetTrackStatus(AliVTrack::kITSrefit);
    //reader->SwitchOnTrackHitSPDSelection();    // Check that the track has at least a hit on the SPD, not much sense to use for hybrid or TPC only tracks
  }
  
  // Calorimeter
  
  reader->SetEMCALClusterListName(kClusterArray);
  if(kClusterArray == "" && !kTender) 
  {
    printf("**************** Standard EMCAL clusters branch analysis **************** \n");
    reader->SwitchOnClusterRecalculation();
    // Check in ConfigureCaloUtils that the recalibration and bad map are ON 
  }
  else 
  {
    printf("**************** Input for analysis is Clusterizer %s **************** \n", kClusterArray.Data());
    reader->SwitchOffClusterRecalculation();
  }  
  
  if(!kNonLinearity) reader->SwitchOffClusterELinearityCorrection();
  
  //if(kCalorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();  
    reader->SwitchOnEMCAL();
  //}
  //if(kCalorimeter == "PHOS") { 
    reader->SwitchOnPHOSCells();  
    reader->SwitchOnPHOS();
  //}
  
  // for case data="deltaAOD", no need to fill the EMCAL/PHOS cluster lists
  if(kData.Contains("delta"))
  {
    reader->SwitchOffEMCAL();
    reader->SwitchOffPHOS();
    reader->SwitchOffEMCALCells(); 
    reader->SwitchOffPHOSCells(); 
  }
  
  //-----------------
  // Event selection
  //-----------------
  
  //reader->RejectFastClusterEvents()  ;

  // Event triggered selection settings
  reader->SwitchOnTriggerPatchMatching();
  reader->SwitchOnBadTriggerEventsRemoval(); // only if SwitchOnTriggerPatchMatching();
  reader->SwitchOnUnMatchedTriggerEventsRemoval(); // only if SwitchOnBadTriggerEventsRemoval();
  //reader->SwitchOffTriggerClusterTimeRecal() ;

  reader->SetTriggerPatchTimeWindow(8,9); // L0
  if     (kRunNumber < 146861) reader->SetEventTriggerL0Threshold(3.);
  else if(kRunNumber < 154000) reader->SetEventTriggerL0Threshold(4.);
  else if(kRunNumber < 165000) reader->SetEventTriggerL0Threshold(5.5);

  //redefine for other periods, triggers

  //if(!kUseKinematics) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
  
  // For mixing with AliAnaParticleHadronCorrelation switch it off
  if(kMix)
  {
    reader->SwitchOffEventTriggerAtSE();
    UInt_t mask =  SetTriggerMaskFromName();
    reader->SetEventTriggerMask(mask); // Only for mixing and SwitchOffEventTriggerAtSE();
    //reader->SetMixEventTriggerMask(AliVEvent::kMB); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
    reader->SetMixEventTriggerMask(AliVEvent::kAnyINT); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
    
    printf("---Trigger selection done in AliCaloTrackReader!!!\n");
  }
  else
    reader->SwitchOnEventTriggerAtSE();  
  
  reader->SetZvertexCut(10.);               // Open cut
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex

  if(kEventSelection)
  {
    reader->SwitchOnEventPileUpRejection(); // remove pileup by default
    reader->SwitchOnV0ANDSelection() ;      // and besides v0 AND
  }
  else 
  {
    reader->SwitchOffPileUpEventRejection();// remove pileup by default
    reader->SwitchOffV0ANDSelection() ;     // and besides v0 AND
  }
    
  if(kCollisions=="PbPb") 
  {
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
    reader->SetCentralityBin(kMinCen,kMaxCen); // Accept all events, if not select range
    
    // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
    reader->SetEventPlaneMethod("V0");
  }
  
  if(kPrint) reader->Print("");
  
  return reader;
}

///
/// Configure the class handling the calorimeter clusters specific methods
///
AliCalorimeterUtils* ConfigureCaloUtils()
{
  
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(kDebug);
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  // Search of local maxima in cluster
  if(kCollisions=="pp")
  {
    cu->SetLocalMaximaCutE(0.1);
    cu->SetLocalMaximaCutEDiff(0.03);
    if(kName.Contains("150"))
    {
      printf("Reclusterize with 150 threshold, set PbPb settings\n");
      cu->SetLocalMaximaCutE(0.2);
      cu->SetLocalMaximaCutEDiff(0.03);
    }
  }
  else 
  {
    cu->SetLocalMaximaCutE(0.2);
    cu->SetLocalMaximaCutEDiff(0.03);
  }
  
  cu->SwitchOffClusterPlot();

  if(kRecalTM) cu->SwitchOnRecalculateClusterTrackMatching(); // Done in clusterization
  else         cu->SwitchOffRecalculateClusterTrackMatching();
  
  cu->SwitchOnBadChannelsRemoval() ;
  
  //EMCAL settings

  if(!kSimulation)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  
  if(!kSimulation)
  {
    cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
    if(kClusterArray == "" && !kTender) cu->SwitchOnRunDepCorrection();
  }
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          kSimulation,                             
                          kExotic,
                          kTRUE,//kNonLinearity,
                          kCalibE, 
                          kBadMap,
                          kCalibT);
  //recou->SetExoticCellDiffTimeCut(50.);

  
  if( kNonLinearity ) 
  { 
    printf("ConfigureCaloUtils() - Apply non linearity to EMCAL\n");
    cu->SwitchOnCorrectClusterLinearity();
  }
    
  printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
  printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
  
  
  if(kCalorimeter=="PHOS")
  {
    if (kYears <  2014) cu->SetNumberOfSuperModulesUsed(3);
    else                cu->SetNumberOfSuperModulesUsed(4);
  }
  else
  {
    if      (kYears == 2010) cu->SetNumberOfSuperModulesUsed(4); //EMCAL first year
    else if (kYears <  2014) cu->SetNumberOfSuperModulesUsed(10);
    else                     cu->SetNumberOfSuperModulesUsed(20);
  }
  
  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  if(kPrint) cu->Print("");
  
  return cu;
}

///
/// Configure the task doing the first photon cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection ...
///
AliAnaPhoton* ConfigurePhotonAnalysis()
{
  AliAnaPhoton *ana = new AliAnaPhoton();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // cluster selection cuts
  
  ana->SwitchOffFiducialCut();

  ana->SetCalorimeter(kCalorimeter);
  
  if(kCalorimeter == "PHOS")
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
    ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off 
    // restrict to less than 100 ns when time calibration is on 
    ana->SetMinDistanceToBadChannel(2, 4, 6); 
    
    // NLM cut, used in all, exclude clusters with more than 2 maxima
    // Not needed if M02 cut is already strong or clusterizer V2
    ana->SetNLMCut(1, 2) ;
  }
  
  if(kTM)
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
  caloPID->SetEMCALLambda0CutMax(0.27);
  caloPID->SetEMCALLambda0CutMin(0.10);
  
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
    
  //PHOS
  caloPID->SetPHOSDispersionCut(2.5);
  caloPID->SetPHOSRCut(2.);
  if(kInputData=="AOD") caloPID->SetPHOSRCut(2000.); // Open cut since dX, dZ not stored
      
  ana->SwitchOffFillShowerShapeHistograms();  // Filled before photon shower shape selection
  if(!kSimulation)ana->SwitchOnFillPileUpHistograms();
  //if(!kSimulation) ana->SwitchOnFillEMCALBCHistograms();

  // Input / output delta AOD settings
  
  if(!kData.Contains("delta")) 
  {
    ana->SetOutputAODName(Form("Photon%s",kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    //ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
  }
  else ana->SetInputAODName(Form("Photon%s",kName.Data()));
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaPhoton_TM%d_",kTM));
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms(20);
  ana->FillNPrimaryHistograms(20);
  
  ana->SwitchOnRealCaloAcceptance(); // primary particle acceptance histograms
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the checks on clusters triggered events
/// For filling all histograms meaninfully, in the reader, time cut must be off
/// and bad triggered events not rejected, and of course analyze triggered events.
///
AliAnaEMCALTriggerClusters* ConfigureEMCALTriggerClusterAnalysis()
{
  AliAnaEMCALTriggerClusters *ana = new AliAnaEMCALTriggerClusters();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // cluster selection cuts
  
  ana->SwitchOffFiducialCut();
  ana->SetNCellCut(1);// At least 2 cells
  ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
  ana->SetMaxEnergy(1000);
  ana->SetM02(1, 2) ;
  ana->SwitchOnTrackMatchRejection() ;
  
  ana->AddToHistogramsName("EMCTriggerClusters_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the checks on clusters pile-up events
/// For filling all histograms meaninfully, in the reader, time cut must be off
/// and bad triggered events in different BC not rejected.
///
AliAnaClusterPileUp* ConfigureClusterPileUpAnalysis()
{
  AliAnaClusterPileUp *ana = new AliAnaClusterPileUp();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // cluster selection cuts
  
  ana->SwitchOffFiducialCut();
  ana->SetNCellCut(1);// At least 2 cells
  ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
  ana->SetMaxEnergy(1000);
  
  ana->AddToHistogramsName("ClusterPileUp_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing electron (or charged hadron) cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection, dEdx, E/p ...
///
AliAnaElectron* ConfigureElectronAnalysis()
{
  AliAnaElectron *ana = new AliAnaElectron();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  ana->FillAODWithElectrons();
  //ana->FillAODWithHadrons();
  //ana->FillAODWithAny();
  
  if(kCalorimeter == "PHOS")
  {
    ana->SetNCellCut(2);// At least 2 cells
    ana->SetMinPt(0.3);
    ana->SetMinDistanceToBadChannel(2, 4, 5);
  }
  else 
  {//EMCAL
    ana->SetNCellCut(1);// At least 2 cells
    ana->SetMinPt(0.5); // no effect minium EMCAL cut.
    ana->SetMaxPt(100); 
    //ana->SetTimeCut(400,900);// Time window of [400-900] ns
    ana->SetMinDistanceToBadChannel(2, 4, 6);
  }
  
  //Electron selection cuts with tracks
  ana->SetEOverP(0.85, 1.2);

  // TO DO, find a more suitable way to set this
  if     (kRunNumber < 146861) ana->SetdEdxCut(72, 90);
  else if(kRunNumber < 154000) ana->SetdEdxCut(54, 70);
  else                         ana->SetdEdxCut(74, 90);
  
  if(kSimulation)  ana->SetdEdxCut(80, 100);
  
  ana->SetCalorimeter(kCalorimeter);
  
  ana->SwitchOnCaloPID();
  
  AliCaloPID* caloPID = ana->GetCaloPID();
  
  caloPID->SetEMCALLambda0CutMax(0.27);
  caloPID->SetEMCALLambda0CutMin(0.10);

  ana->SwitchOffFillShowerShapeHistograms();  
  ana->SwitchOffFillWeightHistograms()  ;
  ana->SwitchOffFiducialCut();
  
  if(!kData.Contains("delta")) 
  {
    ana->SetOutputAODName(Form("Electron%s",kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else ana->SetInputAODName(Form("Electron%s",kName.Data()));
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaElectron_TM%d_",kTM));
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana ;
}

///
/// Configure the task doing random trigger generation
///
AliAnaRandomTrigger* ConfigureRandomTriggerAnalysis(TString detector = "")
{
  AliAnaRandomTrigger *ana = new AliAnaRandomTrigger();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  if(detector=="") detector = kCalorimeter;
  ana->SetTriggerDetector(detector);

  // selection cuts
  ana->SetMinPt(4.); 
  ana->SetMaxPt(51.);   
  
  if     (detector=="EMCAL")
  {
    ana->SetEtaCut(-0.71,0.71);
    ana->SetPhiCut(100*TMath::DegToRad(), 160*TMath::DegToRad());
  }
  else if(detector=="PHOS")
  {
    ana->SetEtaCut(-0.13,0.13);
    ana->SetPhiCut(260*TMath::DegToRad(), 320*TMath::DegToRad());
  }
  else if(detector=="CTS")
  {
    ana->SetEtaCut(-0.9,0.9);
    ana->SetPhiCut(0, TMath::TwoPi());
  }
  
  // AOD branch
  if(!kData.Contains("delta")) 
  {
    ana->SetOutputAODName(Form("RandomTrigger%s%s",detector.Data(),kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else 
    ana->SetInputAODName(Form("RandomTrigger%s%s",detector.Data(),kName.Data()));
  
  printf("Set RandomTrigger%s%s\n",detector.Data(),kName.Data());
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaRandomTrigger%s_",detector.Data()));
  
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  if(detector=="CTS")
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing cluster identification as merged pi0/eta decays
///
AliAnaInsideClusterInvariantMass* ConfigureInClusterIMAnalysis(Bool_t useSS = kTRUE, Bool_t useAsy = kFALSE)
{
  AliAnaInsideClusterInvariantMass *ana = new AliAnaInsideClusterInvariantMass();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // selection cuts
  
  ana->SetMinEnergy(6);
  ana->SetMaxEnergy(200.);
  ana->SetMinNCells(6); // check same as in calopid
  
  ana->SetCalorimeter(kCalorimeter);
  ana->SwitchOnSplitClusterDistToBad();
  
  ana->SwitchOffFillSSWeightHistograms() ;
  ana->SetNWeightForShowerShape(0);
  //ana->SetWeightForShowerShape(0, 4.6);
  
  ana->SwitchOnFillNCellHistograms();
  ana->SwitchOffFillEbinHistograms();
  if(!useSS && !useAsy) ana->SwitchOnFillEbinHistograms();
  
  if(!kTM)
  {
    ana->SwitchOnFillTMHistograms();
    ana->SwitchOnFillTMResidualHistograms();
  }
  else
  {
    ana->SwitchOffFillTMHistograms();
    ana->SwitchOffFillTMResidualHistograms();
  }
  
  //printf("Set correction slope for SS weight \n");
  //ana->SetWCorrectionParameter(0.07);
  //ana->SetNECellCutForShowerShape(0);
  //ana->SetECellCutForShowerShape(0, 0.07);
  //ana->SetECellCutForShowerShape(1, 0.1);
  //ana->SetECellCutForShowerShape(2, 0.2);
  
  if(kSimulation)
  {
    ana->SwitchOnFillMCPrimaryHistograms() ;
    ana->SwitchOffFillMCOverlapHistograms() ; // Off when possible
    if(!useSS && !useAsy) ana->SwitchOnFillMCOverlapHistograms() ;
  }
  
  AliCaloPID* caloPID = ana->GetCaloPID();
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
  caloPID->SetClusterSplittingM02Cut(0,100000); // use parametrized cut, not fixed
  
  caloPID->SetPi0MassRange(0.11, 0.18);
  caloPID->SetEtaMassRange(0.40, 0.60);
  caloPID->SetPhotonMassRange(0.00, 0.08);
  
  caloPID->SetSplitWidthSigma(3.); // cut at 3 sigma of the mean pi0 peak.
  
  caloPID->SetClusterSplittingMinNCells(6);
  
  if(kCollisions=="PbPb" || kName.Contains("150"))
  {
    caloPID->SetClusterSplittingMinNCells(4);
    ana->SetMinNCells(4);
    caloPID->SetPi0MassShiftHighECell(0.005);
    if(kCollisions=="PbPb") ana->SwitchOnFillHighMultHistograms();
  }
  
  ana->AddToHistogramsName("AnaInClusterIM_");
  
  if(useAsy)
  {
    caloPID->SwitchOnSplitAsymmetryCut() ;
  }
  else
  {
    printf("InClusterIM: Do not apply Asy cut on merged pi0 in cluster analysis \n");
    caloPID->SwitchOffSplitAsymmetryCut() ;
    ana->AddToHistogramsName("AnaInClusterIM_OpenAsy_");
  }
  
  if(!useSS)
  {
    printf("InClusterIM: Do not apply SS cut on merged pi0 in cluster analysis \n");
    caloPID->SwitchOffSplitShowerShapeCut() ;
    ana->AddToHistogramsName("AnaInClusterIM_OpenSS_");
  }
  else  caloPID->SwitchOnSplitShowerShapeCut() ;
  
  if(!useAsy && !useSS)
  {
    printf("InClusterIM: Do not apply SS and Asy cut on merged pi0 in cluster analysis \n");
    ana->AddToHistogramsName("AnaInClusterIM_OpenSS_OpenAsy_");
  }
  
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing charged track selection
///
AliAnaChargedParticles* ConfigureChargedAnalysis()
{
  
  AliAnaChargedParticles *ana = new AliAnaChargedParticles();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // selection cuts
  
  ana->SetMinPt(0.5);
  ana->SwitchOnFiducialCut();
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ; //more restrictive cut in reader and after in isolation
  
  ana->SwitchOnFillVertexBC0Histograms() ;
  if(!kSimulation) ana->SwitchOnFillPileUpHistograms();
  
  // Input / output delta AOD settings
  
  if(!kData.Contains("delta"))
  {
    ana->SetOutputAODName(Form("Hadron%s",kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    //ana->SetOutputAODClassName("AliAODPWG4Particle"); // use if no correlation done
  }
  else
    ana->SetInputAODName(Form("Hadron%s",kName.Data()));
  printf("Set Hadron%s\n",kName.Data());
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaHadrons_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
  ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the 2 cluster invariant mass analysis
///
AliAnaPi0* ConfigurePi0Analysis()
{
  AliAnaPi0 *ana = new AliAnaPi0();

  ana->SetDebug(kDebug);//10 for lots of messages
    
  // Input delta AOD settings
  ana->SetInputAODName(Form("Photon%s",kName.Data()));
  
  // Calorimeter settings
  ana->SetCalorimeter(kCalorimeter);
  
  // Acceptance plots
  //  ana->SwitchOnFiducialCut(); // Needed to fill acceptance plots with predefined calorimeter acceptances
  //  ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.7, 100, 180) ; 
  //  ana->GetFiducialCut()->DoEMCALFiducialCut(kTRUE);
  
  ana->SwitchOffFiducialCut();
  ana->SwitchOnRealCaloAcceptance();

  // settings for pp collision mixing
  ana->SwitchOnOwnMix(); //Off when mixing done with general mixing frame
  
  // Cuts
  if(kCalorimeter=="EMCAL") ana->SetPairTimeCut(40);
  
  ana->SetNAsymCuts(1); // no asymmetry cut, previous studies showed small effect.
  // In EMCAL assymetry cut prevents combination of assymetric decays which is the main source of pi0 at high E.
    
  if     (kCollisions=="pp"  )
  {
    ana->SetNCentrBin(1);
    ana->SetNZvertBin(10);
    ana->SetNRPBin(1);
    ana->SetNMaxEvMix(100);
  }
  else if(kCollisions=="PbPb")
  {
    ana->SetNCentrBin(5);
    ana->SetNZvertBin(3);
    ana->SetNRPBin(1);
    ana->SetNMaxEvMix(5);
  }
    
  ana->SwitchOffSMCombinations();
  ana->SwitchOffMultipleCutAnalysis();
    
  // Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("AnaPi0_TM%d_",kTM));
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the pi0 even by event selection via the split method.
/// Here the pairs, clusters, are added to an AOD branch to be used by other analysis
/// unlike in ConfigurePi0Analysis.
///
AliAnaPi0EbE* ConfigurePi0EbEAnalysis(TString particle,
                                      Int_t analysis, Bool_t useSS = kTRUE, Bool_t useAsy = kTRUE)
{
  AliAnaPi0EbE *ana = new AliAnaPi0EbE();
  ana->SetDebug(kDebug);//10 for lots of messages
  
  ana->SetAnalysisType(analysis);
  TString opt = "";
  if(analysis == AliAnaPi0EbE::kIMCaloTracks) opt = "Conv";
  if(analysis == AliAnaPi0EbE::kSSCalo)       opt = "SS";
  
  if(analysis == AliAnaPi0EbE::kIMCalo && kCalorimeter=="EMCAL" && !kSimulation) ana->SetPairTimeCut(100);
  if(analysis == AliAnaPi0EbE::kIMCaloTracks) ana->SetInputAODGammaConvName("PhotonsCTS");

  // Common settings for all 3 type of analysis
  
  ana->SwitchOnSelectedClusterHistoFill();

  ana->SetCalorimeter(kCalorimeter);
  
  //Set Histograms name tag, bins and ranges
  ana->AddToHistogramsName(Form("Ana%s%sEbE_TM%d_",particle.Data(),opt.Data(),kTM));
  
  // Specific settings for different type of analysis
  
  ana->SwitchOffFillWeightHistograms();
  if(!kSimulation) ana->SwitchOnFillPileUpHistograms();
  
  if(kTM)
  {
    //printf("--->>>REMOVE MATCHED Pi0\n");
    ana->SwitchOnTrackMatchRejection() ;
    ana->SwitchOffTMHistoFill() ;
  }
  else
  {
    //printf("---->>>ACCEPT MATCHED Pi0\n");
    ana->SwitchOffTrackMatchRejection() ;
    ana->SwitchOnTMHistoFill() ;
  }
  
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  //ana->SwitchOnFillEMCALBCHistograms();
  
  if(kPrint) ana->Print("");
  
  ConfigureMC(ana);

  if(!kInputDataType.Contains("delta"))
  {
    ana->SetOutputAODName(Form("%s%s%s",particle.Data(), opt.Data(), kName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
    
  }
  else
    ana->SetInputAODName(Form("%s%s%s",particle.Data(),opt.Data(),kName.Data()));
  
  if(analysis!=AliAnaPi0EbE::kSSCalo)
  {
    // Input / output delta AOD settings
    
    ana->SetInputAODName(Form("Photon%s",kName.Data()));
    
    AliNeutralMesonSelection *nms = ana->GetNeutralMesonSelection();
    nms->SetParticle(particle);
    
    // Tighten a bit mass cut with respect to default window
    if(particle=="Pi0") nms->SetInvMassCutRange(0.120,0.150);
    if(particle=="Eta") nms->SetInvMassCutRange(0.520,0.580);
    
    //if(!particle.Contains("SideBand")) nms->SwitchOnAngleSelection();
    //else nms->SwitchOnAngleSelection();
    
    nms->SwitchOffAngleSelection();
    if(particle.Contains("Pi0SideBand")) // For pi0, do not consider left band
      nms->SetSideBandCutRanges(-1,0,0.180,0.220);
    
    nms->KeepNeutralMesonSelectionHistos(kTRUE);
    //nms->SetAngleMaxParam(2,0.2);
    nms->SetHistoERangeAndNBins(0, 20, 80) ;
    //nms->SetHistoIMRangeAndNBins(0, 1, 400);
  }
  else
  { // cluster splitting settings
    ana->SetMinEnergy(6);
    ana->SetMaxEnergy(200.);
    
    ana->SetNLMMinEnergy(0, 10);
    ana->SetNLMMinEnergy(1, 6);
    ana->SetNLMMinEnergy(2, 6);
    
    ana->SetMinDistanceToBadChannel(2, 4, 6); // only use the first one
    ana->SwitchOnSplitClusterDistToBad();
    
    ana->SetTimeCut(-1e10,1e10); // Open time cut

    // NLM cut, used in all, exclude clusters with more than 2 maxima
    ana->SetNLMCut(1, 2) ;
    
    AliCaloPID* caloPID = ana->GetCaloPID();
    
    caloPID->SetSplitWidthSigma(3.); // cut at 3 sigma of the mean pi0 peak.
    
    if(!useSS)
    {
      printf("Do not apply SS cut on merged pi0 analysis \n");
      caloPID->SwitchOffSplitShowerShapeCut() ;
      ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_TM%d_",particle.Data(),opt.Data(),kTM));
      ana->SetOutputAODName(Form("%s%s%s_OpenSS",particle.Data(), opt.Data(), kName.Data()));
      caloPID->SetClusterSplittingM02Cut(0.1,10000); 
    }
    else
    {
      caloPID->SetClusterSplittingM02Cut(0.3,5); // Do the selection in the analysis class and not in the PID method to fill SS histograms
      caloPID->SwitchOnSplitShowerShapeCut() ;
    }
    
    if(useAsy) caloPID->SwitchOnSplitAsymmetryCut() ;
    else
    {
      caloPID->SwitchOffSplitAsymmetryCut() ;
      if(!useSS)
      {
        ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenSS_OpenAsy_TM%d_",particle.Data(),opt.Data(),kTM));
        ana->SetOutputAODName(Form("%s%s%s_OpenSS_OpenAsy",particle.Data(), opt.Data(), kName.Data()));
      }
      else
      {
        ana->AddToHistogramsName(Form("Ana%s%sEbE_OpenAsy_TM%d_",particle.Data(),opt.Data(),kTM));
        ana->SetOutputAODName(Form("%s%s%s_OpenAsy",particle.Data(), opt.Data(), kName.Data()));
      }
    }
    
    //For Pi0 only if  SwitchOnSimpleSplitMassCut()
    caloPID->SetPi0MassRange(0.10, 0.18);
    caloPID->SetEtaMassRange(0.40, 0.60);
    caloPID->SetPhotonMassRange(0.00, 0.08);
    
    caloPID->SetClusterSplittingMinNCells(6);
    
    //caloPID->SetSplitEnergyFractionMinimum(0, 0.95);
    //caloPID->SetSplitEnergyFractionMinimum(1, 0.95);
    //caloPID->SetSplitEnergyFractionMinimum(2, 0.8);
    
    if(kCollisions=="PbPb" || kName.Contains("150"))
    {
      caloPID->SetClusterSplittingMinNCells(4);
      caloPID->SetPi0MassShiftHighECell(0.005);
    }
  }
  
  return  ana;
  
}

///
/// Configure the task doing the identification of merged clusters as pi0 or eta
/// same as one of the options of ConfigurePi0EbEAnalysis, but here no AOD with
/// selected particles is created
///
AliAnaInsideClusterInvariantMass* ConfigureInClusterIMAnalysis(Float_t l0min, Float_t l0max)
{
  AliAnaInsideClusterInvariantMass *ana = new AliAnaInsideClusterInvariantMass();
  ana->SetDebug(kDebug); //10 for lots of messages

  // selection cuts
  
  ana->SetMinEnergy(5);
  ana->SetMaxEnergy(200.);
  ana->SetMinNCells(3);
  ana->SetM02Cut(l0min,l0max);
  ana->SetCalorimeter(kCalorimeter);
    
  //ana->AddToHistogramsName(Form("AnaInClusterIM_%1.2f_%1.2f_",l0min,l0max));
  ana->AddToHistogramsName("AnaInClusterIM_");
    
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  AliCaloPID* caloPID = ana->GetCaloPID();
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
  caloPID->SetClusterSplittingM02Cut(0,100); // Do the selection in the analysis class and not in the PID method to fill SS histograms
    
  caloPID->SetPi0MassRange(0.10, 0.18);
  caloPID->SetEtaMassRange(0.40, 0.60);
  caloPID->SetPhotonMassRange(0.00, 0.08);
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the trigger particle isolation
///
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString particle="Photon", 
                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                    Float_t cone = 0.3,
                                                    Float_t pth  = 0.3,
                                                    Bool_t multi      = kFALSE)
{
  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
  //ana->SetDebug(kDebug);
  ana->SetDebug(kDebug);
  
  ana->SwitchOnFiducialCut();
  //Avoid borders of EMCal
  if(kCalorimeter=="EMCAL")
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;

  // Same Eta as EMCal, cut in phi if EMCAL was triggering
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
   // if(kTrig.Contains("EMC"))
   //   ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
   // else
      ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;    
  }
  
  ana->SetMinPt(5);
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s",particle.Data(),kName.Data()));
  ana->SetAODObjArrayName(Form("IC%s_%s",particle.Data(),kName.Data())); 
  
  ana->SetCalorimeter(kCalorimeter);
  
  if(!kTM)  ana->SwitchOnTMHistoFill();
  else      ana->SwitchOffTMHistoFill();
  
  ana->SwitchOffSSHistoFill();
  if(!kSimulation) ana->SwitchOnFillPileUpHistograms();

  //Do settings for main isolation cut class
  AliIsolationCut * ic =  ana->GetIsolationCut();	
  ic->SetDebug(kDebug);
  
  if(cone >0 && pth > 0)
  {
    ic->SetPtThreshold(pth);
    ic->SetConeSize(cone);
  }
  else
  {
    if(kCollisions=="pp") 
    {
      ic->SetPtThreshold(0.5);
      ic->SetConeSize(0.4);
    }
    if(kCollisions=="PbPb")
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
    ana->SetAODObjArrayName(Form("MultiIC%sTM%d",particle.Data(),kTM));
     
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
  
  if(!multi)ana->AddToHistogramsName(Form("AnaIsol%s_TM%d_",particle.Data(),kTM));
  else      ana->AddToHistogramsName(Form("AnaMultiIsol%s_TM%d_",particle.Data(),kTM));

  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  ana->SwitchOnRealCaloAcceptance(); // primary particle acceptance histograms
  ConfigureMC(ana);
  
  if(kPrint) ic ->Print("");
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the task doing the trigger particle hadron correlation
///
AliAnaParticleHadronCorrelation* ConfigureHadronCorrelationAnalysis(TString particle,
                                                                    Bool_t bIsolated,
                                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                                    Float_t cone = 0.3,
                                                                    Float_t pth  = 0.3)
{
  AliAnaParticleHadronCorrelation *ana = new AliAnaParticleHadronCorrelation();
  ana->SetDebug(kDebug);
  
  ana->SwitchOnAbsoluteLeading();  // Select trigger leading particle of all the selected tracks
  ana->SwitchOffNearSideLeading(); // Select trigger leading particle of all the particles at +-90 degrees, default
  
  //ana->SwitchOnLeadHadronSelection();
  //ana->SetLeadHadronPhiCut(TMath::DegToRad()*100., TMath::DegToRad()*260.);
  //ana->SetLeadHadronPtCut(0.5, 100);
  
  ana->SetTriggerPtRange(5,100);
  ana->SetAssociatedPtRange(0.2,100);
  //ana->SetDeltaPhiCutRange( TMath::Pi()/2,3*TMath::Pi()/2 ); //[90 deg, 270 deg]
  ana->SetDeltaPhiCutRange  (TMath::DegToRad()*120.,TMath::DegToRad()*240.);
  ana->SetUeDeltaPhiCutRange(TMath::DegToRad()*60. ,TMath::DegToRad()*120.);
  ana->SwitchOnFillEtaGapHistograms();

  ana->SetNAssocPtBins(9);
  ana->SetAssocPtBinLimit(0, 0.2) ;
  ana->SetAssocPtBinLimit(1, 0.5) ;
  ana->SetAssocPtBinLimit(2, 1)   ;
  ana->SetAssocPtBinLimit(3, 2)   ;
  ana->SetAssocPtBinLimit(4, 3)   ;
  ana->SetAssocPtBinLimit(5, 4)   ;
  ana->SetAssocPtBinLimit(6, 6)   ;
  ana->SetAssocPtBinLimit(7, 10)  ;
  ana->SetAssocPtBinLimit(8, 30)  ;
  ana->SetAssocPtBinLimit(9, 200) ;
  //ana->SwitchOnFillPtImbalancePerPtABinHistograms();

  ana->SelectIsolated(bIsolated); // do correlation with isolated photons
  
  if(bIsolated)
  {
    //Do settings for main isolation cut class
    AliIsolationCut * ic =  ana->GetIsolationCut();	
    ic->SetDebug(kDebug);
    
    if(cone >0 && pth > 0)
    {
      ic->SetPtThreshold(pth);
      ic->SetConeSize(cone);
    }
    else
    {
      if(kCollisions=="pp") 
      {
        ic->SetPtThreshold(0.5);
        ic->SetConeSize(0.4);
      }
      if(kCollisions=="PbPb")
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
    
  }  
  
  // Mixing with own pool
  if(kMix)
  {
    ana->SwitchOnOwnMix();
    ana->SwitchOnFillNeutralInMixedEvent();
  }
  else
    ana->SwitchOffOwnMix();
  
  ana->SetNZvertBin(20);
  ana->SwitchOffCorrelationVzBin() ;
  ana->SwitchOffFillHighMultiplicityHistograms();

  if(kCollisions=="pp")
  {
    ana->SetNMaxEvMix(100);    
    ana->SwitchOnTrackMultBins();
    ana->SetNTrackMultBin(10);  // same as SetNCentrBin(10);
    ana->SetNRPBin(1);
  }
  else 
  {
    ana->SetNMaxEvMix(10);
    ana->SwitchOffTrackMultBins(); // centrality bins
    ana->SetNCentrBin(3); 
    ana->SetNRPBin(3);
    if(kName.Contains("60_90"))
    {
      printf("*** Set mixing for peripheral\n");
      ana->SetNMaxEvMix(50);    
      ana->SetNCentrBin(2); 
    }    
  }  
  
  ana->SwitchOnFiducialCut();
  
  //Avoid borders of EMCal, same as for isolation
  if(kCalorimeter=="EMCAL")
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;
  
  // Same Eta as EMCal, cut in phi if EMCAL was triggering
  if(particle=="Hadron" || particle.Contains("CTS"))
  {
    //if(kTrig.Contains("EMC"))
    //  ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
    //else
      ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;    
  }
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s",particle.Data(),kName.Data()));
  ana->SetAODObjArrayName(Form("%sHadronCorrIso%d_%s",particle.Data(),bIsolated,kName.Data())); 
  
  // Fill extra plots on tagged decay photons
  // If trigger is pi0/eta found with invariant mass, get the decays
  // If trigger is photon, check if it was tagged as decay previously
  if(particle!="Hadron" )
  {
    if(particle.Contains("Pi0") || particle.Contains("Eta"))
    {
      ana->SwitchOffPi0TriggerDecayCorr();
      ana->SwitchOffDecayTriggerDecayCorr();
    }
    else
    {
      ana->SwitchOffPi0TriggerDecayCorr();
      ana->SwitchOnDecayTriggerDecayCorr(); // Make sure pi0 decay tagging runs before this task
    }
  }
  else
  {
    ana->SwitchOffPi0TriggerDecayCorr();
    ana->SwitchOffDecayTriggerDecayCorr(); 
  }
  
  if(particle=="Photon")
  {
    printf("**** SET M02 limits *** \n");
    ana->SetM02Cut(0.1,0.27);
  }
  
  // if triggering on PHOS and EMCAL is on
  //if(kCalorimeter=="PHOS") ana->SwitchOnNeutralCorr();
  ana->SwitchOffNeutralCorr(); // Do only correlation with TPC
  
  ana->SwitchOffHMPIDCorrelation();
  
  ana->SwitchOffFillBradHistograms();
  
  // Underlying event
  ana->SwitchOnSeveralUECalculation();
  ana->SetUeDeltaPhiCutRange(TMath::Pi()/3, 2*TMath::Pi()/3);
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName(Form("Ana%sHadronCorr_Iso%d_TM%d_",particle.Data(),bIsolated,kTM));
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }  
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");  
  
  return ana;
}

///
/// Configure the task doing standard calorimeter QA
///
AliAnaCalorimeterQA* ConfigureQAAnalysis()
{
  AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
  ana->SetDebug(kDebug); //10 for lots of messages
  ana->SetCalorimeter(kCalorimeter);
  
  ana->SetTimeCut(-1e10,1e10); // Open time cut
  
  // Study inter detector correlation (PHOS, EMCAL, Tracks, V0)
  if(kCalorimeter=="PHOS"  && kTrig=="PHOS")
    ana->SwitchOnCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
  if(kCalorimeter=="EMCAL" && kClusterArray=="")
    ana->SwitchOnCorrelation(); // make sure you switch in the reader PHOS and EMCAL cells and clusters if option is ON
  else 
    ana->SwitchOffCorrelation();
  
  // Study exotic clusters PHOS and EMCAL
  if(kClusterArray=="") ana->SwitchOnStudyBadClusters() ; 
  else                  ana->SwitchOffStudyBadClusters() ;
  
  
  ana->SwitchOffFiducialCut();
  ana->SwitchOffFillAllTH3Histogram();
  ana->SwitchOffFillAllPositionHistogram();
  ana->SwitchOffFillAllPositionHistogram2();
  if(!kExotic)ana->SwitchOnStudyBadClusters();
  else        ana->SwitchOffStudyBadClusters();
  ana->SwitchOffStudyClustersAsymmetry();
  ana->SwitchOffStudyWeight();
  ana->SwitchOnFillAllTrackMatchingHistogram();
  ana->SwitchOnFillAllCellTimeHisto() ;
  
  ana->AddToHistogramsName("QA_"); //Begining of histograms name
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  ConfigureMC(ana);
  
  if(kPrint) ana->Print("");	
  
  return ana;
}

///
/// Configure the task doing analysis at the generator level
/// of high pT photon or pi0 and correlations.
///
AliAnaGeneratorKine* ConfigureGenKineAnalysis()
{
  AliAnaGeneratorKine *ana = new AliAnaGeneratorKine();
  ana->SetDebug(kDebug); //10 for lots of messages
  
  // Trigger detector, acceptance and pT cut
  ana->SetTriggerDetector("EMCAL");
  ana->SetMinPt(10); // Trigger photon, pi0 minimum pT
  ana->GetFiducialCutForTrigger()->SetSimpleEMCALFiducialCut(0.6, 85, 175);
  
  // Particles associated to trigger or isolation cone acceptance and pT cut
  ana->SetCalorimeter("EMCAL");
  ana->SetMinChargedPt(0.2);
  ana->SetMinNeutralPt(0.3);
  ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.65, 81, 179);
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360);
  
  // Isolation paramters
  AliIsolationCut * ic =  ana->GetIsolationCut();
  ic->SetDebug(kDebug);
  ic->SetPtThreshold(0.5);
  ic->SetConeSize(0.5);
  ic->SetSumPtThreshold(1.0) ;
  ic->SetICMethod(AliIsolationCut::kPtThresIC); // kSumPtIC
  
  ana->AddToHistogramsName("AnaGenKine_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges()); // see method below
  
  if(kPrint) ana->Print("");
  
  return ana;
}

///
/// Configure the selection of MC events
///
void ConfigureMC(AliAnaCaloTrackCorrBaseClass* ana)
{
  if(kSimulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else            ana->SwitchOffDataMC() ;

  //Set here generator name, default pythia
  //ana->GetMCAnalysisUtils()->SetMCGenerator("");
}  

///
/// Set common histograms binning and ranges
///
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges)
{
  // Set common bins for all analysis and MC histograms filling
    
  histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
  
  if(kCalorimeter=="EMCAL")
  {
    if ( kYears == 2010 )
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
      histoRanges->SetHistoXRangeAndNBins(-230,90,120); // QA
      histoRanges->SetHistoYRangeAndNBins(370,450,40);  // QA
    }
    else if ( kYears < 2014 )
    {           
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 182*TMath::DegToRad(), 108) ;
      histoRanges->SetHistoXRangeAndNBins(-460,90,200); // QA
      histoRanges->SetHistoYRangeAndNBins(100,450,100); // QA
    }
    else // Run2
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 329*TMath::DegToRad(), 250) ;
      histoRanges->SetHistoXRangeAndNBins(-460,460,230); // QA
      histoRanges->SetHistoYRangeAndNBins(-450,450,225); // QA
    }

    histoRanges->SetHistoEtaRangeAndNBins(-0.72, 0.72, 144) ;
  }
  else
  {
    histoRanges->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 320*TMath::DegToRad(), 70) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
  }
  
  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 4.9, 500);
  
  // Invariant mass histoRangeslysis
  histoRanges->SetHistoMassRangeAndNBins(0., 1., 200) ;
  histoRanges->SetHistoAsymmetryRangeAndNBins(0., 1. , 100) ;
  
  // check if time calibration is on
  histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);
  
  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistodRRangeAndNBins(0.,0.15,150);//QA

  // QA, electron, charged
  histoRanges->SetHistoPOverERangeAndNBins(0,2.,200);
  histoRanges->SetHistodEdxRangeAndNBins(0.,200.,200);
  
  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  histoRanges->SetHistoZRangeAndNBins(-350,350,175);
  histoRanges->SetHistoRRangeAndNBins(430,460,30);
  histoRanges->SetHistoV0SignalRangeAndNBins(0,5000,500);
  histoRanges->SetHistoV0MultiplicityRangeAndNBins(0,5000,500);
  
  // QA, correlation
  if(kCollisions=="PbPb")
  {
    histoRanges->SetHistoNClusterCellRangeAndNBins(0,100,100);
    histoRanges->SetHistoNClustersRangeAndNBins(0,500,50);
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,2000,200);
  }
  else
  {
    histoRanges->SetHistoNClusterCellRangeAndNBins(0,50,50);
    histoRanges->SetHistoNClustersRangeAndNBins(0,50,50);
    histoRanges->SetHistoTrackMultiplicityRangeAndNBins(0,200,200);
  }
  
  // xE, zT
  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,200);
  histoRanges->SetHistoHBPRangeAndNBins  (0.,10.,200);
  
  // Isolation
  histoRanges->SetHistoPtInConeRangeAndNBins(0, 50 , 250);
  histoRanges->SetHistoPtSumRangeAndNBins   (0, 100, 250);
  
}

///
/// Set the trigger requested for the analysis,
/// depending on a string given
///
UInt_t SetTriggerMaskFromName()
{
  if(kTrig=="EMC7")
  {
    printf("CaloTrackCorr trigger EMC7\n");
    return AliVEvent::kEMC7;
  }
  else if (kTrig=="INT7")
  {
    printf("CaloTrackCorr trigger INT7\n");
    return AliVEvent::kINT7;
  }
  else if(kTrig=="EMC1")
  {
    printf("CaloTrackCorr trigger EMC1\n");
    return AliVEvent::kEMC1;
  }
  else if(kTrig=="MB")
  {
    printf("CaloTrackCorr trigger MB\n");
    return AliVEvent::kMB;
  }  
  else if(kTrig=="PHOS")
  {
    printf("CaloTrackCorr trigger PHOS\n");
    return AliVEvent::kPHI7;
  }  
  else if(kTrig=="PHOSPb")
  {
    printf("CaloTrackCorr trigger PHOSPb\n");
    return AliVEvent::kPHOSPb;
  }
  else if(kTrig=="AnyINT")
  {
    printf("CaloTrackCorr trigger AnyINT\n");
    return AliVEvent::kAnyINT;
  }  
  else if(kTrig=="INT")
  {
    printf("CaloTrackCorr trigger AnyINT\n");
    return AliVEvent::kAny;
  }
  else if(kTrig=="EMCEGA")
  {
    printf("CaloTrackCorr trigger EMC Gamma\n");
    return AliVEvent::kEMCEGA;
  } 
  else if(kTrig=="EMCEJE")
  {
    printf("CaloTrackCorr trigger EMC Jet\n");
    return AliVEvent::kEMCEJE;
  }
  else if(kTrig=="Central")
  {
    printf("CaloTrackCorr trigger Central\n");
    return (AliVEvent::kCentral  | AliVEvent::kMB);
  }
  else if(kTrig=="CentralEGA")
  {
    printf("CaloTrackCorr trigger Central+EMCEGA\n");
    return (AliVEvent::kCentral | AliVEvent::kEMCEGA);
  }
  else if(kTrig=="SemiCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral\n");
    return (AliVEvent::kSemiCentral | AliVEvent::kMB);
  }
  else if(kTrig=="SemiOrCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral Or Central\n");
    return (AliVEvent::kSemiCentral | AliVEvent::kCentral  | AliVEvent::kMB);
  }
}

