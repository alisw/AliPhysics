/// \file AddTaskMergedPi0Selection.C
/// \ingroup CaloTrackCorrMacros
/// \brief Configuration analysis identifying single EMCal clusters as merged decays
///
/// Configuration macro for analysis of merged pi0 spectra analysis
/// It does the tagging of EMCAL clusters as pi0 via shower shape and splitting criteria in AliAnaPi0EbE
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)

/// Global name to be composed of the settings, used to set the AOD branch name
TString kAnaMergedPi0 = "";

///
/// Main method calling all the configuration 
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
///   * A bool identifying the data as simulation
///   * The year the data was taken, used to configure some histograms
///   * A string with the colliding system
///   * A string with the trigger class, abbreviated, defined in method below SetTriggerMaskFromName()
///   * An int to reject EMCal triggered events with bad trigger: 0 no rejection, 1 old runs L1 bit, 2 newer runs L1 bit
///   * A string with the array of clusters not being the default (default is empty string)
///   * A bool indicating if the tender was running before this analysis
///   * A bool indicating that the non linearity correction is needed
///   * An int to select the minimum centrality, -1 means no selection
///   * An int to select the maximum centrality, -1 means no selection
///   * A string to change the name of the histograms output file, default is AnalysisResults.root
///   * A bool to enable the print of the settings per task
///   * An int to define the debug level of all the tasks
///
/// The options that can be passed to the macro are:
/// \param simulation : A bool identifying the data as simulation
/// \param year: The year the data was taken, used to configure some histograms
/// \param col: A string with the colliding system
/// \param trigger : A string with the trigger class, abbreviated, defined in method belowSetTriggerMaskFromName()
/// \param rejectEMCTrig : An int to reject EMCal triggered events with bad trigger: 0 no rejection, 1 old runs L1 bit, 2 newer runs L1 bit
/// \param clustersArray : A string with the array of clusters not being the default (default is empty string)
/// \param tender : A bool indicating if the tender was running before this analysis
/// \param nonLinOn : A bool to set the use of the non linearity correction
/// \param minCen : An int to select the minimum centrality, -1 means no selection
/// \param maxCen : An int to select the maximum centrality, -1 means no selection
/// \param outputfile : A string to change the name of the histograms output file, default is AnalysisResults.root
/// \param printSettings : A bool to enable the print of the settings per task
/// \param debug : An int to define the debug level of all the tasks
///
AliAnalysisTaskCaloTrackCorrelation * AddTaskMergedPi0Selection
(
 Bool_t   simulation    = kFALSE,
 Int_t    year          = 2011,
 TString  col           = "pp", 
 TString  trigger       = "EMC7", 
 Int_t    rejectEMCTrig = 0, 
 TString  clustersArray = "",
 Bool_t   tender        = kFALSE,
 Bool_t   nonLinOn      = kFALSE,
 Int_t    minCen        = -1,
 Int_t    maxCen        = -1,
 TString  outputfile    = "",
 Bool_t   printSettings = kFALSE, 
 Int_t    debug         = 0  // Debug level
 )
{  
  // Get the pointer to the existing analysis manager via the static access method.
  
  printf("Passed settings:\n simulation <%d>, year <%d>,\n col <%s>, trigger <%s>, reject EMC <%d>, clustersArray <%s>, tender <%d>, non linearity <%d>\n, minCen <%d>, maxCen <%d>,\n, outputfile <%s>, printSettings <%d>, debug <%d>\n", 
         simulation,year,col.Data(),trigger.Data(), rejectEMCTrig, clustersArray.Data(),tender, nonLinOn,
         minCen,maxCen,outputfile.Data(),printSettings,debug);
  
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
  
  // Name for containers
  
  kAnaMergedPi0 = Form("MergedPi0_Trig%s_Col_%s_Year%d_Cl%s_Ten%d",
                             trigger.Data(),col.Data(),year,clustersArray.Data(),tender);

  if(col=="PbPb" && maxCen>=0) kAnaMergedPi0+=Form("Cen%d_%d",minCen,maxCen);
    
  printf("<<<< NAME: %s >>>>>\n",kAnaMergedPi0.Data());
    
  // #### Configure analysis ####
    
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
 
  // General frame setting and configuration
  maker->SetReader   ( ConfigureReader   (col,simulation,clustersArray,tender,nonLinOn,trigger,rejectEMCTrig,minCen,maxCen,printSettings,debug) ); 
  maker->SetCaloUtils( ConfigureCaloUtils(col,simulation,clustersArray,tender,nonLinOn,year,                               printSettings,debug) );
                       
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important
  
  //
  // Merged pi0 analysis
  //
  
  // Track matching rejection on
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(1,1,1,1,col,simulation,year,printSettings,debug), n++); // all cuts
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(0,1,1,1,col,simulation,year,printSettings,debug), n++); // no ss cut
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(0,0,1,1,col,simulation,year,printSettings,debug), n++); // no ss and asy cuts
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(1,0,1,1,col,simulation,year,printSettings,debug), n++); // no asy cut
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(1,0,0,1,col,simulation,year,printSettings,debug), n++); // no asy and mass cut
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(1,1,0,1,col,simulation,year,printSettings,debug), n++); // no mass cut

  // Track matching rejection off
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(1,1,1,0,col,simulation,year,printSettings,debug), n++); // all cuts
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(0,1,1,0,col,simulation,year,printSettings,debug), n++); // no ss cut
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(0,0,1,0,col,simulation,year,printSettings,debug), n++); // no ss and asy cuts
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(1,0,1,0,col,simulation,year,printSettings,debug), n++); // no asy cut
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(1,0,0,0,col,simulation,year,printSettings,debug), n++); // no asy and mass cut
  maker->AddAnalysis(ConfigurePi0EbEAnalysis(1,1,0,0,col,simulation,year,printSettings,debug), n++); // no mass cut
  
  // Needed for pi0 primary spectra histogram, activate only for simulations
  if(simulation) maker->AddAnalysis(ConfigurePi0KineAnalysis(col,year,printSettings,debug), n++); 
   
  maker->SetAnaDebug(debug)  ;
  
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;
  
  if( simulation || !trigger.Contains("EMC") ) maker->SwitchOffDataControlHistograms();
  
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

  
  if(printSettings) maker->Print("");
  
  // Create task
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("Pi0SS%s",kAnaMergedPi0.Data()));
  
  task->SetDebugLevel(debug);

  //task->SetBranches("ESD:AliESDRun.,AliESDHeader");
  //task->SetConfigFileName(""); //Don't configure the analysis via configuration file.

  task->SetAnalysisMaker(maker);
  
  mgr->AddTask(task);
  
  //Create containers
  
  if(outputfile.Length()==0) outputfile = AliAnalysisManager::GetCommonFileName(); 
  
  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(kAnaMergedPi0, TList::Class(), 
                                                             AliAnalysisManager::kOutputContainer, 
                                                             Form("%s",outputfile.Data()));
	
  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",kAnaMergedPi0.Data()), TList::Class(), 
                                                             AliAnalysisManager::kParamContainer, 
                                                             "AnalysisParameters.root");
  
  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  //if(!kData.Contains("delta")   && outputAOD) mgr->ConnectOutput (task, 0, mgr->GetCommonOutputContainer());
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);
    
   
  UInt_t mask =  SetTriggerMaskFromName(trigger);
  task->SelectCollisionCandidates(mask);
  
  
  return task;
}

///
/// Configure the class handling the events and cluster/tracks filtering.
///
AliCaloTrackReader * ConfigureReader(TString col,           Bool_t simulation, 
                                     TString clustersArray, Bool_t tender, 
                                     Bool_t nonLinOn,
                                     TString trigger,       Bool_t rejectEMCTrig, 
                                     Int_t   minCen,        Int_t  maxCen,
                                     Bool_t printSettings,  Int_t   debug        )
{
  // Get the data type ESD or AOD
  AliAnalysisManager * mgr = AliAnalysisManager::GetAnalysisManager();
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); 
  
  AliCaloTrackReader * reader = 0;
  if     (inputDataType == "AOD") reader = new AliCaloTrackAODReader();
  else if(inputDataType == "ESD") reader = new AliCaloTrackESDReader();
  else printf("AliCaloTrackReader::ConfigureReader() - Data not known InputData=%s\n",inputDataType.Data());
  
  reader->SetDebug(debug);//10 for lots of messages
  
  //
  // MC settings
  //
  // Check if kine stack is available, independent of request of simulation
  Bool_t useKinematics = kFALSE;
  useKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;
  
  if(simulation)
  {
    if (!useKinematics && inputDataType=="AOD") useKinematics = kTRUE; //AOD primary should be available ...
  }

  if(useKinematics)
  {
    if(inputDataType == "ESD")
    {
      reader->SwitchOnStack();          
      reader->SwitchOffAODMCParticles(); 
    }
    else if(inputDataType == "AOD")
    {
      reader->SwitchOffStack();          
      reader->SwitchOnAODMCParticles(); 
    }
  }  
  
  // In case of Pythia pt Hard bin simulations (jet-jet, gamma-jet)
  // reject some special events that bother the cross section
  if(simulation)
  {
    // Event rejection cuts for jet-jet simulations, do not use in other
    reader->SetPtHardAndJetPtComparison(kTRUE);
    reader->SetPtHardAndJetPtFactor(4);
    
    reader->SetPtHardAndClusterPtComparison(kTRUE);
    reader->SetPtHardAndClusterPtFactor(1.5);
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

  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();
  
  // Shower shape smearing
  // Set it in the train configuration page not here for the moment
//  if(simulation)
//  {
//    reader->SwitchOffShowerShapeSmearing(); // Active only on MC, off by default
//    reader->SetShowerShapeSmearWidth(0.005);  
//  }

  //
  // Tracks
  //
  reader->SwitchOnCTS();

  reader->SwitchOffUseTrackTimeCut();
  reader->SetTrackTimeCut(0,50);
  
  reader->SwitchOnFiducialCut();
  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;

  reader->SwitchOffUseTrackDCACut();
  //reader->SetTrackDCACut(0,0.0105);
  //reader->SetTrackDCACut(1,0.035);
  //reader->SetTrackDCACut(2,1.1);
  
  if(inputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
    
    //AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
    //reader->SetTrackCuts(esdTrackCuts);
    
    AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
    reader->SetTrackCuts(esdTrackCuts);
    AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
    reader->SetTrackComplementaryCuts(esdTrackCuts2);
    
    reader->SwitchOnConstrainTrackToVertex();
  }
  else if(inputDataType=="AOD")
  {
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SwitchOnAODTrackSharedClusterSelection();
    reader->SetTrackStatus(AliVTrack::kITSrefit);
    
    //reader->SwitchOnAODPrimaryTrackSelection(); // Used in preliminary results of QM from Nicolas and Xiangrong?
    //reader->SwitchOnTrackHitSPDSelection();    // Check that the track has at least a hit on the SPD, not much sense to use for hybrid or TPC only tracks
    //reader->SetTrackFilterMask(128);           // Filter bit, not mask, use if off hybrid, TPC only
  }
  
  //
  // Calorimeter
  //
  if(clustersArray == "" && !tender) 
  {
    printf("**************** Standard EMCAL clusters branch analysis **************** \n");
    reader->SwitchOnClusterRecalculation();
    // Check in ConfigureCaloUtils that the recalibration and bad map are ON 
  }
  else 
  {
    printf("**************** Input for analysis is Clusterizer %s **************** \n", clustersArray.Data());
    reader->SetEMCALClusterListName(clustersArray);
    reader->SwitchOffClusterRecalculation();
  }  
  
  // Time cuts
  reader->SwitchOffUseParametrizedTimeCut();
  if(simulation) 
  {
    reader->SwitchOffUseEMCALTimeCut();
    reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
  }
  else
  {
    reader->SwitchOnUseEMCALTimeCut();
    reader->SetEMCALTimeCut(-25,20);
  }

  // CAREFUL
  if(nonLinOn) reader->SwitchOnClusterELinearityCorrection();
  else         reader->SwitchOffClusterELinearityCorrection();
  
  
  reader->SwitchOnEMCALCells();  
  reader->SwitchOnEMCAL();
  
  reader->SwitchOffPHOSCells();
  reader->SwitchOffPHOS();
  
  
  //-----------------
  // Event selection
  //-----------------
  
  //if(!simulation) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
  
  // Event triggered by EMCal selection settings
  reader->SwitchOffTriggerPatchMatching();
  reader->SwitchOffBadTriggerEventsRemoval();
  
  if( rejectEMCTrig > 0 && !simulation && (trigger.Contains("EMC") || trigger.Contains("L")))
  {
    printf("=== Remove bad triggers === \n");
    reader->SwitchOnTriggerPatchMatching();
    reader->SwitchOnBadTriggerEventsRemoval();
    
//    reader->SetTriggerPatchTimeWindow(8,9); // default values
//    if     (kRunNumber < 146861) reader->SetEventTriggerL0Threshold(3.);
//    else if(kRunNumber < 154000) reader->SetEventTriggerL0Threshold(4.);
//    else if(kRunNumber < 165000) reader->SetEventTriggerL0Threshold(5.5);
//    //redefine for other periods, triggers
//    
//    if(kRunNumber < 172000)
//    {
//      reader->SetEventTriggerL1Bit(4,5); // current LHC11 data
//      printf("\t Old L1 Trigger data format!\n");
//    }
//    else
//    {
//      reader->SetEventTriggerL1Bit(6,8); // LHC12-13 data
//      printf("\t Current L1 Trigger data format!\n");
//    }
    
    if(clustersArray != "" || tender)
    {
      printf("Trigger cluster calibration OFF\n");
      reader->SwitchOffTriggerClusterTimeRecal() ;
    }
    
  }

  //reader->RejectFastClusterEvents() ;
  
  reader->SwitchOnEventTriggerAtSE();
  
  reader->SetZvertexCut(10.);               // Open cut
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
  reader->SwitchOnRejectNoTrackEvents();

  reader->SwitchOffV0ANDSelection() ;       // and besides v0 AND
  reader->SwitchOffPileUpEventRejection();  // remove pileup by default off, apply it only for MB not for trigger
  
  if(col=="PbPb") 
  {
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
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
AliCalorimeterUtils* ConfigureCaloUtils(TString col,           Bool_t simulation,
                                        TString clustersArray, Bool_t tender, 
                                        Bool_t  nonLinOn,      Int_t year, 
                                        Bool_t  printSettings, Int_t   debug)
{
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  
  cu->SetDebug(debug);

  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder (2);
  
  cu->SetNumberOfSuperModulesUsed(10);
  
  if     (year == 2010) cu->SetNumberOfSuperModulesUsed(4);
  else if(year <= 2013) cu->SetNumberOfSuperModulesUsed(10);
  else if(year >  2013) cu->SetNumberOfSuperModulesUsed(20);
  else                  cu->SetNumberOfSuperModulesUsed(10);
  
  printf("xxx Number of SM set to <%d> xxx\n",cu->GetNumberOfSuperModulesUsed());
  
  // Search of local maxima in cluster
  if(col=="pp")
  {
    cu->SetLocalMaximaCutE(0.1);
    cu->SetLocalMaximaCutEDiff(0.03);
  }
  else 
  {
    cu->SetLocalMaximaCutE(0.2);
    cu->SetLocalMaximaCutEDiff(0.03);
  }
  
  cu->SwitchOffRecalculateClusterTrackMatching();
  
  cu->SwitchOnBadChannelsRemoval() ;
  
  // EMCAL settings

  if(!simulation)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();
  
  // calibrations
  Bool_t calibEner = kFALSE;
  Bool_t calibTime = kFALSE;
  cu->SwitchOffRecalibration(); 
  cu->SwitchOffRunDepCorrection();
  
  if( !tender )
  {
    cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
    cu->SwitchOnRunDepCorrection();
    
    calibEner = kTRUE;
    calibTime = kTRUE;
  }
  
  if( simulation )
  {
    calibEner = kFALSE;
    calibTime = kFALSE;

    cu->SwitchOffRecalibration(); // Check the reader if it is taken into account during filtering
    cu->SwitchOffRunDepCorrection();
  }
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          simulation,                             
                          kTRUE,      // exotic
                          nonLinOn,   // Non linearity
                          calibEner,  // E calib
                          kTRUE,      // bad map
                          calibTime); // time calib   
  
  if( calibTime ) recou->SetExoticCellDiffTimeCut(50);
  
  if( nonLinOn )  cu->SwitchOnCorrectClusterLinearity();
  
  printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
  printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());

  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  if(printSettings) cu->Print("");
  
  return cu;
}

///
/// Configure the task doing the pi0 even by event selection via the split method
///
AliAnaPi0EbE* ConfigurePi0EbEAnalysis(Bool_t useSS,     Bool_t useAsy,     
                                      Bool_t useMass,   Bool_t tm,
                                      TString col,      Bool_t simulation, Int_t  year,      
                                      Bool_t  printSet, Int_t debug                    )
{  
  printf("Set: TM%d_SS%d_Asy%d_Mass%d_%s \n", tm,useSS,useAsy,useMass,kAnaMergedPi0.Data());
  AliAnaPi0EbE *ana = new AliAnaPi0EbE();
  
  ana->SetAnalysisType(AliAnaPi0EbE::kSSCalo);
  
  ana->SwitchOffAllNLMHistoFill();
  ana->SwitchOffSelectedClusterHistoFill();

  ana->SwitchOffFillWeightHistograms();

  if(tm) ana->SwitchOnTrackMatchRejection() ;
  else   ana->SwitchOffTrackMatchRejection() ;
  ana->SwitchOffTMHistoFill() ;

  ana->SetCalorimeter("EMCAL");
  
  // Branch AOD settings
  ana->SetOutputAODName(Form("Pi0SS_TM%d_SS%d_Asy%d_Mass%d_%s",tm,useSS,useAsy,useMass,kAnaMergedPi0.Data()));
  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
   
  // Set Histograms name tag, bins and ranges
  ana->AddToHistogramsName(Form("AnaPi0SSEbE_TM%d_SS%d_Asy%d_Mass%d_",tm,useSS,useAsy,useMass));

  SetAnalysisCommonParameters(ana,"EMCAL",year,col,simulation,printSet,debug); // see method below
  
  // cluster splitting settings
  ana->SetMinEnergy(6);
  ana->SetMaxEnergy(100.);
  
  ana->SetNLMMinEnergy(0, 10);
  ana->SetNLMMinEnergy(1, 6);
  ana->SetNLMMinEnergy(2, 6);
  
  // NLM cut, used in all, exclude clusters with more than 2 maxima
  ana->SetNLMCut(1, 2) ;
  
  //
  ana->SetMinDistanceToBadChannel(2, 4, 6);
  ana->SwitchOnSplitClusterDistToBad();
  ana->SetTimeCut(-1e10,1e10); // Open time cut
  
  AliCaloPID* caloPID = ana->GetCaloPID();
  
  caloPID->SetSplitWidthSigma(3); // cut at 3 sigma of the mean pi0 peak.
  
  if(!useSS)
  {
    printf("Do not apply SS cut on merged pi0 analysis \n");
    caloPID->SwitchOffSplitShowerShapeCut() ;
    caloPID->SetClusterSplittingM02Cut(0.1,10);
  }
  else
  {
    caloPID->SetClusterSplittingM02Cut(0.3,4); // Do the selection in the analysis class and not in the PID method to fill SS histograms
    caloPID->SwitchOnSplitShowerShapeCut() ;
  }
  
  if(useAsy)
  {
    caloPID->SwitchOnSplitAsymmetryCut() ;
    ana->GetCaloPID()->SetSubClusterEnergyMinimum(0,2);
    ana->GetCaloPID()->SetSubClusterEnergyMinimum(1,0.5);
    ana->GetCaloPID()->SetSubClusterEnergyMinimum(2,0.5);
  }
  else
  {
    caloPID->SwitchOffSplitAsymmetryCut() ;
  }
  
  if(!useMass)
  {
    ana->SwitchOffSplitClusterDistToBad();
    caloPID->SwitchOnSimpleSplitMassCut();
    caloPID->SetPi0MassRange(-1, 10000); // Open the cut    
  }

  
  // For Pi0 only if  SwitchOnSimpleSplitMassCut()
  caloPID->SetPi0MassRange(0.10, 0.18);
  caloPID->SetEtaMassRange(0.50, 0.60);
  caloPID->SetPhotonMassRange(0.00, 0.08);
  
  caloPID->SetClusterSplittingMinNCells(6);
  
  //caloPID->SetSplitEnergyFractionMinimum(0, 0.95);
  //caloPID->SetSplitEnergyFractionMinimum(1, 0.95);
  //caloPID->SetSplitEnergyFractionMinimum(2, 0.8);
  
  if(col=="PbPb")
  {
    caloPID->SetClusterSplittingMinNCells(4);
    //caloPID->SetPi0MassShiftHighECell(0.005);
  }
  
  return  ana;
}

///
/// Configure the task pi0/eta kinematic histograms
///
AliAnaPi0* ConfigurePi0KineAnalysis(TString col, Int_t year, Bool_t printSet, Int_t debug)
{
  AliAnaPi0 *ana = new AliAnaPi0();

  ana->SwitchOnFillOnlyMCAcceptanceHisto();
  
  ana->SwitchOffFillOriginHisto();
  
  ana->SwitchOnRealCaloAcceptance();
  //ana->SwitchOnFiducialCut();
  //ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.65, 81, 179);
  
  // Calorimeter settings
  ana->SetCalorimeter("EMCAL");
  
  //Set Histograms name tag, bins and ranges
  
  ana->AddToHistogramsName("AnaPi0Kine_");
  
  SetAnalysisCommonParameters(ana,"EMCAL",year,col,kTRUE,printSet,debug); // see method below

  return ana;
  
}

///
/// Set common histograms binning 
/// plus other analysis common settings like TRD covered super modules
/// the activation of the MC dedicated histograms and the activation of 
/// the debug mode
///
void SetAnalysisCommonParameters(AliAnaCaloTrackCorrBaseClass* ana, 
                                 TString calorimeter,   Int_t  year, 
                                 TString col,           Bool_t simulation,
                                 Bool_t  printSettings, Int_t  debug)
{
  //
  // Histograms ranges
  //
  AliHistogramRanges* histoRanges = ana->GetHistogramRanges();
  
  histoRanges->SetHistoPtRangeAndNBins(0, 100, 200) ; // Energy and pt histograms
  
  if(calorimeter=="EMCAL")
  {
    if(year==2010)
    {
      histoRanges->SetHistoPhiRangeAndNBins(78*TMath::DegToRad(), 122*TMath::DegToRad(), 78) ;
      histoRanges->SetHistoXRangeAndNBins(-230,90,120); // QA
      histoRanges->SetHistoYRangeAndNBins(370,450,40);  // QA
    }
    else if ( year < 2014 )
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
  else if(calorimeter=="PHOS") 
  {
    histoRanges->SetHistoPhiRangeAndNBins(250*TMath::DegToRad(), 320*TMath::DegToRad(), 70) ;
    histoRanges->SetHistoEtaRangeAndNBins(-0.13, 0.13, 130) ;
  }
  else if(calorimeter=="CTS")
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
    
  histoRanges->SetHistoShowerShapeRangeAndNBins(-0.1, 4.9, 500);
  
  // Invariant mass histo
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
  if(col=="PbPb")
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
  
  //
  // TRD SM
  //
  if     (year == 2011) ana->SetFirstSMCoveredByTRD( 6);
  else if(year == 2012 || 
          year == 2013) ana->SetFirstSMCoveredByTRD( 4);
  else                  ana->SetFirstSMCoveredByTRD(-1);

  //
  // MC histograms?
  //
  if(simulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else           ana->SwitchOffDataMC() ;
  
  //Set here generator name, default pythia
  //ana->GetMCAnalysisUtils()->SetMCGenerator("");

  //
  // Debug
  //
  if(printSettings) ana->Print("");
  
  ana->SetDebug(debug); // 10 for lots of messages
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
  else return AliVEvent::kAny;
}

