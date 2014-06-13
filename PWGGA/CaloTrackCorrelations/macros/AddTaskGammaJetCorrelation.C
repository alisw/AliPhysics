// Configuration macro for analysis of photon-jet analysis
// Adam Matyja based on AddTaskIsoPhoton by Gustavo Conesa & Marie Germain.
// need to be updated to the newest version

//new macro
Bool_t  kSimulation    = kFALSE;
Bool_t  kUseKinematics = kFALSE;
Bool_t  kOutputAOD     = kFALSE;
Bool_t  kEventSelection= kFALSE;
Bool_t  kExotic        = kTRUE;
Bool_t  kNonLinearity  = kFALSE;
TString kCollisions    = "pp";
TString kTrig          = "EMC7" ;
TString kFiredTrig     = "EG1" ;
TString kClusterArray  = "";
TString kData          = ""; // MC or deltaAOD
TString kInputDataType = "ESD";
Bool_t  kTM            = kTRUE;
Bool_t  kRecalTM       = kTRUE;
Int_t   kMinCen        = -1;
Int_t   kMaxCen        = -1;
TString kName          = "";
Bool_t  kCalibE        = kTRUE;
Bool_t  kCalibT        = kTRUE;
Bool_t  kBadMap        = kTRUE;
Bool_t  kTender        = kFALSE;
Bool_t  kMix           = kFALSE;
Int_t   kRunNumber     = -1;

AliAnalysisTaskCaloTrackCorrelation *AddTaskGammaJetCorrelation(//const Float_t  cone          = 0.4,
						//const Float_t  pth           = 0.5,
						const TString  data          = "",
						const TString  calorimeter   = "EMCAL", //done
						const Bool_t   simulation    = kFALSE,
						const Bool_t   eventsel      = kFALSE,
						const Bool_t   exotic        = kTRUE,
						const Bool_t   nonlin        = kFALSE,
						TString        outputfile    = "",
						const TString  col           = "pp", 
						const TString  trigger       = "MB", 
						const TString  firedTrigger  = "EG1",
						const TString  clustersArray = "V1",
						const Bool_t   mix           = kTRUE,
						const Bool_t   recaltm       = kTRUE,
						const Bool_t   tm            = kTRUE,
						const Int_t    minCen        = -1,
						const Int_t    maxCen        = -1,
						const Bool_t   calibE        = kTRUE,
						const Bool_t   badmap        = kTRUE,
						const Bool_t   calibT        = kTRUE,
						const Bool_t   tender        = kFALSE,
						const Bool_t   outputAOD     = kFALSE, 
						const Int_t    debug         = -1,//done
						const Bool_t   printSettings = kFALSE,//done
						const Double_t scaleFactor   = -1,
						const Int_t    runNumber     = -1
						)
{
  // Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
  
  kSimulation    = simulation;
  kCollisions    = col;
  kExotic        = exotic;
  kNonLinearity  = nonlin;
  kTrig          = trigger;
  kFiredTrig     = firedTrigger;
  kClusterArray  = clustersArray;
  kData          = data;
  kOutputAOD     = outputAOD;
  kTM            = tm;
  kRecalTM       = recaltm;
  kMinCen        = minCen;
  kMaxCen        = maxCen;
  kEventSelection= eventsel;
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
  
  //kName = Form("%s_Trig%s_Cl%s_TM%d_R%1.1f_Pt%1.1f",calorimeter.Data(), kTrig.Data(),kClusterArray.Data(),kTM,cone,pth);
  //  kName = Form("%s_Trig%s_Cl%s_TM%d",calorimeter.Data(), kTrig.Data(),kClusterArray.Data(),kTM);
  kName = Form("%s_Trig%s_Fired%s_Cl%s_TM%d",calorimeter.Data(), kTrig.Data(),kFiredTrig.Data(),kClusterArray.Data(),kTM);//<<<---changed here

  
  if(kCollisions=="PbPb" && kMaxCen>=0) kName+=Form("Cen%d_%d",kMinCen,kMaxCen);
    
  printf("<<<< NAME: %s >>>>>\n",kName.Data());
  
  // #### Configure analysis ####
    
  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();
  printf("SCALE FACTOR %e\n",scaleFactor);
  maker->SetScaleFactor(scaleFactor); // for MC, negative (not scaled) by default
  
  // General frame setting and configuration
  maker->SetReader   (ConfigureReader(calorimeter,debug,printSettings)   ); 
  maker->SetCaloUtils(ConfigureCaloUtils(debug,printSettings)); 
  
  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important

  // Isolation settings
  Int_t partInCone = AliIsolationCut::kNeutralAndCharged; // kOnlyCharged;
  Int_t thresType  = AliIsolationCut::kPtThresIC;//  AliIsolationCut::kSumPtFracIC ; 
  Float_t cone = -1;
  Float_t pth  = -1;
  
  maker->AddAnalysis(ConfigurePhotonAnalysis(calorimeter,debug,printSettings), n++); // Photon cluster selection
  maker->AddAnalysis(ConfigureIsolationAnalysis(calorimeter,"Photon", partInCone,thresType, cone, pth,kFALSE,debug,printSettings), n++); // Photon isolation   
  maker->AddAnalysis(ConfigurePhotonJetAnalysis(debug,printSettings), n++);// photon-jet correlation analysis

  maker->SetAnaDebug(debug)  ;
  maker->SwitchOnHistogramsMaker()  ;
  if(kData.Contains("delta")) maker->SwitchOffAODsMaker() ;
  else                        maker->SwitchOnAODsMaker()  ;
  
  if(printSettings) maker->Print("");
  
  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, calorimeter.Data());
  // Create task
  
  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("CaloTrackCorr%s",kName.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(debug);
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

//____________________________________
AliCaloTrackReader * ConfigureReader(TString calorimeter = "EMCAL",Int_t debug = -1,Bool_t printSettings = kFALSE)
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
  
  reader->SetDebug(debug);//10 for lots of messages
  //reader->SetDebug(10);//10 for lots of messages
  //reader->SetDebug(2);//10 for lots of messages

  if(kSimulation)
  {
    // Event rejection cuts for jet-jet simulations, do not use in other
    reader->SetPtHardAndJetPtComparison(kTRUE);
    reader->SetPtHardAndJetPtFactor(4);
    
    reader->SetPtHardAndClusterPtComparison(kTRUE);
    reader->SetPtHardAndClusterPtFactor(1.5);
  }

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
  reader->SetInputNonStandardJetBranchName("clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00");//in pp 7 LHC12a15f MC
  reader->SwitchOffBackgroundJets();
  //reader->SetInputBackgroundJetBranchName("jeteventbackground_clustersAOD_KT04_B0_Filter00768_Cut00150_Skip00");//in pp 7 LHC13e4 MC


  // Time cuts
  if(kSimulation) 
  {
    reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
    reader->SwitchOffUseEMCALTimeCut() ;   
    reader->SwitchOffUseParametrizedTimeCut();
  }
  else
  {
    if(kCalibT)
    { 
      printf("Set time cut parameters for run %d\n",kRunNumber);
      //reader->SetEMCALTimeCut(-20,20); 
      reader->SwitchOnUseParametrizedTimeCut();
      if     (kRunNumber >= 151636 && kRunNumber <= 155384 )
      {
        printf("Set time parameters for LHC11c\n");
        reader->SetEMCALParametrizedMinTimeCut(0,-5  ); reader->SetEMCALParametrizedMinTimeCut(1,-1 ); reader->SetEMCALParametrizedMinTimeCut(2, 1.87); reader->SetEMCALParametrizedMinTimeCut(3, 0.4);   
        reader->SetEMCALParametrizedMaxTimeCut(0, 3.5); reader->SetEMCALParametrizedMaxTimeCut(1, 50); reader->SetEMCALParametrizedMaxTimeCut(2, 0.15); reader->SetEMCALParametrizedMaxTimeCut(3, 1.6);   
      }
      else if(kRunNumber >= 156447 && kRunNumber <= 159635 )
      {
        printf("Set time parameters for LHC11d\n");
        reader->SetEMCALParametrizedMinTimeCut(0,-5);  reader->SetEMCALParametrizedMinTimeCut(1,-1 );  reader->SetEMCALParametrizedMinTimeCut(2, 3.5 ); reader->SetEMCALParametrizedMinTimeCut(3, 1.  );   
        reader->SetEMCALParametrizedMaxTimeCut(0, 5);  reader->SetEMCALParametrizedMaxTimeCut(1, 50);  reader->SetEMCALParametrizedMaxTimeCut(2, 0.45); reader->SetEMCALParametrizedMaxTimeCut(3, 1.25);   
      }
      else 
      {
        printf("*** Fixed time cut 20 ns *** \n");
        reader->SetEMCALTimeCut(-20,20);
      }
      
    }
    else
    {
      reader->SwitchOffUseEMCALTimeCut();
      reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
      reader->SwitchOffUseParametrizedTimeCut(); 
    }
  }  
  
  reader->SwitchOnFiducialCut();
  //reader->SwitchOffFiducialCut();// <<<----changed here 
  //reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;
  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.9, 0, 360) ;// <<<--- changed here

  // Tracks

  //  reader->SwitchOffCTS();//here changed 0n->off
  reader->SwitchOnCTS();//here changed 0n->off

  if(kInputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_ROOT/PWGJE/macros/CreateTrackCutsPWGJE.C");
    AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
    reader->SetTrackCuts(esdTrackCuts);
    reader->SwitchOnConstrainTrackToVertex();
  }
  else if(kInputDataType=="AOD")
  {
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SetTrackStatus(AliVTrack::kITSrefit);    
    //reader->SetTrackFilterMask(128);           // Filter bit, not mask, use if off hybrid
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
  
  if(calorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();  
    reader->SwitchOnEMCAL();
  }
  if(calorimeter == "PHOS") { 
    reader->SwitchOnPHOSCells();  
    reader->SwitchOnPHOS();
  }
  
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
  
  //if(!kUseKinematics) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma
  if(!kUseKinematics) {
    if(kCollisions =="pPb" && kTrig=="EMCEGA") {
      reader->SetFiredTriggerClassName(kFiredTrig);
    }
  }

  // For mixing with AliAnaParticleHadronCorrelation switch it off
  if(kMix)
  {
    reader->SwitchOffEventTriggerAtSE();
    UInt_t mask =  SetTriggerMaskFromName();
    reader->SetEventTriggerMaks(mask); // Only for mixing and SwitchOffEventTriggerAtSE();
    //reader->SetMixEventTriggerMaks(AliVEvent::kMB); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT
    reader->SetMixEventTriggerMaks(AliVEvent::kAnyINT); // Careful, not all productions work with kMB, try kINT7, kINT1, kAnyINT

    printf("---Trigger selection done in AliCaloTrackReader!!!\n");
  }
  else 
    reader->SwitchOnEventTriggerAtSE();
  
  reader->SetZvertexCut(10.);                // Open cut
  //reader->SetZvertexCut(50.);                // Open cut <<<--- changed here
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
  //reader->SwitchOffPrimaryVertexSelection(); // and besides primary vertex <<<--- changed here

  if(kEventSelection)
  {
    reader->SwitchOnPileUpEventRejection();   // remove pileup by default  
    //    reader->SwitchOnEventSelection();         // remove pileup by default
    reader->SwitchOnV0ANDSelection() ;        // and besides v0 AND
  }
  else 
  {
    reader->SwitchOffPileUpEventRejection();  // remove pileup by default   
    //    reader->SwitchOffEventSelection();         // remove pileup by default
    reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND
  }
    
  if(kCollisions=="PbPb") 
  {
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(10);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
    reader->SetCentralityBin(kMinCen,kMaxCen); // Accept all events, if not select range
    
    // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
    reader->SetEventPlaneMethod("V0");
  }
  
  if(printSettings) reader->Print("");
  reader->Print("");
  
  return reader;
  
}

//_______________________________________
AliCalorimeterUtils* ConfigureCaloUtils(Int_t debug = -1,Bool_t print = kFALSE)
{
  
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(debug);
  
  // Remove clusters close to borders, at least max energy cell is 1 cell away 
  //cu->SetNumberOfCellsFromEMCALBorder(0);// <<<----changed here 
  cu->SetNumberOfCellsFromEMCALBorder(1);
  cu->SetNumberOfCellsFromPHOSBorder(2);
  
  // Search of local maxima in cluster
  if(kCollisions=="pp" || kCollisions=="pPb")
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

  if(kRecalTM) cu->SwitchOnRecalculateClusterTrackMatching(); // Done in clusterization
  else         cu->SwitchOffRecalculateClusterTrackMatching();
  
  cu->SwitchOnBadChannelsRemoval() ;
  //cu->SwitchOffBadChannelsRemoval() ;// <<<----changed here 
  
  //EMCAL settings

  if(!kSimulation)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();

  if(!kSimulation)
  {
    cu->SwitchOnRecalibration(); // Check the reader if it is taken into account during filtering
    if(kClusterArray == "" && !kTender) cu->SwitchOnRunDepCorrection(); 
  }

  cu->SwitchOnEMCALOADB();//FIX ME!!!

  gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          kSimulation,                             
                          kExotic,
                          kNonLinearity,
                          kCalibE, 
                          kBadMap,
                          kCalibT);   
  recou->SetExoticCellDiffTimeCut(1e10);
  
  if( kNonLinearity ) 
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
    
  
  printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
  printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
  
    
  // PHOS 
  cu->SwitchOffLoadOwnPHOSGeometryMatrices();
    
  if(printSettings) cu->Print("");
  
  return cu;
  
}

//_____________________________________
AliAnaPhoton* ConfigurePhotonAnalysis(TString calorimeter = "EMCAL",Int_t debug = -1,Bool_t printSettings = kFALSE)
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
    //ana->SetNCellCut(0);// At least 2 cells <<<----changed here 
    ana->SetMinEnergy(0.3); // avoid mip peak at E = 260 MeV
    ana->SetMaxEnergy(1000); 
    ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off 
    // restrict to less than 100 ns when time calibration is on 
    ana->SetMinDistanceToBadChannel(2, 4, 6); 
    // Not useful if M02 cut is already strong
    //ana->SetNLMCut(1, 2) ;
    //ana->SetNLMCut(1, 10) ;//<<<----changed here 
    ana->SetNLMCut(1, 1) ;//<<<----changed here 
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
  //ana->SwitchOffCaloPID(); //<<<----changed here 
  ana->SwitchOnCaloPID(); // do PID selection, unless specified in GetCaloPID, selection not based on bayesian
  AliCaloPID* caloPID = ana->GetCaloPID();
  //Not used in bayesian
  
  //EMCAL
  caloPID->SetEMCALLambda0CutMax(0.50);//0.27 was before//0.50//<<<----changed here
  caloPID->SetEMCALLambda0CutMin(0.10);
  
  caloPID->SetEMCALDEtaCut(0.025);
  caloPID->SetEMCALDPhiCut(0.030);
    
  //PHOS
  caloPID->SetPHOSDispersionCut(2.5);
  caloPID->SetPHOSRCut(2.);
  if(kInputData=="AOD") caloPID->SetPHOSRCut(2000.); // Open cut since dX, dZ not stored
      
  ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection <<<--- changed here
  //ana->SwitchOffFillShowerShapeHistograms();  // Filled before photon shower shape selection
  ana->SwitchOffFillPileUpHistograms();

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
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below
  
  // Number of particle type MC histograms
  ana->FillNOriginHistograms(20);
  ana->FillNPrimaryHistograms(20);
  
  ConfigureMC(ana);
  
  if(printSettings) ana->Print("");
  
  return ana;
  
}

//____________________________________________________________________________________________________
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString calorimeter = "EMCAL",
						    TString particle="Photon", 
                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                    Float_t cone = 0.3,
                                                    Float_t pth  = 0.3,
                                                    Bool_t multi      = kFALSE,
						    Int_t debug = -1,
						    Bool_t printSettings = kFALSE)
{
  
  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
  //ana->SetDebug(debug);
  ana->SetDebug(debug);
    
  ana->SwitchOnFiducialCut();
  //Avoid borders of EMCal
  if(calorimeter=="EMCAL")
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;

  // Same Eta as EMCal, cut in phi if EMCAL was triggering
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    //if(kTrig.Contains("EMC"))
    //  ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 260, 360) ;
    //else
      ana->GetFiducialCut()->SetSimpleCTSFiducialCut  (0.6, 0, 360) ;    
  }
  
  ana->SetMinPt(3);//<<---changed here
  
  // Input / output delta AOD settings
  
  ana->SetInputAODName(Form("%s%s",particle.Data(),kName.Data()));
  ana->SetAODObjArrayName(Form("IC%s_%s",particle.Data(),kName.Data())); 
  
  ana->SetCalorimeter(calorimeter);
  
  if(!kTM)  ana->SwitchOnTMHistoFill();
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
    if(kCollisions=="pp") 
    {
      ic->SetPtThreshold(1.);//<<---changed here was 0.5,1
      ic->SetConeSize(0.3);//<<---changed here was 0.4
    }
    if(kCollisions=="pPb")
      {
	ic->SetPtThreshold(1.0);
	ic->SetConeSize(0.3);
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

  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below
  
  ana->SetHistoPtInConeRangeAndNBins(0, 50 , 250);
  ana->SetHistoPtSumRangeAndNBins   (0, 100, 250);
  
  if(particle=="Hadron"  || particle.Contains("CTS"))
  {
    ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
    ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;
  }
  
  ConfigureMC(ana);
  
  if(printSettings) ic ->Print("");
  if(printSettings) ana->Print("");
  
  return ana;
  
}

//________________________________________________________
void ConfigureMC(AliAnaCaloTrackCorrBaseClass* ana)
{
  if(kSimulation) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else            ana->SwitchOffDataMC() ;

  //Set here generator name, default pythia
  //ana->GetMCAnalysisUtils()->SetMCGenerator("");
}  

//________________________________________________________
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges,TString calorimeter = "EMCAL")
{
  // Set common bins for all analysis and MC histograms filling
    
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
  
}

//_________________
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
    return AliVEvent::kCentral;
  } 
  else if(kTrig=="CentralEGA")
  {
    printf("CaloTrackCorr trigger Central+EMCEGA\n");
    return (AliVEvent::kCentral | AliVEvent::kEMCEGA);
  } 
  else if(kTrig=="SemiCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral\n");
    return AliVEvent::kSemiCentral;
  }
  else if(kTrig=="SemiOrCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral Or Central\n");
    return (AliVEvent::kSemiCentral | AliVEvent::kCentral);
  }  
  else if(kTrig=="SemiOrCentralOrAnyINT")
    {
      printf("CaloTrackCorr trigger SemiCentral Or Central Or AnyINT\n");
      return (AliVEvent::kSemiCentral | AliVEvent::kCentral | AliVEvent::kAnyINT);
    }

}

AliAnaParticleJetFinderCorrelation* ConfigurePhotonJetAnalysis(Int_t debug = -1,Bool_t printSettings = kFALSE){

  AliAnaParticleJetFinderCorrelation *ana = new AliAnaParticleJetFinderCorrelation();
  ana->SetDebug(debug);
  TString particle="Photon";
  // ### Correlation with Jet Finder AOD output
  //AliAnaParticleJetFinderCorrelation *anacorrjet = new AliAnaParticleJetFinderCorrelation();
  //anacorrjet->SetInputAODName(Form("%s%s_Trig%s_Cl%s",particle.Data(),calorimeter.Data(), kTrig.Data(),kClusterArray.Data()));
  ana->SetInputAODName(Form("%s%s",particle.Data(),kName.Data()));

  ana->SwitchOffFiducialCut();

  ana->SetConeSize(0.4); //was 1 - cone to calculate FF
  ana->SelectIsolated(kTRUE); // do correlation with isolated photons <<---changed here
  ana->SetMakeCorrelationInHistoMaker(kFALSE);
  ana->SetPtThresholdInCone(0.150);//<<---- change here
  ana->SetDeltaPhiCutRange(TMath::Pi()/2.,TMath::Pi()*3./2.);//Mostly Open Cuts 
  ana->SetJetConeSize(0.4);//jet cone size / check the reco jet name
  ana->SetJetMinPt(5);//min jet pt
  ana->SetJetAreaFraction(0.8);//min area fraction was 0.6
  ana->SetMinPt(0.3);//min cluster pt repeated from reader
  ana->SetGammaConeSize(0.3);//isolation cone repeated from isolation ana


  ana->SetRatioCutRange(0.01,5.); //Mostly Open Cuts //0.01-5//<<---- change here
  ana->UseJetRefTracks(kTRUE); //Working now
  //Set Histograms bins and ranges
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below 0,100,200
  //ana->SetHistoPtRangeAndNBins(0, 50, 200) ;
  //      ana->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 100) ;
  //      ana->SetHistoEtaRangeAndNBins(-0.7, 0.7, 100) ;

  ana->SwitchOnNonStandardJetFromReader();
  ana->SwitchOnBackgroundJetFromReader();
  //background subtraction for photons
  ana->SwitchOnBackgroundSubtractionGamma();

  ana->SwitchOnSaveGJTree();
  //ana->SwitchOnMostOpposite();
  ana->SwitchOnMostEnergetic();


  if(kUseKinematics) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms
  //if(printSettings) 
    ana->Print("");

  return ana;





}
