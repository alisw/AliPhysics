/// \file AddTaskIsoPhoton.C
/// \ingroup CaloTrackCorrMacros
/// \brief Isolated photon spectra configuration.
///
/// Configuration of the isolated photon analysis analysis
/// based on AddTaskIsoPhoton by Gustavo Conesa & Marie Germain.
///
/// \author Marie Germain <Marie.Germain@cern.ch>, SUBATECH, main author.
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

TString kAnaIsoPhotonName = "";   /// Global name to be composed of the settings, used to set the AOD branch name

Int_t   kDebug         = -1;      /// Global debug level

TString kCalorimeter   = "EMCAL"; /// Global setting of calorimeter of photon
TString kData  = "" ;             /// Global string for data type
Bool_t  kPrint = 0  ;             /// Global bool for print option

///
/// Main method calling all the configuration
/// Creates a CaloTrackCorr task, configures it and adds it to the analysis manager.
///
/// The options that can be passed to the macro are:
/// \param cone : A float setting the isolation cone size
/// \param pth : A float setting the isolation pT threshold (sum of particles in cone or leading particle)
/// \param leading : select leading trigger clusters?
/// \param timecut : activate time cut
/// \param calorimeter : A string with he calorimeter used to measure the trigger particle
/// \param simu : A bool identifying the data as simulation
/// \param exotic : reject exotic clusters
/// \param nonlin : A bool to set the use of the non linearity correction
/// \param trigger : A string with the trigger class, abbreviated, defined in method belowSetTriggerMaskFromName()
/// \param tm : A bool to select neutral clusters as triggers
/// \param minCen : An int to select the minimum centrality, -1 means no selection
/// \param maxCen : An int to select the maximum centrality, -1 means no selection
/// \param deltaphicut : track matching residual cut in azimuth
/// \param deltaetacut : track matching residual cut in pseudo-rapidity
/// \param tmin : minimum cluster time
/// \param tmax : maximum cluster time
/// \param trackTcut : apply time cut on tracks
/// \param disttobad : value of cut on distance to bad channel
/// \param nlmMax : maximum value of shower shape parameter
/// \param qaan : activate detector qa analysis
/// \param primvtx : select primary vertex
/// \param notrackcut : reject events without tracks
/// \param rdmtrigger : do the analysis with random triggers
/// \param tag : name to pass to analysis generated branch and histo container
/// \param debug : An int to define the debug level of all the tasks
/// \param print : A bool to enable the print of the settings per task
/// \param tmInCone : A bool to enable the CPV in cone (to reject charged clusters in Eiso calculation)
/// \param SSsmearing : An integer to enable the shower shape smearing, 0: no smearing, 1: Smearing with Gustavo's settings, 2: Smearing with Astrid's settings
/// \param clustListName: name of list with clusters
///
AliAnalysisTaskCaloTrackCorrelation *AddTaskIsoPhoton(const Float_t  cone          = 0.4,
                                                      const Float_t  pth           = 2.,
                                                      const Bool_t   leading       = kFALSE,
                                                      const Bool_t   timecut       = kFALSE,
                                                      const TString  calorimeter   = "EMCAL",
                                                      const Bool_t   simu          = kFALSE,
                                                      const Bool_t   exotic        = kTRUE,
                                                      const Bool_t   nonlin        = kFALSE,
                                                      const TString  trigger       = "EMC7",
                                                      const Bool_t   tm            = kTRUE,
                                                      const Int_t    minCen        = -1,
                                                      const Int_t    maxCen        = -1,
                                                      const Float_t  deltaphicut   = 0.03,
                                                      const Float_t  deltaetacut   = 0.02,
                                                      const Float_t  tmin          = -30.,
                                                      const Float_t  tmax          = 30.,
                                                      const Bool_t   trackTcut     = kFALSE,
                                                      const Int_t    disttobad     = 2,
                                                      const Int_t    nlmMax        =  20,
                                                      const Bool_t   qaan          = kFALSE,
                                                      const Bool_t   primvtx       = kTRUE,
                                                      const Bool_t   notrackcut    = kTRUE,
                                                      const Bool_t   rdmtrigger    = kFALSE,
                                                      const TString  tag           = "",
                                                      const Int_t    debug         = -1,
                                                      const Bool_t   print         = kFALSE,
                                                      const Bool_t   tmInCone      = kTRUE,
                                                      const Int_t    SSsmearing    = 0,
                                                      const TString  clustListName = ""
                                                      )
{
kDebug = debug;
kCalorimeter  = calorimeter ;
kPrint = print ;

  printf("AddTaskIsoPhoton() - Settings: cone %2.2f, pth %2.2f, timeCut On %d, NLM max cut %d, calorimeter %s, simu %d, exotic %d, non lin %d, trigger %s, TM %d, qa %d, debug %d, centrality %d-%d\n",
                                         cone,    pth,    timecut   ,    nlmMax,      calorimeter.Data(),simu, exotic,    nonlin,     trigger.Data(), tm,    qaan,   debug, minCen, maxCen );

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

  Bool_t useKinematics = kFALSE;
  useKinematics = (mgr->GetMCtruthEventHandler())?kTRUE:kFALSE;

  // Name for containers

 if(tag!="")
 kAnaIsoPhotonName = Form("%s_Trig%s_TM%d_%1.3f_dB%d_R%1.1f_Pt%1.1f_%s",calorimeter.Data(), trigger.Data(),tm,deltaphicut,disttobad,cone,pth,tag.Data());
 else
 kAnaIsoPhotonName = Form("%s_Trig%s_TM%d_%1.3f_dB%d_R%1.1f_Pt%1.1f",calorimeter.Data(), trigger.Data(),tm,deltaphicut,disttobad,cone,pth);

  if(maxCen>=0) kAnaIsoPhotonName+=Form("Cen%d_%d",minCen,maxCen);

  printf("<<<< NAME: %s >>>>>\n",kAnaIsoPhotonName.Data());

  // #### Configure analysis ####

  AliAnaCaloTrackCorrMaker * maker = new AliAnaCaloTrackCorrMaker();

  //maker->SetScaleFactor(scaleFactor); // for MC, negative (not scaled) by default

  // General frame setting and configuration
  maker->SetReader   (ConfigureReader   (mgr->GetInputEventHandler()->GetDataType(),useKinematics,simu,
                                         calorimeter,nonlin, timecut, primvtx, notrackcut,tmin,tmax,trackTcut,minCen, maxCen, debug,print,SSsmearing,clustListName));
  maker->SetCaloUtils(ConfigureCaloUtils(nonlin,exotic,simu,timecut,debug,print));

  // Analysis tasks setting and configuration
  Int_t n = 0;//Analysis number, order is important

  // Isolation settings
  Int_t partInCone = AliIsolationCut::kNeutralAndCharged; // kOnlyCharged;
  //  Int_t thresType  = AliIsolationCut::kPtThresIC;//  AliIsolationCut::kSumPtFracIC ;
  Int_t thresType  = AliIsolationCut::kSumPtIC ;

 if(!rdmtrigger)
{
  // Photon analysis
  maker->AddAnalysis(ConfigurePhotonAnalysis(calorimeter,tm,deltaphicut,deltaetacut,disttobad,nlmMax,simu,debug,print), n++); // Photon cluster selection

  // Isolation analysis
  maker->AddAnalysis(ConfigureIsolationAnalysis(calorimeter,"Photon", partInCone,thresType,cone, pth,tm,leading,kFALSE,simu,debug,print,tmInCone), n++); // Photon isolation
}
else
{
  maker->AddAnalysis(ConfigureRandomTriggerAnalysis(), n++);
  maker->AddAnalysis(ConfigureIsolationAnalysis(calorimeter,Form("RandomTrigger%s",kCalorimeter.Data()), partInCone,thresType,cone, pth,tm,leading,kFALSE,simu,debug,print,tmInCone), n++);// Ghost trigger isolation
}


  // QA histograms on clusters or tracks
  if(qaan)
  {
    maker->AddAnalysis(ConfigureQAAnalysis(calorimeter,simu,debug,print),n++);
    maker->AddAnalysis(ConfigureChargedAnalysis(simu,debug), n++); // charged tracks plots
  }

  maker->SetAnaDebug(debug)  ;
  maker->SwitchOnHistogramsMaker()  ;
  maker->SwitchOnAODsMaker()  ;

  if(print) maker->Print("");

  maker->SwitchOffDataControlHistograms();

  printf("<< End Configuration of %d analysis for calorimeter %s >>\n",n, calorimeter.Data());

  // Create task

  AliAnalysisTaskCaloTrackCorrelation * task = new AliAnalysisTaskCaloTrackCorrelation (Form("CaloTrackCorr%s",kAnaIsoPhotonName.Data()));
  task->SetConfigFileName(""); //Don't configure the analysis via configuration file.
  task->SetDebugLevel(debug);
  task->SetBranches("ESD:AliESDRun.,AliESDHeader");
  task->SetAnalysisMaker(maker);
  mgr->AddTask(task);

  //Create containers

  TString outputfile = AliAnalysisManager::GetCommonFileName();

  AliAnalysisDataContainer *cout_pc   = mgr->CreateContainer(kAnaIsoPhotonName, TList::Class(),
                                                             AliAnalysisManager::kOutputContainer,
                                                             Form("%s",outputfile.Data()));

  AliAnalysisDataContainer *cout_cuts = mgr->CreateContainer(Form("Param_%s",kAnaIsoPhotonName.Data()), TList::Class(),
                                                             AliAnalysisManager::kParamContainer,
                                                             "AnalysisParameters.root");

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  mgr->ConnectInput  (task, 0, mgr->GetCommonInputContainer());
  // AOD output slot will be used in a different way in future
  mgr->ConnectOutput (task, 1, cout_pc);
  mgr->ConnectOutput (task, 2, cout_cuts);

  // Set the trigger selection
  UInt_t mask =  SetTriggerMaskFromName(trigger);
  task->SelectCollisionCandidates(mask);

  return task;
}

///
/// Configure the class handling the events and cluster/tracks filtering.
///
AliCaloTrackReader * ConfigureReader(TString inputDataType = "AOD", Bool_t useKinematics = kFALSE, Bool_t simu = kFALSE,
                                     TString calorimeter = "EMCAL", Bool_t nonlin = kTRUE, Bool_t timecut = kFALSE,
                                     Bool_t primvtx = kFALSE, Bool_t notrackcut = kFALSE, Float_t tmin, Float_t tmax,
                                     Bool_t trackTcut = kFALSE, Float_t minCen = -1, Float_t maxCen = -1,
                                     Int_t debug = -1, Bool_t print = kFALSE, Int_t SSsmearing = 0,TString clustListName ="")
{
  if(simu)
  {
    if (!useKinematics && inputDataType=="AOD") useKinematics = kTRUE; //AOD primary should be available ...
  }

  cout<<"********* ACCESS KINE? "<<useKinematics<< endl;

  AliCaloTrackReader * reader = 0;
  if     (inputDataType=="AOD") reader = new AliCaloTrackAODReader();
  else if(inputDataType=="ESD") reader = new AliCaloTrackESDReader();
  else  printf("AliCaloTrackReader::ConfigureReader() - Data not known inputData=%s\n",inputDataType.Data());

  reader->SetDebug(debug);//10 for lots of messages

  reader->SwitchOffTriggerPatchMatching();
  reader->SwitchOffBadTriggerEventsRemoval();

  reader->SwitchOffWriteDeltaAOD()  ;

  if(SSsmearing != 0)
  {
    reader->SwitchOnShowerShapeSmearing();
    if(SSsmearing == 1) //Gustavo's settings
    { 
      reader->SetSmearingFunction(AliCaloTrackReader::kSmearingLandau);
      reader->SetShowerShapeSmearWidth(0.005);
    }
    else if(SSsmearing == 2) //Astrid's settings
    { 
      reader->SetSmearingFunction(AliCaloTrackReader::kSmearingLandauShift);
      reader->SetShowerShapeSmearWidth(0.035);
    }
  }

  //------------------------
  // Detector input filling
  //------------------------

  if(clustListName!="")
    reader->SetEMCALClusterListName(clustListName);
  //Min cluster/track E
  reader->SetEMCALEMin(0.3);
  reader->SetEMCALEMax(1000);
  reader->SetPHOSEMin(0.3);
  reader->SetPHOSEMax(1000);
  reader->SetCTSPtMin(0.2);
  reader->SetCTSPtMax(1000);

  // Time cuts
  if(simu)
  {
    reader->SwitchOffUseTrackTimeCut();
    reader->SwitchOffUseParametrizedTimeCut();
    reader->SwitchOffUseEMCALTimeCut();
    reader->SetEMCALTimeCut(-1e10,1e10); // Open time cut
  }
  else
  {
    reader->SwitchOffUseParametrizedTimeCut();

    if(timecut)
    {
      printf("Set time cut \n");
      reader->SwitchOnUseEMCALTimeCut();
      //Absolute window
      reader->SetEMCALTimeCut(tmin,tmax); // default is -25ns-20ns
    }
    else
    {
      printf("Off time cuts time cut \n");
      reader->SwitchOffUseEMCALTimeCut();
      //Absolute window
      reader->SetEMCALTimeCut(-1.e6,1.e6);
    }
  }

  reader->SwitchOffFiducialCut();
  reader->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ;

  // Tracks
  reader->SwitchOnCTS();


  reader->SwitchOffRecalculateVertexBC();
  reader->SwitchOffVertexBCEventSelection();

if(trackTcut)
  reader->SwitchOnUseTrackTimeCut();
else
  reader->SwitchOffUseTrackTimeCut();

  reader->SetTrackTimeCut(0,50);

  reader->SwitchOffUseTrackDCACut();
  //reader->SetTrackDCACut(0,0.0105);
  //reader->SetTrackDCACut(1,0.035);
  //reader->SetTrackDCACut(2,1.1);

  if(inputDataType=="ESD")
  {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/macros/CreateTrackCutsPWGJE.C");
    //AliESDtrackCuts * esdTrackCuts = CreateTrackCutsPWGJE(10041004);
    //reader->SetTrackCuts(esdTrackCuts);
    //reader->SwitchOnConstrainTrackToVertex();

//    if(kYears>2010)
//    {
      //Hybrids 2011
      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001008);
      reader->SetTrackCuts(esdTrackCuts);
      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10011008);
      reader->SetTrackComplementaryCuts(esdTrackCuts2);
//    }
//    else
//    {
//      //Hybrids 2010
//      AliESDtrackCuts * esdTrackCuts  = CreateTrackCutsPWGJE(10001006);
//      reader->SetTrackCuts(esdTrackCuts);
//      AliESDtrackCuts * esdTrackCuts2 = CreateTrackCutsPWGJE(10041006);
//      reader->SetTrackComplementaryCuts(esdTrackCuts2);
//    }
  }
  else if(inputDataType=="AOD")
  {
    //reader->SetTrackFilterMask(128);           // Filter bit, not mask, use if off hybrid
    reader->SwitchOnAODHybridTrackSelection(); // Check that the AODs have Hybrids!!!!
    reader->SetTrackStatus(AliVTrack::kITSrefit);
    //reader->SwitchOnTrackHitSPDSelection();    // Check that the track has at least a hit on the SPD, not much sense to use for hybrid or TPC only tracks
  }

  // Calorimeter

  reader->SwitchOffClusterRecalculation();


  // CAREFUL
  if(!nonlin) reader->SwitchOffClusterELinearityCorrection();
  else        reader->SwitchOnClusterELinearityCorrection();

  if(calorimeter == "EMCAL") {
    reader->SwitchOnEMCALCells();
    reader->SwitchOnEMCAL();
  }
  if(calorimeter == "PHOS") { // Should be on if QA is activated with correlation on
    reader->SwitchOffPHOSCells();
    reader->SwitchOffPHOS();
  }

  //-----------------
  // Event selection
  //-----------------

  //if(!kUseKinematics) reader->SetFiredTriggerClassName("CEMC7EGA-B-NOPF-CENTNOTRD"); // L1 Gamma

  // reader->RejectFastClusterEvents() ;

  reader->SwitchOnEventTriggerAtSE();

  reader->SetZvertexCut(10.);               // Open cut
  if(primvtx)
  reader->SwitchOnPrimaryVertexSelection(); // and besides primary vertex
  else
  reader->SwitchOffPrimaryVertexSelection();
  if(notrackcut)
  reader->SwitchOnRejectNoTrackEvents();
  else
  reader->SwitchOffRejectNoTrackEvents();

  reader->SwitchOffPileUpEventRejection();   // remove pileup
  reader->SwitchOffV0ANDSelection() ;        // and besides v0 AND

  if(maxCen > 0 )
  {
    // Centrality
    reader->SetCentralityClass("V0M");
    reader->SetCentralityOpt(100);  // 10 (c= 0-10, 10-20 ...), 20  (c= 0-5, 5-10 ...) or 100 (c= 1, 2, 3 ..)
    reader->SetCentralityBin(minCen,maxCen); // Accept all events, if not select range

    // Event plane (only used in Maker and mixing for AliAnaPi0/AliAnaHadronCorrelation for the moment)
    reader->SetEventPlaneMethod("V0");
  }

  if(print) reader->Print("");

  return reader;
}

///
/// Configure the class handling the calorimeter clusters specific methods
///
AliCalorimeterUtils* ConfigureCaloUtils(Bool_t nonlin = kTRUE, Bool_t exotic = kTRUE ,Bool_t simu = kFALSE, Bool_t timecut = kFALSE, Int_t debug = -1, Bool_t print = kFALSE)
{
  AliCalorimeterUtils *cu = new AliCalorimeterUtils;
  cu->SetDebug(debug);

  // Remove clusters close to borders, at least max energy cell is 1 cell away
  cu->SetNumberOfCellsFromEMCALBorder(0);//this was originally set to one
  cu->SetNumberOfCellsFromPHOSBorder(2);

  cu->SwitchOffRecalculateClusterTrackMatching();

  cu->SwitchOffBadChannelsRemoval() ;

  //EMCAL settings

  if(simu)
    cu->SwitchOnLoadOwnEMCALGeometryMatrices();

  /*  AliEMCALRecoUtils * recou = cu->GetEMCALRecoUtils();

  cu->SwitchOffRecalibration();
  cu->SwitchOffRunDepCorrection();

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  ConfigureEMCALRecoUtils(recou,
                          simu,
                          exotic,
                          nonlin,
                          kFALSE, // e calib
                          kFALSE, // bad map
                          kFALSE); // time calib
  if(timecut) recou->SetExoticCellDiffTimeCut(50.);
  */
  if( nonlin)
  {
    printf("ConfigureCaloUtils() - Apply non linearity to EMCAL\n");
    cu->SwitchOnCorrectClusterLinearity();
  }
  /*
  printf("ConfigureCaloUtils() - EMCAL Recalibration ON? %d %d\n",recou->IsRecalibrationOn(), cu->IsRecalibrationOn());
  printf("ConfigureCaloUtils() - EMCAL BadMap        ON? %d %d\n",recou->IsBadChannelsRemovalSwitchedOn(), cu->IsBadChannelsRemovalSwitchedOn());
  */
  cu->SetNumberOfSuperModulesUsed(10);

  if(print) cu->Print("");

  return cu;
}

///
/// Configure the task doing the first photon cluster selections
/// Basically the track matching, minor shower shape cut, NLM selection ...
///
AliAnaPhoton* ConfigurePhotonAnalysis(TString calorimeter = "EMCAL", Bool_t tm = kFALSE, Float_t deltaphicut = 0.02, Float_t deltaetacut = 0.03,Int_t disttobad=0,Int_t nlmMax = 2, Bool_t simu = kFALSE, Int_t debug = -1, Bool_t print = kFALSE)
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
    //    ana->SetTimeCut(-1e10,1e10); // open cut, usual time window of [425-825] ns if time recalibration is off
    // restrict to less than 100 ns when time calibration is on
    ana->SetMinDistanceToBadChannel(disttobad, 4, 6);

    // NLM cut, used in all, exclude clusters with more than 2 maxima
    // Not needed if M02 cut is already strong or clusterizer V2
    ana->SetNLMCut(1, nlmMax) ;
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
  caloPID->SetEMCALLambda0CutMax(1000.);
  caloPID->SetEMCALLambda0CutMin(0.);

  // caloPID->SetEMCALDEtaCut(0.025);
  // caloPID->SetEMCALDPhiCut(0.030);
  caloPID->SetEMCALDEtaCut(deltaetacut);
  caloPID->SetEMCALDPhiCut(deltaphicut);

  ana->SwitchOnFillShowerShapeHistograms();  // Filled before photon shower shape selection
  if(!simu) ana->SwitchOnFillPileUpHistograms();

  // Input / output delta AOD settings

  ana->SetOutputAODName(Form("Photon%s",kAnaIsoPhotonName.Data()));
  ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");

  //Set Histograms name tag, bins and ranges

  ana->AddToHistogramsName(Form("AnaPhoton_TM%d_",tm));
  SetHistoRangeAndNBins(ana->GetHistogramRanges(), calorimeter); // see method below

  // Number of particle type MC histograms
  ana->FillNOriginHistograms(20);
  ana->FillNPrimaryHistograms(20);

  ConfigureMC(ana,simu);

  if(print) ana->Print("");

  return ana;
}

///
/// Configure the task doing the trigger cluster/random isolation
///
AliAnaParticleIsolation* ConfigureIsolationAnalysis(TString calorimeter = "EMCAL",
                                                    TString particle="Photon",
                                                    Int_t  partInCone = AliIsolationCut::kOnlyCharged,
                                                    Int_t  thresType  = AliIsolationCut::kSumPtFracIC,
                                                    Float_t cone = 0.3,
                                                    Float_t pth  = 0.3,
                                                    Bool_t tm = kFALSE,
                                                    Bool_t leading = kTRUE,
                                                    Bool_t multi = kFALSE, Bool_t simu = kFALSE,
                                                    Int_t debug = -1, Bool_t print = kFALSE,
                                                    Bool_t tmInCone = kTRUE )
{
  AliAnaParticleIsolation *ana = new AliAnaParticleIsolation();
  ana->SetDebug(debug);

  ana->SwitchOnFiducialCut();
  //Avoid borders of EMCal
  if(calorimeter=="EMCAL")
  {
    //ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.6, 86, 174) ;
    ana->GetFiducialCut()->SetSimpleEMCALFiducialCut(0.27, 103, 157) ;
  }

  ana->SetMinPt(5);

  // Input / output delta AOD settings

  ana->SetInputAODName(Form("%s%s",particle.Data(),kAnaIsoPhotonName.Data()));
  ana->SetAODObjArrayName(Form("IC%s_%s",particle.Data(),kAnaIsoPhotonName.Data()));

  ana->SetCalorimeter(calorimeter);

  if(!tm)  ana->SwitchOnTMHistoFill();
  else     ana->SwitchOffTMHistoFill();
  //   ana->SwitchOnTMHistoFill();

  // ana->SwitchOffSSHistoFill();
  // if(!kSimulation) ana->SwitchOnFillPileUpHistograms();
   ana->SwitchOnSSHistoFill();
  if(leading) ana->SwitchOnLeadingOnly();
  else ana->SwitchOffLeadingOnly();
  if(!simu) ana->SwitchOnFillPileUpHistograms();

  //Do settings for main isolation cut class
  AliIsolationCut * ic =  ana->GetIsolationCut();
  ic->SetDebug(debug);

  printf("\t *** Set: R = %2.2f, Threshold %2.2f, Method %d, Neutral/Charged option %d ***\n",cone,pth,thresType,partInCone);

  //Main parameters
  //****
  ic->SetConeSize(cone);

  ic->SetPtFraction    (0.1);
  ic->SetPtThreshold   (pth);
  ic->SetSumPtThreshold(pth);

  ic->SetParticleTypeInCone(partInCone);

  ic->SetICMethod(thresType);

  ic->SetTrackMatchedClusterRejectionInCone(tmInCone);
  //****

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
  caloPID->SetEMCALDEtaCut(0.02);
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

  ConfigureMC(ana,simu);

  if(print) ic ->Print("");
  if(print) ana->Print("");

  return ana;
}

///
/// Configure the task doing standard calorimeter QA
///
AliAnaCalorimeterQA* ConfigureQAAnalysis(TString calorimeter = "EMCAL", Bool_t simu = kFALSE, Int_t debug = -1, Bool_t print = kFALSE)
{
  AliAnaCalorimeterQA *ana = new AliAnaCalorimeterQA();
  ana->SetDebug(debug); //10 for lots of messages
  ana->SetCalorimeter(calorimeter);

  ana->SetTimeCut(-1e10,1e10); // Open time cut
  ana->SwitchOffCorrelation();

  // Study exotic clusters PHOS and EMCAL
  ana->SwitchOffStudyBadClusters() ;

  ana->SwitchOffFiducialCut();
  ana->SwitchOffFillAllTH3Histogram();
  ana->SwitchOffFillAllPositionHistogram();
  ana->SwitchOffFillAllPositionHistogram2();
  ana->SwitchOnFillAllTrackMatchingHistogram();
  ana->SwitchOnFillAllCellTimeHisto() ;

  ana->AddToHistogramsName("QA_"); //Begining of histograms name
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),calorimeter); // see method below

  ConfigureMC(ana,simu);

  if(print) ana->Print("");

  return ana;
}

///
/// Configure the task doing charged track selection
///
AliAnaChargedParticles* ConfigureChargedAnalysis(Bool_t simulation, Int_t debugLevel)
{
  AliAnaChargedParticles *ana = new AliAnaChargedParticles();
  ana->SetDebug(debugLevel); //10 for lots of messages

  // selection cuts

  ana->SetMinPt(0.5);
  ana->SwitchOnFiducialCut();
  ana->GetFiducialCut()->SetSimpleCTSFiducialCut(0.8, 0, 360) ; //more restrictive cut in reader and after in isolation

  ana->SwitchOffFillVertexBC0Histograms() ;
  //if(!simulation) ana->SwitchOnFillPileUpHistograms();

  // Input / output delta AOD settings

  ana->SetOutputAODName(Form("Hadron%s",kAnaIsoPhotonName.Data()));
  ana->SetOutputAODClassName("AliAODPWG4Particle");
  ana->SetInputAODName(Form("Hadron%s",kAnaIsoPhotonName.Data()));

  //Set Histograms name tag, bins and ranges

  ana->AddToHistogramsName("AnaHadrons_");
  SetHistoRangeAndNBins(ana->GetHistogramRanges(),""); // see method below

  ana->GetHistogramRanges()->SetHistoPhiRangeAndNBins(0, TMath::TwoPi(), 200) ;
  ana->GetHistogramRanges()->SetHistoEtaRangeAndNBins(-1.5, 1.5, 300) ;

  if(debugLevel > 0) ana->Print("");

  return ana;
}

///
/// Configure the task doing random trigger generation
///
AliAnaRandomTrigger* ConfigureRandomTriggerAnalysis(TString detector = "")
{
  AliAnaRandomTrigger *ana = new AliAnaRandomTrigger();
  ana->SetDebug(kDebug); //10 for lots of messages

  if(detector=="") detector = kCalorimeter;
  ana->SetDetector(detector);

  // selection cuts
  ana->SetMinPt(4.);
  ana->SetMaxPt(61.);

  if     (detector=="EMCAL")
  {
    ana->SetEtaCut(-0.27,0.27);
    //ana->SetPhiCut(103*TMath::DegToRad(), 157*TMath::DegToRad());
    ana->SetPhiCut(1.8, 2.75);
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
    ana->SetOutputAODName(Form("RandomTrigger%s%s",detector.Data(),kAnaIsoPhotonName.Data()));
    ana->SetOutputAODClassName("AliAODPWG4ParticleCorrelation");
  }
  else
    ana->SetInputAODName(Form("RandomTrigger%s%s",detector.Data(),kAnaIsoPhotonName.Data()));

  printf("Set RandomTrigger%s%s\n",detector.Data(),kAnaIsoPhotonName.Data());

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
/// Configure the selection of MC events
///
void ConfigureMC(AliAnaCaloTrackCorrBaseClass* ana, Bool_t simu = kFALSE)
{
  if(simu) ana->SwitchOnDataMC() ;//Access MC stack and fill more histograms, AOD MC not implemented yet.
  else     ana->SwitchOffDataMC() ;

  //Set here generator name, default pythia
  //ana->GetMCAnalysisUtils()->SetMCGenerator("");
}

///
/// Set common histograms binning and ranges
///
void SetHistoRangeAndNBins (AliHistogramRanges* histoRanges, TString calorimeter = "EMCAL")
{
  histoRanges->SetHistoPtRangeAndNBins(0., 100., 200) ; // Energy and pt histograms

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
  histoRanges->SetHistoTimeRangeAndNBins(-1000.,1000,1000);
  histoRanges->SetHistoDiffTimeRangeAndNBins(-200, 200, 800);

  // track-cluster residuals
  histoRanges->SetHistoTrackResidualEtaRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistoTrackResidualPhiRangeAndNBins(-0.15,0.15,300);
  histoRanges->SetHistodRRangeAndNBins(0.,0.15,150);//QA

  // QA, electron, charged
  histoRanges->SetHistoPOverERangeAndNBins(0,  2.5 ,500);
  histoRanges->SetHistodEdxRangeAndNBins  (0.,250.0,500);

  // QA
  histoRanges->SetHistoFinePtRangeAndNBins(0, 10, 200) ; // bining for fhAmpId
  histoRanges->SetHistoRatioRangeAndNBins(0.,2.,100);
  histoRanges->SetHistoVertexDistRangeAndNBins(0.,500.,500);
  histoRanges->SetHistoNClusterCellRangeAndNBins(0,50,50);
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
/// Configure the task doing the trigger cluster-jet correlation
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
    return (AliVEvent::kCentral  | AliVEvent::kMB);
  }
  else if(trigger=="CentralEGA")
  {
    printf("CaloTrackCorr trigger Central+EMCEGA\n");
    return (AliVEvent::kCentral | AliVEvent::kEMCEGA);
  }
  else if(trigger=="SemiCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral\n");
    return (AliVEvent::kSemiCentral | AliVEvent::kMB);
  }
  else if(trigger=="SemiOrCentral")
  {
    printf("CaloTrackCorr trigger SemiCentral Or Central\n");
    return (AliVEvent::kSemiCentral | AliVEvent::kCentral  | AliVEvent::kMB);
  }
}


