///
/// EMCal Tender configuration macro
///
AliAnalysisTaskSE *AddTaskEMCALTender(
  Bool_t distBC               = kTRUE,    // distance to bad channel
  Bool_t recalibClus          = kTRUE,    // recalibrate cluster energy
  Bool_t recalcClusPos        = kTRUE,    // recalculate cluster position
  Bool_t nonLinearCorr        = kTRUE,    // apply non-linearity
  Bool_t remExoticCell        = kTRUE,    // remove exotic cells
  Bool_t remExoticClus        = kTRUE,    // remove exotic clusters
  Bool_t fidRegion            = kFALSE,   // apply fiducial cuts
  Bool_t calibEnergy          = kTRUE,    // calibrate energy
  Bool_t calibTime            = kTRUE,    // calibrate timing
  Bool_t remBC                = kTRUE,    // remove bad channels
  UInt_t nonLinFunct          = AliEMCALRecoUtils::kBeamTestCorrectedv3, // For MC use kPi0MCv3
  Bool_t reclusterize         = kTRUE,    // reclusterize
  Float_t seedthresh          = 0.100,    // seed threshold
  Float_t cellthresh          = 0.050,    // cell threshold
  UInt_t clusterizer          = AliEMCALRecParam::kClusterizerv2,
  Bool_t trackMatch           = kTRUE,    // track matching
  Bool_t updateCellOnly       = kFALSE,   // only change if you run your own clusterizer task
  Float_t timeMin             =-1000e-9,  // minimum time of physical signal in a cell/digit (s)
  Float_t timeMax             = 1000e-9,  // maximum time of physical signal in a cell/digit (s)
  Float_t timeCut             = 1000e-9,  // maximum time difference between the digits inside EMC cluster (s)
  const char *pass            = 0,        // string defining pass (use none if figured out from path)
  Bool_t  remapMcAod          = kFALSE,   // switch on the remaping for the MC labels in AOD productions,
  TString cdbStorage          = "raw://", // "local://"
  Float_t diffEAggregation    = 0.,       // difference E in aggregation of cells (i.e. stop aggregation if E_{new} > E_{prev} + diffEAggregation)
  Bool_t enableFracEMCRecalc  = 0,        // enables the recalculation of the MC labels including the fractional eneryg on cell level
  Int_t  removeNMCGenerators  = 0,        // set number of accepted MC generators input (only for enableFracEMCRecalc=1)
  Bool_t enableMCGenRemovTrack= 1,        // apply the MC generators rejection also for track matching  
  TString removeMCGen1        = "",       // name of generator input to be accepted
  TString removeMCGen2        = "",       // name of generator input to be accepted
  TString customBCmap         = "",       // location of custom bad channel map (full path including file)
  Bool_t useNewRWTempCalib    = kFALSE,   // switch for usage of new temperature calib parameters (available for run1 and run2)
  TString customSMtemps       = "",       // location of custom SM-wise temperature OADB file (full path including file)
  TString customTempParams    = "",        // location of custom temperature calibration parameters OADB file (full path including file)
  Bool_t useOneHistAllBCS     = kFALSE,    // flag to use on histogram for the all the BCs
  Bool_t load1DBCmap          = kFALSE     // Flag to load 1D bad channel map
) 
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALTender", "No analysis manager to connect to.");
    return NULL;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();

  // Create the task and configure it.
  //===========================================================================

  AliAnalysisTaskSE *ana = 0;
  AliAnalysisDataContainer *coutput1 = 0;
  AliEMCALTenderSupply *EMCALSupply = 0;
  
  #ifdef __CLING__
    // ROOT6 version of the Config macro. JIT cannot handle load and execute macro (compiler error) - need to call via gROOT->ProcessLine(...)
    std::stringstream configbuilder;
    configbuilder << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/ConfigEmcalTenderSupply.C(";
    configbuilder << (distBC ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (recalibClus ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (recalcClusPos ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (nonLinearCorr ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (remExoticCell ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (remExoticClus ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (fidRegion ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (calibEnergy ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (calibTime ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << (remBC ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << nonLinFunct << ", ";
    configbuilder << (reclusterize ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << seedthresh << ", ";
    configbuilder << cellthresh << ", ";
    configbuilder << clusterizer << ", ";
    configbuilder << trackMatch << ", ";
    configbuilder << (updateCellOnly ? "kTRUE" : "kFALSE") << ", ";
    configbuilder << timeMin << ", ";
    configbuilder << timeMax << ", ";
    configbuilder << timeCut << ", ";
    configbuilder << diffEAggregation;
    configbuilder << ")";
    std::string configbuilderstring = configbuilder.str();
    std::cout << "Running config macro " << configbuilderstring << std::endl;
    EMCALSupply = (AliEMCALTenderSupply *)gROOT->ProcessLine(configbuilderstring.c_str());
  #else
    // ROOT5 version, allows loading a macro
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/ConfigEmcalTenderSupply.C");

  EMCALSupply = ConfigEmcalTenderSupply(distBC, recalibClus, recalcClusPos, nonLinearCorr, remExoticCell, remExoticClus, 
                                        fidRegion, calibEnergy, calibTime, remBC, nonLinFunct, reclusterize, seedthresh, 
                                        cellthresh, clusterizer, trackMatch, updateCellOnly, timeMin, timeMax, timeCut, diffEAggregation);
  #endif

  EMCALSupply->SwitchUseMergedBCs(useOneHistAllBCS);

  if(load1DBCmap) EMCALSupply->Load1DBadChannelMap();

  if (pass) 
    EMCALSupply->SetPass(pass);
  if (customBCmap!="")
    EMCALSupply->SetCustomBC(customBCmap);
  if (useNewRWTempCalib)
    EMCALSupply->SwitchUseNewRunDepTempCalib(useNewRWTempCalib);
  if(customSMtemps!="" && customTempParams!="")
    EMCALSupply->SetCustomTimeCalibration(customSMtemps,customTempParams);
  if (evhand->InheritsFrom("AliESDInputHandler")) {
    #ifdef __CLING__
        AliTender* alitender = dynamic_cast<AliTender *>(mgr->GetTopTasks()->FindObject("AliTender"));
    #else
        AliTender* alitender = (AliTender *)mgr->GetTopTasks()->FindObject("AliTender");
    #endif
    if (!alitender)
      alitender = new  AliTender("AliTender");
    
    alitender->AddSupply(EMCALSupply);
    alitender->SetDefaultCDBStorage(cdbStorage); 
    
    ana = alitender;
    
    coutput1 = mgr->CreateContainer("emcal_tender_event", 
                                      AliESDEvent::Class(), 
                                      AliAnalysisManager::kExchangeContainer, 
                                      "default_tender");
  } else if (evhand->InheritsFrom("AliAODInputHandler")) {
    #ifdef __CLING__
        AliEmcalTenderTask* emcaltender = dynamic_cast<AliEmcalTenderTask *>(mgr->GetTopTasks()->FindObject("AliEmcalTenderTask"));
    #else
        AliEmcalTenderTask* emcaltender = (AliEmcalTenderTask *)mgr->GetTopTasks()->FindObject("AliEmcalTenderTask");
    #endif
      
    if (!emcaltender)
        emcaltender = new  AliEmcalTenderTask("AliEmcalTenderTask");
    
    if (remapMcAod)
        EMCALSupply->SwitchOnRemapMCLabelForAODs();
    
     emcaltender->SetEMCALTenderSupply(EMCALSupply);
    
    ana = emcaltender;
    
    coutput1 = mgr->CreateContainer(  "emcal_tender_event",
                                      AliAODEvent::Class(), 
                                      AliAnalysisManager::kExchangeContainer, 
                                      "default_tender");
  } else {
    ::Error("AddTaskEMCALTender", "Input event handler not recognized, AOD/ESD expected. Returning...");
    return NULL;
  }

  if (enableFracEMCRecalc)
  {
    EMCALSupply->SwitchOnUseMCEdepFracLabelForCell();
    
    printf("Enable frac MC recalc, remove generators %d \n",removeNMCGenerators);
    if(removeNMCGenerators > 0)
    {
      printf("\t gen1 <%s>, gen2 <%s>, remove tracks %d\n",removeMCGen1.Data(),removeMCGen2.Data(),enableMCGenRemovTrack);
      AliEMCALRecoUtils* ru = new AliEMCALRecoUtils;
      ru->SetNumberOfMCGeneratorsToAccept(removeNMCGenerators) ;
      ru->SetNameOfMCGeneratorsToAccept(0,removeMCGen1);
      ru->SetNameOfMCGeneratorsToAccept(1,removeMCGen2);
      
      if(!enableMCGenRemovTrack) ru->SwitchOffMCGeneratorToAcceptForTrackMatching() ;
      
      EMCALSupply->SetRecoUtils(ru);
    }
  }

  mgr->AddTask(ana);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  mgr->ConnectInput(ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(ana, 1, coutput1 );

  ::Info("AddTaskEMCALTender", "Tender configured");
   
  return ana;
}
