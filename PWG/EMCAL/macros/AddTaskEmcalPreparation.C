AliAnalysisTaskSE *AddTaskEmcalPreparation(
    const char *perstr  = "LHC11h",
    UInt_t clusterizer  = AliEMCALRecParam::kClusterizerv2,
    const char *pass    = 0 /*should not be needed; will be recovered from path of AOD/ESD; no need to specify by user; */
) {

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    Error("AddTaskEmcalPreparation","No analysis manager found.");
    return NULL;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();
  if (!evhand) {
    Error("AddTaskEmcalPreparation", "This task requires an input event handler");
    return NULL;
  }

  TString period(perstr);
  period.ToLower();

  Bool_t isMC = kFALSE;
  if(period.Sizeof()>7) isMC = kTRUE;
  //  Bool_t isMC = (mgr->GetMCtruthEventHandler() != NULL); //only works for ESD
  if(isMC) Printf("AddTaskEmcalPreparation: Running on MC");
  else     Printf("AddTaskEmcalPreparation: Running on DATA");

  //----------------------- Add tender -------------------------------------------------------
  Bool_t distBC         = kFALSE; //switch for recalculation cluster position from bad channel 
  Bool_t recalibClus    = kFALSE;
  Bool_t recalcClusPos  = kFALSE;
  Bool_t nonLinearCorr  = kFALSE;
  Bool_t remExoticCell  = kFALSE;
  Bool_t remExoticClus  = kFALSE;
  Bool_t fidRegion      = kFALSE;
  Bool_t calibEnergy    = kTRUE;
  Bool_t calibTime      = kTRUE;
  Bool_t remBC          = kTRUE;
  UInt_t nonLinFunct    = 0;
  Bool_t reclusterize   = kFALSE;
  Float_t seedthresh    = 0.1;      // 100 MeV
  Float_t cellthresh    = 0.05;     // 50 MeV 
  UInt_t clusterizerT   = 0;
  Bool_t trackMatch     = kFALSE;
  Bool_t updateCellOnly = kFALSE;
  Float_t timeMin       = -50e-6;   // minimum time of physical signal in a cell/digit
  Float_t timeMax       =  50e-6;   // maximum time of physical signal in a cell/digit
  Float_t timeCut       = 1e6;
  if(period.Contains("lhc11h")) {
    timeMin = -50e-9;
    timeMax = 100e-9;
  }
  if(isMC) {
    calibEnergy = kFALSE;
    calibTime   = kFALSE;
    timeMin = -1.;
    timeMax = 1e6;
  }
#ifdef __CLING__
  // ROOT6 version of the Config macro. JIT cannot handle load and execute macro (compiler error) - need to call via gROOT->ProcessLine(...)
  std::stringstream emcaltenderadd;
  emcaltenderadd << ".x " << gSystem->Getenv("ALICE_PHYSICS") <<"/PWG/EMCAL/macros/AddTaskEMCALTender.C(";
  emcaltenderadd << (distBC ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (recalibClus ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (recalcClusPos ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (nonLinearCorr ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (remExoticCell ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (remExoticClus ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (fidRegion ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (calibEnergy ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (calibTime ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (remBC ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << nonLinFunct << ", ";
  emcaltenderadd << (reclusterize ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << seedthresh << ", ";
  emcaltenderadd << cellthresh << ", ";
  emcaltenderadd << clusterizerT << ", ";
  emcaltenderadd << (trackMatch ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << (updateCellOnly ? "kTRUE" : "kFALSE") << ", ";
  emcaltenderadd << timeMin << ", ";
  emcaltenderadd << timeMax << ", ";
  emcaltenderadd << timeCut;// << ", ";
  //emcaltenderadd << ((pass != NULL) ? pass : 0);
  emcaltenderadd << ")";
  std::string emcaltenderaddstring = emcaltenderadd.str();
  std::cout << "Calling Add macro using command string " << emcaltenderaddstring << std::endl;
  //gROOT->SetMacroPath(Form("%s", Form("%s/PWG/EMCAL/macros", gSystem->Getenv("ALICE_PHYSICS"))));
  AliAnalysisTaskSE *tender = (AliAnalysisTaskSE *)gROOT->ProcessLine(emcaltenderaddstring.c_str());
#else
  // ROOT5 version, allows loading a macro
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");//tendertasks
  AliAnalysisTaskSE *tender = AddTaskEMCALTender(distBC, recalibClus, recalcClusPos, nonLinearCorr, remExoticCell, remExoticClus,
      fidRegion, calibEnergy, calibTime, remBC, nonLinFunct, reclusterize, seedthresh,
      cellthresh, clusterizerT, trackMatch, updateCellOnly, timeMin, timeMax, timeCut,pass);
#endif

  //----------------------- Add clusterizer -------------------------------------------------------
  remExoticCell  = kFALSE;
#ifdef __CLING__
  // ROOT6 version of the Config macro. JIT cannot handle load and execute macro (compiler error) - need to call via gROOT->ProcessLine(...)
  std::stringstream clusterizeradd;
  clusterizeradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskClusterizerFast.C(";
  clusterizeradd << "\"ClusterizerFast\", ";
  clusterizeradd << "\"\", ";
  clusterizeradd << "\"\", ";
  clusterizeradd << clusterizer << ", ";
  clusterizeradd << cellthresh << ", ";
  clusterizeradd << seedthresh << ", ";
  clusterizeradd << timeMin << ", ";
  clusterizeradd << timeMax << ", ";
  clusterizeradd << timeCut << ", ";
  clusterizeradd << (remExoticCell ? "kTRUE" : "kFALSE") << ", ";
  clusterizeradd << (distBC ? "kTRUE" : "kFALSE") << ", ";
  clusterizeradd << AliAnalysisTaskEMCALClusterizeFast::kFEEData;
  clusterizeradd << ")";
  std::string clusterizeraddstring = clusterizeradd.str();
  std::cout << "Calling Add macro using command string " << clusterizeraddstring << std::endl;
  AliAnalysisTaskEMCALClusterizeFast *clusterizerTask = (AliAnalysisTaskEMCALClusterizeFast *)gROOT->ProcessLine(clusterizeraddstring.c_str());
#else
  // ROOT5 version, allows loading a macro
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
  AliAnalysisTaskEMCALClusterizeFast *clusterizerTask = AddTaskClusterizerFast("ClusterizerFast","","",clusterizer,cellthresh,seedthresh,
      timeMin,timeMax,timeCut,remExoticCell,distBC,
      AliAnalysisTaskEMCALClusterizeFast::kFEEData);
#endif

  if(isMC){
    if((period.Contains("lhc10") || period.Contains("lhc11") || period.Contains("lhc12a") 
        || period.Contains("lhc12b") || period.Contains("lhc12c")
        || period.Contains("lhc12d") || period.Contains("lhc12e"))
        && !(period.Contains("lhc12a15d") || period.Contains("lhc12a15f")
            || period.Contains("lhc12a15g") || period.Contains("lhc12a15h")
            || period.Contains("lhc12c4") || period.Contains("lhc12d1")
            || period.Contains("lhc12d2") || period.Contains("lhc12d3")
            || period.Contains("lhc12e1") || period.Contains("lhc12e2") || period.Contains("lhc12e3")))
    {
      clusterizerTask->SetCellMCLabelFromCluster(1);
    }
  }

  //----------------------- Add cluster maker -----------------------------------------------------
  nonLinFunct = AliEMCALRecoUtils::kBeamTestCorrected;
  if(isMC) {
    if(period == "lhc12a15a") 
      nonLinFunct = AliEMCALRecoUtils::kPi0MCv2;
    else
      nonLinFunct = AliEMCALRecoUtils::kPi0MCv3;
  }
  remExoticClus  = kTRUE;
#ifdef __CLING__
  // ROOT6 version of the Config macro. JIT cannot handle load and execute macro (compiler error) - need to call via gROOT->ProcessLine(...)
  std::stringstream clustermakeradd;
  clustermakeradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C(";
  clustermakeradd << nonLinFunct << ", ";
  clustermakeradd << (remExoticClus ? "kTRUE" : "kFALSE") << ", ";
  clustermakeradd << "0, ";
  clustermakeradd << "\"EmcCaloClusters\", ";
  clustermakeradd << "0., ";
  clustermakeradd << "kTRUE";
  clustermakeradd << ")";
  std::string clustermakeraddstring = clustermakeradd.str();
  std::cout << "Calling Add macro using command string " << clustermakeraddstring << std::endl;
  AliEmcalClusterMaker *clusMaker = (AliEmcalClusterMaker *)gROOT->ProcessLine(clustermakeraddstring.c_str());
#else
  // ROOT5 version, allows loading a macro
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); //cluster maker: non-linearity,
  AliEmcalClusterMaker *clusMaker = AddTaskEmcalClusterMaker(nonLinFunct,remExoticClus,0,"EmcCaloClusters",0.,kTRUE);
#endif
  clusMaker->GetClusterContainer(0)->SetClusPtCut(0.);
  clusMaker->GetClusterContainer(0)->SetClusECut(0.);

  return clusMaker;

}
