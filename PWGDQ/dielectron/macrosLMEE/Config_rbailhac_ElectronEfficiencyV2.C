TString names("cutTPC3sigma;cutTPChr;cutTOF;cutpid;cutTPC3sigma_itss;cutTPChr_itss;cutTOF_itss;cutpid_itss");
bool DoPairing = kTRUE;
bool DoULSLS = kTRUE;

bool GetResolutionFromAlien = kTRUE;
std::string resoFilename = "pp_resolution_deltaXvsX_Oton.root";
std::string resoFilenameFromAlien = "/alice/cern.ch/user/r/rbailhac/supportFiles/pp_resolution_deltaXvsX_Oton.root";

bool GetCentralityFromAlien = kFALSE;
std::string centralityFilename = "";
std::string centralityFilenameFromAlien = "";

const Int_t triggerNames = AliVEvent::kINT7;

const int nMCSignal = 0;
const int nCutsetting = 0;

const double centMin = -1.;
const double centMax = -1.;

const double minGenPt = 0.1;
const double maxGenPt = 100;
const double minGenEta = -1.5;
const double maxGenEta =  1.5;


// fiducial cuts
const double minPtCut = 0.2;
const double maxPtCut = 15.0;
const double minEtaCut = -0.8;
const double maxEtaCut = 0.8;



// binning of single leg histograms
bool usePtVector = true;
const Double_t ptBins[] = {
  0.00,0.10,0.11,0.12,0.13,0.14,0.15,0.155,0.16,0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.20,0.205,0.21,0.215,0.22,0.225,0.23,0.235,0.24,0.245,0.25,0.255,
  0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,0.30,0.32,0.34,0.36,0.38,0.40,0.43,0.46,
  0.49,0.52,0.55,0.60,0.65,0.70,0.75,0.80,0.90,1.00,1.10,1.20,
  1.40,1.60,1.80,2.00,2.40,2.80,3.20,3.70,4.50,6.00,8.00,12.0
};
const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;

const double minPtBin = 0;
const double maxPtBin = 8;
const int    stepsPtBin = 800;

const double minEtaBin = -1.0;
const double maxEtaBin =  1.0;
const int    stepsEtaBin = 40;

const double minPhiBin = 0;
const double maxPhiBin =  TMath::TwoPi();
const int    stepsPhiBin = 90;

const double minThetaBin = 0;
const double maxThetaBin =  TMath::TwoPi();
const int    stepsThetaBin = 60;

const double minMassBin = 0;
const double maxMassBin =  4;
const int    stepsMassBin = 400;
const double minPairPtBin = 0;
const double maxPairPtBin =  10;
const int    stepsPairPtBin = 400;

// Binning of resolution histograms
const int    NbinsDeltaMom    = 1200;
const double DeltaMomMin   =-10.0;
const double DeltaMomMax   =  2.0;
const int    NbinsRelMom      = 400;
const double RelMomMin     =  0.0;
const double RelMomMax     =  2.0;
const int    NbinsDeltaEta    = 200;
const double DeltaEtaMin   = -0.4;
const double DeltaEtaMax   =  0.4;
const int    NbinsDeltaTheta  = 200;
const double DeltaThetaMin = -0.4;
const double DeltaThetaMax =  0.4;
const int    NbinsDeltaPhi    = 200;
const double DeltaPhiMin   = -0.4;
const double DeltaPhiMax   =  0.4;

AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition);
AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t cutDefinition, Bool_t isAOD);
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition);
AliAnalysisCuts* SetupEventCuts(Bool_t isAOD);
void AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task);
void AddPairMCSignal(AliAnalysisTaskElectronEfficiencyV2* task);
void SetEtaCorrectionTOFMean(AliAnalysisTaskElectronEfficiencyV2 *task, Int_t corrXdim, Int_t corrYdim);
void SetEtaCorrectionTOFRMS(AliAnalysisTaskElectronEfficiencyV2 *task, Int_t corrXdim, Int_t corrYdim); 

// #########################################################
// #########################################################

//________________________________________________________________
AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t cutDefinition, Bool_t isAOD)
{
  std::cout << "SetupTrackCutsAndSettings()" <<std::endl;


  AliAnalysisFilter *anaFilter = new AliAnalysisFilter(Form("anaFilter_%d",cutDefinition),Form("anaFilter_%d",cutDefinition)); // named constructor seems mandatory!
  // do not change these initial values!
  //Int_t selectedPID=-1;
  //Bool_t isPrefilterCutset=kFALSE;

  AliDielectronV0Cuts *noconv = new AliDielectronV0Cuts("IsGamma","IsGamma");
  // which V0 finder you want to use
  noconv->SetV0finder(AliDielectronV0Cuts::kAll);  // kAll(default), kOffline or kOnTheFly
  // add some pdg codes (they are used then by the KF package and important for gamma conversions)
  noconv->SetPdgCodes(22,11,11); // mother, daughter1 and 2
  // add default PID cuts (defined in AliDielectronPID)
  // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
  //noconv->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
  // add the pair cuts for V0 candidates
  noconv->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kR,                             3.0,  90.00, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kM,                             0.0,   0.10, kFALSE);
  noconv->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  // selection or rejection of V0 tracks
  noconv->SetExcludeTracks(kTRUE); 
  
  anaFilter->AddCuts( SetupTrackCuts(cutDefinition) );
  anaFilter->AddCuts( SetupPIDcuts(cutDefinition) );
  anaFilter->AddCuts( noconv );
  std::cout << "...cuts added!" <<std::endl; 
  
  return anaFilter;
}
//________________________________________________________________
AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition)
{
  std::cout << "SetupTrackCuts()" <<std::endl;
  //AliAnalysisCuts* trackCuts=0x0;

  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA); // 1<<4
  trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
  trackCutsDiel->SetRequireITSRefit(kTRUE);
  trackCutsDiel->SetRequireTPCRefit(kTRUE);


  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  if(cutDefinition!=0 && cutDefinition!=4) {  
    trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.8,   0.8);
    trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.2,   15.);
  }
  else {
    trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -1.1,   1.1);
    trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.1,   100.);
  }
  if(cutDefinition > 3) {
    printf("Add shared cluster cut\n");
    trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS, 1.0, 6.0, kTRUE);
  }
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY, -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,  -3.0,   3.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,      3.0, 999.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,      80.0, 999.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,    -999.,   4.5);
  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,    -999.,   4.0); 
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,    100.0, 999.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,     0.8, 999.);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSFracTPC,     -999., 0.4);
  

  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("Trackcuts","Trackcuts",AliDielectronCutGroup::kCompAND);
  trackCuts->AddCut(trackCutsAOD);
  trackCuts->AddCut(trackCutsDiel);
  

  return trackCuts;

}

//________________________________________________________________
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition)
{
  
  std::cout << "SetupPIDcuts()" <<std::endl;

  AliDielectronPID *mastermind_TPC_3sigma = new AliDielectronPID("mastermind_TPC_3sigma","mastermind_TPC_3sigma");
  AliDielectronPID *mastermind_TPC = new AliDielectronPID("mastermind_TPC","mastermind_TPC");
  AliDielectronPID *mastermind_TOF = new AliDielectronPID("mastermind_TOF","mastermind_TOF");

  // Simple PID
  mastermind_TPC_3sigma->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  
  //TPC electrons: includes electrons and exclude all possible other contributions using the TPC
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4.,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kMuon,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton,    -3. ,3.,0.0, 100., kTRUE, AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  
  //TOF electrons: includes all electrons, exlcludes Pions using the TPC
  mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4. , 0.0 ,100., kTRUE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  mastermind_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0. ,100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  
  

  AliAnalysisCuts* fancyCut=0x0;

  if(cutDefinition==0 || cutDefinition==4){
    printf("Simple 3 sigma in TPC\n");
    AliDielectronCutGroup* mastermind_cg = new AliDielectronCutGroup("mastermind_cg_1","mastermind_cg_1",AliDielectronCutGroup::kCompOR);
    mastermind_cg->AddCut(mastermind_TPC_3sigma);
    fancyCut = mastermind_cg;
  }
  if(cutDefinition==1 || cutDefinition==5){
    printf("TPC hardon rejection\n");
    AliDielectronCutGroup* mastermind_cg = new AliDielectronCutGroup("mastermind_cg_1","mastermind_cg_1",AliDielectronCutGroup::kCompOR);
    mastermind_cg->AddCut(mastermind_TPC);
    fancyCut = mastermind_cg;
  }
  if(cutDefinition==2 || cutDefinition==6){
    printf("TPC and TOF\n");
    AliDielectronCutGroup* mastermind_cg = new AliDielectronCutGroup("mastermind_cg_1","mastermind_cg_1",AliDielectronCutGroup::kCompOR);
    mastermind_cg->AddCut(mastermind_TOF);
    fancyCut = mastermind_cg;
  }
  if(cutDefinition==3 || cutDefinition==7){
    printf("Combined\n");
    AliDielectronCutGroup* mastermind_cg = new AliDielectronCutGroup("mastermind_cg_1","mastermind_cg_1",AliDielectronCutGroup::kCompOR);
    mastermind_cg->AddCut(mastermind_TPC);
    mastermind_cg->AddCut(mastermind_TOF);
    fancyCut = mastermind_cg;
  }
    
    return fancyCut;
}

// #########################################################
// #########################################################

AliAnalysisCuts* SetupEventCuts(Bool_t isAOD)
{
  // event cuts are identical for all analysis 'cutDefinition's that run together!
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  //not sure if this can be used, probably not:
  // Patrick:
  //  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  // Mahmut:
  //  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxTracksOrSPD);
  //  eventCuts->SetRequireV0and();
  return eventCuts;
}

// #########################################################
// #########################################################
void AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){
  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  task->AddSingleLegMCSignal(eleFinalState);

  AliDielectronSignalMC eleFinalStateFromPion("eleFinalStateFromPion","eleFinalStateFromPion");
  eleFinalStateFromPion.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromPion.SetMotherPDGs(111, 111); // open charm mesons and baryons together
  eleFinalStateFromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  task->AddSingleLegMCSignal(eleFinalStateFromPion);
  //
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddSingleLegMCSignal(eleFinalStateFromD);
  //
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  task->AddSingleLegMCSignal(eleFinalStateFromB);

  // AliDielectronSignalMC eleSecondary("eleSecondary","eleSecondary");
  // eleSecondary.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  // eleSecondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
  // eleSecondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  // task->AddSingleLegMCSignal(eleSecondary);
  //
  // AliDielectronSignalMC eleDontCare("eleDontCare","eleDontCare");
  // eleDontCare.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  // eleDontCare.SetCheckBothChargesLegs(kTRUE,kTRUE);
  // eleDontCare.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // task->AddSingleLegMCSignal(eleDontCare);
}


// #########################################################
// #########################################################
void AddPairMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){

    AliDielectronSignalMC pair_sameMother("sameMother","sameMother");
    pair_sameMother.SetLegPDGs(11,-11);
    pair_sameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother.SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. should have no effect on final state ele.
    task->AddPairMCSignal(pair_sameMother);

    AliDielectronSignalMC pair_sameMother_pion("sameMother_pion","sameMother_pion");
    pair_sameMother_pion.SetLegPDGs(11,-11);
    pair_sameMother_pion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_pion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_pion.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_pion.SetMotherPDGs(111,111); //
    task->AddPairMCSignal(pair_sameMother_pion);

    // AliDielectronSignalMC pair_sameMother_eta("sameMother_eta","sameMother_eta");
    // pair_sameMother_eta.SetLegPDGs(11,-11);
    // pair_sameMother_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
    // pair_sameMother_eta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    // //mother
    // pair_sameMother_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    // pair_sameMother_eta.SetMotherPDGs(221,221); //
    // task->AddPairMCSignal(pair_sameMother_eta);

    // AliDielectronSignalMC pair_conversion("pair_conversion","pair_conversion");
    // pair_conversion.SetLegPDGs(11,-11);
    // pair_conversion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    // pair_conversion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    // //mother
    // pair_conversion.SetMothersRelation(AliDielectronSignalMC::kSame);
    // task->AddPairMCSignal(pair_conversion);
}
void SetEtaCorrectionTOFMean(AliAnalysisTaskElectronEfficiencyV2 *task, Int_t corrXdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the TOF, can be one/two/three-dimensional
  //
  std::cout << "Set eta correction mean for TOF\n";
  std::string file_name = "CorrTOF_MC_07_11_2018.root";

  TFile* _file = TFile::Open(file_name.c_str());
  std::cout << _file << std::endl;
  if (_file == 0x0){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/r/rbailhac/PWGDQ/dielectron/macrosLMEE/CorrTOF_MC_07_11_2018.root  .");
    std::cout << "Copy TOF correction from Alien" << std::endl;
    _file = TFile::Open("CorrTOF_MC_07_11_2018.root");
    if(_file == 0x0) {
      printf("Did not find the file for mean\n");
      return;
    }
    else printf("Correction loaded\n");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }

  TH2D* mean = dynamic_cast<TH2D*>(_file->Get("Corrmean"));
  if(mean) {
    printf("Set the mean correction for TOF\n");
    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, mean,  corrXdim, corrYdim);
  }
  

}
void SetEtaCorrectionTOFRMS(AliAnalysisTaskElectronEfficiencyV2 *task, Int_t corrXdim, Int_t corrYdim) {
  //
  // eta correction for the centroid and width of electron sigmas in the TOF, can be one/two/three-dimensional
  //
  std::cout << "Set eta correction width for TOF\n";
  std::string file_name = "CorrTOF_MC_07_11_2018.root";

  TFile* _file = TFile::Open(file_name.c_str());
  std::cout << _file << std::endl;
  if (_file == 0x0){
    gSystem->Exec("alien_cp alien:///alice/cern.ch/user/r/rbailhac/PWGDQ/dielectron/macrosLMEE/CorrTOF_MC_07_11_2018.root .");
    std::cout << "Copy TOF correction from Alien" << std::endl;
    _file = TFile::Open("CorrTOF_MC_07_11_2018.root");
    if(_file == 0x0) {
      printf("Did not find the file for sigma\n");
      return;
    }
    else printf("correction loaded\n");
  }
  else {
    std::cout << "Correction loaded" << std::endl;
  }
  
  TH2D* width= dynamic_cast<TH2D*>(_file->Get("Corrsigma"));
  if(width){
    printf("Set the width correction for TOF\n");
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, width, corrXdim, corrYdim);
  }
  
  

}
