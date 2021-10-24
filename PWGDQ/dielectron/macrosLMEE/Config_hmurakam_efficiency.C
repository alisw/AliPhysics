TString names("NewPid_Phiv");
bool DoPairing = kTRUE;
bool DoULSLS = kTRUE;

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
/* const Double_t ptBins[] = { */
/*   0.00,0.10,0.11,0.12,0.13,0.14,0.15,0.155,0.16,0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.20,0.205,0.21,0.215,0.22,0.225,0.23,0.235,0.24,0.245,0.25,0.255, */
/*   0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,0.30,0.32,0.34,0.36,0.38,0.40,0.43,0.46, */
/*   0.49,0.52,0.55,0.60,0.65,0.70,0.75,0.80,0.90,1.00,1.10,1.20, */
/*   1.40,1.60,1.80,2.00,2.40,2.80,3.20,3.70,4.50,6.00,8.00,12.0 */
/* }; */
const Double_t ptBins[] = {
  0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,10.0,20.0
};
const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;

// Ivan's ptee bins
//const Double_t PteeBins[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 20.0};
//const Int_t nBinsPtee = ( sizeof(PteeBins) / sizeof(PteeBins[0]) )-1;

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
const double maxMassBin =  4.;//4.
const int    stepsMassBin = 400;//400

const double minPhiVBin = 0;
const double maxPhiVBin =  TMath::Pi();
const int    stepsPhiVBin = 180;

// Binning of resolution histograms
const int    NbinsDeltaMom    = 1;
const double DeltaMomMin   =-10.0;
const double DeltaMomMax   =  2.0;
const int    NbinsRelMom      = 1;
const double RelMomMin     =  0.0;
const double RelMomMax     =  2.0;
const int    NbinsDeltaEta    = 1;
const double DeltaEtaMin   = -0.4;
const double DeltaEtaMax   =  0.4;
const int    NbinsDeltaTheta  = 1;
const double DeltaThetaMin = -0.4;
const double DeltaThetaMax =  0.4;
const int    NbinsDeltaPhi    = 1;
const double DeltaPhiMin   = -0.4;
const double DeltaPhiMax   =  0.4;

AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t cutDefinition, Bool_t isAOD);
AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition);
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition);
AliAnalysisCuts* SetupEventCuts(Bool_t isAOD);

//________________________________________________________________
AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t cutDefinition, Bool_t isAOD)
{
  std::cout << "SetupTrackCutsAndSettings()" <<std::endl;


  AliAnalysisFilter *anaFilter = new AliAnalysisFilter(Form("anaFilter_%d",cutDefinition),Form("anaFilter_%d",cutDefinition)); // named constructor seems mandatory!
  // do not change these initial values!
  //Int_t selectedPID=-1;
  //Bool_t isPrefilterCutset=kFALSE;

  // AliDielectronV0Cuts *noconv = new AliDielectronV0Cuts("IsGamma","IsGamma");
  // // which V0 finder you want to use
  // noconv->SetV0finder(AliDielectronV0Cuts::kAll);  // kAll(default), kOffline or kOnTheFly
  // // add some pdg codes (they are used then by the KF package and important for gamma conversions)
  // noconv->SetPdgCodes(22,11,11); // mother, daughter1 and 2
  // // add default PID cuts (defined in AliDielectronPID)
  // // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
  // //noconv->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
  // // add the pair cuts for V0 candidates
  // noconv->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.00, kFALSE);
  // noconv->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.00, kFALSE);
  // noconv->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
  // noconv->AddCut(AliDielectronVarManager::kR,                             3.0,  90.00, kFALSE);
  // noconv->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
  // noconv->AddCut(AliDielectronVarManager::kM,                             0.0,   0.10, kFALSE);
  // noconv->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
  // // selection or rejection of V0 tracks
  // noconv->SetExcludeTracks(kTRUE); 
  
  anaFilter->AddCuts( SetupTrackCuts(cutDefinition) );
  anaFilter->AddCuts( SetupPIDcuts(cutDefinition) );
  // anaFilter->AddCuts( noconv );
  std::cout << "...cuts added!" <<std::endl; 

  anaFilter->Print();
  
  return anaFilter;
}
//________________________________________________________________
AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition)
{
  std::cout << "SetupTrackCuts()" <<std::endl;
  //AliAnalysisCuts* trackCuts=0x0;

  AliDielectronTrackCuts *trackCutsDiel = new AliDielectronTrackCuts("trackCutsDiel","trackCutsDiel");
  trackCutsDiel->SetAODFilterBit(AliDielectronTrackCuts::kGlobalNoDCA);
  trackCutsDiel->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst);
  trackCutsDiel->SetRequireITSRefit(kTRUE);
  trackCutsDiel->SetRequireTPCRefit(kTRUE);

  AliDielectronVarCuts* trackCutsAOD =new AliDielectronVarCuts("trackCutsAOD","trackCutsAOD");
  // pT and eta
  trackCutsAOD->AddCut(AliDielectronVarManager::kPt, 0.2,   1e30);
  trackCutsAOD->AddCut(AliDielectronVarManager::kEta, -0.8,   0.8);

  //primary selection
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParXY,    -1.0,   1.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kImpactParZ,     -3.0,   3.0);

  //TPC
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCr,     100.0, 160.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsTPC,        80.0, 160.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kTPCchi2Cl,       0.0,   4.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNFclsTPCfCross,  0.8,   1.5);
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSFracTPC,    0.0,   0.4);
  //ITS
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsITS,         3.0, 100.0);
  trackCutsAOD->AddCut(AliDielectronVarManager::kITSchi2Cl,       0.0,   4.5);

  //  if(cutDefinition==1){
  printf("Add shared cluster cut\n");
  trackCutsAOD->AddCut(AliDielectronVarManager::kNclsSITS, 1.0, 6.0, kTRUE);//accept no shared cls hit (default)
  //  }

  // AliDielectronVarCuts *pairCutsPhiV(0x0);//Does not work
  // if(cutDefinition==1 || cutDefinition==3){
  //   pairCutsPhiV = new AliDielectronVarCuts("PhiV","PhiV");//mass and Phiv together
  //   pairCutsPhiV->AddCut(AliDielectronVarManager::kPhivPair, 2.0 , 3.2,kTRUE);
  //   pairCutsPhiV->AddCut(AliDielectronVarManager::kM, 0.0 , 0.040,kTRUE);
  //   pairCutsPhiV->SetCutType(1);
  //   //    die->GetPairFilter().AddCuts(pairCutsPhiV);              //PhiV cut
  // }
  
  AliDielectronCutGroup* trackCuts = new AliDielectronCutGroup("Trackcuts","Trackcuts",AliDielectronCutGroup::kCompAND);
  trackCuts->AddCut(trackCutsAOD);
  trackCuts->AddCut(trackCutsDiel);
  //  trackCuts->AddCut(pairCutsPhiV);//Does not work

  trackCuts->Print();

  return trackCuts;

}

//________________________________________________________________
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition)
{

  std::cout << "SetupPIDcuts()" <<std::endl;

  //PID1
  AliDielectronPID *pidTPCTOFreq = new AliDielectronPID("pidTPCTOFreq","pidTPCTOFreq");//<--> same as recoverTOF
  pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,     -3. , 3., 0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCTOFreq->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. , 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCTOFreq->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,      -3., 3., 0.4, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kP);// TOF required
  //PID2
  AliDielectronPID *pidTPCHadRej = new AliDielectronPID("pidTPCHadRej","pidTPCHadRej");//pure TPC pid
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,  -3., 3., 0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kPion,    -100., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kKaon,      -4., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRej->AddCut(AliDielectronPID::kTPC, AliPID::kProton,    -4., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  //PID3
  AliDielectronPID *pidTPCHadRejTOFif = new AliDielectronPID("pidTPCHadRejTOFif","pidTPCHadRejTOFif");//pure TPC pid + TOF cut if avilable
  pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,  -3., 3., 0.0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kPion,    -100., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kKaon,      -4., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTPC, AliPID::kProton,    -4., 4., 0.0, 1e30, kTRUE,  AliDielectronPID::kRequire, AliDielectronVarManager::kPt);//TPC
  pidTPCHadRejTOFif->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,  -3., 3., 0.4, 1e30, kFALSE,  AliDielectronPID::kIfAvailable, AliDielectronVarManager::kP);// TOF if available

  AliAnalysisCuts* fancyCut=0x0;
  // combine 2 cut sets with OR option
  if(cutDefinition==0){//New pid
    AliDielectronCutGroup* combinedPIDcuts = new AliDielectronCutGroup("combinedPIDcuts","combinedPIDcuts",AliDielectronCutGroup::kCompOR);
    combinedPIDcuts->AddCut(pidTPCTOFreq);//PID1
    combinedPIDcuts->AddCut(pidTPCHadRej);//PID2
    fancyCut = combinedPIDcuts;
  }else if (cutDefinition==1){// Old pid
    AliDielectronCutGroup* combinedPIDcuts2 = new AliDielectronCutGroup("combinedPIDcuts2","combinedPIDcuts2", AliDielectronCutGroup::kCompOR);
    combinedPIDcuts2->AddCut(pidTPCTOFreq);//PID1
    combinedPIDcuts2->AddCut(pidTPCHadRejTOFif);//PID3
    fancyCut = combinedPIDcuts2;
  }

  return fancyCut;

}

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

  // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  // ordering is according to MCSignals of single legs
  printf("Init the DielectronsPairNotFromSameMother vector\n");
  std::vector<bool> DielectronsPairNotFromSameMother;
  DielectronsPairNotFromSameMother.push_back(false);
  DielectronsPairNotFromSameMother.push_back(false);
  DielectronsPairNotFromSameMother.push_back(false);
  DielectronsPairNotFromSameMother.push_back(false);
  task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);

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

    //    AliDielectronSignalMC pair_sameMother_pion("sameMother_pion","sameMother_pion");
    //    pair_sameMother_pion.SetLegPDGs(11,-11);
    //    pair_sameMother_pion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    //    pair_sameMother_pion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    //    pair_sameMother_pion.SetMothersRelation(AliDielectronSignalMC::kSame);
    //    pair_sameMother_pion.SetMotherPDGs(111,111); //
    //task->AddPairMCSignal(pair_sameMother_pion);

    AliDielectronSignalMC eleFromJPsi("eleFromJPsi", "eleFromJPsi");
    eleFromJPsi.SetLegPDGs(11,-11);
    eleFromJPsi.SetCheckBothChargesLegs(kTRUE,kTRUE);
    eleFromJPsi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    eleFromJPsi.SetMotherPDGs(443, 443);
    eleFromJPsi.SetMothersRelation(AliDielectronSignalMC::kSame);
    eleFromJPsi.SetCheckBothChargesMothers(kTRUE,kTRUE);
    task->AddPairMCSignal(pair_sameMother);
    task->AddPairMCSignal(eleFromJPsi);

    // D mesons
    AliDielectronSignalMC diEleOpenCharmNeutralCharged("Dmesons","di-electrons from open charm no B grandmother");  // dielectrons originating from open charm hadrons
    diEleOpenCharmNeutralCharged.SetLegPDGs(11,-11);
    diEleOpenCharmNeutralCharged.SetMotherPDGs(402,402);
    diEleOpenCharmNeutralCharged.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    diEleOpenCharmNeutralCharged.SetMothersRelation(AliDielectronSignalMC::kDifferent);
    //diEleOpenCharmNeutralCharged.SetCheckStackForPDG(kTRUE);
    //diEleOpenCharmNeutralCharged.SetPDGforStack(503);
    diEleOpenCharmNeutralCharged.SetCheckBothChargesLegs(kTRUE,kTRUE);
    diEleOpenCharmNeutralCharged.SetCheckBothChargesMothers(kTRUE,kTRUE);
    //task->AddPairMCSignal(diEleOpenCharmNeutralCharged);

    // AliDielectronSignalMC pair_sameMother_eta("sameMother_eta","sameMother_eta");
    // pair_sameMother_eta.SetLegPDGs(11,-11);
    // pair_sameMother_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
    // pair_sameMother_eta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    // //mother
    // pair_sameMother_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    // pair_sameMother_eta.SetMotherPDGs(221,221); //
    // task->AddPairMCSignal(pair_sameMother_eta);

    AliDielectronSignalMC pair_conversion("pair_conversion","pair_conversion");
    pair_conversion.SetLegPDGs(11,-11);
    pair_conversion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_conversion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_conversion.SetMothersRelation(AliDielectronSignalMC::kSame);
    task->AddPairMCSignal(pair_conversion);
}
//______________________________________________________________________________________
void SetTOFSigmaEleCorrection(AliAnalysisTaskElectronEfficiencyV2 *task, Int_t corrXdim, Int_t corrYdim, TString year) {
  //
  // eta correction for the centroid and width of electron sigmas in the TOF
  //

  //
  printf("starting SetTOFSigmaEleCorrection() MC\n");
  printf(" correction Xdim = %s\n", AliDielectronVarManager::GetValueName(corrXdim));
  printf(" correction Ydim = %s\n", AliDielectronVarManager::GetValueName(corrYdim));

  if (corrYdim!=AliDielectronVarManager::kEta) { printf(" correction only available for Ydim = eta!\n"); return; }

  if (corrXdim!=AliDielectronVarManager::kP)
    {
      printf(" no correction available for Xdim = %s!\n", AliDielectronVarManager::GetValueName(corrXdim));
      printf(" no correction applied!\n");
      return;
    }

  TH2F* histMean2DTOF;
  TH2F* histWidth2DTOF;

  TString rootFile = "alien://alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/calibLMEE/calMaps_TOF_mc.root";
  TGrid::Connect("alien://");
  gSystem->Exec(Form("alien_cp %s .",rootFile.Data()));
  TString fileName= "calMaps_TOF_mc.root";
  TFile* f = TFile::Open(fileName.Data());
  //  printf("Load correction map:%s mean:%s and width:%s",fileName.Data(),Form("m%s",year.Data(), Form("w%s",year.Data());
  TString meanName =Form("m%s",year.Data());
  TString widthName =Form("w%s",year.Data());
  f->GetObject(meanName.Data(),histMean2DTOF);
  f->GetObject(widthName.Data(),histWidth2DTOF);
  printf("%s and %s\n",meanName.Data(),widthName.Data());

  for (Int_t i = 0; i <= histMean2DTOF->GetNbinsX()+1; i++){
    for (Int_t k = 0; k <= histMean2DTOF->GetNbinsY()+1; k++){
      if ( (i == 0) || (k == 0) || (i > histMean2DTOF->GetNbinsX()) || (k > histMean2DTOF->GetNbinsY())) { // under/overflows
	histMean2DTOF->SetBinContent(i, k, 0.0 );
	histWidth2DTOF->SetBinContent(i, k, 1.0 );
      }
    }
  }

  /* die->SetCentroidCorrFunctionTOF(histMean2DTOF, corrXdim, corrYdim); */
  /* die->SetWidthCorrFunctionTOF(histWidth2DTOF, corrXdim, corrYdim); */
  //    printf("no default TOF PID correction! Corrections were not applied!\n");

  task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, histWidth2DTOF, corrXdim, corrYdim);
  task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, histMean2DTOF,  corrXdim, corrYdim);

}
