
TString names=("cut5_pt200");

// TString generatorNameForMCSignal  = "pizero_0;eta_1;etaprime_2;rho_3;omega_4;phi_5;jpsi_6"; // neue MC
// TString generatorNameForMCSignal  = "Hijing"; // neue MC
// TString generatorNameForMCSignal  = "pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7";
TString generatorNameForMCSignal  = "Hijing_0";
// TString generatorNameForMCSignal  = "";

bool SetGeneratedSmearingHistos = false;

bool DoPairing = true;
bool DoFourPairing = true;

bool GetResolutionFromAlien = kTRUE;
std::string resoFilename = "";
std::string resoFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/resolution_PbPb2015_0080_deltaXvsP.root";

bool DoCocktailWeighting = kFALSE;
bool GetCocktailFromAlien = kTRUE;
std::string CocktailFilename = "Cocktail_PbPb0080_5TeV_MotherPt.root";
std::string CocktailFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/Cocktail_PbPb0080_5TeV_MotherPt.root";

bool GetCentralityFromAlien = kFALSE;
std::string centralityFilename = "";
std::string centralityFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/centralityLHC16g1.root";

const Int_t triggerNames = AliVEvent::kMB;

const int nMCSignal = 0;
const int nCutsetting = 0;


const double minGenPt = 0.05;
const double maxGenPt = 100;
const double minGenEta = -1.5;
const double maxGenEta =  1.5;


// const double minPtCut = 0.2;
// const double maxPtCut = 8.0;
// const double minEtaCut = -0.8;
// const double maxEtaCut = 0.8;

const double minPtCut = 0.;
const double maxPtCut = 100.;
const double minEtaCut = -100;
const double maxEtaCut = 100;


// binning of single leg histograms
bool usePtVector = true;
double ptBins[] = {0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,8.0
  };
const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;

const double minPtBin = 0;
const double maxPtBin = 8;
const int    stepsPtBin = 320;

const double minEtaBin = -1.0;
const double maxEtaBin =  1.0;
const int    stepsEtaBin = 20;

const double minPhiBin = 0;
const double maxPhiBin =  6.3;
const int    stepsPhiBin = 20;

const double minThetaBin = 0;
const double maxThetaBin =  TMath::TwoPi();
const int    stepsThetaBin = 60;

const double minMassBin = 0;
const double maxMassBin =  5;
const int    stepsMassBin = 250;
const double minPairPtBin = 0;
const double maxPairPtBin =  8;
const int    stepsPairPtBin = 160;


Int_t centrality = -1; // DUMMY!!!!!!! will be overwritten

void GetCentrality(const int centrality, double& CentMin, double& CentMax){
  std::cout << "GetCentrality with centrality " << centrality << std::endl;
  if      (centrality == 0) {CentMin = 0;  CentMax = 10;}
  else if (centrality == 1) {CentMin = 10; CentMax = 50;}
  else if (centrality == 2) {CentMin = 50; CentMax = 80;}
  else if (centrality == 3) {CentMin = 0;  CentMax = 80;}
  else if (centrality == 4) {CentMin = 0;  CentMax = 10;}
  else if (centrality == 5) {CentMin = 10; CentMax = 20;}
  else if (centrality == 6) {CentMin = 20; CentMax = 40;}
  else if (centrality == 7) {CentMin = 40; CentMax = 50;}
  else if (centrality == 8) {CentMin = 50; CentMax = 60;}
  else if (centrality == 9) {CentMin = 60; CentMax = 80;}

  else if (centrality == 10) {CentMin = 0; CentMax = 10;}
  else if (centrality == 11) {CentMin = 10; CentMax = 20;}
  else if (centrality == 12) {CentMin = 20; CentMax = 30;}
  else if (centrality == 13) {CentMin = 30; CentMax = 40;}
  else if (centrality == 14) {CentMin = 40; CentMax = 50;}
  else if (centrality == 15) {CentMin = 50; CentMax = 60;}
  else if (centrality == 16) {CentMin = 60; CentMax = 70;}
  else if (centrality == 17) {CentMin = 70; CentMax = 80;}
  else if (centrality == 18) {CentMin = 80; CentMax = 90;}

  else if (centrality == 20) {CentMin = 0; CentMax = 1;}
  else if (centrality == 21) {CentMin = 1; CentMax = 2;}
  else if (centrality == 22) {CentMin = 2; CentMax = 3;}
  else if (centrality == 23) {CentMin = 3; CentMax = 4;}
  else if (centrality == 24) {CentMin = 4; CentMax = 5;}
  else if (centrality == 25) {CentMin = 5; CentMax = 6;}
  else if (centrality == 26) {CentMin = 6; CentMax = 7;}
  else if (centrality == 27) {CentMin = 7; CentMax = 8;}
  else if (centrality == 28) {CentMin = 8; CentMax = 9;}
  else if (centrality == 29) {CentMin = 9; CentMax = 10;}

  else if (centrality == 30) {CentMin = 0; CentMax = 5;}
  else if (centrality == 31) {CentMin = 5; CentMax = 10;}
  else if (centrality == 32) {CentMin = 10; CentMax = 15;}
  else if (centrality == 33) {CentMin = 15; CentMax = 20;}

  else                      {CentMin = 0;  CentMax = 0;}
  return;
}



// #########################################################
// #########################################################
AliAnalysisCuts* SetupEventCuts(Bool_t isAOD)
{
  std::cout << "Setup Event Cuts" << std::endl;
  // event cuts are identical for all analysis 'cutInstance's that run together!
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  if (isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  eventCuts->Print();
  return eventCuts;
}



// #########################################################
// #########################################################
std::vector<bool> AddSingleLegMCSignal(AliAnalysisTaskEtaReconstruction* task){
  /*AliDielectronSignalMC partFinalState("partFinalState","partFinalState");
  partFinalState.SetLegPDGs(0,1);//dummy second leg (never MCtrue)\n"
  partFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  partFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  */
  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC eleFinalStateFromSameMotherMeson("eleFinalStateFromSameMotherMeson","eleFinalStateFromSameMotherMeson");
  eleFinalStateFromSameMotherMeson.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromSameMotherMeson.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromSameMotherMeson.SetMotherPDGs(600, 600); // open charm mesons and baryons together
  eleFinalStateFromSameMotherMeson.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromEta("eleFinalStateFromEta","eleFinalStateFromEta");
  eleFinalStateFromEta.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromEta.SetMotherPDGs(221, 221); // open charm mesons and baryons together
  eleFinalStateFromEta.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromPion("eleFinalStateFromPion","eleFinalStateFromPion");
  eleFinalStateFromPion.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalStateFromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromPion.SetMotherPDGs(111, 111); // open charm mesons and baryons together
  eleFinalStateFromPion.SetCheckBothChargesMothers(kTRUE,kTRUE);

  AliDielectronSignalMC PhotonFinalState("PhotonFinalState","PhotonFinalState");
  PhotonFinalState.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC PhotonFinalStateFromSameMotherMeson("PhotonFinalStateFromSameMotherMeson","PhotonFinalStateFromSameMotherMeson");
  PhotonFinalStateFromSameMotherMeson.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromSameMotherMeson.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromSameMotherMeson.SetMotherPDGs(600, 600); // open charm mesons and baryons together
  PhotonFinalStateFromSameMotherMeson.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //
  AliDielectronSignalMC PhotonFinalStateFromD("PhotonFinalStateFromD","PhotonFinalStateFromD");
  PhotonFinalStateFromD.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PhotonFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  PhotonFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC PhotonFinalStateFromB("PhotonFinalStateFromB","PhotonFinalStateFromB");
  PhotonFinalStateFromB.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PhotonFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  PhotonFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC PhotonFinalStateFromEta("PhotonFinalStateFromEta","PhotonFinalStateFromEta");
  PhotonFinalStateFromEta.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PhotonFinalStateFromEta.SetMotherPDGs(221, 221); // open charm mesons and baryons together
  PhotonFinalStateFromEta.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC PhotonFinalStateFromPion("PhotonFinalStateFromPion","PhotonFinalStateFromPion");
  PhotonFinalStateFromPion.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonFinalStateFromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonFinalStateFromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  PhotonFinalStateFromPion.SetMotherPDGs(111, 111); // open charm mesons and baryons together
  PhotonFinalStateFromPion.SetCheckBothChargesMothers(kTRUE,kTRUE);

  AliDielectronSignalMC eleSecondary("eleSecondary","eleSecondary");
  eleSecondary.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleSecondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);

  AliDielectronSignalMC eleDontCare("eleDontCare","eleDontCare");
  eleDontCare.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleDontCare.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleDontCare.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);

  AliDielectronSignalMC eleSecondaryFromPhoton("eleSecondaryFromPhoton","eleSecondaryFromPhoton");
  eleSecondaryFromPhoton.SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleSecondaryFromPhoton.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleSecondaryFromPhoton.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  eleSecondaryFromPhoton.SetMotherPDGs(22, 1); // open charm mesons and baryons together
  eleSecondaryFromPhoton.SetCheckBothChargesMothers(kTRUE,kTRUE);
  // eleSecondaryFromPhoton.SetMotherSource(AliDielectronSignalMC::kSecondary);

  AliDielectronSignalMC PhotonSecondary("PhotonSecondary","PhotonSecondary");
  PhotonSecondary.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  PhotonSecondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
  PhotonSecondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);

  AliDielectronSignalMC PhotonDontCare("PhotonDontCare","PhotonDontCare");
  eleDontCare.SetLegPDGs(22,1);//dummy second leg (never MCtrue)
  eleDontCare.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleDontCare.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);


  // task->AddSingleLegMCSignal(eleFinalState);
  // task->AddSingleLegMCSignal(partFinalState);
  // task->AddSingleLegMCSignal(eleFinalStateFromPion);
  // task->AddSingleLegMCSignal(eleFinalStateFromD);
  // task->AddSingleLegMCSignal(eleFinalStateFromB);
  // task->AddSingleLegMCSignal(eleFinalStateFromSameMotherMeson);
  task->AddSingleLegMCSignal(eleFinalStateFromEta);
  // //
  // task->AddSingleLegMCSignal(PhotonFinalState);
  // task->AddSingleLegMCSignal(PhotonFinalStateFromPion);
  // task->AddSingleLegMCSignal(PhotonFinalStateFromD);
  // task->AddSingleLegMCSignal(PhotonFinalStateFromB);
  // task->AddSingleLegMCSignal(PhotonFinalStateFromSameMotherMeson);
  // task->AddSingleLegMCSignal(PhotonFinalStateFromEta);
  //
  // task->AddSingleLegMCSignal(eleDontCare);
  // task->AddSingleLegMCSignal(eleSecondary);
  task->AddSingleLegMCSignal(eleSecondaryFromPhoton);
  // task->AddSingleLegMCSignal(PhotonSecondary);
  // task->AddSingleLegMCSignal(PhotonDontCare);


  // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  // ordering is according to MCSignals of single legs
 std::vector<bool> DielectronsPairNotFromSameMother;
 DielectronsPairNotFromSameMother.push_back(true);
 DielectronsPairNotFromSameMother.push_back(true);
 // DielectronsPairNotFromSameMother.push_back(true);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 // DielectronsPairNotFromSameMother.push_back(false);
 return DielectronsPairNotFromSameMother;
}


// #########################################################
// #########################################################
void AddPairMCSignal(AliAnalysisTaskEtaReconstruction* task){

    AliDielectronSignalMC pair_sameMother_photon_finalstate("pair_sameMother_photon_finalstate","pair_sameMother_photon_finalstate");
    pair_sameMother_photon_finalstate.SetLegPDGs(11,-11);
    pair_sameMother_photon_finalstate.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_finalstate.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_photon_finalstate.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_finalstate.SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. should have no effect on final state ele.

    AliDielectronSignalMC pair_sameMother_photon_secondary("pair_sameMother_photon_secondary","pair_sameMother_photon_secondary");
    pair_sameMother_photon_secondary.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_sameMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_sameMother_photon_secondaryfromMaterial("pair_sameMother_photon_secondaryfromMaterial","pair_sameMother_photon_secondaryfromMaterial");
    pair_sameMother_photon_secondaryfromMaterial.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondaryfromMaterial.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondaryfromMaterial.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
    //mother
    pair_sameMother_photon_secondaryfromMaterial.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondaryfromMaterial.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_sameMother_photon_secondaryfromWD("pair_sameMother_photon_secondaryfromWD","pair_sameMother_photon_secondaryfromWD");
    pair_sameMother_photon_secondaryfromWD.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondaryfromWD.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondaryfromWD.SetLegSources(AliDielectronSignalMC::kSecondaryFromWeakDecay, AliDielectronSignalMC::kSecondaryFromWeakDecay);
    //mother
    pair_sameMother_photon_secondaryfromWD.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondaryfromWD.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_UndifinedMother_photon_secondary("pair_UndefinedMother_photon_secondary","pair_UndefinedMother_photon_secondary");
    pair_UndifinedMother_photon_secondary.SetLegPDGs(11,-11);
    pair_UndifinedMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_UndifinedMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_UndifinedMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kUndefined);
    pair_UndifinedMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_DifferentMother_photon_secondary("pair_DifferentMother_photon_secondary","pair_DifferentMother_photon_secondary");
    pair_DifferentMother_photon_secondary.SetLegPDGs(11,-11);
    pair_DifferentMother_photon_secondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_DifferentMother_photon_secondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_DifferentMother_photon_secondary.SetMothersRelation(AliDielectronSignalMC::kDifferent);
    pair_DifferentMother_photon_secondary.SetMotherPDGs(22,22);

    AliDielectronSignalMC pair_sameMother_photon_secondary_eta("pair_sameMother_photon_secondary_eta","pair_sameMother_photon_secondary_eta");
    pair_sameMother_photon_secondary_eta.SetLegPDGs(11,-11);
    pair_sameMother_photon_secondary_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_photon_secondary_eta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_sameMother_photon_secondary_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondary_eta.SetMotherPDGs(22,22);
    //grand-mother
    pair_sameMother_photon_secondary_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_photon_secondary_eta.SetGrandMotherPDGs(221,221);

    AliDielectronSignalMC pair_sameMother_eta("sameMother_eta","sameMother_eta");
    pair_sameMother_eta.SetLegPDGs(11,-11);
    pair_sameMother_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_eta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_eta.SetMotherPDGs(221,221); //

    AliDielectronSignalMC pair_conversion("pair_conversion","pair_conversion");
    pair_conversion.SetLegPDGs(11,-11);
    pair_conversion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_conversion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_conversion.SetMothersRelation(AliDielectronSignalMC::kSame);

    AliDielectronSignalMC pair_random("pair_random","pair_random");
    pair_random.SetLegPDGs(11,-11);
    pair_random.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_random.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);

    AliDielectronSignalMC pair_sameMother_pion("pair_sameMother_pion","pair_sameMother_pion");
    pair_sameMother_pion.SetLegPDGs(11,-11);
    pair_sameMother_pion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_pion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_pion.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_pion.SetMotherPDGs(111,111); //

    AliDielectronSignalMC pair_sameMother_CharmedMesonsWithSameMother("CharmedMesonsWithSameMother","CharmedMesonsWithSameMother");
    pair_sameMother_CharmedMesonsWithSameMother.SetLegPDGs(11,-11);
    pair_sameMother_CharmedMesonsWithSameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_CharmedMesonsWithSameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_CharmedMesonsWithSameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_CharmedMesonsWithSameMother.SetMotherPDGs(402, 402); //

    AliDielectronSignalMC pair_sameMother_BeautyMesonsWithSameMother("BeautyMesonsWithSameMother","BeautyMesonsWithSameMother");
    pair_sameMother_BeautyMesonsWithSameMother.SetLegPDGs(11,-11);
    pair_sameMother_BeautyMesonsWithSameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_BeautyMesonsWithSameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_BeautyMesonsWithSameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_BeautyMesonsWithSameMother.SetMotherPDGs(502, 502); //



    task->AddPairMCSignal(pair_sameMother_photon_finalstate);
    task->AddPairMCSignal(pair_sameMother_photon_secondary);
    task->AddPairMCSignal(pair_sameMother_photon_secondaryfromMaterial);
    task->AddPairMCSignal(pair_sameMother_photon_secondaryfromWD);
    task->AddPairMCSignal(pair_UndifinedMother_photon_secondary);
    task->AddPairMCSignal(pair_DifferentMother_photon_secondary);
    task->AddPairMCSignal(pair_sameMother_photon_secondary_eta);
    task->AddPairMCSignal(pair_sameMother_eta);
    task->AddPairMCSignal(pair_conversion);
    task->AddPairMCSignal(pair_random);
    task->AddPairMCSignal(pair_sameMother_pion);
    // task->AddPairMCSignal(pair_sameMother_CharmedMesonsWithSameMother);
    // task->AddPairMCSignal(pair_sameMother_BeautyMesonsWithSameMother);
}


  void AddFourPairULSMCSignal(AliAnalysisTaskEtaReconstruction* task){
    //#########################################
    //########### Unlike Signe Signals ########
    //#########################################
    // AliDielectronSignalMC for 4 particles
    // Seperate it into 2 pair MCSignals
    //
    // First Pair
      AliDielectronSignalMC ULSFourElePair1_FromEta("ULSFourElePair1_FromEta","ULSFourElePair1_FromEta");
      ULSFourElePair1_FromEta.SetLegPDGs(11,-11);
      ULSFourElePair1_FromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
      ULSFourElePair1_FromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      //mother
      ULSFourElePair1_FromEta.SetMothersRelation(AliDielectronSignalMC::kSame);
      ULSFourElePair1_FromEta.SetMotherPDGs(221, 221); //

    // Second Pair
      AliDielectronSignalMC ULSFourElePair2_FromEta("ULSFourElePair2_FromEta","ULSFourElePair2_FromEta");
      ULSFourElePair2_FromEta.SetLegPDGs(11,-11);
      ULSFourElePair2_FromEta.SetCheckBothChargesLegs(kTRUE,kTRUE);
      ULSFourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
      // FourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      // FourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
      // mother
      ULSFourElePair2_FromEta.SetMothersRelation(AliDielectronSignalMC::kSame);
      ULSFourElePair2_FromEta.SetMotherPDGs(22, 22); //
      // grand-mother
      ULSFourElePair2_FromEta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
      ULSFourElePair2_FromEta.SetGrandMotherPDGs(221, 221); //

//________________________________________________________
      // First Pair
        AliDielectronSignalMC ULSFourElePair1_FromPion("ULSFourElePair1_FromPion","ULSFourElePair1_FromPion");
        ULSFourElePair1_FromPion.SetLegPDGs(11,-11);
        ULSFourElePair1_FromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
        ULSFourElePair1_FromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        //mother
        ULSFourElePair1_FromPion.SetMothersRelation(AliDielectronSignalMC::kSame);
        ULSFourElePair1_FromPion.SetMotherPDGs(111, 111); //

      // Second Pair
        AliDielectronSignalMC ULSFourElePair2_FromPion("ULSFourElePair2_FromPion","ULSFourElePair2_FromPion");
        ULSFourElePair2_FromPion.SetLegPDGs(11,-11);
        ULSFourElePair2_FromPion.SetCheckBothChargesLegs(kTRUE,kTRUE);
        ULSFourElePair2_FromPion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
        // FourElePair2_FromPion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
        // FourElePair2_FromPion.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
        // mother
        ULSFourElePair2_FromPion.SetMothersRelation(AliDielectronSignalMC::kSame);
        ULSFourElePair2_FromPion.SetMotherPDGs(22, 22); //
        // grand-mother
        ULSFourElePair2_FromPion.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
        ULSFourElePair2_FromPion.SetGrandMotherPDGs(111, 111); //

//________________________________________________________
      task->AddFourPairULSMCSignal(ULSFourElePair1_FromEta);
      task->AddFourPairULSMCSignal(ULSFourElePair2_FromEta);
      // task->AddFourPairULSMCSignal(ULSFourElePair1_FromPion);
      // task->AddFourPairULSMCSignal(ULSFourElePair2_FromPion);

}

void AddFourPairLSMCSignal(AliAnalysisTaskEtaReconstruction* task){
    //#########################################
    //########### Like Signe Signals ##########
    //#########################################
    // AliDielectronSignalMC for 4 particles
    // Seperate it into 2 pair MCSignals
    //
    // First Pair
      AliDielectronSignalMC LSFourElePair1_FromEta("LSFourElePair1_FromEta","LSFourElePair1_FromEta");
      LSFourElePair1_FromEta.SetLegPDGs(11,11);
      LSFourElePair1_FromEta.SetCheckBothChargesLegs(kFALSE,kFALSE);
      LSFourElePair1_FromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      //mother
      // LSFourElePair1_FromEta.SetMothersRelation(AliDielectronSignalMC::kSame);
       LSFourElePair1_FromEta.SetMotherPDGs(221, 221); //

    // Second Pair
      AliDielectronSignalMC LSFourElePair2_FromEta("LSFourElePair2_FromEta","LSFourElePair2_FromEta");
      LSFourElePair2_FromEta.SetLegPDGs(-11,-11);
      LSFourElePair2_FromEta.SetCheckBothChargesLegs(kFALSE,kFALSE);
      LSFourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
      // FourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
      // FourElePair2_FromEta.SetLegSources(AliDielectronSignalMC::kSecondaryFromMaterial, AliDielectronSignalMC::kSecondaryFromMaterial);
      // mother
      // LSFourElePair2_FromEta.SetMothersRelation(AliDielectronSignalMC::kSame);
      LSFourElePair2_FromEta.SetMotherPDGs(22, 22); //
      // grand-mother
      // LSFourElePair2_FromEta.SetGrandMothersRelation(AliDielectronSignalMC::kSame);
       LSFourElePair2_FromEta.SetGrandMotherPDGs(221, 221); //

//________________________________________________________
        // task->AddFourPairLSMCSignal(LSFourElePair1_FromEta);
        // task->AddFourPairLSMCSignal(LSFourElePair2_FromEta);

}
