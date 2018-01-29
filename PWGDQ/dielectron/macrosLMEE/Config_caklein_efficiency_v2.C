
TString names=("kPbPb2015_Pt400_looseTOFif;kPbPb2015_Pt400_tightTOFreq;NewPID;kPbPb2015_Pt400_looseTOFif_0SITSCl;kPbPb2015_Pt400_tightTOFreq_0SITSCl;NewPID_0SITSCl;kPbPb2015_Pt400_looseTOFif_ITSMAP;kPbPb2015_Pt400_tightTOFreq_ITSMAP;NewPID_ITSMAP");

TString generatorName = "";
// TString generatorName = "Hijing_0";
// TString generatorName = "Jpsi2ee_1";

bool DoPairing = kTRUE;
bool DoULSLS = true;

bool GetResolutionFromAlien = kTRUE;
std::string resoFilename = "resolution_PbPb2015_0080_deltaXvsP.root";
std::string resoFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/resolution_PbPb2015_0080_deltaXvsP.root";

bool GetCocktailFromAlien = kTRUE;
std::string CocktailFilename = "cocktail_PbPb0080_5TeV.root";
std::string CocktailFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/cocktail_PbPb0080_5TeV.root";

bool GetCentralityFromAlien = kFALSE;
std::string centralityFilename = "";
std::string centralityFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/centralityLHC16g1.root";

const Int_t triggerNames = AliVEvent::kMB;

const int nMCSignal = 0;
const int nCutsetting = 0;

const double centMin = 0.;
const double centMax = 80.;

const double minGenPt = 0.1;
const double maxGenPt = 100;
const double minGenEta = -1.5;
const double maxGenEta =  1.5;


const double minPtCut = 0.4;
const double maxPtCut = 8.0;
const double minEtaCut = -0.8;
const double maxEtaCut = 0.8;



// binning of single leg histograms
bool usePtVector = true;
double ptBins[] = {0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,8.0
  };
const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;

const double minPtBin = 0;
const double maxPtBin = 8;
const int    stepsPtBin = 800;

const double minEtaBin = -1.0;
const double maxEtaBin =  1.0;
const int    stepsEtaBin = 20;

const double minPhiBin = 0;
const double maxPhiBin =  6.3;
const int    stepsPhiBin = 30;

const double minThetaBin = 0;
const double maxThetaBin =  TMath::TwoPi();
const int    stepsThetaBin = 60;

const double minMassBin = 0;
const double maxMassBin =  5;
const int    stepsMassBin = 250;
const double minPairPtBin = 0;
const double maxPairPtBin =  8;
const int    stepsPairPtBin = 160;

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


// #########################################################
// #########################################################

AliAnalysisFilter* SetupTrackCutsAndSettings(TString cutDefinition, Bool_t isAOD)
{
  std::cout << "TEST" << std::endl;
  std::cout << "SetupTrackCutsAndSettings( cutInstance = " << cutDefinition << " )" <<std::endl;
  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!

  AnalysisCut AnaCut;


  LMEECutLib* LMcutlib = new LMEECutLib();
  if (cutDefinition == "kPbPb2015_Pt100_ResolutionCuts"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt100_ResolutionCuts);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kResolutionTrackCuts);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "NewPID"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_tightTOFreq"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_tightTOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_looseTOFif"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_looseTOFif);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "NewPID_0SITSCl"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst_0SITSCl);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_tightTOFreq_0SITSCl"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_tightTOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst_0SITSCl);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_looseTOFif_0SITSCl"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_looseTOFif);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst_0SITSCl);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "NewPID_ITSMAP"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_PID_cutoff_pion_kaon_proton);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst_ITSMAP);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_tightTOFreq_ITSMAP"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_tightTOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst_ITSMAP);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_looseTOFif_ITSMAP"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_looseTOFif);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst_ITSMAP);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }

  if (!isAOD) anaFilter->AddCuts( LMcutlib->GetESDTrackCutsAna(AnaCut) );
  anaFilter->AddCuts( LMcutlib->GetPIDCutsAna(AnaCut) );
  anaFilter->SetName(cutDefinition);
  return anaFilter;
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
