//TString generatorNameForMCSignal  = "pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7";
//TString generatorNameForMCSignal  = "Hijing_0;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;Pythia CC_8;Pythia BB_8;Pythia B_8";
//TString generatorNameForMCSignal  = "Pythia CC_1;Pythia BB_1;Pythia B_1;Jpsi2ee_1;B2JPsi2ee_1";
// TString generatorNameForMCSignal  = "Hijing_0;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7";
 
//TString generatorNameForULSSignal = "Hijing_0;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7;Pythia CC_8;Pythia BB_8;Pythia B_8";
//TString generatorNameForULSSignal = "Hijing_0;pizero_1;eta_2;etaprime_3;rho_4;omega_5;phi_6;jpsi_7";




void Config_slehner_Efficiency(AliAnalysisTaskElectronEfficiencyV2 *task,  Bool_t useAODFilterCuts,  TString TMVAweight){
  Int_t trackCut=0;
  Int_t PIDCut=0;
  Int_t MVACut=0;
  
  for(int glcut = 0; glcut <=30; ++glcut){
//  for(int glcut = 0; glcut <=0; ++glcut){
    ////////DEFINE THE CUTS AS FUNCTION OF GLCUT//////
    if(glcut>0 && glcut<21) continue;    
    PIDCut=glcut-10;
    trackCut=glcut;
    if(glcut==0) trackCut=-1;
    if(glcut==0) PIDCut=0;
    
    std::cout << "Config_slehner_Efficiency: CutTr: "<<trackCut<<" CutPID: "<<PIDCut<<" MVA Cut: "<<0<<" added"<< std::endl;
    AliAnalysisFilter* filter = SetupTrackCutsAndSettings(trackCut, PIDCut, MVACut, useAODFilterCuts,TMVAweight);
    task->AddTrackCuts(filter); 
    }
}

Bool_t setGens=kTRUE;  //decides if generator to be used are set (e.g. for LHC18b5a) or not (e.g. LHC16g1)

Bool_t SetTPCCorrection = kTRUE;
Bool_t SetITSCorrection = kTRUE;
Bool_t SetTOFCorrection = kFALSE;

Bool_t SetGeneratedSmearingHistos = kFALSE;

Bool_t DoPairing    = kTRUE;
Bool_t DoULSLS      = kTRUE;
Bool_t DeactivateLS = kFALSE;

// Leave blank to not use resolution files
std::string resoFilename = "resolution_PbPb2015_0080_deltaXvsP.root";
std::string resoFilenameFromAlien = "/alice/cern.ch/user/s/selehner/reso/resolution_PbPb2015_0080_deltaXvsP.root";

Bool_t DoCocktailWeighting  = kFALSE;
Bool_t GetCocktailFromAlien = kFALSE;
std::string CocktailFilename = "";
std::string CocktailFilenameFromAlien = "/alice/cern.ch/user/s/slehner/.root";
// std::string CocktailFilenameFromAlien = "/alice/cern.ch/user/c/cklein/data/cocktail_PbPb0080_5TeV.root";

Bool_t GetCentralityFromAlien = kFALSE;
std::string centralityFilename = "";
std::string centralityFilenameFromAlien = "/alice/cern.ch/user/s/slehner/.root";

const Int_t triggerNames = AliVEvent::kINT7;

const Int_t nMCSignal   = 0;
const Int_t nCutsetting = 0;

const Double_t minGenPt  = 0.15;
const Double_t maxGenPt  = 10;
const Double_t minGenEta = -1.2;
const Double_t maxGenEta = 1.2;

const Double_t minPtCut  = 0.2;
const Double_t maxPtCut  = 8.0;
const Double_t minEtaCut = -0.8;
const Double_t maxEtaCut = 0.8;



// binning of single leg histograms
Bool_t usePtVector = kTRUE;
Double_t ptBins[] = {0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,8.0
  };
const Int_t nBinsPt =  ( sizeof(ptBins) / sizeof(ptBins[0]) )-1;

const Double_t minPtBin   = 0;
const Double_t maxPtBin   = 10;
const Int_t    stepsPtBin = 400;

const Double_t minEtaBin   = -1.0;
const Double_t maxEtaBin   = 1.0;
const Int_t    stepsEtaBin = 20;

const Double_t minPhiBin   = 0;
const Double_t maxPhiBin   = 6.3;
const Int_t    stepsPhiBin = 20;

const Double_t minThetaBin   = 0;
const Double_t maxThetaBin   = TMath::TwoPi();
const Int_t    stepsThetaBin = 60;

const Double_t minMassBin     = 0;
const Double_t maxMassBin     = 10;
const Int_t    stepsMassBin   = 200;
const Double_t minPairPtBin   = 0;
const Double_t maxPairPtBin   = 10;
const Int_t    stepsPairPtBin = 20;

//varying bin size
//lmee mass spectrum
//  double mbinsarr[] = { 0.00, 0.02 ,0.04 ,0.08 ,0.14 ,0.22 ,0.38 ,0.54 ,1.1 ,1.7 ,2.5 ,2.9 ,3.0 ,3.1 ,3.3 ,3.5 ,4.0 ,5.0}; //Carsten's binning
//  double ptbinsarr[]= {0.0,0.4,0.6,1,2.5,8};
  
//low ptee
double mbinsarr[]= { 0.0,0.1,0.4,0.5 ,0.6 ,0.7 ,1.1, 1.5,2.0 ,2.7,2.9,3.1 ,5.0}; // for low ptee
double ptbinsarr[]= {0.0, 0.025, 0.05, 0.075, 0.1,0.125, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1, 2.0, 5.0, 8.0};// for low ptee
//  TVectorD* centbins= AliDielectronHelper::MakeLinBinning(10,0,100);

// Binning of resolution histograms
const Int_t    NbinsDeltaMom   = 2000;
const Double_t DeltaMomMin     = -10.0;
const Double_t DeltaMomMax     = 10.0;
const Int_t    NbinsRelMom     = 400;
const Double_t RelMomMin       = 0.0;
const Double_t RelMomMax       = 2.0;
const Int_t    NbinsDeltaEta   = 200;
const Double_t DeltaEtaMin     = -0.4;
const Double_t DeltaEtaMax     = 0.4;
const Int_t    NbinsDeltaTheta = 200;
const Double_t DeltaThetaMin   = -0.4;
const Double_t DeltaThetaMax   = 0.4;
const Int_t    NbinsDeltaPhi   = 200;
const Double_t DeltaPhiMin     = -0.4;
const Double_t DeltaPhiMax     = 0.4;

void GetCentrality(const Int_t centrality, Double_t& CentMin, Double_t& CentMax){
  std::cout << "GetCentrality with centrality " << centrality << std::endl;
  if     (centrality == 0) {CentMin = 0;  CentMax = 90;}
  else if(centrality == 1) {CentMin = 0;  CentMax = 20;}
  else if(centrality == 2) {CentMin = 20; CentMax = 60;}
  else if(centrality == 3) {CentMin = 60; CentMax = 90;}
  else                      {CentMin = 0; CentMax = 0;}
  return;
}

void setPIDCorrections(AliAnalysisTaskElectronEfficiencyV2* task){
  std::cout << task << std::endl;

  std::cout << "starting SetPIDCorrections()\n";
  if(SetTPCCorrection == kTRUE){
    std::cout << "Loading TPC correction" << std::endl;
    std::string file_name = "recalib_mc_tpc_nsigmaele.root";
    TFile* _file = TFile::Open(file_name.c_str());

    if(!_file){
      TString path="alien:///alice/cern.ch/user/s/selehner/recal/recalib_mc_tpc_nsigmaele.root";
      gSystem->Exec(TString::Format("alien_cp %s .",path.Data()));
      std::cout << "Copy TPC correction from Alien" << std::endl;
      _file = TFile::Open(file_name.c_str());
      if(!_file) std::cout << "Could not get alien:///alice/cern.ch/user/s/selehner/recal/recalib_mc_tpc_nsigmaele.root" << std::endl;      
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTPC, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTPC, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
  if(SetITSCorrection == kTRUE){
    std::cout << "Loading ITS correction" << std::endl;
    std::string file_name = "recalib_mc_its_nsigmaele.root";
    TFile* _file = TFile::Open(file_name.c_str());

    if(!_file){
      TString path="alien:///alice/cern.ch/user/s/selehner/recal/recalib_mc_its_nsigmaele.root";
      gSystem->Exec(TString::Format("alien_cp %s .",path.Data()));
      std::cout << "Copy ITS correction from Alien" << std::endl;
      _file = TFile::Open(file_name.c_str());
      if(!_file) std::cout << "Could not get alien:///alice/cern.ch/user/s/selehner/recal/recalib_mc_its_nsigmaele.root" << std::endl;
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kITS, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kITS, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
  if(SetTOFCorrection == kTRUE){
    std::cout << "Loading TOF correction" << std::endl;
    std::string file_name = "recalib_mc_tof_nsigmaele.root";
    TFile* _file = TFile::Open(file_name.c_str());

    if(!_file){
      TString path="alien:///alice/cern.ch/user/s/selehner/recal/recalib_mc_tof_nsigmaele.root";
      gSystem->Exec(TString::Format("alien_cp %s .",path.Data()));
      std::cout << "Copy TOF correction from Alien" << std::endl;
      _file = TFile::Open(file_name.c_str());
      if(!_file) std::cout << "Could not get alien:///alice/cern.ch/user/s/selehner/recal/recalib_mc_tof_nsigmaele.root" << std::endl;
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
}

// #########################################################
// #########################################################

AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t selTr, Int_t selPID, Int_t MVACut, Bool_t useAODFilterCuts,TString TMVAweight)
{
  std::cout<<"SetupTrackCutsAndSettings: TrackCut: "<<selTr<<", PIDCut; "<<selPID<<" MVACut: "<<MVACut*0.2<<std::endl;
  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!

  LMEECutLib* LMcutlib = new LMEECutLib();

  anaFilter->AddCuts(LMcutlib->GetTrackCuts(selTr, selPID, MVACut, useAODFilterCuts,TMVAweight));     // Setting MVA cut for efficiency to 0 - no efficiency correction for MVA cut here
  anaFilter->SetName(TString::Format("CutTr%d_PID%d_MVA%d",selTr, selPID,MVACut,TMVAweight));
  anaFilter->Print();
  return anaFilter;
}

// #########################################################
// #########################################################
/* AliAnalysisCuts* SetupEventCuts(Bool_t isAOD) */
/* { */
/*   std::cout << "Setup Event Cuts" << std::endl; */
/*   // event cuts are identical for all analysis 'cutInstance's that run together! */
/*   AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0"); */
/*   eventCuts->SetRequireVertex(); */
/*   eventCuts->SetMinVtxContributors(1); */
/*   eventCuts->SetVertexZ(-10.,10.); */
/*   if(isAOD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD */
/*   eventCuts->Print(); */
/*   return eventCuts; */
/* } */



// #########################################################
// #########################################################
std::vector<Bool_t> AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){
  AliDielectronSignalMC partFinalState("partFinalState","partFinalState");
  partFinalState.SetLegPDGs(0,1);//dummy second leg (never MCkTRUE)\n"
  // partFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  partFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCkTRUE)\n"
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);    //in slight O(1%) disagreement with ML tree maker non-conversion electrons
  eleFinalState.SetMotherPDGs(22,22,kTRUE,kTRUE); // this line leads to results in agreement with ML tree maker when counting non-conversion electrons
  
  AliDielectronSignalMC eleFinalStateFromSameMotherMeson("eleFinalStateFromSameMotherMeson","eleFinalStateFromSameMotherMeson");
  eleFinalStateFromSameMotherMeson.SetLegPDGs(11,1);//dummy second leg (never MCkTRUE)\n"
  eleFinalStateFromSameMotherMeson.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromSameMotherMeson.SetMotherPDGs(600, 600); // open charm mesons and baryons together
  eleFinalStateFromSameMotherMeson.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);//dummy second leg (never MCkTRUE)\n"
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);//dummy second leg (never MCkTRUE)\n"
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);

  // AliDielectronSignalMC eleSecondary("eleSecondary","eleSecondary");
  // eleSecondary.SetLegPDGs(11,1);//dummy second leg (never MCkTRUE)\n"
  // eleSecondary.SetCheckBothChargesLegs(kTRUE,kTRUE);
  // eleSecondary.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
  // task->AddSingleLegMCSignal(eleSecondary);
  //
//   AliDielectronSignalMC eleDontCare("eleDontCare","eleDontCare");
//   eleDontCare.SetLegPDGs(11,1);//dummy second leg (never MCkTRUE)\n"
//   eleDontCare.SetCheckBothChargesLegs(kTRUE,kTRUE);
//   eleDontCare.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
//   task->AddSingleLegMCSignal(eleDontCare);

  AliDielectronSignalMC partConv("partConv","partConv");
  partConv.SetLegPDGs(0,1);//dummy second leg (never MCkTRUE)\n"
  partConv.SetCheckBothChargesLegs(kTRUE,kTRUE);
  partConv.SetMotherPDGs(22,22,kFALSE,kFALSE); 

  AliDielectronSignalMC eleConv("eleConv","eleConv");
  eleConv.SetLegPDGs(11,1);//dummy second leg (never MCkTRUE)\n"
  eleConv.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleConv.SetMotherPDGs(22,22,kFALSE,kFALSE); 

//  task->AddSingleLegMCSignal(partConv);
//  task->AddSingleLegMCSignal(eleConv);
  task->AddSingleLegMCSignal(eleFinalState);
  // task->AddSingleLegMCSignal(eleFinalStateFromPion);
  task->AddSingleLegMCSignal(eleFinalStateFromD);
  task->AddSingleLegMCSignal(eleFinalStateFromB);
//  task->AddSingleLegMCSignal(eleFinalStateFromSameMotherMeson);

  // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  // ordering is according to MCSignals of single legs
 std::vector<Bool_t> DielectronsPairNotFromSameMother;
 DielectronsPairNotFromSameMother.push_back(kFALSE);
 DielectronsPairNotFromSameMother.push_back(kTRUE);
 DielectronsPairNotFromSameMother.push_back(kTRUE);
// DielectronsPairNotFromSameMother.push_back(kFALSE);
 return DielectronsPairNotFromSameMother;
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

    AliDielectronSignalMC pair_sameMother_pion("sameMother_pion","sameMother_pion");
    pair_sameMother_pion.SetLegPDGs(11,-11);
    pair_sameMother_pion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_pion.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_pion.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_pion.SetMotherPDGs(111,111); //

    AliDielectronSignalMC pair_sameMother_eta("sameMother_eta","sameMother_eta");
    pair_sameMother_eta.SetLegPDGs(11,-11);
    pair_sameMother_eta.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_eta.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_eta.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_eta.SetMotherPDGs(221,221); //

    AliDielectronSignalMC pair_sameMother_etaP("sameMother_etaP","sameMother_etaP");
    pair_sameMother_etaP.SetLegPDGs(11,-11);
    pair_sameMother_etaP.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_etaP.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_etaP.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_etaP.SetMotherPDGs(331,331); //

    AliDielectronSignalMC pair_sameMother_rho("sameMother_rho","sameMother_rho");
    pair_sameMother_rho.SetLegPDGs(11,-11);
    pair_sameMother_rho.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_rho.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_rho.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_rho.SetMotherPDGs(113, 113); //

    AliDielectronSignalMC pair_sameMother_omega("sameMother_omega","sameMother_omega");
    pair_sameMother_omega.SetLegPDGs(11,-11);
    pair_sameMother_omega.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_omega.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_omega.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_omega.SetMotherPDGs(223, 223); //

    AliDielectronSignalMC pair_sameMother_phi("sameMother_phi","sameMother_phi");
    pair_sameMother_phi.SetLegPDGs(11,-11);
    pair_sameMother_phi.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_phi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_phi.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_phi.SetMotherPDGs(333, 333); //

    AliDielectronSignalMC pair_sameMother_jpsi("sameMother_jpsi","sameMother_jpsi");
    pair_sameMother_jpsi.SetLegPDGs(11,-11);
    pair_sameMother_jpsi.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_jpsi.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_jpsi.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_jpsi.SetMotherPDGs(443, 443); //

    AliDielectronSignalMC pair_sameMother_CharmedMesonsWithSameMother("CharmedMesonsWithSameMother","CharmedMesonsWithSameMother");
    pair_sameMother_CharmedMesonsWithSameMother.SetLegPDGs(11,-11);
    pair_sameMother_CharmedMesonsWithSameMother.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_CharmedMesonsWithSameMother.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
    //mother
    pair_sameMother_CharmedMesonsWithSameMother.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_CharmedMesonsWithSameMother.SetMotherPDGs(402, 402); //

    AliDielectronSignalMC pair_conversion("pair_conversion","pair_conversion");
    pair_conversion.SetLegPDGs(11,-11);
    pair_conversion.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_conversion.SetLegSources(AliDielectronSignalMC::kSecondary, AliDielectronSignalMC::kSecondary);
    //mother
    pair_conversion.SetMothersRelation(AliDielectronSignalMC::kSame);



  task->AddPairMCSignal(pair_sameMother);
//  task->AddPairMCSignal(pair_sameMother_pion);
//  task->AddPairMCSignal(pair_sameMother_eta);
//  task->AddPairMCSignal(pair_sameMother_etaP);
//  task->AddPairMCSignal(pair_sameMother_rho);
//  task->AddPairMCSignal(pair_sameMother_omega);
//  task->AddPairMCSignal(pair_sameMother_phi);
//  task->AddPairMCSignal(pair_sameMother_jpsi);
//  task->AddPairMCSignal(pair_sameMother_CharmedMesonsWithSameMother);
//  task->AddPairMCSignal(pair_conversion);
}
