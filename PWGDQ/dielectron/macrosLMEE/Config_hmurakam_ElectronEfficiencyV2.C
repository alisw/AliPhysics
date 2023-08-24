void DoAdditionalWork(AliAnalysisTaskElectronEfficiencyV2* task);
void DoAdditionalWork2(AliAnalysisTaskElectronEfficiencyV2* task);
AliAnalysisFilter* SetupCuts(Int_t cutDefinition);
AliAnalysisCuts *SetupPIDcuts(Int_t cutDefinition);
AliDielectronVarCuts *SetupTrackCuts(Int_t cutDefinition);
AliDielectronEventCuts *GetEventCuts();
std::vector<bool> AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task);
void AddPairMCSignal(AliAnalysisTaskElectronEfficiencyV2* task);

AliAnalysisTaskElectronEfficiencyV2* Config_Devel (TString name, int wagonnr);

TString names("cut0;cut1");

TString year = "16";

Int_t  isAOD = 1;
Bool_t DoCentralityCorrection = kFALSE;
Bool_t cutlibPreloaded = kFALSE;
Int_t  centrality = -1; // -1 for else, so min bias

TString generatorNameForMCSignal  = "";//GP
TString generatorNameForULSSignal = "";//GP
//TString generatorNameForMCSignal  = "Pythia CC_1;Pythia BB_1;Pythia B_1";//2016 HF MC
//TString generatorNameForULSSignal = "Pythia CC_1;Pythia BB_1;Pythia B_1";//2016 HF MC
//TString generatorNameForMCSignal  = "Pythia CC_0;Pythia BB_0;Pythia B_0";//201[7,8] HF MC
//TString generatorNameForULSSignal = "Pythia CC_0;Pythia BB_0;Pythia B_0";//201[7,8] HF MC
//TString generatorNameForMCSignal  = "Jpsi2ee_1;B2Jpsi2ee_1";//Jpsi MC
//TString generatorNameForULSSignal = "Jpsi2ee_1;B2Jpsi2ee_1";//Jpsi MC

Bool_t SetTOFCorrection = kTRUE;
bool SetGeneratedSmearingHistos = true;

bool DoPairing = true;
bool DoULSLS   = true;

bool DeactivateLS   = true;

//const Int_t triggerNames = AliVEvent::kMB;
const Int_t triggerNames = AliVEvent::kINT7;

Bool_t usePhiV           = kFALSE;
Double_t maxMee          = 0.04;
Double_t minphiv         = 2.0;

const int nMCSignal = 0;
const int nCutsetting = 0;

const double minGenPt = 0.05;
const double maxGenPt = 100;
const double minGenEta = -1.5;
const double maxGenEta =  1.5;


// const double minPtCut = 0.2;
// const double minPtCut = 0.;
// const double maxPtCut = 100.0;
// const double minEtaCut = -100;
// const double maxEtaCut =  100;
const double minPtCut = 0.2;
// const double minPtCut = 0.2;
const double maxPtCut = 10.0;
const double minEtaCut = -0.8;
const double maxEtaCut = 0.8;

// binning of single leg histograms
bool usePtVector = true;
double ptBins[] = {0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
                   1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,4.00,5.0,6.0,7.0,8.0, 9.0,10.};
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

// Binning of resolution histograms
const int    NbinsDeltaMom    = 1000;
const double DeltaMomMin   =-10.0;
const double DeltaMomMax   = 10.0;
const int    NbinsRelMom      = 200;
const double RelMomMin     =  0.0;
const double RelMomMax     =  2.0;
const int    NbinsDeltaEta    = 100;
const double DeltaEtaMin   = -0.4;
const double DeltaEtaMax   =  0.4;
const int    NbinsDeltaTheta  = 100;
const double DeltaThetaMin = -0.4;
const double DeltaThetaMax =  0.4;
const int    NbinsDeltaPhi    = 100;
const double DeltaPhiMin   = -0.4;
const double DeltaPhiMax   =  0.4;


//Wrapper to return task with right name (same as config file)
//forward declarations
//get rif of cutlib

AliAnalysisTaskElectronEfficiencyV2* Config_Devel(TString name, int wagonnr)
{
  AliAnalysisTaskElectronEfficiencyV2* task = new AliAnalysisTaskElectronEfficiencyV2(Form("%s%d",name.Data(), wagonnr));

  // #########################################################
  // #########################################################
  // Possibility to set generator. If nothing set all generators are taken into account
  // task->SetGeneratorName(generatorName);
  task->SetGeneratorMCSignalName(generatorNameForMCSignal);
  task->SetGeneratorULSSignalName(generatorNameForULSSignal);

  // #########################################################
  // #########################################################
  // Event selection. Is the same for all the different cutsettings
  task->SetEnablePhysicsSelection(kTRUE);
  task->SetTriggerMask(triggerNames);
  task->SetEventFilter(GetEventCuts()); //returns eventCuts from Config.

  double centMin = 0.;
  double centMax = 100.;

  std::cout << "CentMin = " << centMin << "  CentMax = " << centMax << std::endl;
  task->SetCentrality(centMin, centMax);

  // #########################################################
  // #########################################################
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  // Do not set here your analysis pt-cuts
  task->SetMinPtGen(minGenPt);
  task->SetMaxPtGen(maxGenPt);
  task->SetMinEtaGen(minGenEta);
  task->SetMaxEtaGen(maxGenEta);


  // #########################################################
  // #########################################################
  // Set minimum and maximum values of generated tracks. Only used to save computing power.
  task->SetKinematicCuts(minPtCut, maxPtCut, minEtaCut, maxEtaCut);

  // #########################################################
  // #########################################################
  // Set Binning
  if (usePtVector == true) {
    std::vector<double> ptBinsVec;
    for (unsigned int i = 0; i < nBinsPt+1; ++i){
      ptBinsVec.push_back(ptBins[i]);
    }
    task->SetPtBins(ptBinsVec);
  }
  else task->SetPtBinsLinear   (minPtBin,  maxPtBin, stepsPtBin);
  task->SetEtaBinsLinear  (minEtaBin, maxEtaBin, stepsEtaBin);
  task->SetPhiBinsLinear  (minPhiBin, maxPhiBin, stepsPhiBin);
  task->SetThetaBinsLinear(minThetaBin, maxThetaBin, stepsThetaBin);
  // Make arbitrary bins matching data!
  std::vector<double> massBins = {0.00,0.02,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58,0.60,0.62,0.64,0.66,0.68,0.70,0.72,0.74,0.76,0.78,0.80,0.82,0.84,0.86,0.88,0.90,0.92,0.94,0.96,0.98,1.00,1.02,1.04,1.06,1.08,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.20,2.30,2.40,2.50,2.60,2.70,2.80,2.90,3.00,3.02,3.04,3.06,3.08,3.10,3.12,3.30,3.50,4.00,5.00};
  task->SetMassBins (massBins);
  // task->SetMassBinsLinear (minMassBin, maxMassBin, stepsMassBin);
  std::vector<double> pteeBins= {0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.100,2.200,2.300,2.400,2.500,2.600,2.700,2.800,2.900,3.000,3.100,3.200,3.300,3.400,3.500,3.600,3.700,3.800,3.900,4.000,4.100,4.200,4.300,4.400,4.500,5.000,5.500,6.000,6.500,7.000,8.000};
  task->SetPairPtBins(pteeBins);
  // task->SetPairPtBinsLinear(minPairPtBin, maxPairPtBin, stepsPairPtBin);

  // #########################################################
  // #########################################################
  task->SetSmearGenerated                (SetGeneratedSmearingHistos);
  task->SetResolutionDeltaPtBinsLinear   (DeltaMomMin, DeltaMomMax, NbinsDeltaMom);
  task->SetResolutionRelPtBinsLinear     (RelMomMin, RelMomMax, NbinsRelMom);
  task->SetResolutionEtaBinsLinear       (DeltaEtaMin, DeltaEtaMax, NbinsDeltaEta);
  task->SetResolutionPhiBinsLinear       (DeltaPhiMin, DeltaPhiMax, NbinsDeltaPhi);
  task->SetResolutionThetaBinsLinear     (DeltaThetaMin, DeltaThetaMax, NbinsDeltaTheta);

  // #########################################################
  // #########################################################
  // Set centrality correction. If resoFilename = "" no correction is applied

  // #########################################################
  // #########################################################
  // Set MCSignal and Cutsetting to fill the support histograms
  task->SetSupportHistoMCSignalAndCutsetting(nMCSignal, nCutsetting);

  // #########################################################
  // #########################################################
  // Pairing related config
  task->SetDoPairing(DoPairing);
  task->SetULSandLS(DoULSLS);
  task->SetDeactivateLS(DeactivateLS);
  task->SetPhiVBinsLinear(0, TMath::Pi(), 180);
  task->SetFillPhiV(kFALSE);

  //Set Phiv Cut
  task->SetPhiVCut(usePhiV,maxMee,minphiv);
  
  // #########################################################
  // #########################################################
  // Add MCSignals. Can be set to see differences of:
  // e.g. secondaries and primaries. or primaries from charm and resonances
  AddSingleLegMCSignal(task);
  AddPairMCSignal(task);
  std::vector<bool> DielectronsPairNotFromSameMother = AddSingleLegMCSignal(task);
  task->AddMCSignalsWhereDielectronPairNotFromSameMother(DielectronsPairNotFromSameMother);

  // #########################################################
  // #########################################################
  // Adding cutsettings
  TObjArray*  arrNames=names.Tokenize(";");
  const Int_t nDie=arrNames->GetEntriesFast();

  for (int iCut = 0; iCut < nDie; ++iCut){
    TString cutDefinition(arrNames->At(iCut)->GetName());
    AliAnalysisFilter* filter = SetupCuts(iCut);
    filter->SetName(cutDefinition);
    task->AddTrackCuts(filter);
    DoAdditionalWork(task);
  }

  return task;
}
// #########################################################
// #########################################################
// Set mean and width correction for ITS, TPC and TOF


void DoAdditionalWork(AliAnalysisTaskElectronEfficiencyV2* task){
  
  //Load PID post calibration
  std::cout << task << std::endl;

  std::cout << "starting DoAdditionalWork()\n";
  
  if (SetTOFCorrection == true){
    std::cout << "Loading TOF correction" << std::endl;
    std::string file_name = "calMaps_TOF_mc.root";    
    TFile* _file = TFile::Open(file_name.c_str());    
    if (!_file){
      gSystem->Exec(("alien_cp alien:///alice/cern.ch/user/h/hmurakam/PWGDQ/dielectron/calibLMEE/" + file_name + " file:.").c_str());
      
      std::cout << "Copy TOF correction from Alien" << std::endl;
      _file = TFile::Open(file_name.c_str());
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }
    
    TH2F* histMean2DTOF  = dynamic_cast<TH2F*>(_file->Get(Form("m%s",year.Data())));
    TH2F* histWidth2DTOF = dynamic_cast<TH2F*>(_file->Get(Form("w%s",year.Data())));
    printf("%s and %s\n",Form("m%s",year.Data()),Form("w%s",year.Data()));
    
    for (Int_t i = 0; i <= histMean2DTOF->GetNbinsX()+1; i++){
      for (Int_t k = 0; k <= histMean2DTOF->GetNbinsY()+1; k++){
	if ( (i == 0) || (k == 0) || (i > histMean2DTOF->GetNbinsX()) || (k > histMean2DTOF->GetNbinsY())) { // under/overflows
	  histMean2DTOF->SetBinContent(i, k, 0.0 );
	  histWidth2DTOF->SetBinContent(i, k, 1.0 );
	}
      }
    }
    
    cout<<"Adding mean & width TOF PID correction" <<endl;
    task->SetWidthCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, histWidth2DTOF, AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, histMean2DTOF, AliDielectronVarManager::kPIn, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
}

void DoAdditionalWork2(AliAnalysisTaskElectronEfficiencyV2* task){
  //Load PID post calibration
  std::cout << task << std::endl;

  std::cout << "starting DoAdditionalWork()\n";

  if (SetTOFCorrection == true){
    std::cout << "Loading TOF correction" << std::endl;
    std::string file_name = "outputTOF_MC.root";
    TFile* _file = TFile::Open(file_name.c_str());

    if (!_file){
      gSystem->Exec(("alien_cp alien:///alice/cern.ch/user/h/hscheid/supportFiles/PIDrecalibration/TOF_MC/" + file_name + " file:.").c_str());
      std::cout << "Copy TOF correction from Alien" << std::endl;
      _file = TFile::Open(file_name.c_str());
    }
    else {
      std::cout << "Correction loaded" << std::endl;
    }

    TH3D* mean = dynamic_cast<TH3D*>(_file->Get("sum_mean_correction"));
    TH3D* width= dynamic_cast<TH3D*>(_file->Get("sum_width_correction"));

    task->SetCentroidCorrFunction(AliAnalysisTaskElectronEfficiencyV2::kTOF, mean,  AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
    task->SetWidthCorrFunction   (AliAnalysisTaskElectronEfficiencyV2::kTOF, width, AliDielectronVarManager::kP, AliDielectronVarManager::kEta, AliDielectronVarManager::kRefMultTPConly);
  }
}


//______________________________________________________________________________________
AliAnalysisFilter* SetupCuts(Int_t cutDefinition)
{

  AliAnalysisFilter *die = new AliAnalysisFilter("cuts","cuts");
  // Setup specified cuts
  // AOD track filter (needs to be first cut to speed up)
  AliDielectronTrackCuts *trkFilter = new AliDielectronTrackCuts("TrkFilter","TrkFilter");
  trkFilter->SetAODFilterBit(AliDielectronTrackCuts::kTPCqual);
  trkFilter->SetRequireITSRefit(kTRUE);
  trkFilter->SetRequireTPCRefit(kTRUE);
  // trkFilter->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD,AliDielectronTrackCuts::kFirst); // removed for track study

  //Add track cuts to dielectron object
  die->AddCuts(trkFilter);
  // shared cluster cut -> maybe add to track cuts later, for variation
  AliDielectronVarCuts *nSharedClsITS = new AliDielectronVarCuts("nSharedClsITS","nSharedClsITS");
  nSharedClsITS->AddCut(AliDielectronVarManager::kNclsSITS, 1.0,   6.0, kTRUE);

  if (cutDefinition != 0){  //for resolution maps: dont use any analysis (Track/PID) cuts
    die->AddCuts(SetupPIDcuts(cutDefinition));
    die->AddCuts(SetupPIDcuts(cutDefinition)); // removed for track study
    die->AddCuts(SetupTrackCuts(cutDefinition)); // removed for track study
    die->AddCuts(nSharedClsITS);
  }  

  return die;

}
//______________________________________________________________________________________
//-----------------------------------pid------------------------------------------------

AliAnalysisCuts *SetupPIDcuts(Int_t cutDefinition){

  std::cout << ">>>>>>>>>>>> Setup PID cuts! <<<<<<<<<<<<" << '\n';
  //New PID standard cut: full PID setting with recovery of ITS and TOF tracks
  //accepts an electron if it's identified in one of the 3 detector samples

  AliDielectronPID *recover_TPC = new AliDielectronPID("recover_TPC","recover_TPC");
  AliDielectronPID *recover_TOF = new AliDielectronPID("recover_TOF","recover_TOF");

  //TPC electrons: includes electrons and exclude all possible other contributions using the TPC
  //possible elemination of contamination using ITS and TOF
  recover_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3. , 0. ,100., kFALSE);
  recover_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 3.5 , 0.00 ,100., kTRUE);
  recover_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kKaon, -3, 3. , 0.200 ,100., kTRUE);
  recover_TPC->AddCut(AliDielectronPID::kTPC,AliPID::kProton, -3, 3. , 0.200 ,100., kTRUE);
  
  //TOF electrons: includes all electrons, exlcludes Pions using the TPC
  //possible elemination of contamination using ITS
  recover_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kElectron, -3., 3. , 0.0 ,100., kFALSE);
  recover_TOF->AddCut(AliDielectronPID::kTPC,AliPID::kPion, -100, 4. , 0.0 ,100., kTRUE);
  recover_TOF->AddCut(AliDielectronPID::kTOF,AliPID::kElectron, -3. , 3. , 0.0 ,100., kFALSE, AliDielectronPID::kRequire);
  // if the ITS identify a track as an electron include it back in
  
  // AliDielectronCutGroup* recover_cg = new AliDielectronCutGroup("recover_cg","recover_cg",AliDielectronCutGroup::kCompAND);
  AliDielectronCutGroup* recover_cg = new AliDielectronCutGroup("recover_cg","recover_cg",AliDielectronCutGroup::kCompOR);
  recover_cg->AddCut(recover_TPC);
  recover_cg->AddCut(recover_TOF);

  AliAnalysisCuts* returnCut = NULL;
  returnCut = reinterpret_cast<AliAnalysisCuts*> (recover_cg);

  return returnCut;
}


//______________________________________________________________________________________
//-----------------------------------track cuts-----------------------------------------

//Reimplementation of SetupESDtrackCuts. Uses AliDielectronVarCuts to also be used in AOD analysis
AliDielectronVarCuts *SetupTrackCuts(Int_t cutDefinition){

  AliDielectronVarCuts *fTrackCuts = new AliDielectronVarCuts("fTrackCuts","fTrackCuts");
  //global
  fTrackCuts->AddCut(AliDielectronVarManager::kPt,  0.2, 100.); //SetPtRange( 0.2 , 100. );
  fTrackCuts->AddCut(AliDielectronVarManager::kEta,-0.8, 0.8);  //SetEtaRange( -0.8 , 0.8 );
  fTrackCuts->AddCut(AliDielectronVarManager::kImpactParXY,   -1., 1.); //SetMaxDCAToVertexZ(3.);
  fTrackCuts->AddCut(AliDielectronVarManager::kImpactParZ,    -3., 3.); //SetMaxDCAToVertexXY(1.);
  fTrackCuts->AddCut(AliDielectronVarManager::kNclsTPC,       100., 160.);  //SetMinNClustersTPC(100);
  fTrackCuts->AddCut(AliDielectronVarManager::kNFclsTPCr,      100., 160.); //SetMinNClustersTPC(100);
  fTrackCuts->AddCut(AliDielectronVarManager::kNFclsTPCfCross, 0.5, 1.1);   //SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
  fTrackCuts->AddCut(AliDielectronVarManager::kTPCchi2Cl,      0.0, 4.);    //SetMaxChi2PerClusterTPC(4);
  fTrackCuts->AddCut(AliDielectronVarManager::kNclsITS,        3. , 10.);   //SetMinNClustersITS(3);
  fTrackCuts->AddCut(AliDielectronVarManager::kITSchi2Cl,      0.0, 5.5);   //SetMaxChi2PerClusterITS(5.5);
 
  AliDielectronVarCuts *returnCut = NULL;
  returnCut = fTrackCuts;
  return returnCut;


}



AliDielectronEventCuts *GetEventCuts(){

  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  eventCuts->SetRequireVertex();
  eventCuts->SetVertexZ(-10.,10.);
  eventCuts->SetMinVtxContributors(1);

  return eventCuts;
}


// #########################################################
// #########################################################
std::vector<bool> AddSingleLegMCSignal(AliAnalysisTaskElectronEfficiencyV2* task){
  AliDielectronSignalMC partFinalState("partFinalState","partFinalState");
  partFinalState.SetLegPDGs(0,1);//dummy second leg (never MCtrue)\n"
  // partFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  partFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC eleFinalState("eleFinalState","eleFinalState");
  eleFinalState.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalState.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);

  AliDielectronSignalMC eleFinalStateFromSameMotherMeson("eleFinalStateFromSameMotherMeson","eleFinalStateFromSameMotherMeson");
  eleFinalStateFromSameMotherMeson.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromSameMotherMeson.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromSameMotherMeson.SetMotherPDGs(600, 600); // open charm mesons and baryons together
  eleFinalStateFromSameMotherMeson.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  //
  AliDielectronSignalMC eleFinalStateFromD("eleFinalStateFromD","eleFinalStateFromD");
  eleFinalStateFromD.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromD.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromD.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromD.SetMotherPDGs(402, 402); // open charm mesons and baryons together
  eleFinalStateFromD.SetCheckBothChargesMothers(kTRUE,kTRUE);
  //
  AliDielectronSignalMC eleFinalStateFromB("eleFinalStateFromB","eleFinalStateFromB");
  eleFinalStateFromB.SetLegPDGs(11,1);//dummy second leg (never MCtrue)\n"
  eleFinalStateFromB.SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalStateFromB.SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalStateFromB.SetMotherPDGs(502, 502); // open charm mesons and baryons together
  eleFinalStateFromB.SetCheckBothChargesMothers(kTRUE,kTRUE);

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

  // task->AddSingleLegMCSignal(partFinalState);
  task->AddSingleLegMCSignal(eleFinalState);
  // task->AddSingleLegMCSignal(eleFinalStateFromPion);
  task->AddSingleLegMCSignal(eleFinalStateFromD);
  task->AddSingleLegMCSignal(eleFinalStateFromB);
  // task->AddSingleLegMCSignal(eleFinalStateFromSameMotherMeson);

  // this is used to get electrons from charmed mesons in a environment where GEANT is doing the decay of D mesons, like in LHC18b5a
  // ordering is according to MCSignals of single legs
 std::vector<bool> DielectronsPairNotFromSameMother;
 DielectronsPairNotFromSameMother.push_back(false);
 DielectronsPairNotFromSameMother.push_back(false);//true
 DielectronsPairNotFromSameMother.push_back(false);//true
 // DielectronsPairNotFromSameMother.push_back(false);
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

    AliDielectronSignalMC pair_sameMother_pion_anySource("sameMother_pion_anySource","sameMother_pion_anySource");
    pair_sameMother_pion_anySource.SetLegPDGs(11,-11);
    pair_sameMother_pion_anySource.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_pion_anySource.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
    //mother
    pair_sameMother_pion_anySource.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_pion_anySource.SetMotherPDGs(111,111); //

    AliDielectronSignalMC pair_sameMother_pion_fromK("sameMother_pion_fromK0","sameMother_pion_fromK0");
    pair_sameMother_pion_fromK.SetLegPDGs(11,-11);
    pair_sameMother_pion_fromK.SetCheckBothChargesLegs(kTRUE,kTRUE);
    pair_sameMother_pion_fromK.SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
    //mother
    pair_sameMother_pion_fromK.SetMothersRelation(AliDielectronSignalMC::kSame);
    pair_sameMother_pion_fromK.SetMotherPDGs(111,111); //
    // grand mother
    pair_sameMother_pion_fromK.SetGrandMotherPDGs(310,310);

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
    //    task->AddPairMCSignal(pair_sameMother_pion);
    //    task->AddPairMCSignal(pair_sameMother_pion_anySource);
    //    task->AddPairMCSignal(pair_sameMother_pion_fromK);
    // task->AddPairMCSignal(pair_sameMother_eta);
    // task->AddPairMCSignal(pair_sameMother_etaP);
    // task->AddPairMCSignal(pair_sameMother_rho);
    // task->AddPairMCSignal(pair_sameMother_omega);
    // task->AddPairMCSignal(pair_sameMother_phi);
    //    task->AddPairMCSignal(pair_sameMother_jpsi);
    // task->AddPairMCSignal(pair_sameMother_CharmedMesonsWithSameMother);
    // task->AddPairMCSignal(pair_conversion);
}
