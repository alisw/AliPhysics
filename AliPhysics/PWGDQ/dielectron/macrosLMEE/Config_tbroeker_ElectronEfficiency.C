TString names("cut01;cut02;cut03;cut04;cut05;cut06;cut07;cut08;cut09;cut10;cut11;cut12;cut13;cut14;cut15;cut16;cut17;cut18;cut19;cut20;cut21;cut22;cut23");

TObjArray*  arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();
Int_t       selectedCentrality=-1; // not yet implemented
// the following settings must be initialized each time SetupTrackCutsAndSettings() is called.
Int_t       selectedPID;
Bool_t      isPrefilterCutset;
Double_t    rejCutMee;
Double_t    rejCutTheta;
Double_t    rejCutPhiV;
// -----

//________________________________________________________________
// binning of 3D output histograms
// eta bins
const Double_t EtaMin   = -1.;
const Double_t EtaMax   =  1.;
const Int_t    nBinsEta = 20; //flexible to rebin
// phi bins
const Double_t PhiMin   = 0.;
const Double_t PhiMax   = 6.2832;
const Int_t    nBinsPhi = 40; //flexible to rebin
const Double_t PtBins[] = {
  0.00,0.10,0.12,0.14,0.16,0.18,0.20,0.21,0.22,0.23,0.24,0.25,
  0.26,0.27,0.28,0.29,0.30,0.32,0.34,0.36,0.38,0.40,0.43,0.46,
  0.49,0.52,0.55,0.60,0.65,0.70,0.75,0.80,0.90,1.00,1.10,1.20,
  1.40,1.60,1.80,2.00,2.40,2.80,3.20,3.70,4.50,6.00,8.00,12.0
};

Bool_t bUseRelPResolution = kTRUE;
Bool_t bUseEtaResolution = kFALSE;
Bool_t CalcEfficiencyRec = kFALSE;
Bool_t CalcEfficiencyPoslabel = kFALSE;
Bool_t CalcResolution = kTRUE;
Bool_t MakeResolutionSparse = kTRUE;
Bool_t doPairing = kFALSE;
// resolution binnings
Int_t    NbinsMom        = 2000;
Double_t MomMin          = 0.;
Double_t MomMax          = 10.;
Int_t    NbinsDeltaMom   = 1001;
Double_t DeltaMomMin     = -9.005;
Double_t DeltaMomMax     =  1.005;
Int_t    NbinsRelMom     =  1201;
Double_t RelMomMin       = -0.0005;
Double_t RelMomMax       =  1.2005;
Int_t    NbinsDeltaEta   = 1001;
Double_t DeltaEtaMin     = -0.5005;
Double_t DeltaEtaMax     =  0.5005;
Int_t    NbinsDeltaTheta = 1001;
Double_t DeltaThetaMin   = -0.5005;
Double_t DeltaThetaMax   =  0.5005;
Int_t    NbinsDeltaPhi   = 601;
Double_t DeltaPhiMin     = -0.3005;
Double_t DeltaPhiMax     =  0.3005;
Int_t    NbinsDeltaAngle = 401;
Double_t DeltaAngleMin   = -0.2005;
Double_t DeltaAngleMax   =  0.6005;

// mee bins
  const Double_t MeeBins[] = { 0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
                               0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,
                               0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,
                               0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,
                               0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,
                               0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,
                               0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,
                               0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,
                               0.80,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,
                               0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,
                               1.00,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,
                               1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,
                               2.20,2.40,2.60,2.80,3.00,3.01,3.02,3.03,3.04,3.05,
                               3.06,3.07,3.08,3.09,3.10,3.30,3.50,4.00 };
  const Int_t nBinsMee = ( sizeof(MeeBins) / sizeof(MeeBins[0]) )-1;                               

// ptee bins
  const Double_t PteeBins[] = { 0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,
                                0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,
                                0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,
                                0.30,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,
                                0.40,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,
                                0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,
                                1.00,1.05,1.10,1.15,1.20,1.25,1.30,1.35,1.40,1.45,
                                1.50,1.55,1.60,1.65,1.70,1.75,1.80,1.85,1.90,1.95,
                                2.00,2.05,2.10,2.15,2.20,2.25,2.30,2.35,2.40,2.45,
                                2.50,2.60,2.70,2.80,2.90,3.00,3.10,3.20,3.30,3.40,
                                3.50,3.60,3.70,3.80,3.90,4.00,4.10,4.20,4.30,4.40,
                                4.50,5.00,5.50,6.00,6.50,7.00 };
  const Int_t nBinsPtee = ( sizeof(PteeBins) / sizeof(PteeBins[0]) )-1;

// run dependency (currently only "TPC_dEdx_P_run")
// run string must be sorted in increasing order!
//AOD_115_goodPID //TString sRuns("167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593");
//MC LHC12a17g_fix at GSI (5.3.2014) (one more run than 17h)
const TString sRuns("167915, 167920, 167985, 167987, 168069, 168076, 168105, 168107, 168108, 168115, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169965, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593");
//
// ^^^^^^^^^^ [/end binning histograms] ^^^^^^^^^^

//________________________________________________________________
// specify if track tree shall be filled and written to file (only recommended for small checks!)
const Bool_t    writeTree = kFALSE;
// specify for which "cutDefinition" the support histos should be filled!
const Int_t     supportedCutInstance = 1;
//
//________________________________________________________________
// settings which are identical for all configs that run together
// event cuts
const Bool_t    reqVertex = kTRUE;
const Double_t  vertexZcut = 10.;
const Double_t  CentMin =  -2.;
const Double_t  CentMax = 102.;
// MC cuts
const Double_t  EtaMinGEN = -1.;    // make sure to be within 3D histogram binning (EtaMin, EtaMax, PtBins[]).
const Double_t  EtaMaxGEN =  1.;
const Double_t  PtMinGEN  =  0.100; // 100 MeV as absolute lower limit for any setting.
const Double_t  PtMaxGEN  =  20.;    // 8 GeV is current upper limit of PtBins[]. Dont want overflow bin filled.

const Bool_t    CutInjectedSignals = kFALSE;
const UInt_t    NminEleInEventForRej = 2;
// ^^^^^^^^^^ [/end common settings] ^^^^^^^^^^

//________________________________________________________________
//
// Strategy:  One cutDefinition for analysis tracking&PID efficiency (as usual).
//            Optional, separate cutDefinition for prefilter efficiencies: it also produces the usual tracking&PID efficiency
//            (but of course for the specified prefilter track sample, so mainly for convenience and curiosity),
//            and then additionally the pair rejection efficiency, using random rejection of pions with the selected electrons.
//________________________________________________________________


// TODO: implement this:
//
//AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
//noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
//die->GetTrackFilter().AddCuts(noconv);

//________________________________________________________________
void SetupMCSignals(AliAnalysisTaskElectronEfficiency* task)
{
  AliDielectronSignalMC* eleFinalState = new AliDielectronSignalMC("eleFinalState","eleFinalState");
  eleFinalState->SetFillPureMCStep(kFALSE);
  eleFinalState->SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalState->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
  eleFinalState->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);//equiv. to IsPrimary();
  task->AddSignalMC(eleFinalState);
}

//________________________________________________________________
AliAnalysisCuts* SetupEventCuts()
{
  // event cuts are identical for all analysis 'cutDefinition's that run together!
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
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

//________________________________________________________________
AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t cutDefinition)
{
  std::cout << "SetupTrackCutsAndSettings()" <<std::endl;
  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!
  // do not change these initial values!
  selectedPID=-1;
  isPrefilterCutset=kFALSE;
  rejCutMee=-1;
  rejCutTheta=-1;
  rejCutPhiV=3.2; // relevant are values below pi, so initialization to 3.2 means disabling.
  
  // -----
  // produce analysis filter by using functions in this config:
  // -----
  
  if (cutDefinition==100) {
    // kinematic cuts for the legs during pair efficiency determination:
    AliDielectronVarCuts *kineCuts = new AliDielectronVarCuts("kineCuts","kineCuts");  
    kineCuts->AddCut(AliDielectronVarManager::kEta, -0.8, 0.8);
    kineCuts->AddCut(AliDielectronVarManager::kPt,   0.2, 1000.);
    anaFilter->AddCuts( kineCuts );
    return anaFilter; // return here because we dont want any other cuts.
  }
  if(cutDefinition == 101){
    rejCutMee=-1;
    rejCutTheta=-1.;/*50.e-3*/;
    rejCutPhiV=3.2;
    return 0x0;
  }
  
  AliDielectronVarCuts *sharedClustersCut(0x0);
  if(cutDefinition > -1000){ 
    Printf("cutting on number of shared clusters in ITS! Max: 0");
    sharedClustersCut  = new AliDielectronVarCuts("sharedClustersCut","sharedClustersCut");
    sharedClustersCut->AddCut(AliDielectronVarManager::kNclsSITS, -0.5, 0.5);
  }
  
  AliDielectronV0Cuts *noconv(0x0);
  if(cutDefinition > -1000){
    std::cout << "Adding V0 conversion cut!" <<std::endl;
    noconv = new AliDielectronV0Cuts("IsGamma","IsGamma");
    // which V0 finder you want to use
    noconv->SetV0finder(AliDielectronV0Cuts::kAll);  // kAll(default), kOffline or kOnTheFly
    // add some pdg codes (they are used then by the KF package and important for gamma conversions)
    noconv->SetPdgCodes(22,11,11); // mother, daughter1 and 2
    // add default PID cuts (defined in AliDielectronPID)
    // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
    noconv->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
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
  }
  
  anaFilter->AddCuts( SetupTrackCuts(cutDefinition) );
  anaFilter->AddCuts( SetupPIDcuts(cutDefinition) );
  if(noconv)            anaFilter->AddCuts( noconv );
  if(sharedClustersCut) anaFilter->AddCuts( sharedClustersCut );
  std::cout << "...cuts added!" <<std::endl; 
  
  return anaFilter;
}


//________________________________________________________________
Int_t SetupPrefilterPairCuts(Int_t cutDefinition)
{
  std::cout << "SetupPrefilterPairCuts()" <<std::endl;
  
  if (selectedPID > -1) { // case with LMcutlib
    LMEECutLib* LMcutlib = new LMEECutLib();
    
    AliDielectronCutGroup* cgPairCutsPre = (AliDielectronCutGroup*) LMcutlib->GetPairCutsPre(selectedPID);
    if(!cgPairCutsPre) {
      std::cout << "WARNING: no Prefilter PairCuts given (bad cutgroup)!" << std::endl;
      return -1;
    }
    cgPairCutsPre->Print();
    Double_t dummy=-1;
    //
    // using some new (Feb 2015) functions:
    //   AliDielectronCutGroup::GetNCuts() and ::GetCut()
    //   AliDielectronVarCuts::IsCutOnVariableX() and ::GetCutLimits()
    //
    // this is NOT failsafe!!! needs a cutgroup as top level object and one or more varcuts inside, which have internal cuts on the needed variables.
    for (Int_t iCutGroupCut=0; iCutGroupCut<cgPairCutsPre->GetNCuts(); iCutGroupCut++) {
      AliDielectronVarCuts* varcuti = (AliDielectronVarCuts*) cgPairCutsPre->GetCut(iCutGroupCut);
      
      for (Int_t iCut=0; iCut<varcuti->GetNCuts(); iCut++) {
        if ( varcuti->IsCutOnVariableX(iCut, AliDielectronVarManager::kM) ) {
          if (rejCutMee>-1) { std::cout << "WARNING: rejCutMee was defined two times!" << std::endl; return -1; } // should take the stronger cut in that case...
          // fill rejCutMee
          varcuti->GetCutLimits(iCut, dummy, rejCutMee);
        }
        if ( varcuti->IsCutOnVariableX(iCut, AliDielectronVarManager::kOpeningAngle) ) {
          if (rejCutTheta>-1) { std::cout << "WARNING: rejCutTheta was defined two times!" << std::endl; return -1; }
          // fill rejCutTheta
          varcuti->GetCutLimits(iCut, dummy, rejCutTheta);
        }
        if ( varcuti->IsCutOnVariableX(iCut, AliDielectronVarManager::kPhivPair) ) {
          if (rejCutPhiV<3.14159) { std::cout << "WARNING: rejCutPhiV was defined two times!" << std::endl; return -1; }
          // fill rejCutPhiV
          varcuti->GetCutLimits(iCut, dummy, rejCutPhiV);
        }
      }
    }
  }
  else { // case without LMcutlib...
    switch (cutDefinition) {
      case 0:
        //rejCutMee=...;
        //rejCutTheta=...;
        //break;
      default:
        std::cout << "WARNING: no Prefilter PairCuts given!" << std::endl;
        return -1;
    }
  }
  
  return 1;
}    


//________________________________________________________________
AliAnalysisCuts* SetupTrackCuts(Int_t cutDefinition)
{
  std::cout << "SetupTrackCuts()" <<std::endl;
  //AliAnalysisCuts* trackCuts=0x0;

  AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts();

  //global
//  fesdTrackCuts->SetPtRange( 0.2 , 100. );
//  fesdTrackCuts->SetEtaRange( -0.8 , 0.8 );
  fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fesdTrackCuts->SetDCAToVertex2D(kFALSE);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.);
  fesdTrackCuts->SetMaxDCAToVertexXY(1.);
 
  fesdTrackCuts->SetRequireTPCRefit(kTRUE);
  fesdTrackCuts->SetRequireITSRefit(kTRUE);

  // resolution cuts
  if(cutDefinition == -1){
//    fesdTrackCuts->SetPtRange( 0.2 , 100. );
//    fesdTrackCuts->SetEtaRange( -0.8 , 0.8 );
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 0){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 1){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(130);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 2){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 3){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(120);
    fesdTrackCuts->SetMinNCrossedRowsTPC(130);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 4){
    fesdTrackCuts->SetMinNClustersITS(6);
    fesdTrackCuts->SetMaxChi2PerClusterITS(2.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 5){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 6){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 7){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(120);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 8){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(130);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(2.5);
  }
  if(cutDefinition == 9){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(120);
    fesdTrackCuts->SetMinNCrossedRowsTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 10){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 11){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(2.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 12){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(120);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 13){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(120);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 14){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(2.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(120);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 15){
    fesdTrackCuts->SetMinNClustersITS(6);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 16){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 17){
    fesdTrackCuts->SetMinNClustersITS(6);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.7);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 18){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);    
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);
  }
  if(cutDefinition == 19){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 20){
    fesdTrackCuts->SetMinNClustersITS(3);
    fesdTrackCuts->SetMaxChi2PerClusterITS(2.5);
    fesdTrackCuts->SetMinNClustersTPC(100);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
  }
  if(cutDefinition == 21){
    fesdTrackCuts->SetMinNClustersITS(4);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3);  
  }
  if(cutDefinition == 22){
    fesdTrackCuts->SetMinNClustersITS(5);
    fesdTrackCuts->SetMaxChi2PerClusterITS(3.5);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersTPC(120);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(3); 
  }
  return fesdTrackCuts;
}

//________________________________________________________________
AliAnalysisCuts* SetupPIDcuts(Int_t cutDefinition)
{
  std::cout << "SetupPIDcuts()" <<std::endl;
  AliAnalysisCuts* pidCuts=0x0;
  
  AliDielectronPID *pid = new AliDielectronPID();
  
  // resolution cuts
  if(cutDefinition == -1){
//    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
//    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -5.,5. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 0){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }	
  if(cutDefinition == 1){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 2){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 3){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -4. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3.5,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 4){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3.5,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 5){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }	
  if(cutDefinition == 6){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 7){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 8){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3.5,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -4. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 9){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3.5,0. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 10){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 11){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 12){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 13){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 14){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,2. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 15){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 16){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -2. ,2. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 17){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 18){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 19){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 20){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,3. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 21){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -4. ,1.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5,3. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  if(cutDefinition == 22){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,    -100. ,4. ,0.0, 100., kTRUE ,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3. ,3. ,0.0, 100., kFALSE,AliDielectronPID::kIfAvailable,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kITS,AliPID::kElectron,  -3. ,0.5,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -0.5,4. ,0.0, 100., kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  pidCuts = pid;   
  return pidCuts;
}



