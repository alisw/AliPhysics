TString names=("pt200_cut1");

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
const Int_t    nBinsEta = 10; //flexible to rebin
// phi bins
const Double_t PhiMin   = 0.;
const Double_t PhiMax   = 6.2832;
const Int_t    nBinsPhi = 60; //flexible to rebin
const Double_t PtBins[] = {
  0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,10.0,20.0
};

// Bool_t bUseRelPResolution = kTRUE; //not used
Bool_t bUseEtaResolution = kTRUE; // use eta or theta resolution?
Bool_t CalcEfficiencyRec = kTRUE;
Bool_t CalcEfficiencyPoslabel = kFALSE;
Bool_t CalcResolution = kTRUE;
Bool_t MakeResolutionSparse = kFALSE;
Bool_t doPairing = kTRUE;

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


//// mee bins
//const Double_t MeeMin    = 0.;
//const Double_t MeeMax    = 5.;
//const Int_t    nBinsMee  = 500;

const Double_t MeeBins[] = { 0., 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.65, 0.688, 0.725,
    0.75, 0.775, 0.8, 0.85, 0.95, 0.975, 1.0, 1.025, 1.05, 1.125, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 2.85,
    2.95, 3.05, 3.1, 3.15, 3.3, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0 };

const Int_t nBinsMee = ( sizeof(MeeBins) / sizeof(MeeBins[0]) )-1;

//// ptee bins
//const Double_t PteeMin   = 0.;
//const Double_t PteeMax   = 6.;
//const Int_t    nBinsPtee = 600;

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
    4.50,5.00,5.50,6.00,6.50,7.00,10.00,15.00,20.00,100.00};
const Int_t nBinsPtee = ( sizeof(PteeBins) / sizeof(PteeBins[0]) )-1;

// in increasing order
const TString sRuns("258919, 258923, 258962, 258964, 259088, 259090, 259091, 259096, 259099, 259117, 259118, 259162, 259164, 259204, 259261, 259263, 259264, 259269, 259270, 259271, 259272, 259273, 259274, 259302, 259303, 259305, 259307, 259334, 259336, 259339, 259340, 259341, 259342, 259378, 259382, 259388, 259389, 259394, 259395, 259396, 259473, 259477, 259649, 259650, 259668, 259697, 259700, 259703, 259704, 259705, 259711, 259713, 259747, 259748, 259750, 259751, 259752, 259756, 259781, 259788, 259789, 259822, 259841, 259842, 259860, 259866, 259867, 259868, 259888");

//
// ^^^^^^^^^^ [/end binning histograms] ^^^^^^^^^^

//________________________________________________________________
// specify if track tree shall be filled and written to file (only recommended for small checks!)
const Bool_t    writeTree = kFALSE;
// specify for which "cutInstance" the support histos should be filled!
const Int_t     supportedCutInstance = 0;
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
const Double_t  PtMaxGEN  =  50.;

const Bool_t    CutInjectedSignals = kFALSE;
const UInt_t    NminEleInEventForRej = 2;
// ^^^^^^^^^^ [/end common settings] ^^^^^^^^^^

//________________________________________________________________
//
// Strategy:  One cutInstance for analysis tracking&PID efficiency (as usual).
//            Optional, separate cutInstance for prefilter efficiencies: it also produces the usual tracking&PID efficiency
//            (but of course for the specified prefilter track sample, so mainly for convenience and curiosity),
//            and then additionally the pair rejection efficiency, using random rejection of pions with the selected electrons.
//________________________________________________________________


// TODO: implement this:
//
//AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
//noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
//die->GetTrackFilter().AddCuts(noconv);

//________________________________________________________________
//void SetupMCSignals(AliAnalysisTaskElectronEfficiency* task)
//{
//  AliDielectronSignalMC* eleFinalState = new AliDielectronSignalMC("eleFinalState","eleFinalState");
//  eleFinalState->SetFillPureMCStep(kFALSE);
//  eleFinalState->SetLegPDGs(11,1);//dummy second leg (never MCtrue)
//  eleFinalState->SetCheckBothChargesLegs(kTRUE,kTRUE);
//  eleFinalState->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);
//  eleFinalState->SetMotherSources(AliDielectronSignalMC::kPrimary, AliDielectronSignalMC::kPrimary);//equiv. to IsPrimary();
//  task->AddSignalMC(eleFinalState);
//}

//________________________________________________________________
AliAnalysisCuts* SetupEventCuts()
{
  // event cuts are identical for all analysis 'cutInstance's that run together!
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  return eventCuts;
}

//________________________________________________________________
AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t cutInstance)
{
  std::cout << "SetupTrackCutsAndSettings()" <<std::endl;
  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!
  // do not change these initial values!
  selectedPID=-1;
  isPrefilterCutset=kFALSE;
  rejCutMee=-1;
  rejCutTheta=-1;
  rejCutPhiV=3.2; // relevant are values below pi, so initialization to 3.2 means disabling.
    
    if (cutInstance==100) {
        // kinematic cuts for the legs during pair efficiency determination:
        AliDielectronVarCuts *kineCuts = new AliDielectronVarCuts("kineCuts","kineCuts");
        kineCuts->AddCut(AliDielectronVarManager::kEta, -0.8, 0.8);
        kineCuts->AddCut(AliDielectronVarManager::kPt,   0.2, 1e30);
        anaFilter->AddCuts( kineCuts );
        return anaFilter; // return here because we dont want any other cuts.
    }
    
    if(cutInstance == 101){
        rejCutMee=-1;
        rejCutTheta=-1.;/*50.e-3*/;
        rejCutPhiV=3.2;
        return 0x0;
    }
  
  // -----
  // produce analysis filter by using functions in this config:
  // -----
  
  anaFilter->AddCuts( SetupTrackCuts(cutInstance) );
  anaFilter->AddCuts( SetupPIDcuts(cutInstance) );
  
    
//  AliDielectronV0Cuts *noconv = new AliDielectronV0Cuts("IsGamma","IsGamma");
//  // which V0 finder you want to use
//  noconv->SetV0finder(AliDielectronV0Cuts::kAll);  // kAll(default), kOffline or kOnTheFly
//  // add some pdg codes (they are used then by the KF package and important for gamma conversions)
//  noconv->SetPdgCodes(22,11,11); // mother, daughter1 and 2
//  // add default PID cuts (defined in AliDielectronPID)
//  // requirement can be set to at least one(kAny) of the tracks or to both(kBoth)
//  //noconv->SetDefaultPID(16, AliDielectronV0Cuts::kAny);
//  // add the pair cuts for V0 candidates
//  noconv->AddCut(AliDielectronVarManager::kCosPointingAngle, TMath::Cos(0.02),   1.00, kFALSE);
//  noconv->AddCut(AliDielectronVarManager::kChi2NDF,                       0.0,  10.00, kFALSE);
//  noconv->AddCut(AliDielectronVarManager::kLegDist,                       0.0,   0.25, kFALSE);
//  noconv->AddCut(AliDielectronVarManager::kR,                             3.0,  90.00, kFALSE);
//  noconv->AddCut(AliDielectronVarManager::kPsiPair,                       0.0,   0.05, kFALSE);
//  noconv->AddCut(AliDielectronVarManager::kM,                             0.0,   0.10, kFALSE);
//  noconv->AddCut(AliDielectronVarManager::kArmPt,                         0.0,   0.05, kFALSE);
//  // selection or rejection of V0 tracks
//  if(cutInstance < 20)  
//  anaFilter->AddCuts( noconv );
  
  std::cout << "...cuts added!" <<std::endl; 
  
//  std::cout << "__________ anaFilter->GetCuts()->Print() __________ cutInstance = " << cutInstance <<std::endl;
//  anaFilter->GetCuts()->Print();
//  std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<std::endl;
//  std::cout << "__________ anaFilter->GetCuts()->Dump() __________ cutInstance = " << cutInstance <<std::endl;
//  anaFilter->GetCuts()->Dump();
//  std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<std::endl;
  return anaFilter;
}


// prefilter cuts are not used at the moment
//________________________________________________________________
Int_t SetupPrefilterPairCuts(Int_t cutInstance)
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
    switch (cutInstance) {
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
AliAnalysisCuts* SetupTrackCuts(Int_t cutInstance)
{
  std::cout << "SetupTrackCuts()" <<std::endl;
  //AliAnalysisCuts* trackCuts=0x0;
    
    AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts();
    
    //primary selection
    fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fesdTrackCuts->SetDCAToVertex2D(kFALSE);
    fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    fesdTrackCuts->SetMaxDCAToVertexZ(3.0);
    fesdTrackCuts->SetMaxDCAToVertexXY(1.0);
    
    //TPC
    fesdTrackCuts->SetRequireTPCRefit(kTRUE);
    fesdTrackCuts->SetMinNCrossedRowsTPC(100);
    fesdTrackCuts->SetMinNClustersTPC(80);
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4.0);
    fesdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
    
    //ITS
    fesdTrackCuts->SetRequireITSRefit(kTRUE);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fesdTrackCuts->SetMinNClustersITS(3);
    fesdTrackCuts->SetMaxChi2PerClusterITS(4.5);
    
  if(cutInstance == 0){
      
      // pT and eta
      fesdTrackCuts->SetPtRange(0.2, 1e30);
      fesdTrackCuts->SetEtaRange(-0.8, 0.8);
      
  }

  return fesdTrackCuts;
}

//________________________________________________________________
AliAnalysisCuts* SetupPIDcuts(Int_t cutInstance)
{
  std::cout << "SetupPIDcuts()" <<std::endl;
  AliAnalysisCuts* pidCuts=0x0;
  
  AliDielectronPID *pid = new AliDielectronPID();
  
  // loose PID for resolution cuts to increase statistics
  // suppose that we have same responce here for all particles
  if(cutInstance == -1){
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -5.,5. ,0.0, 1e30, kFALSE,AliDielectronPID::kRequire    ,AliDielectronVarManager::kPt);
  }
  
  if(cutInstance == 0){
      
      //TPC
      pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron,   -3. ,  3.,   0.0, 1e30,  kFALSE,  AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
      pid->AddCut(AliDielectronPID::kTPC, AliPID::kPion,       -100. ,  4.,   0.0, 1e30,  kTRUE,   AliDielectronPID::kRequire,     AliDielectronVarManager::kPt);
      
      // TOF
      pid->AddCut(AliDielectronPID::kTOF, AliPID::kElectron,     -3. ,  3.,   0.4, 1e30,    kFALSE,  AliDielectronPID::kRequire, AliDielectronVarManager::kP);
      
  }
    
  pidCuts = pid;   
  return pidCuts;
}



