TString names=("cut_resolution");
//TString names=("cut01_SPD1_PID1;cut02_SPD2_PID2;cut03_SPD3_PID3;cut04_SPD4_PID4;cut05_SPD5_PID1;cut06_SPD6_PID6;cut07_SPD7_PID7;cut08_SPD8_PID8;cut09_SPD9_PID9;cut10_SPD1_PID10;cut11_SPD11_PID11;cut12_SPD11_PID12;cut13_SPD6_PID13;cut14_SPDorSDD14_PID14;cut15_SPDorSDD15_PID13;cut16_SPDorSDD14_PID16;cut17_SPDorSDD17_PID16;cut18_SPDorSDD15_PID1;cut19_SPDorSDD17_PID6;cut20_SPDorSDD17_PID10");
TObjArray*  arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();
//________________________________________________________________
//
// Strategy:  One cutInstance for analysis tracking&PID efficiency (as usual).
//            Optional, separate cutInstance for prefilter efficiencies: it also produces the usual tracking&PID efficiency
//            (but of course for the specified prefilter track sample, so mainly for convenience and curiosity),
//            and then additionally the pair rejection efficiency, using random rejection of "testparticles" (pions) with the selected electrons.
//            By now this is better done by the AliAnalysisTaskRandomRejection, which runs over real data.
//________________________________________________________________
//
//________________________________________________________________
//
// WARNING:   The post-PID-corrections are identical for all cutInstances in this Config.
//            --> Define only cutInstances which have the same kinematic cuts!
//________________________________________________________________
//
//________________________________________________________________
// main task settings
// fill resolutions for one cutInstance (step 1).
const Bool_t CalcResolution   = kTRUE;
Bool_t CalcEfficiencyGen      = kTRUE; // can be turned off to save memory (especially in step 2).
// use previously extracted resolutions (step 2).
TString resolutionfile = "resolution_PbPb2011_CENTRALITY_deltaXvsP.root";
Bool_t CalcEfficiencyRec      = kFALSE;  // use given resolution file to smear the kinematics.
Bool_t bUseRelPResolution     = kTRUE;  // specify if the file contains a relative or an absolute momentum resolution array.
Bool_t bUseEtaResolution      = kFALSE; // kFALSE means using theta instead of eta.
// determine efficiency from only positive label tracks (in addition to using all labels).
Bool_t CalcEfficiencyPoslabel = kFALSE;
// determine pair efficiency for all cutInstances. (Consider high combinatorics if not only MC-true electrons are selected.)
const Bool_t doPairing = kFALSE;
// specify for which "cutInstance" the support histos should be filled!
const Int_t     supportedCutInstance = 0;
// specify if track tree shall be filled and written to file (only recommended for small checks!)
const Bool_t    writeTree = kFALSE;
// activate UsePhysicsSelection and SetTriggerMask for MC (may be needed for new MC productions according to Mahmut)
//const Bool_t    forcePhysSelAndTrigMask = kFALSE; // default kFALSE
// ^^^^^^^^^^ [/end main task settings] ^^^^^^^^^^
//
//________________________________________________________________
// binning of output histograms
// eta bins
const Double_t EtaMin   = -1.;
const Double_t EtaMax   =  1.;
const Int_t    nBinsEta = 100; //flexible to rebin
// phi bins
const Double_t PhiMin   = 0.;
const Double_t PhiMax   = 6.2832;
const Int_t    nBinsPhi = 60; //flexible to rebin
const Double_t PtBins[] = {
  0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,8.0
};
// resolution binning
Int_t NbinsDeltaMom    =1200;
Double_t DeltaMomMin   =-10.0;
Double_t DeltaMomMax   =  2.0;
Int_t NbinsRelMom      = 400;
Double_t RelMomMin     =  0.0;
Double_t RelMomMax     =  2.0;
Int_t NbinsDeltaEta    = 200;
Double_t DeltaEtaMin   = -0.4;
Double_t DeltaEtaMax   =  0.4;
Int_t NbinsDeltaTheta  = 200;
Double_t DeltaThetaMin = -0.4;
Double_t DeltaThetaMax =  0.4;
Int_t NbinsDeltaPhi    = 200;
Double_t DeltaPhiMin   = -0.4;
Double_t DeltaPhiMax   =  0.4;

// mee bins
const Double_t MeeMin    = 0.;
const Double_t MeeMax    = 5.;
const Int_t    nBinsMee  = 250;
// ptee bins
const Double_t PteeMin   = 0.;
const Double_t PteeMax   = 5.;
const Int_t    nBinsPtee = 100;
// run dependency (currently only "TPC_dEdx_P_run")
// run string must be sorted in increasing order!
//AOD_115_goodPID //TString sRuns("167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593");
//MC LHC12a17g_fix at GSI (5.3.2014) (one more run than 17h)
const TString sRuns("167915, 167920, 167985, 167987, 168069, 168076, 168105, 168107, 168108, 168115, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 169923, 169965, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593");
// ^^^^^^^^^^ [/end binning histograms] ^^^^^^^^^^
//
//________________________________________________________________
// settings which are identical for all configs that run together
// event cuts done via 'SetupEventCuts()'
// centrality cuts done in AddTask
// MC cuts
const Double_t  EtaMinGEN = -1.;    // make sure to be within 3D histogram binning (EtaMin, EtaMax, PtBins[]).
const Double_t  EtaMaxGEN =  1.;
const Double_t  PtMinGEN  =  0.100; // 100 MeV as absolute lower limit for any setting.
const Double_t  PtMaxGEN  =  8.;    // 8 GeV is current upper limit of PtBins[]. Dont want overflow bin filled.

const UInt_t    NminEleInEventForRej = 2;
// ^^^^^^^^^^ [/end common settings] ^^^^^^^^^^
//________________________________________________________________
// the following settings must be initialized each time SetupTrackCutsAndSettings() is called. (do not give values here!)
Int_t       selectedPairCutsPre;
Bool_t      isPrefilterCutset;
AliAnalysisFilter *anaFilterExtra;
Double_t    rejCutMee;
Double_t    rejCutTheta;
Double_t    rejCutPhiV;
// ^^^^^^^^^^ [/end automatic settings] ^^^^^^^^^^

// TODO: implement this:
//
//AliDielectronTrackCuts *noconv=new AliDielectronTrackCuts("noConv","conversion tagging");
//noconv->SetV0DaughterCut(AliPID::kElectron,kTRUE);
//die->GetTrackFilter().AddCuts(noconv);


//________________________________________________________________
AliAnalysisCuts* SetupEventCuts(Bool_t isESD=kTRUE)
{
  // event cuts are identical for all analysis 'cutInstance's that run together!
  //
  AliDielectronEventCuts *eventCuts=new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");
  eventCuts->SetRequireVertex();
  eventCuts->SetMinVtxContributors(1);
  eventCuts->SetVertexZ(-10.,10.);
  //not sure if this can be used, probably not:
  // Patrick:
  if (!isESD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  // Mahmut:
  //  eventCuts->SetVertexType(AliDielectronEventCuts::kVtxTracksOrSPD);
  //  eventCuts->SetRequireV0and();
  return eventCuts;
}


//________________________________________________________________
void SetupITSSigmaEleCorrection(AliAnalysisTaskElectronEfficiency* task)
{
  LMEECutLib* LMcutlib = new LMEECutLib();
  //LMcutlib->SetITSSigmaEleCorrectionMC(task, AliDielectronVarManager::kNacc, AliDielectronVarManager::kEta);
  LMcutlib->SetITSSigmaEleCorrectionMC(task, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
  return;
}

//________________________________________________________________
void SetupTPCSigmaEleCorrection(AliAnalysisTaskElectronEfficiency* task)
{
  LMEECutLib* LMcutlib = new LMEECutLib();
  //LMcutlib->SetTPCSigmaEleCorrectionMC(task, AliDielectronVarManager::kNacc, AliDielectronVarManager::kEta);
  LMcutlib->SetTPCSigmaEleCorrectionMC(task, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
  return;
}

//________________________________________________________________
void SetupMCSignals(AliAnalysisTaskElectronEfficiency* task)
{
  LMEECutLib* LMcutlib = new LMEECutLib();
  LMcutlib->AddMCSignals(task, -1);
}


//________________________________________________________________
AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t cutInstance, Bool_t isESD=kTRUE)
{
  std::cout << "SetupTrackCutsAndSettings( cutInstance = " << cutInstance << " )" <<std::endl;
  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!
  // do not change these initial values!
  selectedPairCutsPre=-1;
  isPrefilterCutset=kFALSE;
  anaFilterExtra = new AliAnalysisFilter("anaFilterExtra","anaFilter");
  rejCutMee=-1;
  rejCutTheta=-1;
  rejCutPhiV=3.2; // relevant are values below pi, so initialization to 3.2 means disabling.

  // -----
  // produce analysis filter by using functions in this config:
  // -----
  int nCutsUsingConfigFunctions = 0; // set it manually!
  //
  if (cutInstance == -1) { // set negative when not needed.
    anaFilter->AddCuts( SetupTrackCuts(cutInstance) );
    anaFilter->AddCuts( SetupPIDcuts(cutInstance) );
    std::cout << "...cuts added!" <<std::endl;
    //nCutsUsingConfigFunctions++; this is logically not possible, set it manually above!
  }
  // -----
  // produce analysis filter by using LMEECutLib:
  // -----
  else // "else" is important to not add additional cuts to previous cutInstance by accident.
  {
    LMEECutLib* LMcutlib = new LMEECutLib();
    LMcutlib->SetIsESDTask(isESD);

    // --------------------------------------------------
    // common settings:
    // --------------------------------------------------
    //    LMcutlib->selectedCentrality  = LMEECutLib::kPbPb2011_00to10; // use AddTask for centrality cuts!
    // prefilter settings:
    //    LMcutlib->selectedPIDPre      = LMEECutLib::kPbPb2011PID_TPCITSif_2;
    //    LMcutlib->selectedQualityPre  = LMEECutLib::kPbPb2011TRK_FilterBit0;
    //    LMcutlib->selectedKineCutsPre = LMEECutLib::kKineCut_pt50_eta090;
    //    LMcutlib->selectedPairCutsPre = LMEECutLib::kPairCut_mee40_theta80;
    // ana settings:
    LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_pt200_eta080;
    //LMcutlib->selectedPairCutsAna = LMEECutLib::kPairCut_theta50;
    // ESD cuts:
    //LMcutlib->SetUseITScutsESD(kTRUE);
    //
    // --------------------------------------------------
    // specific settings for each cutset:
    // --------------------------------------------------
    /*** step 1 ***/
    if (cutInstance==0+nCutsUsingConfigFunctions) {
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_p50inf_eta150; // maximum acceptance
      LMcutlib->selectedPIDAna      = LMEECutLib::kPbPb2011PID_TPCITS_3;  // ITS+TPC, 100% efficiency
      LMcutlib->selectedQualityAna  = LMEECutLib::kCut16;                 // combined tracks SPDorSDD
    }
    else if (cutInstance==99) {
      // cuts for the resolution extraction and usage
      LMcutlib->selectedKineCutsAna = LMEECutLib::kKineCut_p50inf_eta150; // maximum acceptance
      LMcutlib->selectedPIDAna      = LMEECutLib::kPbPb2011PID_TPCITS_3;  // ITS+TPC, 100% efficiency
      LMcutlib->selectedQualityAna  = LMEECutLib::kCut16;                 // combined tracks SPDorSDD
    }
    /*** step 2 ***/
//    if (cutInstance==0+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut01;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut01;
//    }
//    else if (cutInstance==1+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut02;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut02;
//    }
//    else if (cutInstance==2+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut03;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut03;
//    }
//    else if (cutInstance==3+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut04;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut04;
//    }
//    else if (cutInstance==4+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut05;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut05;
//    }
//    else if (cutInstance==5+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut06;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut06;
//    }
//    else if (cutInstance==6+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut07;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut07;
//    }
//    else if (cutInstance==7+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut08;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut08;
//    }
//    else if (cutInstance==8+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut09;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut09;
//    }
//    else if (cutInstance==9+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut10;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut10;
//    }
//    else if (cutInstance==10+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut11;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut11;
//    }
//    else if (cutInstance==11+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut12;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut12;
//    }
//    else if (cutInstance==12+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut13;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut13;
//    }
//    else if (cutInstance==13+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut14;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut14;
//    }
//    else if (cutInstance==14+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut15;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut15;
//    }
//    else if (cutInstance==15+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut16;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut16;
//    }
//    else if (cutInstance==16+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut17;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut17;
//    }
//    else if (cutInstance==17+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut18;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut18;
//    }
//    else if (cutInstance==18+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut19;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut19;
//    }
//    else if (cutInstance==19+nCutsUsingConfigFunctions) {
//      LMcutlib->selectedPIDAna      = LMEECutLib::kCut20;
//      LMcutlib->selectedQualityAna  = LMEECutLib::kCut20;
//    }
    else if (cutInstance==101) {
      // pair cuts during pair efficiency determination:
      SetupPairCutsAna( LMcutlib->selectedPairCutsAna ); // class member is public...
      return 0x0; // return here because we dont want any other cuts.
    }
    else {
      cout << " =============================== " << endl;
      cout << " ==== INVALID CONFIGURATION ==== " << endl;
      cout << " cutInstance = " << cutInstance << endl;
      cout << " =============================== " << endl;
      return 0x0;
    }

    if (!isPrefilterCutset) {
      anaFilter->AddCuts( LMcutlib->GetTrackCutsAna() );
    }
    else { // cutInstance for prefilter efficiency determination:
      anaFilter->AddCuts( LMcutlib->GetTrackCutsPre() );
      // cuts for final analysis electrons
      anaFilterExtra->AddCuts( LMcutlib->GetTrackCutsAna() );
    }

    // export selectedPairCutsPre so that function SetupPrefilterPairCuts() can use it. (not really nice but...)
    selectedPairCutsPre = LMcutlib->selectedPairCutsPre;
  }
  // -----


  //  std::cout << "__________ anaFilter->GetCuts()->Print() __________ cutInstance = " << cutInstance <<std::endl;
  //  anaFilter->GetCuts()->Print();
  //  std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<std::endl;
  //  std::cout << "__________ anaFilter->GetCuts()->Dump() __________ cutInstance = " << cutInstance <<std::endl;
  //  anaFilter->GetCuts()->Dump();
  //  std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<std::endl;
  return anaFilter;
}


//________________________________________________________________
Int_t SetupPairCutsAna(Int_t selectedPairCutsAna)
{
  std::cout << "SetupPairCutsAna()" <<std::endl;
  selectedPairCutsPre = selectedPairCutsAna;
  return SetupPrefilterPairCuts(-1);
}

//________________________________________________________________
Int_t SetupPrefilterPairCuts(Int_t cutInstance)
{
  std::cout << "SetupPrefilterPairCuts()" <<std::endl;

  if (selectedPairCutsPre > -1) { // case with LMcutlib
    LMEECutLib* LMcutlib = new LMEECutLib();
    LMcutlib->selectedPairCutsPre = selectedPairCutsPre;
    //cout << " LMcutlib->selectedPairCutsPre = " << LMcutlib->selectedPairCutsPre << endl;
    Double_t dummy=-1;

    // this function is NOT failsafe!!!
    // the version commented out by /* ... */ needs a cutgroup as top level object and one or more varcuts inside, which have internal cuts on the needed variables.

    /*    AliDielectronCutGroup* cgPairCutsPre = (AliDielectronCutGroup*) LMcutlib->GetPairCutsPre();
     if(!cgPairCutsPre) {
     std::cout << "WARNING: no Prefilter PairCuts given (bad cutgroup)!" << std::endl;
     return -1;
     }
     cgPairCutsPre->Print();  */
    AliDielectronVarCuts* varcuti = (AliDielectronVarCuts*) LMcutlib->GetPairCutsPre();
    if(!varcuti) {
      std::cout << "WARNING: no Prefilter PairCuts given!" << std::endl;
      return -1;
    }
    varcuti->Print();
    //
    // using some new (Feb 2015) functions:
    //   AliDielectronCutGroup::GetNCuts() and ::GetCut()
    //   AliDielectronVarCuts::IsCutOnVariableX() and ::GetCutLimits()
    //

    /*  for (Int_t iCutGroupCut=0; iCutGroupCut<cgPairCutsPre->GetNCuts(); iCutGroupCut++) {
     AliDielectronVarCuts* varcuti = (AliDielectronVarCuts*) cgPairCutsPre->GetCut(iCutGroupCut);  */

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
        varcuti->GetCutLimits(iCut, rejCutPhiV, dummy);
      }
    }
    /*  }  */
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

  std::cout << "SetupPrefilterPairCuts() done!" <<std::endl;
  std::cout << "  rejCutMee   = " << rejCutMee <<std::endl;
  std::cout << "  rejCutTheta = " << rejCutTheta <<std::endl;
  std::cout << "  rejCutPhiV  = " << rejCutPhiV <<std::endl;

  return 1;
}


//________________________________________________________________
AliAnalysisCuts* SetupTrackCuts(Int_t cutInstance)
{
  std::cout << "SetupTrackCuts()" <<std::endl;
  //AliAnalysisCuts* trackCuts=0x0;

  if(cutInstance == 0) {
    // reproduce AOD filter bit 4:
    AliESDtrackCuts *esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    //esdTrackCuts->SetMaxDCAToVertexXY(2.4);
    //esdTrackCuts->SetMaxDCAToVertexZ(3.2);
    esdTrackCuts->SetDCAToVertex2D(kTRUE);
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    // additional or modified cuts: like in LMEECutLib::GetTrackCutsAna()
    esdTrackCuts->SetMaxDCAToVertexXY(1.);
    esdTrackCuts->SetMaxDCAToVertexZ(3.);
    esdTrackCuts->SetMinNClustersITS(4);
    esdTrackCuts->SetMinNCrossedRowsTPC(100); //default is 70 (in GetStandardITSTPCTrackCuts2011())
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8); //default is 0.8
    //
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);

    // kinematic cuts:
    esdTrackCuts->SetPtRange(   0.4 , 3.5 );
    esdTrackCuts->SetEtaRange( -0.9 , 0.9 );
    //    esdTrackCuts->SetPtRange(   0.399 , 3.499 );
    //    esdTrackCuts->SetEtaRange( -0.899 , 0.901 ); // on purpose to check if something gets mixed up...
  }

  return esdTrackCuts;
  //  trackCuts = fesdTrackCuts;
  //  trackCuts->Print();
  //  return trackCuts;
}

//________________________________________________________________
AliAnalysisCuts* SetupPIDcuts(Int_t cutInstance)
{
  std::cout << "SetupPIDcuts()" <<std::endl;
  AliAnalysisCuts* pidCuts=0x0;

  if(cutInstance == 0) {
    AliDielectronPID *pid = new AliDielectronPID("pidXtraPIn","pidXtraPIn");
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kPion,      -3.  ,3. ,0.0, 100., kTRUE);
    pid->AddCut(AliDielectronPID::kTPC,AliPID::kElectron,  -1.5 ,3. ,0.0, 100., kFALSE);
    pid->AddCut(AliDielectronPID::kTOF,AliPID::kElectron,  -3.  ,3. ,0.0, 1.7 , kFALSE);
  }

  pidCuts = pid;
  //pidCuts->Print();
  return pidCuts;
}
