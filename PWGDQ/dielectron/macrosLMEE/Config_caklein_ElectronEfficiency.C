// TString names=("resolutionCuts");
// TString names=("kPbPb2015_Pt400_tightTOFreq_MB;kPbPb2015_Pt400_tightTOFif_MB;kPbPb2015_Pt400_looseTOFif_MB");
TString names=("cut0;cut1;cut2;cut3;cut4;cut5;cut6;cut7;cut8;cut9");
// TString names=("cut10;cut11;cut12;cut13;cut14;cut15;cut16;cut17;cut18;cut19");

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
const Bool_t CalcResolution   = kFALSE;
// use previously extracted resolutions (step 2).
TString resolutionfile = "resolution_PbPb2015_CENTRALITY_deltaXvsP.root";
Bool_t CalcEfficiencyRec      = kTRUE;  // use given resolution file to smear the kinematics.
Bool_t bUseRelPResolution     = kTRUE;  // specify if the file contains a relative or an absolute momentum resolution array.
Bool_t bUseEtaResolution      = kFALSE; // kFALSE means using theta instead of eta.
// determine efficiency from only positive label tracks (in addition to using all labels).
Bool_t CalcEfficiencyPoslabel = kFALSE;
// determine pair efficiency for all cutInstances. (Consider high combinatorics if not only MC-true electrons are selected.)
const Bool_t doPairing = kFALSE;
// specify for which "cutInstance" the support histos should be filled!
const Int_t     supportedCutInstance = 1;
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
const Int_t    nBinsEta = 40; //flexible to rebin
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
Int_t NbinsDeltaMom    = 1200;
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
const Int_t    nBinsMee  = 500;
// ptee bins
const Double_t PteeMin   = 0.;
const Double_t PteeMax   = 8.;
const Int_t    nBinsPtee = 800;

// run string must be sorted in increasing order!
//AOD_115_goodPID //TString sRuns("167987, 167988, 168310, 168311, 168322, 168325, 168341, 168342, 168361, 168362, 168458, 168460, 168464, 168467, 168511, 168512, 168514, 168777, 168826, 168988, 168992, 169035, 169040, 169044, 169045, 169091, 169094, 169099, 169138, 169144, 169145, 169148, 169156, 169160, 169167, 169238, 169411, 169415, 169417, 169418, 169419, 169420, 169475, 169498, 169504, 169506, 169512, 169515, 169550, 169553, 169554, 169555, 169557, 169586, 169587, 169588, 169590, 169591, 169835, 169837, 169838, 169846, 169855, 169858, 169859, 170027, 170040, 170081, 170083, 170084, 170085, 170088, 170089, 170091, 170155, 170159, 170163, 170193, 170203, 170204, 170207, 170228, 170230, 170268, 170269, 170270, 170306, 170308, 170309, 170311, 170312, 170313, 170315, 170387, 170388, 170572, 170593");
const TString sRuns("246087");
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
  if (!isESD) eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD
  return eventCuts;
}


//________________________________________________________________
void SetupITSSigmaEleCorrection(AliAnalysisTaskElectronEfficiency* task)
{
  // LMEECutLib* LMcutlib = new LMEECutLib();
  // //LMcutlib->SetITSSigmaEleCorrectionMC(task, AliDielectronVarManager::kNacc, AliDielectronVarManager::kEta);
  // LMcutlib->SetITSSigmaEleCorrectionMC(task, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
  return;
}

//________________________________________________________________
void SetupTPCSigmaEleCorrection(AliAnalysisTaskElectronEfficiency* task)
{
  // LMEECutLib* LMcutlib = new LMEECutLib();
  // //LMcutlib->SetTPCSigmaEleCorrectionMC(task, AliDielectronVarManager::kNacc, AliDielectronVarManager::kEta);
  // LMcutlib->SetTPCSigmaEleCorrectionMC(task, AliDielectronVarManager::kP, AliDielectronVarManager::kEta);
  return;
}

//________________________________________________________________
void SetupMCSignals(AliAnalysisTaskElectronEfficiency* task){
  AliDielectronSignalMC* eleFinalState = new AliDielectronSignalMC("eleFinalState","eleFinalState");
  eleFinalState->SetFillPureMCStep(kFALSE);
  eleFinalState->SetLegPDGs(11,1);//dummy second leg (never MCtrue)
  eleFinalState->SetCheckBothChargesLegs(kTRUE,kTRUE);
  eleFinalState->SetLegSources(AliDielectronSignalMC::kFinalStateFromBGEvent, AliDielectronSignalMC::kFinalStateFromBGEvent);
  eleFinalState->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);//equiv. to IsPrimary();
  task->AddSignalMC(eleFinalState);
}

//________________________________________________________________
AliAnalysisFilter* SetupTrackCutsAndSettings(TString cutDefinition, Bool_t isESD=kTRUE)
{
  std::cout << "SetupTrackCutsAndSettings( cutInstance = " << cutDefinition << " )" <<std::endl;
  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!
  // do not change these initial values!
  selectedPairCutsPre=-1;
  isPrefilterCutset=kFALSE;
  anaFilterExtra = new AliAnalysisFilter("anaFilterExtra","anaFilter");
  rejCutMee=-1;
  rejCutTheta=-1;
  rejCutPhiV=3.2; // relevant are values below pi, so initialization to 3.2 means disabling.

  AnalysisCut AnaCut;


  LMEECutLib* LMcutlib = new LMEECutLib();
  if (cutDefinition == "kPbPb2015_Pt400_tightTOFreq_MB"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_tightTOFreq);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_looseTOFif_MB"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_looseTOFif);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetMixing(LMEECutLib::kEventMixing_1);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "kPbPb2015_Pt400_tightTOFif_MB"){
    AnaCut.SetPIDAna(LMEECutLib::kPbPb2015_Pt400_tightTOFif);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kSPDfirst);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut0"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_0);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_0);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut1"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_1);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_1);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut2"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_2);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_2);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut3"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_3);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_3);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut4"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_4);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_4);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut5"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_5);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_5);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut6"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_6);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_6);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut7"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_7);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_7);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut8"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_8);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_8);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut9"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_9);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_9);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut10"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_10);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_10);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut11"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_11);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_11);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut12"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_12);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_12);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut13"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_13);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_13);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut14"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_14);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_14);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut15"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_15);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_15);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut16"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_16);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_16);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut17"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_17);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_17);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut18"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_18);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_18);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else if (cutDefinition == "cut19"){
    AnaCut.SetPIDAna(LMEECutLib::kPIDcut_19);
    AnaCut.SetTrackSelectionAna(LMEECutLib::kTRACKcut_19);
    AnaCut.SetCentrality(LMEECutLib::kPbPb_00to80);
    AnaCut.SetStandardCut();
  }
  else {
    cout << " =============================== " << endl;
    cout << " ==== INVALID CONFIGURATION ==== " << endl;
    cout << " cutInstance = " << cutDefinition << endl;
    cout << " =============================== " << endl;
    return 0x0;
  }

  if (!isPrefilterCutset) {
    anaFilter->AddCuts( LMcutlib->GetPIDCutsAna(AnaCut) );
  }
  else { // cutInstance for prefilter efficiency determination:
    anafilter->AddCuts( LMcutlib->GetESDTrackCutsAna(AnaCut) );
    anaFilter->AddCuts( LMcutlib->GetTrackSelectionPre(AnaCut) );
    anaFilterExtra->AddCuts( LMcutlib->GetPIDCutsAna(AnaCut) );
  }

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
