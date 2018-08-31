
TString names=("TTree_cuts");
TObjArray*  arrNames=names.Tokenize(";");
const Int_t nDie=arrNames->GetEntriesFast();

// main task settings
// fill resolutions for one cutInstance (step 1).
const Bool_t CalcResolution   = kTRUE;
// use previously extracted resolutions (step 2).
TString resolutionfile = "";
Bool_t CalcEfficiencyRec      = kFALSE;  // use given resolution file to smear the kinematics.
Bool_t bUseRelPResolution     = kTRUE;  // specify if the file contains a relative or an absolute momentum resolution array.
Bool_t bUseEtaResolution      = kTRUE; // kFALSE means using theta instead of eta.
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
// Needs to be set for correct usage of direct pair efficiency
const Double_t EtaMinCut = -0.8;
const Double_t EtaMaxCut = 0.8;
const Double_t PtMinCut  = 0.2;
const Double_t PtMaxCut  = 8.0;
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
const Double_t PteeMax   = 10.;
const Int_t    nBinsPtee = 1000;

// run string must be sorted in increasing order!
const TString sRuns("265309, 265332, 265334, 265336, 265338,
                    265339, 265342, 265343, 265344, 265377, 
                    265378, 265381, 265383, 265384, 265385, 
                    265387, 265388, 265419, 265420, 265421,
                    265422, 265424, 265425, 265426, 265427, 
                    265435, 265499, 265500, 265501, 265212,
										265525");

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
const Double_t  PtMaxGEN  =  10.;    // 8 GeV is current upper limit of PtBins[]. Dont want overflow bin filled.

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
AliAnalysisCuts* SetupEventCuts()
{
  // event cuts are identical for all analysis 'cutInstance's that run together!
	AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");

  /* eventCuts->SetVertexType(AliDielectronEventCuts::kVtxSPD); // AOD */
	/* eventCuts->SetRequireVertex(); */
	/* eventCuts->SetMinVtxContributors(1); */
	/* eventCuts->SetVertexZ(-10.,10.); */

	return eventCuts;}


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
  eleFinalState->SetLegSources(AliDielectronSignalMC::kFinalState, AliDielectronSignalMC::kFinalState);


  // eleFinalState->SetLegSources(AliDielectronSignalMC::kDontCare, AliDielectronSignalMC::kDontCare);
  // eleFinalState->SetMotherPDGs(22,22,kTRUE,kTRUE); // exclude conversion electrons. has no effect for final state ele.

  // eleFinalState->SetMotherSources(AliDielectronSignalMC::kDirect, AliDielectronSignalMC::kDirect);//equiv. to IsPrimary();
  task->AddSignalMC(eleFinalState);
}

//________________________________________________________________
AliAnalysisFilter* SetupTrackCutsAndSettings(TString cutDefinition, Bool_t hasITS = kTRUE)
{

	std::cout << "SetupTrackCutsAndSettings()" <<std::endl;
	AliAnalysisFilter *AnaCut = new AliAnalysisFilter("AnaCut","AnaCut"); // named constructor seems mandatory!
  isPrefilterCutset=kFALSE;
  rejCutMee=-1;
  rejCutTheta=-1;
  rejCutPhiV=3.2; // relevant are values below pi, so initialization to 3.2 means disabling.
	// produce analysis filter by using functions in this config:
  LMEECutLib* LMcutlib = new LMEECutLib(hasITS);
  if(cutDefinition == "TTree_Cuts"){
    AnaCut->AddCuts(LMcutlib->GetCentralityCuts(LMEECutLib::kTTreeCuts));
		Printf("Cent cuts loaded");
    AnaCut->AddCuts(LMcutlib->GetTrackCuts(LMEECutLib::kTTreeCutsESD, LMEECutLib::kTTreeCuts));
		Printf("Track cuts loaded");
  }
	else{
		std::cout << "Cut definition not defined...oida" << std::endl;
		return 0x0;
	}

	std::cout << "...cuts added!" <<std::endl; 
	
	return AnaCut;
}

