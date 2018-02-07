TString names=("track_cuts");

TObjArray*  arrNames           = names.Tokenize(";");
const Int_t nDie               = arrNames->GetEntriesFast();
Int_t       selectedCentrality = -1; // not yet implemented
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
const Int_t    nBinsEta = 40; //flexible to rebin
// phi bins
const Double_t PhiMin   = 0.;
const Double_t PhiMax   = 6.2832;
const Int_t    nBinsPhi = 60; //flexible to rebin

const Double_t PtBins[] = {0.0, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
						   0.75, 0.80, 0.85, 0.90, 1.00, 1.10, 1.15, 1.25, 1.35, 1.55, 1.80,
						   2.05, 2.30, 2.60, 2.90, 3.30, 3.60, 4.00, 5.00, 6.50, 8.00, 10.0};
//Ivan binning
/*const Double_t PtBins[] = {
  0.000,0.050,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700,0.750,0.800,0.850,0.900,0.950,
  1.000,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00,2.10,2.30,2.50,3.00,3.50,
  4.00,5.0,6.0,7.0,10.0
};*/

// Bool_t bUseRelPResolution = kTRUE; //not used
Bool_t bUseEtaResolution = kTRUE; // use eta or theta resolution?
Bool_t CalcEfficiencyRec = kTRUE;
Bool_t CalcEfficiencyPoslabel = kFALSE;
Bool_t CalcResolution = kFALSE;
Bool_t MakeResolutionSparse = kFALSE;
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
const TString sRuns("265309, 265334, 265335, 265338, 265339, 
                     265342, 265343, 265344, 265377, 265378, 
                     265381, 265383, 265384, 265385, 265387, 
                     265388, 265419, 265420, 265421, 265422, 
                     265424, 265427, 265435, 265499, 265500, 
                     265501, 265521, 265525");

//
// ^^^^^^^^^^ [/end binning histograms] ^^^^^^^^^^

//________________________________________________________________
// specify if track tree shall be filled and written to file (only recommended for small checks!)
const Bool_t    writeTree = kTRUE;
// specify for which "cutInstance" the support histos should be filled!
const Int_t     supportedCutInstance = 0;
//
//________________________________________________________________
// settings which are identical for all configs that run together
// event cuts
const Bool_t    reqVertex = kTRUE;
const Double_t  vertexZcut = 10.;
//Set centrality in AddTask arguments
/* const Double_t  CentMin =  -2.; */
/* const Double_t  CentMax = 102.; */
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


//________________________________________________________________
AliAnalysisCuts* SetupEventCuts()
{
	  // event cuts are identical for all analysis 'cutInstance's that run together!
	  AliDielectronEventCuts *eventCuts = new AliDielectronEventCuts("eventCuts","Vertex Track && |vtxZ|<10 && ncontrib>0");

	  eventCuts->SetRequireVertex();
	  eventCuts->SetMinVtxContributors(1);
	  eventCuts->SetVertexZ(-10.,10.);

	  return eventCuts;
}

//________________________________________________________________
AliAnalysisFilter* SetupTrackCutsAndSettings(Int_t cutInstance, Bool_t hasITS = kTRUE)
{
	  std::cout << "SetupTrackCutsAndSettings()" <<std::endl;
	  AliAnalysisFilter *anaFilter = new AliAnalysisFilter("anaFilter","anaFilter"); // named constructor seems mandatory!
	  // do not change these initial values!
	  selectedPID=-1;
	  isPrefilterCutset=kFALSE;
	  rejCutMee=-1;
	  rejCutTheta=-1;
	  rejCutPhiV=3.2; // relevant are values below pi, so initialization to 3.2 means disabling.
    
   
	// produce analysis filter by using functions in this config:
	// -----
	anaFilter->AddCuts( SetupTrackCuts(cutInstance, hasITS) );
	anaFilter->AddCuts( SetupPIDcuts(cutInstance) );
	std::cout << "...cuts added!" <<std::endl; 
	  
	//std::cout << "__________ anaFilter->GetCuts()->Print() __________ cutInstance = " << cutInstance << std::endl;
	//anaFilter->GetCuts()->Print();
	//std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<std::endl;
	//std::cout << "__________ anaFilter->GetCuts()->Dump() __________ cutInstance = " << cutInstance <<std::endl;
	//anaFilter->GetCuts()->Dump();
	//std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" <<std::endl;
	
	return anaFilter;
}


// prefilter cuts are not used at the moment
//________________________________________________________________
Int_t SetupPrefilterPairCuts(Int_t cutInstance)
{
  std::cout << "SetupPrefilterPairCuts()" <<std::endl;
  
  //PID cuts turned off. Only want track reconstruction efficiency
  return 1;
}    


//________________________________________________________________
AliAnalysisCuts* SetupTrackCuts(Int_t cutInstance, Bool_t hasITS = kTRUE)
{
  std::cout << "SetupTrackCuts()" <<std::endl;
  //AliAnalysisCuts* trackCuts=0x0;
    
    AliESDtrackCuts *fesdTrackCuts = new AliESDtrackCuts();

    //Cuts implemented in TreeMaker
		//FilterBit 4 used to filter AODs
		//Set via GetStandardITSTPCTrackCuts 2011(kFALSE, 1)
    fesdTrackCuts->SetMinNCrossedRowsTPC(70);
    //fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fesdTrackCuts->SetMaxChi2PerClusterTPC(4);
    fesdTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    fesdTrackCuts->SetRequireTPCRefit(kTRUE);
    fesdTrackCuts->SetRequireITSRefit(kTRUE);
    fesdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
		fesdTrackCuts->SetRequireSigmaToVertex(kFALSE);
    fesdTrackCuts->SetMaxChi2PerClusterITS(36);

		//Manually set by fitler bit 4
		////Override setting in TrackCuts2011
    fesdTrackCuts->SetMaxDCAToVertexXY(2.4); 
    fesdTrackCuts->SetMaxDCAToVertexZ(3.2); 
    fesdTrackCuts->SetDCAToVertex2D(kTRUE);

		//Manually implemented track cuts in TreeMaker
		//General
		fesdTrackCuts->SetPtRange(0.2, 10);
    fesdTrackCuts->SetEtaRange(-0.8, 0.8);
		fesdTrackCuts->SetMaxDCAToVertexZ(3.0); 
    fesdTrackCuts->SetMaxDCAToVertexXY(1.0);

    //TPC
    fesdTrackCuts->SetMinNClustersTPC(70);
    //fesdTrackCuts->SetMinNCrossedRowsTPC(60); FilterBit4 stronger cut
    fesdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.3);
    //fesdTrackCuts->SetMaxChi2PerClusterTPC(4.5);
    fesdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
    
    //ITS
		if(hasITS){
			fesdTrackCuts->SetMinNClustersITS(4);
		}else{
			fesdTrackCuts->SetMinNClustersITS(2);
		}

  	return fesdTrackCuts;
}

//________________________________________________________________
AliAnalysisCuts* SetupPIDcuts(Int_t cutInstance)
{
  	std::cout << "SetupPIDcuts()" <<std::endl;
  	AliAnalysisCuts* pidCuts=0x0;
  
  	AliDielectronPID *pid = new AliDielectronPID();
  
	//The only PID cut applied when creating Trees
  	pid->AddCut(AliDielectronPID::kTPC, AliPID::kElectron, -4, 4, 0, 1e30, kFALSE, AliDielectronPID::kRequire, AliDielectronVarManager::kPt);
    
  	pidCuts = pid;   
  	return pidCuts;
}



