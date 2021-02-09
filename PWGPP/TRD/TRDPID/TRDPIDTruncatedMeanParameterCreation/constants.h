
#ifndef CONSTANTS_H 
#define CONSTANTS_H 



//======== cuts ================================
bool fEtaCorrection                         = kFALSE; // switch on/off corrections
bool NclsCorrection                         = kFALSE;
bool CentCorrection                         = kFALSE; // only required for PbPb
Int_t Iteration                             = 0;      // test for iterative procedure
Double_t GaussFrac                          = 0.5;    // gaussian fraction to be considered in fit
bool CorrectToMean                          = kTRUE;  // switch if corrections maps correct data to the mean in each binning (kTRUE, default) or to the fit (kFALSE)
bool FixedResStat                           = kFALSE; // switch to choose if all bins in the resolution plot should have equal statistics
Int_t ResolutionMax                         = 10000;  // value related to switch
bool MomentumCuts                           = kFALSE; // to switch all momentum cuts on or off (defaul == off)

// centrality binning
//const Double_t CentBins[15] ={0,3,6, 9,12, 15, 20,25, 30,35, 42,50,70,90,100};
//const Int_t NCentBins = 14;
const Double_t CentBins[7]={0,100}; //10,30, 50, 70, 90, 100};
const Int_t NCentBins=1; //6

// eta binning
const Double_t EtaBins[6]={-1, -0.56, -0.15, 0.15, 0.56, 1};
const Int_t NEtaBins=5;

const Int_t     nch                         = 6;
const Int_t     nchMax                      = 6;  // no longer in use. Originally used to separate the different numbers of channels. Now: Nch >=4 is used. But not alter this value (output etc.)!
Int_t            nchMin                     = 4;
const Int_t     ncls[3]                     = {17, 17, 17}; // default min average #cls/tracklet
const Int_t     minNumberOfEntries          = 250; //250; //450;
const Double_t  maxFitSliceFitRMS           = 0.9; // 0.8; //0.9 //0.7;
const Int_t     minNumberOfEntriesRes       = 500; //5000; //3000 //10000;
const Double_t  maxFitSliceFitRMSRes        = 0.4; // 0.4; //0.5  0.2; Ppbpb pPB


// nSigma cuts (check for each period)
// For  LHC 15oHIR
const Double_t  cutNSigmaTOF[6]             = {-1.4,     0.6,      -0.9,     1.1,      -0.9,    1.1};         // {e_min, e_max, pi_min, pi_max, p_min, p_max}
const Double_t  cutNSigmaTPC[6]             = {-1.2,      0.8,      -0.9,     1.1,      -1.,    1.1};         // {e_min, e_max, pi_min, pi_max, p_min, p_max}

// For LHC 13bc
//const Double_t  cutNSigmaTOF[6]             = {-1.3,     0.7,      -1,     1,      -1.1,    0.9};         // {e_min, e_max, pi_min, pi_max, p_min, p_max}
//const Double_t  cutNSigmaTPC[6]             = {-1.1,      0.9,      -1,     1,      -1,    1.};         // {e_min, e_max, pi_min, pi_max, p_min, p_max}



// not used anymore
const Double_t  exclusionCutNSigmaTOF[4]    = {0.5,     1.5,    -1,     0.5};                       // {e_min, e_max, pi_min, pi_max}
const Double_t  exclusionCutNSigmaTPC[4]    = {-0.5,    0.5,    -0.5,   0.5};                       // {e_min, e_max, pi_min, pi_max}
const Double_t  yMaxLayer[6]                = {46.8,    49.05,  51.25,  53.45,  55.7,   57.9};      // half of active width per layer (cm)
const Double_t  etaUpperBound[5]            = {0.851,   0.527,  0.145,  -0.157, -0.536};
const Double_t  etaLowerBound[5]            = {0.536,   0.157,  -0.145, -0.527, -0.851};

// slice fit
const Double_t  maxChi = 20;
Int_t     ChiVetoCount = 0;

// options to be set before running the macro
// pid option
// 2 == elec + pion
// 3 == elec + pion + proton
const Int_t pidOpt = 3;



//======= set path to input tree ===============

TString RunPeriodArray[4]={"LHC13", "LHC13", "LHC15n", "LHC18d"};
TString RunPeriod("LHC18d");
TString inputPath("/media/ikp/DATA/Data/TRDPID_NTuple");
TString inputSpecArray[4]={"/LHC13bc", "/LHC13bcNew", "/LHC15nCorrect", "/LHC18d"};
TString inputSpec("/LHC18d");
TString inputFile("/TRDPIDTree_tree.root"); ///TRDPIDTree_tree.root");
TString input = inputPath + inputSpec + inputFile;

//====== set ouput path =======================


TString EtaSpecStr("_etaTPCtglCorr");
TString NclsCorrStr("_NclsCorr");
TString CentCorrStr("_CentCorr");
//TString FileSpec("elAddnSigmaTPC_noLocalEtaY");
TString FileSpecArray[4]={Form("LHC13Pass2_from%ito%iTracklets", nchMin, nchMax), Form("LHC13Pass4_from%ito%iTracklets", nchMin, nchMax), Form("LHC15n_from%ito%iTracklets", nchMin, nchMax), Form("LHC18d_from%ito%iTracklets", nchMin, nchMax)};
TString FileSpec(Form("LHC18d_from%ito%iTracklets", nchMin, nchMax) );
//TString FileSpec(Form("LHC13bcPass4_from%ito%iTracklets", nchMin, nchMax) );



//======== general definitions ==================
#define EPSILON 1e-12
typedef Double_t (*FFunc)(const Double_t *xx, const Double_t *par);
Int_t fgFitNData    = 0;
Double_t * fgFitX   = 0x0;
Double_t * fgFitY   = 0x0;
Double_t * fgFitEx  = 0x0;
Double_t * fgFitEy  = 0x0;
FFunc fgFitFunc     = 0x0;
Double_t fTol       = 1;
Int_t fMaxIter      = (Int_t)1e8;


//======== particle masses ======================
const Double_t mProton    = 0.938272;
const Double_t mPion      = 0.139570;
const Double_t mElectron  = 5.109989e-04;

#endif
