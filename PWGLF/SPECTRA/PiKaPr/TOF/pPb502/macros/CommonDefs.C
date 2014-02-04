
#ifndef COMMONDEFS
#define COMMONDEFS

#define PROTONLEAD 1
#define LEADLEAD 0
#define PROTONPROTON 0

/* analysis dependent stuff */

/***** p-Pb *****/
#if PROTONLEAD

Double_t tofReso = 80.;
Double_t tofTail = 75.;

Double_t scaletexpreso[5] = {1., 1., 1., 1., 1.};
Double_t scaletimezerosigma = 1.;
Double_t forcetimezeroineff = 0.;
Double_t timezero_spread = 200.;

Int_t acceptEventType = 0;
Int_t centralityEstimator = 1;
const Int_t NcentralityBins = 7;
Double_t centralityBin[NcentralityBins + 1] = {0., 5., 10., 20., 40., 60., 80., 100.};

Float_t rapidityShift = -0.465;
Float_t rapidityMinCut = 0.0;
Float_t rapidityMaxCut = 0.5;

/* TZERO corrections */
Double_t TZEROFILL_shift = 0.;
Double_t TZEROA_shift = 0.;
Double_t TZEROC_shift = 0.;
Double_t TZEROTOF_shift = 0.;

Double_t TZEROvertexCorr = 0.5;

Double_t TZEROA_sigma = 1000.;//247.187;
Double_t TZEROC_sigma = 1000.;//235.113;
Double_t TZEROTOF_sigma = 1000.;//234.423;

Double_t TZEROTOF_resoScaleFactor = 1.25;


Int_t multcentColor[7] = {
  kRed,
  kPink+1,
  kOrange+1,
  kYellow+1,
  kGreen+1,
  kAzure+1,
  kViolet+1,
};

/***** Pb-Pb *****/
#elif LEADLEAD

Float_t rapidityShift = 0.;
Float_t rapidityMinCut = -0.5;
Float_t rapidityMaxCut = 0.5;
Float_t rapidityCut = 0.5;

Int_t acceptEventType = 0;
Int_t centralityEstimator = AliAnalysisEvent::kCentEst_V0M;
const Int_t NcentralityBins = 10;
Double_t centralityBin[NcentralityBins + 1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90.};

/* TZERO corrections */
Double_t TZEROFILL_shift = 0.;
Double_t TZEROA_shift = 0.;
Double_t TZEROC_shift = 0.;
Double_t TZEROTOF_shift = 0.;

Double_t TZEROvertexCorr = 0.5;

Double_t TZEROA_sigma = 1000.;//247.187;
Double_t TZEROC_sigma = 1000.;//235.113;
Double_t TZEROTOF_sigma = 1000.;//234.423;

Double_t TZEROTOF_resoScaleFactor = 1.25;


Int_t multcentColor[10] = {
  kRed,
  kOrange+1,
  kOrange,
  kYellow,
  kYellow+1,
  kGreen,
  kGreen+1,
  kCyan+1,
  kBlue,
  kMagenta,
  //  kMagenta+1  
};

/***** p-p *****/
#elif PROTONPROTON /* LHC10d, run 126088 */

Float_t rapidityShift = 0.;
Float_t rapidityMinCut = -0.5;
Float_t rapidityMaxCut = 0.5;
Float_t rapidityCut = 0.5;

Int_t acceptEventType = 1;
Int_t centralityEstimator = 999;
const Int_t NcentralityBins = 6;
Double_t centralityBin[NcentralityBins + 1] = {0., 5., 10., 15., 20., 25., 100};

/* TZERO corrections */
Double_t TZEROFILL_shift = 0.;//-30.6174;
Double_t TZEROA_shift = -6.5e6;// + -13.8783;
Double_t TZEROC_shift = -6.5e6;// + -21.7581;
Double_t TZEROTOF_shift = 0.;//-45.1779;

Double_t TZEROvertexCorr = 0.5;

Double_t TZEROA_sigma = 1000.;//247.187;
Double_t TZEROC_sigma = 1000.;//235.113;
Double_t TZEROTOF_sigma = 1000.;//234.423;

Double_t TZEROA_resolution[NcentralityBins] = {86.4697,
					       81.6549,
					       79.0132,
					       76.5572,
					       74.8259,
					       73.3834};

Double_t TZEROC_resolution[NcentralityBins] = {58.0423,
					       52.9241,
					       50.8822,
					       49.05,
					       48.6712,
					       46.5776};

Double_t TZEROAND_resolution[NcentralityBins] = {48.383,
						 51.9983,
						 51.2704,
						 50.1974,
						 51.0198,
						 48.6588};
						  
Double_t TOFTZEROADIFF[NcentralityBins] = {-20.2311,
					   -19.539,
					   -18.8685,
					   -19.5327,
					   -18.4323,
					   -18.0336};

Double_t TOFTZEROCDIFF[NcentralityBins] = {-16.0008,
					   -11.9249,
					   -9.29299,
					   -8.46063,
					   -8.1128,
					   -4.89226};

Double_t TOFTZEROTOFDIFF[NcentralityBins] = {-11.9418,
					     -12.4526,
					     -14.0921,
					     -18.1132,
					     -17.5909,
					     -18.9787};

Double_t TZEROTOF_resoScaleFactor = 1.25;

Double_t TZEROFILL_sigma = 1000.;

Int_t multcentColor[10] = {
  kRed,
  kOrange+1,
  kOrange,
  kYellow,
  kYellow+1,
  kGreen,
  kGreen+1,
  kCyan+1,
  kBlue,
  kMagenta,
  //  kMagenta+1  
};
#endif

/* commom binning and similar business */

const Char_t *t0FillOnlineFileName = "T0FillOnline.139465.extended.root";
Double_t t0Fill_offset = -1.26416e+04;

//const Char_t *enabledChannelsFileName = "enabledChannels.root";
const Char_t *enabledChannelsFileName = NULL;

const Int_t NptBins = 46;
Double_t ptBin[NptBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};

const Int_t NpBins = 46;
Double_t pBin[NpBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};

const Int_t NmtBins = 46;
Double_t mtBin[NmtBins + 1] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};

Int_t particleColor[5] = {1, 1, 4, 8, 2};
Int_t chargeMarker[2] = {20, 25};

const Char_t *partChargeName[5][2] = {"e^{+}", "e^{-}", "#mu^{+}", "#mu^{-}", "#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"};

const Double_t kEpsilon = 0.001;

enum ECharge_t {
  kPositive,
  kNegative,
  kNCharges
};
const Char_t *chargeName[kNCharges] = {
  "positive",
  "negative"
};

const Int_t NetaBins = 10;
Double_t etaMin = -1.;
Double_t etaMax = 1.;
Double_t etaStep = (etaMax - etaMin) / NetaBins;
Double_t etaBin[NetaBins + 1]; /* computed at run-time */

const Int_t NyBins = 20;
Double_t yMin = -1.;
Double_t yMax = 1.;
Double_t yStep = (yMax - yMin) / NyBins;
Double_t yBin[NyBins + 1]; /* computed at run-time */

const Int_t NphiBins = 10;
Double_t phiMin = 0.;
Double_t phiMax = 2. * TMath::Pi();
Double_t phiStep = (phiMax - phiMin) / NphiBins;
Double_t phiBin[NphiBins + 1]; /* computed at run-time */

const Int_t NptsubBins = 4;
Double_t ptsubBin[NptsubBins + 1] = {0.2, 0.5, 1.0, 1.5, 5.0};
Int_t ptsubBinMin[NptsubBins] = {0, 6, 16, 21};
Int_t ptsubBinMax[NptsubBins] = {5, 15, 20, 45};

const Int_t NdcaBins = 2000;
Double_t dcaBin[NdcaBins + 1];
Double_t dcaMin = -5., dcaMax = 5., dcaStep = (dcaMax - dcaMin) / NdcaBins;



#endif
