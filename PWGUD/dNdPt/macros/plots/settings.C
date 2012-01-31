//------------------------------------------------------------------------------
// settings.C
//
// setting for all data files and folder
// plot options are not set exclusively here
//------------------------------------------------------------------------------


//
// store output plots? (gif+eps) directory?
//
Bool_t SAVE_FIGURES = kTRUE;
// Bool_t SAVE_FIGURES = kFALSE;
TString outputDir = "/u/jacek/CMS_ALICE_UA1_Michael/figures/";

//
// draw ALICE logo and/or PRELIMINARY tag?
//
Bool_t SHOW_LOGO = kTRUE;
Bool_t SHOW_PRELIM = kTRUE;
//Bool_t SHOW_LOGO = kFALSE;
//Bool_t SHOW_PRELIM = kFALSE;

//
// graphics and plot options
//
Int_t textFont = 43;
Int_t titleFont = 43;
Int_t labelFont = 43;
Int_t font = 43;

Float_t titleFontSize = 24;
Float_t markerSize = 1.25;
Float_t legendTextSize = 20;
Float_t labelSize = 24;
Float_t titleSize = 24;

Float_t leftMargin = 0.16; 

const Color_t colorAlice = kGray+2;
const Color_t colorAliceFit = kBlack;
const Color_t colorPhojet = kMagenta+2;
const Color_t colorPythia109 = kCyan+3;
const Color_t colorPythia306 = kRed+2;
const Color_t colorPythia320 = kBlue+1;
const Color_t colorAtlas = kBlue+1;
const Color_t colorCms = kRed+2;
const Color_t colorAliceErrors = kGray; //kSpring-8;
const Color_t colorUa1 = kRed+2;


//
// various settings
//

// limit for shifting pt bins
const Double_t fit_limit = 1e-14;
const Double_t fit_limit2 = 1e-14;


// ua1 inelastic cross section (mb); used to calculate yield/
const Double_t sigmaInelUa1 = 43.5;  
// yields for (h+ + h-)/2
const Double_t avgToHadr = 2.0; 

// eta range used for model comparison
// -0.8 < eta < 0.8 
const Double_t etaRange = 1.6;

// pi
const Double_t M_PI = 3.1415926535897932384626433832795028841971693993751;


//
// set filenames and number of data points
//

// compatibility setting
const Int_t binsAlice = 46;

// ALICE INEL data file (for comparison to models)
Char_t* filenameInelAlice = "/u/jacek/alice/dNdPt/output/data_points/ALICE_Yield_INEL_900GeV.txt";
const Int_t binsInelAlice = 46;                 // number of bins (= lines in data textfile)

// ALICE NSD data file (for comparison to ATLAS and CMS)
Char_t* filenameNsdAlice = "/u/jacek/alice/dNdPt/output/data_points/ALICE_Yield_NSD_900GeV.txt";
const Int_t binsNsdAlice = 46;                 // number of bins (= lines in data textfile)

// ALICE NSD Invariant Yield data file (for comparison to UA1)
Char_t* filenameYieldAlice = "/u/jacek/alice/dNdPt/output/data_points/ALICE_InvYield_NSD_900GeV.txt";
const Int_t binsYieldAlice = 46;                 // number of bins (= lines in data textfile)

// ATLAS NSD data file
Char_t* filenameAtlas = "/u/mknichel/dNdPt/data/atlas_data.txt";
const Int_t binsAtlas = 33;                 // number of bins (= lines in data textfile)

// CMS NSD data file
Char_t* filenameCms = "/u/mknichel/dNdPt/data/cms_data.txt";
const Int_t binsCms = 24;                   // number of bins (= lines in data textfile)

// UA1 invariant Cross section data file
Char_t* filenameUa1 = "/u/mknichel/dNdPt/data/UA1data_0.9TeV_pt_CrossSec_errCrossSec.txt";
const Int_t binsUa1 = 52;              // number of bins (= lines in data textfile)

// photjet root file
Char_t* filenamePhojet = "/u/jacek/alice/dNdPt/output/fastsim/out_fastsim_900Phojet_EtaRange0.8.root";   
const Int_t binsPhojet = 54;   // number of bins

// Pythia D6T (109) root file
Char_t* filenamePythia109 = "/u/jacek/alice/dNdPt/output/fastsim/out_fastsim_900Pythia109_EtaRange0.8.root";
const Int_t binsPythia109 = 54;   // number of bins

// Pythia ATLAS-CSC (306) root file
Char_t* filenamePythia306 = "/u/jacek/alice/dNdPt/output/fastsim/out_fastsim_900Pythia306_EtaRange0.8.root";
const Int_t binsPythia306 = 54;   // number of bins

// Pythia Perugia0 (320) root file
Char_t* filenamePythia320 = "/u/jacek/alice/dNdPt/output/fastsim/out_fastsim_900Pythia320_EtaRange0.8.root";
const Int_t binsPythia320 = 54;   // number of bins

// ALICE Logo file
Char_t* filenameAliceLogo = "/u/mknichel/dNdPt/LogoALICE-DEF-transp.png";
TImage* logo = TImage::Open(filenameAliceLogo);

//
// fit functions
// these have to be set also in shiftPtxxx.C and makePlotsxxx.C
//
TF1* fitNsd = new TF1("fitNsd","[0]*(x*x/sqrt(0.14*0.14+x*x))*(1+x/[2])^((-1) * [1])",0.0,10.0); 
TF1* fitInel = new TF1("fitInel","[0]*(x*x/sqrt(0.14*0.14+x*x))*(1+x/[2])^((-1) * [1])",0.0,10.0); 
TF1* fitYield = new TF1("fitYield","[0]*(x*x/sqrt(0.14*0.14+x*x))*(1+x/[2])^((-1) * [1])",0.0,10.0); 
