#ifndef TMCALCONSTANTS_H
#define TMCALCONSTANTS_H

//System constants
static const int fgkMaxModules = 12;
//static const char fgkDeadMapName[] = "Deadmap.txt"; 
//A static const char[] apparently cannot be initialized
//in class, so it's either a #define or a constant in the constructor...
#define fgkDeadMapName "Deadmap"
#define fgkReferenceObject "RunCalib"

// Path to data and reference probably different for PHOS and EMCAL..
// PHOS
#define fgkSourceFileTemplatePhos "/data1/phos/run%.4u.root"
// #define fgkSourceFileTemplate "/data1/phos/Period_LHC07a.Run_%.9u.Host_001.Seq_10.root"
#define fgkReferenceFilePhos "/home/phos/talho/PHOS/ana/RunCalib.root"
// EMCAL
#define fgkSourceFileTemplateEmCal "/local/data/Run_%.9u.Seq_1A.Stream_0.root"
#define fgkReferenceFileEmCal "/home/emcaldaq/ref/RunCalib.root"
 
//Analysis constants
static const int fgkDeadThreshold = 5; //The threshold for declaring a cell dead
static const unsigned int fgkDefaultRunNo = 7519; //The default run number

//GUI constants
static const int fgkHistosPerRow = 5;
static const int fgk2DHistoWidth = 188;
static const int fgk2DHistoHeight = 172;
static const int fgk2DHistoMarginX = 5;
static const int fgk2DHistoMarginY = 5;
static const int fgk2DSingleHistoWidth = 600;
static const int fgk2DSingleHistoHeight = 500;
static const int fgk2DSingleHistoMarginX = 5;
static const int fgk2DSingleHistoMarginY = 5;
static const double fgkMaxMaxRatioExp = 10.0;//These are the extremes of user control, that's why 'maxmax' and 'minmax'
static const double fgkMinMaxRatioExp = -10.0;
static const double fgkRatioStepScale = 10.0; //The number of steps per unit exponent. Determines the granularity of the ratio sliders.

// having a custom DeadMap palette costs some time (nicer, but slower), so make this
// usage optional with a define switch
//#define fgkCustomDeadMapPalette 1 // comment out to de-activate

#ifdef fgkCustomDeadMapPalette
//Dead map colors. These are indices to the
//standard root palette.
static const Int_t fgkLiveTowerColor = 10;
static const Int_t fgkDeadTowerColor = 1;
static const Int_t fgkResurrectedTowerColor = 4;
static const Int_t fgkRecentlyDeceasedTowerColor = 2;
#else
// when we don't use the custom palette we get stuck with these colors for palette1
static const Int_t fgkLiveTowerColor = 0; // kWhite
static const Int_t fgkDeadTowerColor = 432; // kCyan
static const Int_t fgkResurrectedTowerColor = 827; // kSpring + 7
static const Int_t fgkRecentlyDeceasedTowerColor = 632; // kRed
#endif

static const int fgkDeadCountTextMargin = 30;

#endif // TMCALCONSTANTS_H
