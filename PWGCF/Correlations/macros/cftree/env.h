#ifndef env_h
#define env_h 1

#include <TStyle.h>
#include <TPad.h>
#include <TMath.h>

//variables used in calc_sum_of_ratios
const Float_t pi = TMath::Pi();
Int_t color[]={1,2,3,4,6,7,8,9,11,12,13,14,15,1,2,6,4,8,3,9,11,12,13,14,15};
Int_t line[]={1,1,1,1,1,1,2,3,4,6,7,8,9,10};
Int_t fillstyle[]={3001,3004,3005,3006,3007,3021,3022,3020,3001,3004,3005,3006,3007,3021,3022,3020,3001,3004,3005,3006,3007,3021,3022,3020,3001,3004,3005,3006,3007,3021,3022,3020};
Int_t marker[]={2,3,4,5,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};
Int_t openmarker[]={kOpenCircle,kOpenSquare,kOpenTriangleUp,kOpenTriangleDown,kOpenStar,kOpenDiamond,kOpenCross,kOpenCircle,kOpenSquare,kOpenTriangleUp,kOpenTriangleDown,kOpenStar,kOpenDiamond,kOpenCross,kOpenCircle,kOpenSquare,kOpenTriangleUp,kOpenTriangleDown,kOpenStar,kOpenDiamond,kOpenCross,kOpenCircle,kOpenSquare,kOpenTriangleUp,kOpenTriangleDown,kOpenStar,kOpenDiamond,kOpenCross};
Int_t fullmarker[]={kFullCircle,kFullSquare,kFullTriangleUp,kFullTriangleDown,kFullStar,kFullDiamond,kFullCross,kFullCircle,kFullSquare,kFullTriangleUp,kFullTriangleDown,kFullStar,kFullDiamond,kFullCross,kFullCircle,kFullSquare,kFullTriangleUp,kFullTriangleDown,kFullStar,kFullDiamond,kFullCross,kFullCircle,kFullSquare,kFullTriangleUp,kFullTriangleDown,kFullStar,kFullDiamond,kFullCross};
Int_t colors[] = {kViolet+2,kBlue,kGreen+3,kOrange+7,kRed,kCyan,kViolet+2,kBlue,kGreen+3,kOrange+7,kRed,kCyan};

void PadFor2DCorr(){
  gPad->SetPad(0, 0, 1, 1);
  gPad->SetLeftMargin(0.17);
  gPad->SetTopMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetTheta(55);
  gPad->SetPhi(45);
}

//trigger selection
enum                          {k_CINT7_B=0, k_CMSL7_B, k_CMSH7_B, k_CVHMV0M_B, k_CVHMSH2_B,  k_CVHMV0MMSL_B, k_CVHMSH2MSL_B,kNTriggers};
TString TriggerName[]=        {  "CINT7_B", "CMSL7_B", "CMSH7_B", "CVHMV0M_B", "CVHMSH2_B",  "CVHMV0MMSL_B", "CVHMSH2MSL_B"};
UInt_t alltriggerclasses13TeV=(1<<k_CINT7_B | 1<<k_CMSL7_B | 1<<k_CMSH7_B | 1<<k_CVHMV0M_B | 1<<k_CVHMSH2_B | 1<<k_CVHMV0MMSL_B | 1<<k_CVHMSH2MSL_B);

// flag for type of analysis   0        1       2           3           4           5          6            7                 8                9                   10
enum                  {kTrkTrk=0, kTklTkl, kTrkTrkITS, kTrkTrkGen, kTklTklGen, kTklTklMC, kTrkTrkITSGen};
TString AnaTypeName[]={" TrkTrk","TklTkl","TrkTrkITS","TrkTrkGen","TklTklGen","TklTklMC","TrkTrkITSGen"};
// flag for track type
enum                    {kTrk=0, kTkl, kTrkITS, kTrkGen, kTklGen, kTklMC,kTrkITSGen};
//track type
TString sTrackType[11] = {"trk", "tkl", "trkits",    "mc","tklgen","tklmc","mcits"};
// filter bit selection
enum kFB{fb56,fb89,fb4,fb1,fbITSsa,kNfb};
TString fbname[] = {"fb56","fb89","fb4","fb1","fbITSsa"};
// flag for energy
enum {k7TeV=0, k13TeV};

//variable of the 6D histos
enum {kTrackDeta=0,kTrackPtAs=1,kTrackPtTr=2,kTrackCent=3,kTrackDphi=4,kTrackZvtx=5,kTrackNvar=6};
enum {kEventPtTr=0,kEventCent=1,kEventZvtx=2,kEventNvar=3};


const Float_t etaMin_tkl = -1.2;
const Float_t etaMax_tkl = +1.2;
const Float_t etaMin_mu = -4.;
const Float_t etaMax_mu = -2.5;
const Float_t etaMin_trk = -1.1;
const Float_t etaMax_trk =  1.1;
const Float_t etaMin_trkits = -.9;
const Float_t etaMax_trkits =  .9;

// Track histograms:                       deta, pt_ass, pt_tr,    cent,    dphi, zvtx
// tracklet-tracklet
const Int_t   ntbins_tkltkl[kTrackNvar] = {  48,      5,     5,      13,      72,    8};
const Double_t xtmin_tkltkl[kTrackNvar] = {-2.4,      0,     0,       0, -0.5*pi,   -4};
const Double_t xtmax_tkltkl[kTrackNvar] = {+2.4,      5,     5,      10, +1.5*pi,   +4};
Double_t bins_tkl[] = {0., 1., 2., 3., 4., 5.};
// tracklet-tracklet Gen
const Int_t   ntbins_tkltklGen[kTrackNvar] = {  48,      6,     6,      13,      72,    8};
const Double_t xtmin_tkltklGen[kTrackNvar] = {-2.4,    0.1,   0.1,       0, -0.5*pi,   -4};
const Double_t xtmax_tkltklGen[kTrackNvar] = {+2.4,      8,     8,      10, +1.5*pi,   +4};
Double_t bins_tklGen[] = {0.1,0.3,0.7,1.,2.,5, 8};
// track-track
const Int_t   ntbins_trktrk[kTrackNvar] = {   44,  5,    5,       10,      40,   9};
const Double_t xtmin_trktrk[kTrackNvar] = {  -2.2,  0.7,  0.7,    0.,  -1.*pi,  -4};
const Double_t xtmax_trktrk[kTrackNvar] = {  +2.2,  8,    8,    100.,      pi,  +5};
Double_t bins_trk[] = {0.7,1,1.5,2,3,8};
//integrated pt bin**********************************
//const Int_t   ntbins_trktrk[kTrackNvar] = {   44,   1,   1,     10,      40,   9};
//const Double_t xtmin_trktrk[kTrackNvar] = {  -2.2,  1,   1,     0.,  -1.*pi,  -4};
//const Double_t xtmax_trktrk[kTrackNvar] = {  +2.2,  4,   4,   100.,      pi,  +5};
//Double_t bins_trk[] = {1,4}; 
// track-track - its
const Int_t   ntbins_trktrkits[kTrackNvar] = {   36,    5,    5,      10,      40,   9};
const Double_t xtmin_trktrkits[kTrackNvar] = { -1.8,  0.7,  0.7,      0.,  -1.*pi,  -4};
const Double_t xtmax_trktrkits[kTrackNvar] = { +1.8,    8,    8,     100.,     pi,  +5};
Double_t bins_trkits[] = {0.7,1,1.5,2,3,8};
//for the calc_sum_of_ratios we may want a different binning
const Int_t   ntbins_trktrkits_calc[kTrackNvar] = {   36,    6,    6,      13,       72,   10};
const Double_t xtmin_trktrkits_calc[kTrackNvar] = { -1.8,  0.1,  0.1,      0.,  -0.5*pi,  -10};
const Double_t xtmax_trktrkits_calc[kTrackNvar] = { +1.8,    8,    8,     10.,  +1.5*pi,  +10};
Double_t bins_trkits_calc[] = {0.1, 0.3, 0.5, 0.7, 1., 3., 5.};

const Int_t nc = ntbins_trktrk[kTrackCent];
Double_t binsc[] =    {0, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
Double_t binscmin[] = {0, 1, 5, 10, 15, 20, 30, 40, 50, 70};
Double_t binscmax[] =    {1, 5, 10, 15, 20, 30, 40, 50, 70, 100};

//systematic error
TString systCheck[] = {"","_vertex02","_mixed_norm","_LowMultScaleAS","_NoLMScaling","_CMSSub","_ATLAS_IP","_ATLAS_IP_02_to_05","_ATLAS_IP_0_to_02","_ATLAS_EP","_0_to_02","_02_to_05","_hist_integral","_parabolic_fit","_constant_fit","_baselineGausFit","_baselineHighMult","_baselineHighMultParab","_no_v3","_fit_opt","_exclusion18","_exclusion10","_v2_analytical","_no_subtraction","_subtract7TeV",   "_mixed_integrated",   "_no_mixed"};
Int_t isystCheck[] =  {         0,          1,            2,                3,             4,        5,          6,                   7,                   8,          9,        10,         11,              12,              13,             14,                15,                 16,                      17,      18,        19,            20,            21,              22,               23,             24,                    31,            32};
enum                  {kDefault=0, k_vertex, k_mixed_norm, k_LowMultScaleAS, k_NoLMScaling, k_CMSSub, k_ATLAS_IP, k_ATLAS_IP_02_to_05,  k_ATLAS_IP_0_to_02,  k_ATLAS_EP, k_0_to_02, k_02_to_05, k_hist_integral, k_parabolic_fit, k_constant_fit, k_baselineGausFit, k_baselineHighMult, k_baselineHighMultParab, k_no_v3, k_fit_opt, k_exclusion18, k_exclusion10, k_v2_analytical, k_no_subtraction, k_subtract7TeV, k_mixed_integrated=31, k_no_mixed=32};
const Int_t nsyst=24;
Int_t gStudySystematic = kDefault;
//Int_t gStudySystematic = k_mixed_integrated ;
//Int_t gStudySystematic = k_no_v3 ;
//Int_t gStudySystematic = k_LowMultScaleAS;
//Int_t gStudySystematic = k_NoLMScaling;
//Int_t gStudySystematic = k_no_subtraction;
//Int_t gStudySystematic = k_CMSSub;
//Int_t gStudySystematic = k_ATLAS_IP;
//Int_t gStudySystematic = k_ATLAS_EP;
//Int_t gStudySystematic = k_subtract7TeV;
//Int_t gStudySystematic = k_exclusion18;

//specific correlation.C settings
Bool_t useMCforTracklets=0;
Bool_t cutPtMCforTracklets=0;
Bool_t atLeastOnePair = 0;

//specific calc_sum_of_ratios.C settings
Int_t debug=0; //plots for debugging
Int_t SaveMixed=0; //save mixed event distribution
Int_t SaveSame=0; //same same event distribution
const Double_t gCutStatErrorSM=-1;
Char_t drawOption[]="surf1";//"colz";

//subtraction
const char* kCorrFuncTitle     = "#frac{1}{#it{N}_{trig}} #frac{d^{2}#it{N}^{uncorr}_{assoc}}{d#Delta#etad#Delta#varphi} (rad^{-1})";
const char* kProjYieldTitlePhi = "1/#it{N}_{trig} d#it{N}^{uncorr}_{assoc}/d#Delta#varphi per #Delta#eta (rad^{-1})";
const char* kProjYieldTitleEta = "1/#it{N}_{trig} d#it{N}^{uncorr}_{assoc}/d#Delta#eta per #Delta#varphi (rad^{-1})";
const char* kTitlePhi          = "#Delta#varphi (rad)";
const char* kTitleEta          = "#Delta#eta";
const char * fitOption = "I0Q";
const char* kSystemEnergy = "p-p #sqrt{s} = 13 TeV";
Float_t fontSize = 0.04;
const Char_t* format = ".png";

//syst_mu_tkl
TString centMethods[] = { "V0M","TKL"};
enum                    {kV0M=0, kTKL,kNcentMethods};
const Int_t nm = 2;//kNcentMethods;
const Int_t firstbinc=3;



#endif

