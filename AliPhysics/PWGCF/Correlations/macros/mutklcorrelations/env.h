#ifndef env_h
#define env_h 1
#include "TPad.h"
const Int_t np = 2;
const Int_t nbins_tkl      = 5; Double_t bins_tkl[]      = {0,1,2,3,4,5};
const Int_t nbins_trk      = 1; Double_t bins_trk[]      = {0.5,4};
const Int_t nbins_muon_tkl = 5; Double_t bins_muon_tkl[] = {0.5,1,1.5,2,3,4};
const Int_t nbins_muon_trk = 3; Double_t bins_muon_trk[] = {0.5,1,2,4};
const Int_t nbins_muon_ampt_gen = 8; Double_t bins_muon_ampt_gen[] = {0., 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0};
TString period[] = {"bcde","f"};
TString centMethods[] = { "LEE","LEW","TIM","TIC",              "V0S","V0D","V0M","ZDC","CL1","U0S","ET1","ET2","ET3","V0L",             "AMA","AMP","AMK","AMM","AMR","DPM"};
enum                    {kLEE=0, kLEW, kTIM, kTIC, kNtrk = kTIC, kV0S, kV0D, kV0M, kZDC, kCL1, kU0S, kET1, kET2, kET3, kV0L,kNtkl = kV0L, kAMA, kAMP, kAMK, kAMM, kAMR, kDPM,kNcentMethods};
//const Int_t nm = 14;//kNcentMethods;
const Int_t nm = kNcentMethods;
TString systCheck[] = {""        ,"_vertex05","_mixed_norm","_LowMultScale","_60_to_70","_hist_integral","_parabolic_fit","_constant_fit","_baselineGausFit","_baselineHighMult","_baselineHighMultParab","_no_v3","_fit_opt","_exclusion10","_exclusion08","_v2_analytical","_no_subtraction","_vertex01","_vertex02","_vertex03","_vertex04",   "_mixed_integrated",  "_no_mixed"};
Int_t isystCheck[] =  {         0,          1,            2,              3,          4,              5 ,               6,              7,                 8,                  9,                      10,      11,        12,            13,            14,              15,               16,         17,         18,         19,         20,                    30,           31};
enum                  {kDefault=0, k_vertex05, k_mixed_norm, k_LowMultScale, k_60_to_70, k_hist_integral, k_parabolic_fit, k_constant_fit, k_baselineGausFit, k_baselineHighMult, k_baselineHighMultParab, k_no_v3, k_fit_opt, k_exclusion10, k_exclusion08, k_v2_analytical, k_no_subtraction, k_vertex01, k_vertex02, k_vertex03, k_vertex04, k_mixed_integrated=30,k_no_mixed=31};
const Int_t nsyst=17; //no mixed integrated
const char* kCorrFuncTitle     = "#frac{1}{#it{N}_{trig}} #frac{d^{2}#it{N}^{uncorr}_{assoc}}{d#Delta#etad#Delta#varphi} (rad^{-1})";
const char* kProjYieldTitlePhi = "1/#it{N}_{trig} d#it{N}^{uncorr}_{assoc}/d#Delta#varphi per #Delta#eta (rad^{-1})";
const char* kProjYieldTitleEta = "1/#it{N}_{trig} d#it{N}^{uncorr}_{assoc}/d#Delta#eta per #Delta#varphi (rad^{-1})";
const char* kTitlePhi          = "#Delta#varphi (rad)";
const char* kTitleEta          = "#Delta#eta";
//const char* kSystemEnergy = "p-p #sqrt{s} = 8 TeV";
const char* kSystemEnergy = "p-Pb #sqrt{s_{NN}} = 5.02 TeV";
//const char* kSystemEnergy = "Pb-p #sqrt{s_{NN}} = 5.02 TeV";
Int_t color[]={1,2,3,4,6,7,8,9,11,12,13,14,15,1,2,3,4,5,6};
Int_t marker[]={kFullSquare,kFullTriangleUp,kFullTriangleDown,kOpenCircle,kOpenSquare,kFullCircle,kOpenTriangleUp,kOpenTriangleDown,kFullTriangleUp,kFullTriangleDown,kOpenCircle,kOpenSquare,kFullCircle,kOpenTriangleUp,kOpenTriangleDown,kOpenCircle,kOpenSquare};
Float_t fontSize = 0.04;
const char * fitOption = "I0Q";
const Char_t* format = ".png";
//Int_t gStudySystematic = k_LowMultScale;
//Int_t gStudySystematic = k_constant_fit;
//Int_t gStudySystematic = k_parabolic_fit;
//Int_t gStudySystematic = k_hist_integral;
//Int_t gStudySystematic = k_no_subtraction;
Int_t gStudySystematic = kDefault;
const Float_t pi = TMath::Pi();
enum {kTrackDeta=0,kTrackPtAs=1,kTrackPtTr=2,kTrackCent=3,kTrackDphi=4,kTrackZvtx=5,kTrackNvar=6};
enum {kEventPtTr=0,kEventCent=1,kEventZvtx=2,kEventNvar=3};

void PadFor2DCorr(){
  gPad->SetPad(0, 0, 1, 1);
  gPad->SetLeftMargin(0.17);
  gPad->SetTopMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.05);
  gPad->SetTheta(55);
  gPad->SetPhi(45);
}

#endif
