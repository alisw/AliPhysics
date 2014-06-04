#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TClonesArray.h"
#include "AliParticleYield.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include <fstream>
#include "TLatex.h"
#include "TLegend.h"
#include "TList.h"
#include "TF1.h"
#include "AliPWGHistoTools.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TPaveText.h"
#include <map>
#include <fstream>
#include <istream>
#include "TMarker.h"
#include "TObjString.h"
#include "TLegendEntry.h"

#endif

// Plots ratios for QM and saves input files for thermal models

#if 1  // DUMMY IFDEF USED TO HIDE PREAMBLE in EMACS

enum MyParticles { kPDGPi = 211, kPDGK = 321, kPDGProton = 2212, kPDGKS0 = 310, kPDGLambda=3122, kPDGXi=3312,kPDGOmega=3334,kPDGPhi=333,kPDGKStar=313,kPDGDeuteron=1000010020,kPDGHE3 = 1000020030, kPDGHyperTriton = 1010010030, kPDGSigmaStarPlus=3224,kPDGSigmaStarMinus=3114,kPDGXiStar=3324};

typedef enum {kStatError, kSystError, kTotalError} myerror_t;

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle, Int_t icolor, Int_t imarker = kOpenSquare, Int_t errorsType = kTotalError, Float_t shift = 0) ;
TH1F * GetHistoYields(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle, Int_t icolor, Int_t imarker = kOpenSquare, Int_t errorsType = kTotalError, Float_t shift = 0) ;
void   PrepareThermalModelsInputFiles(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, Bool_t separateCharges=0) ;
void SetStyle(Bool_t graypalette=0) ;
TLegend * NewLegendQM(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Bool_t isYield = 0) ;
void DrawRatio(TString what, Bool_t isYield = kFALSE, Double_t shiftloc=0.);

void DrawFrame(Bool_t yields = 0) ;
void DrawExtrapolatedSymbolsAndLegendPbPb0010() ;
void DrawExtrapolatedSymbolsAndLegendpPb0005() ;
void DrawMarkerKStarNoFit(Bool_t plotLegend = 0) ;
void DrawMarkerNucleiNoFit() ;
void DrawExtrapNotInFitpPb0005(Bool_t drawExtrap = 1) ;

void DrawExtrapolatedSymbolsYieldsPbPb0010(Double_t x1=0.144578, Double_t y1=0.408249, Double_t x2=0.351406, Double_t y2=0.542403, Bool_t plotExtraploatedLegend=1);
Float_t shiftRatioDataModel = 0;
Double_t GetGraphRatioAndStdDev(TGraphErrors * gModel, TGraphErrors * &gRatio, TGraphErrors *&gStdDev) ;
TString particlesToExcludeFromChi2 = "";// The above method recomputes the chi2. This string is used to store values to be excluded from this calculation. Each PDG cocde has to be enclosed in [..]. This is obviously not efficient as it involves many string conversions. But efficiency is not an issue here.


void LoadArrays() ;
//void AddLabel(Float_t x, Float_t y, TString text);
void myLatexDraw(TLatex *currentLatex, Float_t currentSize=0.5, Int_t currentColor=1);
void myPaveSetup(float rRatio=0, float rRange3=0, float rRange5=0,
		 int rFillColor=0);
void myPadSetUp(TPad *currentPad);
TGraphErrors*  PlotThermusYields(const char * filename, Int_t color, Int_t lineStyle,
                                 const char * tag);
TGraphErrors * PlotGSIYields(const char * fileName, Int_t color=kBlack, Int_t lineStyle = kSolid,
                             const char * tag ="", Bool_t isPbPb = 1);
TGraphErrors*  PlotFlorenceYields(const char * filename, Int_t color, Int_t lineStyle,
                                  const char * tag) ;

void AddLineToThermalLegend(TObject * obj, TString line, const char * optFirst = "L");

void SaveCanvas(const char * name) ;


// Ratios to be draw. Remember to change the labels in DrawFrame if you change this
const Int_t nratio = 10;
Int_t num  [nratio]            = {kPDGK  , kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron , kPDGHE3      , kPDGHyperTriton , kPDGPhi , kPDGKStar};
Int_t denum[nratio]            = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGProton   , kPDGDeuteron , kPDGPi          , kPDGK   , kPDGK};
Int_t isSum[nratio]            = {1      , 1          ,  1         ,   1     , 1         , 1            , 0            , 1               , 1       , 1      };

// const Int_t nratio = 10;
// Int_t num  [nratio]            = {kPDGK  , kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron , kPDGHE3      , kPDGHyperTriton , kPDGPhi , kPDGKStar};
// Int_t denum[nratio]            = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGProton   , kPDGDeuteron , kPDGPi          , kPDGK   , kPDGK};
// Int_t isSum[nratio]            = {1      , 1          ,  1         ,   1     , 1         , 1            , 0            , 1               , 1       , 1      };

const char * ratiosLabels[]          = {"#frac{K^{+}+K^{-}}{#pi^{+}+#pi^{-}}", 
                                        "#frac{p+#bar{p}}{#pi^{+}+#pi^{-}}", 
                                        "#frac{2#Lambda}{K_{S}^{0}}", 
                                        "#frac{#Xi^{-}+#Xi^{+}}{#pi^{+}+#pi^{-}}",
                                        "#frac{#Omega^{-}+#Omega^{+}}{#pi^{+}+#pi^{-}}",
                                        "#frac{d}{p+#bar{p}}",
                                        "#frac{{}^{3}He }{d}",
                                        "#frac{{}^{3}_{#Lambda}H+{}^{3}_{#Lambda}#bar{H} }{#pi^{+}+#pi^{-}}",
                                        "#frac{#phi}{K^{+}+K^{-}}",
                                        "#frac{K*+#bar{K}*}{K^{+}+K^{-}}",};
//static const Double_t scale[]  = {1      , 3          ,  0.5       ,  30     ,  250      , 50           , 100          , 4e5             , 2       , 1      };
static const Double_t scale[]  = {1      , 3          ,  0.5       ,  80     ,  1000      , 50           , 100          , 4e5             , 2       , 1      };
//static const Double_t scale[]  = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,};
const Int_t npart = 12;
Int_t particleYields  [npart] = {kPDGPi ,kPDGK   ,kPDGKS0, kPDGKStar, kPDGPhi, kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron, kPDGHyperTriton, kPDGHE3    };
Int_t isSumYields[npart]      = {1      ,1       ,0      , 1        , 0      , 1           ,1            ,1        ,1          ,0           , 1              , 0          };
//Int_t isSumInputFiles[npart]  = {1      ,1       ,0      , 1        , 0      , 1           ,1            ,1        ,1          ,0           , 1              , 0          };
const char * yieldsLabels[]          = {"#frac{#pi^{+}+#pi^{-}}{2}",
                                        "#frac{K^{+}+K^{-}}{2}",
                                        "K_{S}^{0}",
                                        "#frac{K*+K*}{2}",
                                        "#phi",
                                        "#frac{p+#bar{p}}{2}",
                                        "#Lambda", 
                                        "#frac{#Xi^{-}+#Xi^{+}}{2}",
                                        "#frac{#Omega^{-}+#Omega^{+}}{2}",
                                        "d",
                                        "#frac{{}^{3}_{#Lambda}H+{}^{3}_{#Lambda}#bar{H}}{2}",
                                        "{}^{3}He",
};

//
Int_t markerNoFit  = 28;
Int_t markerExtrap = 27;
Double_t maxy = 0.5;


// Data arrays;
TClonesArray *arrPbPb=0, *arrpp7=0, *arrpPb=0, * arrpp276=0, * arrpp900=0, * arrThermus=0;
TClonesArray *arrSTARPbPb=0, *arrPHENIXPbPb=0, *arrBRAHMSPbPb=0;
TClonesArray *arrSTARpp  =0, *arrPHENIXpp=0;

//const Double_t *scaleRatios = 0;
Double_t *correlatedUnc = 0;
//Double_t correlatedUncLocalPbPb[14] = {0.0424 , 0.0424     ,  0.041     ,  0      , 0         , 0.0424       , 0.0424       , 0               , 0.05    , 0.05   };
Double_t correlatedUncLocalPbPbOnlyKStarPhi[14] = {0. , 0.     ,  0.     ,  0      , 0         , 0.       , 0.       , 0               , 0.05    , 0.05   };
Double_t correlatedUncLocalPP  [14] = {0.0424 , 0.0424     ,  0.041     ,  0      , 0         , 0.0424       , 0.0424       , 0               , 0.0424    , 0.0424   };
Double_t correlatedUncZero[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

TCanvas *myCan = 0;
TPad    *myPadStdDev =0;
TPad    *myPadRatio  =0;
TPad    *myPadHisto  =0;
TPad    *myPadLabel  =0;
TLegend * legThermal = 0;
#endif

Bool_t saveCanvas = 1; // if true, the canvas is saved and copied to the analysis note folder

TClonesArray * PlotRatiosForQM14() {
#if !(!defined (__CINT__) || (defined(__MAKECINT__)))
  LoadLibs();
#endif

  //
  LoadArrays();

  // Uncomment stuff in this section to save the inputs for thermal models
  //#define SAVE_INPUT_THERMAL_MODEL
#ifdef SAVE_INPUT_THERMAL_MODEL
    PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/1);
  // PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/1);
  // PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A0005", /*separateCharges*/1);
  // PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A2040", /*separateCharges*/1);
  // PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A6080", /*separateCharges*/1);
  // PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M6080", /* separateCharges*/1);
  //  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M2040", /*separateCharges*/1);

  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A0005", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A2040", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A6080", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M6080", /*separateCharges*/0);  
  // PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M2040", /*separateCharges*/0);

  return 0;
#endif


  SetStyle();

  DrawRatio("allpp");  
  //  DrawRatio("allppWithRHIC");  
  //  DrawRatio("PbPbWithPP7TeV");
  //DrawRatio("allsyst");
  //DrawRatio("PbPb6080andpPb0005");
  //  DrawRatio("pp_vsRHIC");
  //  DrawRatio("PbPb_vsRHIC");
  //DrawRatio("aliceall");


  // Yields and FITS
  //  maxy=2000;

  // DrawRatio("fit_ReferenceFit_PbPb0010", 1);
  //  DrawRatio("fit_ReferenceFit_GSIONLY_PbPb0010", 1);
  //  DrawRatio("fit_ReferenceFit_GSITHERMUS_PbPb0010",1);
  // DrawRatio("fitSHARE_NoPionsNoProtons_PbPb0010",1);
  //  DrawRatio("fitGSI_NoPionsNoProtons_PbPb0010", 1);
  // DrawRatio("fitShare_All_PbPb0010", 1);

  //  DrawRatio("fitShareWithWithoutNuclei_PbPb0010", 1);
  // maxy=200;
  // DrawRatio("fitGSI_PbPb6080", 1);
  //   maxy=150;
  //   DrawRatio("fitGSI_PbPb2040", 1);
  //  maxy = 60;
  // DrawRatio("fitThermus_GammaSFree_pPb0005");
  //  DrawRatio("fitShare_pPb0005");
  //  DrawRatio("fitShare_pPb0005_NoOmega", 1);
  //  maxy=20;
  //  DrawRatio("fitThermus_GammaSFree_pPb2040");
  //  maxy=9;
  //  DrawRatio("fitThermus_GammaSFree_pPb6080");
  //  maxy=9;
  //    DrawRatio("fitGSI_pp");
  //  DrawRatio("fitFlorence_pp");

  //  NewLegendQM();
  return arrPbPb;
}

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle, Int_t icolor, Int_t imarker, Int_t errorType, Float_t shift) {
  // FIXME: THIS SHOULD BE REVIEWED TO MAKE SURE THE PLOTS ARE LABELLED CORRECTLY

  
  //  scaleRatios = scale;
  TH1F * h = new TH1F(Form("hRatio_%d_%0.0f_%s_%d",system,energy,centrality.Data(),errorType), histotitle, nratio, 1+shift, nratio+1+shift);

  TClonesArray * arrRatios = new TClonesArray ("AliParticleYield");// We save to disk the precomputed ratios
  Int_t iratioArr = 0;// Index used only to add particles to the array above

  //  Double_t isSum = -1; // if this is -1, then the sum criterion is ignored
  for(Int_t iratio = 1; iratio <= nratio; iratio++){
    std::cout << "******** " << num[iratio-1] << "/" <<  denum[iratio-1]<< " ("<<isSum[iratio-1]<<")*******" << std::endl ;

    AliParticleYield * ratio = AliParticleYield::FindRatio(arr, num[iratio-1], denum[iratio-1], system, energy, centrality,isSum[iratio-1]);
    if(ratio) 
      {
        ratio = new AliParticleYield(*ratio); // We need to clone it to avoid a mess if we need to use this particle again later (e.g. double scalings)
        std::cout << "  " ;        
        ratio->Print("short");
      }
    
    if(!ratio) {
      // If the ratio is not found, try to build it!
      std::cout << "  Looking for " <<  num[iratio-1] << " ("<<isSum[iratio-1]<<")"<<std::endl;
      AliParticleYield * part1 = AliParticleYield::FindParticle(arr, num[iratio-1], system, energy, centrality,  isSum[iratio-1]);
      if(part1) {
        part1 = new AliParticleYield(*part1); // We need to clone it to avoid a mess if we need to use this particle again later
        if(isSum[iratio-1] && part1->IsTypeAverage()) {
          std::cout << "Sum requested, found average, scaling x2" << std::endl;        
          part1->Scale(2.);
        } 
      }
      // Try with the !sum, if part 1 is not found
      if(!part1) {
        std::cout << "  Looking for " <<  num[iratio-1] << " ("<<!isSum[iratio-1]<<")"<<std::endl;
        part1 = AliParticleYield::FindParticle(arr, num[iratio-1], system, energy, centrality,!isSum[iratio-1]);
        if(part1) 
          {
            part1 = new AliParticleYield(*part1); // We need to clone it to avoid a mess if we need to use this particle again later
            // If the sum was requested, try to recover it!
            if(isSum[iratio-1]) { 
              std::cout << "  Looking for " <<  -num[iratio-1] <<std::endl;
              // If it's lambda and ALICE, use 2L instead of L + Lbar // FIXME: do the same for deuterons?
              if((num[iratio-1] == kPDGLambda || num[iratio-1] == kPDGDeuteron) && energy > 300) {
                std::cout << "   It's Lambda or deuteron ALICE: Scaling x2 " << std::endl;
                part1->Print();
                part1->Scale(2.);
              } else {
                AliParticleYield * part1_bar = AliParticleYield::FindParticle(arr, -num[iratio-1], system, energy, centrality,0);
                if(part1 && part1_bar) {
                  std::cout << "Adding " << part1_bar->GetPartName() << " " << part1->GetPartName() << std::endl;            
                  part1 = AliParticleYield::Add(part1, part1_bar);
                  
                } else if (TDatabasePDG::Instance()->GetParticle(-num[iratio-1])){ // Before scaling x2 check it it makes sense (antiparticle exists) 
                  std::cout << "WARNING: Sum requested but not found, scaling x2 " << part1->GetName() << std::endl;
                  part1->Scale(2);
                }
              }
            } else if(part1) { // if we are here, it means the sum was *not* requested (isSum=0), but we found something with (!isSum) = 1
              // if the not sum was requested, but the sum is found, divide by 2 so that it is comparable
              std::cout << "WARNING: Using sum/2 for " << part1->GetName() << std::endl;
          
              part1->Scale(0.5);
            }
        }
 
      }


      std::cout << "  Looking for " <<  denum[iratio-1] << " ("<<isSum[iratio-1]<<")"<<std::endl;
      AliParticleYield * part2 = AliParticleYield::FindParticle(arr, denum[iratio-1], system, energy, centrality,isSum[iratio-1]);
      if(part2) {
        part2 = new AliParticleYield(*part2); // We need to clone it to avoid a mess if we need to use this particle again later
        if(isSum[iratio-1] && part2->IsTypeAverage()) {
          std::cout << "Sum requested, found average, scaling x2" << std::endl;        
          part2->Scale(2.);
        } 
      }
      if(!part2) {// Try with the !sum, if part 2 is not found
        std::cout << "  Looking for " <<  denum[iratio-1] << " ("<<!isSum[iratio-1]<<")"<<std::endl;
        part2 = AliParticleYield::FindParticle(arr, denum[iratio-1], system, energy, centrality,!isSum[iratio-1]);
        if(part2) 
          {
            part2 = new AliParticleYield(*part2); // We need to clone it to avoid a mess if we need to use this particle again later
            if(isSum[iratio-1]) { 
              std::cout << "  Looking for " <<  -denum[iratio-1] << std::endl;
              AliParticleYield * part2_bar = AliParticleYield::FindParticle(arr, -denum[iratio-1], system, energy, centrality,0);
              if(part2 && part2_bar){
                std::cout << "Adding " << part2_bar->GetPartName() << " " << part2->GetPartName() << std::endl;            
                part2 = AliParticleYield::Add(part2, part2_bar);
              } else if (TDatabasePDG::Instance()->GetParticle(-denum[iratio-1])){ // Before scaling x2 check it it makes sense (antiparticle exists) 
                std::cout << "WARNING: Sum requested but not found, scaling x2 " << part2->GetName() << std::endl;
                part2->Scale(2);
              }
            } else if(part2){
              // if the not sum was requested, but the sum is found, divide by 2 so that it is comparable
              std::cout << "WARNING: Using sum/2 for " << part2->GetName() << std::endl;
              part2->Scale(0.5);
            } 
          }

      }
      ratio = AliParticleYield::Divide(part1, part2, correlatedUnc[iratio-1], "YQ"); // Assume by that the systematics of part1 and part2 are uncorrelated.

      if(ratio) {
        std::cout << "" << std::endl;
        std::cout << "WARNING: building ratio " << num[iratio-1] <<"/"<<denum[iratio-1]<<": Check uncertainties!!" << std::endl;
        ratio->Print("");
        std::cout << "" << std::endl;
      }
    }
    if(ratio){
      ratio->Scale(scale[iratio-1]);
      h->SetBinContent(iratio, ratio->GetYield());
      if(errorType == kTotalError) {
        h->SetBinError  (iratio, ratio->GetTotalError(0/* 0 = no normalization error */));
      } else if (errorType == kStatError) {
        h->SetBinError  (iratio, ratio->GetStatError());
      } else if (errorType == kSystError) {
        h->SetBinError  (iratio, ratio->GetSystError());
      } else {
        std::cout <<  "ERROR: Unknown Error Type " << errorType << std::endl;
      }

      //      h->GetXaxis()->SetBinLabel(iratio, Form("#splitline{%s}{%s}",Form("#times%2.2f",  scale[iratio-1]), ratio->GetLatexName()));
      h->GetXaxis()->SetBinLabel(iratio, ratio->GetLatexName());
      new ((*arrRatios)[iratioArr++]) AliParticleYield(*ratio);
    }
    else {
      h->GetXaxis()->SetBinLabel(iratio, Form("#frac{%d}{%d}",num[iratio-1], denum[iratio-1]));
      
    }
    std::cout << "*** END OF " << num[iratio-1] << "/" <<  denum[iratio-1]<< " *******" << std::endl ;

  }
  
  h->GetYaxis()->SetRangeUser(0, maxy);
  h->SetLineColor(icolor);
  h->SetMarkerColor(icolor);
  h->SetMarkerStyle(imarker);
  // the "if" avoids saving twice the same ratios
  if(errorType == kSystError) AliParticleYield::SaveAsASCIIFile(arrRatios, TString("ratios_")+h->GetName());
  return h;



}

void   PrepareThermalModelsInputFiles(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, Bool_t separateCharges)  {
  // If "Separate charges" is true, tries to dump both charges are dumped
  TClonesArray * arrOut = new TClonesArray("AliParticleYield");
  TClonesArray * arrOutGSI = new TClonesArray("AliParticleYield"); // We add dummy lines to the GSI output file if needed!
;

  Int_t ipartOut = 0; // Index for the array
  Int_t ipartOutGSI = 0; // Index for the array
  
  for(Int_t ipart = 0; ipart < npart; ipart++){
    if(!separateCharges) {
      AliParticleYield * part = AliParticleYield::FindParticle(arr, particleYields[ipart], system, energy, centrality,  isSumYields[ipart]);
      if(!part && isSumYields[ipart]) {
        //Could not find the particle, but the sum was requested: build the sum!
        part = AliParticleYield::FindParticle(arr, particleYields[ipart], system, energy, centrality,  0);
        AliParticleYield * part2 = AliParticleYield::FindParticle(arr, -particleYields[ipart], system, energy, centrality,  0);
        if(part2 && part) part = AliParticleYield::Add(part, part2);        
        else if(part) part->Scale(2.); // If we only found a particle, we can scale it by a factor 2.
        else part = 0;
      }
      // We want to save the average of particle and antiparticle in this case
      if(part) {
        if(isSumYields[ipart] && !part->IsTypeAverage()) part->Scale(0.5); // If it's not already an average, but just a sum, divide by 2
        new((*arrOut   )[ipartOut++]) AliParticleYield(*part);
        new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(*part);
      } else { // Add dummy particle to the GSI list
        new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(particleYields[ipart], system, energy, -10, -10, -10, -10, -10, -10, 5, 256, "DUMMY", 1, "ALICE");        
      }
    }
    else {
      // ignore isSumYields and try to find both particleYields
      Bool_t notFound = 0;
      AliParticleYield * part = AliParticleYield::FindParticle(arr, particleYields[ipart], system, energy, centrality,  0);
      if(part) {
	new((*arrOut)[ipartOut++]) AliParticleYield(*part);
	new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(*part);
      }
      else {
	// std::cout << "ADDING DUMMY part " << particleYields[ipart] << std::endl;	
	// new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(particleYields[ipart], system, energy, -10, -10, -10, -10, -10, -10, 5, 256, "DUMMY", 1, "ALICE");        
	notFound=1;
      }
      // Try to find antiparticle (-pdg code)
      part = AliParticleYield::FindParticle(arr, -particleYields[ipart], system, energy, centrality,  0);
      if(part) {
        new((*arrOut)[ipartOut++]) AliParticleYield(*part);
        new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(*part);
      }
      else if (notFound) {
        // If neither charge was found, check if we at least have the sum 
        part = AliParticleYield::FindParticle(arr, abs(particleYields[ipart]), system, energy, centrality,  1);	
        if (part) { 
          if(!part->IsTypeAverage()) part->Scale(0.5); // If it's a sum (not an average) divide by 2
          new((*arrOut)[ipartOut++]) AliParticleYield(*part);
          new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(*part);
	  if(TDatabasePDG::Instance()->GetParticle(-particleYields[ipart]) && 
	     (particleYields[ipart] != kPDGLambda) && 
	     (particleYields[ipart] != kPDGKStar)
	     ){// if only the sum was found, add a dummy entry to the
	       // GSI file, so that anton always has the same # of
	       // lines. However, we don't do this for the Lambda and
	       // KStar (always one of the charges or the average)
	    new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(-particleYields[ipart], system, energy, -10, -10, -10, -10, -10, -10, 5, 256, "DUMMY", 1, "ALICE");        
	  }
        }
        else {
	  std::cout << "ADDING DUMMY sum " << particleYields[ipart] << std::endl;
          new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(particleYields[ipart], system, energy, -10, -10, -10, -10, -10, -10, 5, 256, "DUMMY", 1, "ALICE");        
	  if (particleYields[ipart]==kPDGHyperTriton) {
	    // If is the 3LH, add another one for the antiparticle
	    AliParticleYield(-particleYields[ipart], system, energy, -10, -10, -10, -10, -10, -10, 5, 256, "DUMMY", 1, "ALICE");        
	  }
        }
      }
      
    }
  }
  std::cout << "Particles for thermal model fits:" << std::endl; 
  arrOut->Print("short");
  //  arrOut->Print("");
  std::cout << "" << std::endl;
  // Write GSI input file
  TIter it(arrOutGSI);
  AliParticleYield * part = 0;
  ofstream fout(Form("gsi_System_%d_Energy_%0.0f_Centr_%s_BothCharges_%d", system, energy, centrality.Data(), separateCharges));
  while ((part = (AliParticleYield*) it.Next())){
    fout << part->GetYield() << " " << part->GetTotalError() << std::endl;
  }
  fout.close();
  // Write thermus file
  AliParticleYield::WriteThermusFile(arrOut, Form("thermus_System_%d_Energy_%0.0f_Centr_%s_BothCharges_%d", system, energy, centrality.Data(), separateCharges));
}


TH1F * GetHistoYields(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle, 
                      Int_t icolor, Int_t imarker, Int_t errorsType, Float_t shift) {

  TH1F * h = new TH1F(Form("hPart_%d_%0.0f_%s",system,energy,centrality.Data()), histotitle, npart, 1+shift, npart+1+shift);

  for(Int_t ipart = 1; ipart <= npart; ipart++){
    std::cout << "----- Searching " << particleYields[ipart-1] << " -------"  << std::endl;
    
    AliParticleYield * part = AliParticleYield::FindParticle(arr, particleYields[ipart-1], system, energy, centrality,isSumYields[ipart-1]);
    if(part) {
      std::cout << "found" << std::endl;
      part->Print();
    }
    if(!part && isSumYields[ipart-1]) {
      //Could not find the particle, but the sum was requested: build the sum!
      part = AliParticleYield::FindParticle(arr, particleYields[ipart-1], system, energy, centrality,  0);
      AliParticleYield * part2 = AliParticleYield::FindParticle(arr, -particleYields[ipart-1], system, energy, centrality,  0);
      if(part2 && part) {
	std::cout << " Building sum" << std::endl;
	part->Print();
	part2->Print();
	part = AliParticleYield::Add(part, part2);        
      }
      else if(part) {
	std::cout << "Scaling part" << std::endl;
	part->Print();
	part = new AliParticleYield(*part); // Always clone before scaling
	part->Scale(2.); // If we only found a particle, we can scale it by a factor 2.
      }
      else part = 0;
    }    
    if(!part){
      std::cout << "Cannot find " << particleYields[ipart-1] << std::endl;
      continue;
    }
    if(isSumYields[ipart-1] && !part->IsTypeAverage()) {
      std::cout << " scaling /2" << std::endl;
      part = new AliParticleYield(*part); // Always clone before scaling
      part->Scale(0.5); // take average
    }
    std::cout << " Plotting " << particleYields[ipart-1] << std::endl;
    part->Print();
    //    part->Scale(scale[ipart-1]);
    h->SetBinContent(ipart, part->GetYield());
    if(errorsType == kTotalError) {
      h->SetBinError  (ipart, part->GetTotalError(0/* 0 = no normalization error */));
    } else if (errorsType == kSystError) {
      h->SetBinError  (ipart, part->GetSystError());
    } else if (errorsType == kStatError) {
      h->SetBinError  (ipart, part->GetStatError());
    }
    h->GetXaxis()->SetBinLabel(ipart, part->GetLatexName());

  }
  h->SetMarkerStyle(imarker);
  h->SetMarkerColor(icolor);
  h->SetLineColor(icolor);
  h->SetMarkerSize(1.4);
  
  return h;
}

void LoadArrays() {
  arrPbPb = AliParticleYield::ReadFromASCIIFile("PbPb_2760_Cascades.txt");
  arrPbPb->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_DeuHelium3.txt"));
  arrPbPb->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_Hypertriton.txt"));
  arrPbPb->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_Kstar892.txt"));
  arrPbPb->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_LambdaK0.txt"));
  arrPbPb->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_PiKaPr.txt"));
  arrPbPb->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_phi1020.txt"));
  arrPbPb->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_AveragedNumbers.txt"));

  arrpp7   = AliParticleYield::ReadFromASCIIFile("pp_7000.txt");

  arrpp276 = AliParticleYield::ReadFromASCIIFile("pp_2760.txt");
  arrpp900 = AliParticleYield::ReadFromASCIIFile("pp_900.txt");

  arrpPb   = AliParticleYield::ReadFromASCIIFile("pPb_5020_MultiStrange.txt");
  arrpPb->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("pPb_5020_PiKaPrLamndaK0.txt"));
  arrpPb->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("pPb_5020_deuteron.txt"));
  arrpPb->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("pPb_5020_AveragedNumbers.txt"));
  arrpPb->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("pPb_5020_phi.txt"));
  arrpPb->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("pPb_5020_Kstar.txt"));
  arrThermus = AliParticleYield::ReadFromASCIIFile("PbPb_2760_Thermus_Boris_20140407.txt");
  
  // RHIC data
  arrSTARPbPb   = AliParticleYield::ReadFromASCIIFile("PbPb_200_STAR-AntonQM12.txt");
  arrPHENIXPbPb = AliParticleYield::ReadFromASCIIFile("PbPb_200_PHENIX-AntonQM12.txt");
  arrBRAHMSPbPb = AliParticleYield::ReadFromASCIIFile("PbPb_200_BRAHMS-AntonQM12.txt");
  arrSTARpp     = AliParticleYield::ReadFromASCIIFile("pp_200_STAR.txt");
  arrPHENIXpp   = AliParticleYield::ReadFromASCIIFile("pp_200_PHENIX.txt");

}

void SetStyle(Bool_t graypalette) {
  std::cout << "Setting style!" << std::endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetDrawBorder(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(1);
  gStyle->SetLabelSize(0.045,"yz");
  gStyle->SetLabelSize(0.06,"x");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.25,"y");
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetTickLength(0.012,"X");  gStyle->SetTickLength(0.012,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);

  gStyle->SetErrorX(0);
  gStyle->SetEndErrorSize(5);
}

TLegend * NewLegendQM(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Bool_t isYield) {

  const char * style = "p";
  const char ** labels=0;
  //  Bool_t beautify=kFALSE;
  Bool_t useTitle=kTRUE;

  TLegend * l = new TLegend(x1, y1, x2, y2);
  l->SetFillColor(kWhite);
  l->SetTextFont(43);
  l->SetTextSize(25);
  l->SetBorderSize(1);
  l->SetLineWidth(1);
  l->SetMargin(0.1);
  // const Int_t markers[] = {20,24,21,25,23,28,33,20,24,21,25,23,28,33};
  // const Int_t colors[]  = {1,2,3,4,6,7,8,9,10,11,1,2,3,4,6,7,8,9,10};

  TList * list = gPad->GetListOfPrimitives(); 
  TIterator * iter = list->MakeIterator();
  TObject * obj = 0;
  Int_t ilabel = -1;
  while ((obj = (TObject*) iter->Next())){
    if (obj->InheritsFrom("TH1") || obj->InheritsFrom("TGraph") || obj->InheritsFrom("TF1")) {
      if( (TString(obj->GetName()) == "hframe" ) ) continue; 
      ilabel++;
      if (labels != 0)
	l->AddEntry(obj, labels[ilabel], style);
      else{
	if (useTitle)  {
          if(TString(obj->GetTitle()).Contains("NoLegend")) continue;
	  TString title = obj->GetTitle();
	  if(title.Contains("p-Pb")) {
	    TObjArray * tokens = title.Tokenize(",");
	    if(tokens) {
	      TString system = ((TObjString*) tokens->At(0))->String().Strip(TString::kBoth, ' ');
	      TString centr  = ((TObjString*) tokens->At(1))->String().Strip(TString::kBoth, ' ');
	      l->AddEntry(obj, system, style);
	      l->AddEntry(obj, centr, "0");
	      delete tokens;
	    }

	  } else {
	    l->AddEntry(obj, obj->GetTitle(), style);
	  }

        }
	else 
	  l->AddEntry(obj, obj->GetName(), style);	  
      }
      // if(beautify) {
	
      //   if(!obj->InheritsFrom("TF1")){
      //     ((TH1*)obj)->SetLineColor(colors[ilabel]);
      //     ((TH1*)obj)->SetMarkerColor(colors[ilabel]);
      //     ((TH1*)obj)->SetMarkerStyle(markers[ilabel]);
      //   } else {
      //     ((TF1*)obj)->SetLineColor(colors[ilabel]);
      //   }
      
    }
  }
  if(isYield) {
    // Add some details on excluded stuff
    l->SetTextSize(22);
    l->SetBorderSize(0);
    l->GetEntry()->SetOption("0");
    l->SetMargin(0.01);
  }
    l->Draw();
  return l;
}



void DrawRatio(TString what, Bool_t isYield, Double_t shift) {
  // This is used to simplify the code above
  // In order to draw syst error bars, we need to convert to graphs the syst errors histos
  // if isYield == true plots yields rather than ratios

  // Sample colors
  //  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2  , kWhite};

  TClonesArray * array = 0;
  Int_t system,  color, marker;
  Float_t energy = 0, shiftloc = shift;
  TString centrality, label;
  // FIXME: move this in the different sections below
  correlatedUnc = correlatedUncZero;
  std::cout << "Plotting " << what.Data() << std::endl;
  

  if (what == "frame" ) {
    correlatedUnc = correlatedUncZero;
    DrawFrame(isYield);
    // TH1 * h = GetHistoRatios(arrPbPb,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "NoLegend", kWhite);
    // h->Draw();
    // h->GetYaxis()->SetDecimals(1);
    // h->GetYaxis()->SetNdivisions(505);
    //    h->GetXaxis()->CenterLabels(1);
    //    Int_t nratio = h->GetNbinsX();
    if(!isYield) {
      for(Int_t iratio = 0; iratio < nratio; iratio++){
        Double_t exp = TMath::Floor(TMath::Log10(TMath::Abs(scale[iratio])));  
        Double_t man = scale[iratio] / TMath::Power(10, exp);
        if(exp > 2) {
          //        TLatex * scaleLabel = new TLatex(iratio+1+0.2,maxy*1.01, Form("#times %0.0f 10^{%0.0f}", man, exp));
          TLatex * scaleLabel = new TLatex(iratio+1+0.2,0.005, Form("#times %0.0f 10^{%0.0f}", man, exp));
          scaleLabel->SetTextFont(43);
          scaleLabel->SetTextSize(20);
          scaleLabel->Draw();
        } else {
          Double_t shiftloc2 = scale[iratio] < 50 ? 0.3 : 0.2;
          TLatex * scaleLabel = new TLatex(iratio+1+shiftloc2,  0.005, Form("#times %g", scale[iratio]));
          scaleLabel->SetTextFont(43);
          scaleLabel->SetTextSize(20);
          scaleLabel->Draw();
        }      
    
      }
    }
    if(isYield) {
      TPaveText *pt = new TPaveText(0.691767, 0.86069, 0.893574, 0.944865,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetTextAlign(12);
      pt->SetTextFont(43);
      pt->SetTextSize(23);
      pt->AddText("ALICE Preliminary");
      pt->Draw();
    }
    else {
      TPaveText *pt = new TPaveText(    0.176, 0.842881, 0.378514, 0.929595,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillColor(0);
      pt->SetLineWidth(1);
      pt->SetTextFont(43);
      pt->SetTextSize(23);
      pt->AddText("ALICE Preliminary");
      pt->Draw();
    }

    gPad->SetGridx();

  }
  else if (what == "PbPb_0010") {
    array = arrPbPb;
    system = 2; energy = 2760.; centrality = "V0M0010";
    label = "Pb-Pb #sqrt{s}_{NN} = 2.76 TeV, 0-10%";
    color = kRed+1;
    marker = kFullCircle;
    if(!shift)    shiftloc =  0;
    correlatedUnc = correlatedUncLocalPbPbOnlyKStarPhi;

  }

  else if (what == "PbPb_6080") {
    array = arrPbPb;
    system = 2; energy = 2760.; centrality = "V0M6080";
    label = "Pb-Pb #sqrt{s}_{NN} = 2.76 TeV, 60-80%";
    color  = kBlue+1;
    marker = kFullCircle;
    if(!shift)    shiftloc =  0.0;
    correlatedUnc = correlatedUncLocalPbPbOnlyKStarPhi;
  }
  else if (what == "PbPb_2040") {
    array = arrPbPb;
    system = 2; energy = 2760.; centrality = "V0M2040";
    label = "Pb-Pb #sqrt{s}_{NN} = 2.76 TeV, 20-40%";
    color  = kBlue+1;
    marker = kFullCircle;
    if(!shift)    shiftloc =  0.0;
    correlatedUnc = correlatedUncLocalPbPbOnlyKStarPhi;
  }
  else if (what == "pp7") {
    array = arrpp7;
    system = 0; energy = 7000.; centrality = "";
    label = "pp #sqrt{s} = 7 TeV";
    color  = kMagenta+1;
    marker = kFullCircle;
    if(!shift)    shiftloc =  0.2;
    //    correlatedUnc = correlatedUncLocalPP;
  }
  else if (what == "pp900") {
    array = arrpp900;
    system = 0; energy = 900.; centrality = "";
    label = "pp #sqrt{s} = 0.9 TeV";
    color  = kCyan+2;
    marker = kFullCircle;
    if(!shift)    shiftloc =  -0.2;
    //    correlatedUnc = correlatedUncLocalPP;  
  }
  else if (what == "pp276") {
    array = arrpp276;
    system = 0; energy = 2760.; centrality = "";
    label = "pp #sqrt{s} = 2.76 TeV";
    color  = kYellow+2;
    marker = kFullCircle;
    if(!shift)    shiftloc = 0;
    //    correlatedUnc = correlatedUncLocalPP;
  }
  else if (what == "pPb0005") {
    array = arrpPb;
    system = 1; energy = 5020.; centrality = "V0A0005";
    //    label = "p-Pb #sqrt{s}_{NN} = 5.02 TeV, V0A 0-5%";
    label = "p-Pb #sqrt{s}_{NN} = 5.02 TeV, V0A Multiplicity (Pb-Side) 0-5%";
    color  = kBlack;
    marker = kFullCircle;
    if(!shift)    shiftloc = -0.2;
    //    correlatedUnc = correlatedUncLocalPP;
  } 
  else if (what == "pPb2040") {
    array = arrpPb;
    system = 1; energy = 5020.; centrality = "V0A2040";
    label = "p-Pb #sqrt{s}_{NN} = 5.02 TeV, V0A Multiplicity (Pb-Side) 20-40%";
    color  = kBlue+1;
    marker = kFullCircle;
    if(!shift)    shiftloc = 0.;
    //    correlatedUnc = correlatedUncLocalPP;
  } 
  else if (what == "pPb6080") {
    array = arrpPb;
    system = 1; energy = 5020.; centrality = "V0A6080";
    label = "p-Pb #sqrt{s}_{NN} = 5.02 TeV, V0A Multiplicity (Pb-Side) 60-80%";
    color  = kGreen+3;
    marker = kFullCircle;
    if(!shift)    shiftloc = 0.;
    //    correlatedUnc = correlatedUncLocalPP;
  } 
  else if (what == "PbPbSTAR") {
    array = arrSTARPbPb;
    system = 2; energy = 200.; centrality = "0005";
    label = "STAR, Au-Au #sqrt{s}_{NN} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenStar;
    if(!shift)    shiftloc = +0.2;
    //    correlatedUnc = correlatedUncZero;
  }
  else if (what == "PbPbPHENIX") {
    array = arrPHENIXPbPb;
    system = 2; energy = 200.; centrality = "0005";
    label = "PHENIX, Au-Au #sqrt{s}_{NN} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenSquare;
    if(!shift)    shiftloc = -0.15;
    //    correlatedUnc = correlatedUncZero;
  }
  else if (what == "PbPbBRAHMS") {
    array = arrBRAHMSPbPb;
    system = 2; energy = 200.; centrality = "0010";
    label = "BRAHMS, Au-Au #sqrt{s}_{NN} = 0.2 TeV, 0-10%";
    color  = kBlack;
    marker = kOpenCross;
    if(!shift)    shiftloc = -0.3;
    //    correlatedUnc = correlatedUncZero;
  }
  else if (what == "ppSTAR") {
    array = arrSTARpp;
    system = 0; energy = 200.; centrality = "";
    label = "STAR, pp #sqrt{s} = 0.2 TeV";
    color  = kBlack;
    marker = kOpenStar;
    if(!shift)    shiftloc = -0.15;
    //    correlatedUnc = correlatedUncZero;
  }
  else if (what == "ppPHENIX") {
    array = arrPHENIXpp;
    system = 0; energy = 200.; centrality = "";
    label = "PHENIX, pp #sqrt{s} = 0.2 TeV";
    color  = kBlack;
    marker = kOpenSquare;
    if(!shift)    shiftloc = -0.2;
    //    correlatedUnc = correlatedUncZero;
  } 
  // From here on, it's meta names, to draw several series of ratios
  else if (what == "allpp"){
    DrawRatio("frame");    
    DrawRatio("pp900");
    DrawRatio("pp276");
    DrawRatio("pp7");
    array =0;
    NewLegendQM(0.62249, 0.635734, 0.910643, 0.94673);

    SaveCanvas("Ratios_pponly");
  }
  else if (what == "allppWithRHIC"){
    DrawRatio("frame");    
    DrawRatio("pp900", 0, 0.05);
    DrawRatio("pp276", 0, 0.20);
    DrawRatio("pp7", 0, 0.35);
    DrawRatio("ppPHENIX", 0, -0.35);
    DrawRatio("ppSTAR", 0, -0.15);

    array =0;
    
    NewLegendQM(0.588353, 0.636857, 0.910643, 0.948352);

    SaveCanvas("Ratios_pponly_withRHIC");
  }
  else if (what == "PbPbWithPP7TeV"){
    
    DrawRatio("frame");
    DrawRatio("PbPb_0010", 0, 0.1);
    DrawRatio("pp7", 0, -0.1);
    array =0;
    

    NewLegendQM(0.538153, 0.749397, 0.893574, 0.950362);
    DrawExtrapolatedSymbolsAndLegendPbPb0010();

    SaveCanvas("Ratios_withpp7tev");

  } else if(what == "allsyst") {

    DrawRatio("frame");
    DrawRatio("pp7", 0, -0.2);
    DrawRatio("pPb0005", 0, 0.00001);
    DrawRatio("PbPb_0010", 0, 0.2);
    array =0;


    NewLegendQM(0.462851, 0.631722, 0.89257, 0.936697);
    DrawExtrapolatedSymbolsAndLegendPbPb0010();
    DrawExtrapolatedSymbolsAndLegendpPb0005();

    SaveCanvas("Ratios_allsystems");

  }else if(what =="PbPb6080andpPb0005") {
    DrawRatio("frame");
    DrawRatio("pPb0005",0, -0.1);
    DrawRatio("PbPb_6080", 0, 0.1);
    array=0;

    NewLegendQM(    0.46988, 0.730036, 0.910643, 0.948736);
    DrawExtrapolatedSymbolsAndLegendpPb0005();
    SaveCanvas("Ratios_6080vspPb");
    
    
  }else if(what =="pp_vsRHIC") {
    DrawRatio("frame");
    DrawRatio("pp7");
    DrawRatio("ppPHENIX", 0, 0.00001);
    DrawRatio("ppSTAR");
    array=0;

    NewLegendQM(    0.554217, 0.677869, 0.910643, 0.948736);

    SaveCanvas("Ratios_vsRHIC_pp");

        
  } else if (what =="PbPb_vsRHIC") {
    DrawRatio("frame");
    DrawRatio("PbPb_0010");
    DrawRatio("PbPbSTAR");
    DrawRatio("PbPbPHENIX");
    DrawRatio("PbPbBRAHMS");
    array = 0;



    NewLegendQM(    0.434739, 0.591593, 0.939759, 0.936697);
    DrawExtrapolatedSymbolsAndLegendPbPb0010();
 
    SaveCanvas("Ratios_vsRHIC_PbPb");
  } else if (what == "aliceall") {
    DrawRatio("frame");
    DrawRatio("PbPb_0010");
    DrawRatio("PbPb_6080");
    DrawRatio("pPb0005");
    DrawRatio("pp7");
    DrawRatio("pp276");
    DrawRatio("pp900");
  }

  // FROM HERE: IT's yields
  else if( what == "fit_ReferenceFit_PbPb0010") {

    DrawRatio("frame",1);  
    DrawRatio("PbPb_0010",1);  
    particlesToExcludeFromChi2="[313]"; // do not consider K* 
    legThermal->SetNColumns(4);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #chi^{2}/NDF", "0");    
    // FIXME: sistemare valori rapporti   
    shiftRatioDataModel=-0.2;
    PlotThermusYields("lhc2760_final_0005_single_gc_output_gs1_wevc_nkst.txt" , kBlack     , kSolid      , "THERMUS 2.3 , 155 #pm 2 , 5924 #pm 543 , 23.6/9");
    shiftRatioDataModel =0;
    PlotGSIYields("data+therm_fit2_s2760_0-10qm14.dat"                        , kOrange-1  , kDashed     , "GSI      , 156 #pm 2  , 5330 #pm 505 , 17.4/9");
    shiftRatioDataModel = 0.2;
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_PbPb_0010_with_nuclei.txt"     , kBlue+1    , kDashDotted , "SHARE 3     , 156 #pm 3 , 4476 #pm 696 , 14.1/9");


    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);
    DrawMarkerKStarNoFit() ;
    DrawExtrapolatedSymbolsYieldsPbPb0010();

    SaveCanvas("Fit_PbPb0010_Reference");

    array =0;
  }  else if( what == "fit_ReferenceFit_GSIONLY_PbPb0010") {

    DrawRatio("frame",1);  
    DrawRatio("PbPb_0010",1);  
    particlesToExcludeFromChi2="[313]"; // do not consider K* 
    legThermal->SetNColumns(4);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #chi^{2}/NDF", "0");    
    // FIXME: sistemare valori rapporti   
    shiftRatioDataModel =0;
    PlotGSIYields("data+therm_fit2_s2760_0-10qm14.dat"                        , kBlack  , kSolid     , "GSI      , 156 #pm 2  , 5330 #pm 505 , 17.4/9");

    myPadHisto->cd();
    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);
    DrawMarkerKStarNoFit() ;
    DrawExtrapolatedSymbolsYieldsPbPb0010();

    SaveCanvas("Fit_PbPb0010_Reference_GSI");

    array =0;
  } else if( what == "fit_ReferenceFit_GSITHERMUS_PbPb0010") {

    DrawRatio("frame",1);  
    DrawRatio("PbPb_0010",1);  
    particlesToExcludeFromChi2="[313]"; // do not consider K* 
    legThermal->SetNColumns(4);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #chi^{2}/NDF", "0");    
    // FIXME: sistemare valori rapporti   
    shiftRatioDataModel =-0.1;
    PlotThermusYields("lhc2760_final_0005_single_gc_output_gs1_wevc_nkst.txt" , kBlack     , kSolid      , "THERMUS 2.3 , 155 #pm 2 , 5924 #pm 543 , 23.6/9");
    shiftRatioDataModel =0.1;
    PlotGSIYields("data+therm_fit2_s2760_0-10qm14.dat"                        , kOrange-1  , kDashed     , "GSI      , 156 #pm 2  , 5330 #pm 505 , 17.4/9");


    myPadHisto->cd();
    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);
    DrawMarkerKStarNoFit() ;
    DrawExtrapolatedSymbolsYieldsPbPb0010();

    SaveCanvas("Fit_PbPb0010_Reference_GSITHERMUS");
    //    SaveCanvas("");

    array =0;
  }else if( what == "fitSHARE_NoPionsNoProtons_PbPb0010") {
    array =0;

    DrawRatio("frame",1);  
    DrawRatio("PbPb_0010",1);  
    legThermal->SetNColumns(4);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #chi^{2}/NDF", "0");

    particlesToExcludeFromChi2 = "[313]";
    shiftRatioDataModel=-0.2;
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_PbPb_0010_with_nuclei.txt"         ,kBlack    , kSolid      ,"SHARE 3          , 156 #pm 3 , 4476 #pm 696 , 14.1/9 ");
    particlesToExcludeFromChi2 = "[313][2212]";
    shiftRatioDataModel=0;
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_PbPb_0010_with_nuclei_protons.txt" ,kOrange-1 , kDashed     ,"SHARE 3 (no p)   , 156 #pm 3 , 4520 #pm 623  , 8.2/8");
    particlesToExcludeFromChi2 = "[313][211]";
    shiftRatioDataModel=0.2;
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_PbPb_0010_with_nuclei_pions.txt"   ,kBlue+1   , kDashDotted ,"SHARE 3 (no #pi) , 157 #pm 3 , 4103 #pm 677 , 12.2/8");

    myPadHisto->cd();
    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);

    // Add markers for additional particles not in fit
    DrawMarkerKStarNoFit() ;
    DrawExtrapolatedSymbolsYieldsPbPb0010();
    
    SaveCanvas("Fit_PbPb0010_SHARE_NoPiNoP");


  } else if( what == "fitGSI_NoPionsNoProtons_PbPb0010") {
    array =0;

    DrawRatio("frame",1);  
    DrawRatio("PbPb_0010",1);  
    legThermal->SetNColumns(4);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #chi^{2}/NDF", "0");

    particlesToExcludeFromChi2 = "[313]";
    shiftRatioDataModel=-0.2;
    PlotGSIYields("data+therm_fit2_s2760_0-10qm14.dat"     , kBlack    , kSolid      , "GSI          , 156 #pm 1.5 , 5330 #pm 505  , 18.1/9" );
    particlesToExcludeFromChi2 = "[313][2212]";
    shiftRatioDataModel=0.;
    PlotGSIYields("data+therm_fit2_s2760_0-10qm14NOp.dat"  , kOrange-1 , kDashed     , "GSI (no p)   , 156 #pm 2   , 5590 #pm 330  ,  7.7/8" );
    particlesToExcludeFromChi2 = "[313][211]";
    shiftRatioDataModel=0.2;
    PlotGSIYields("data+therm_fit2_s2760_0-10qm14NOpi.dat" ,  kBlue+1  , kDashDotted , "GSI (no #pi) , 157 #pm 2   ,  4990 #pm 630 , 16.5/8");

    myPadHisto->cd();
    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);

    // Add markers for additional particles not in fit
    DrawMarkerKStarNoFit() ;
    DrawExtrapolatedSymbolsYieldsPbPb0010();

    SaveCanvas("Fit_PbPb0010_GSI_NoPiNoP");
   

  } else if (what == "fitShareWithWithoutNuclei_PbPb0010") {
    DrawRatio("frame",1);  
    DrawRatio("PbPb_0010",1);  
    legThermal->SetNColumns(4);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #chi^{2}/NDF", "0");

    shiftRatioDataModel = -0.1;
    particlesToExcludeFromChi2="[313]"; // do not consider K* and nuclei 
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_PbPb_0010_with_nuclei.txt"  , kBlack  , kSolid      , "SHARE 3             , 156 #pm 3 , 4476 #pm 696 , 14.1/9");
    shiftRatioDataModel = 0.1;
    particlesToExcludeFromChi2="[313][1000010020][1000020030][1010010030]"; // do not consider K* and nuclei 
    PlotThermusYields("NEW_fit_gamma_q_fixed_PbPb_0010_without_nuclei.txt" , kCyan+2 , kDashDotted , "SHARE 3  (no nuclei), 156 #pm 4 , 4364 #pm 848 , 9.6/6");

    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);

    DrawMarkerKStarNoFit() ;
    DrawExtrapolatedSymbolsYieldsPbPb0010();

    SaveCanvas("Fit_PbPb0010_SHARE_WithWoNuclei");

    
  } else if( what == "fitThermus_GammaSFree_PbPb0010") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitShare_GammaSGammaQFree_PbPb0010") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitShare_All_PbPb0010") {
    DrawRatio("frame",1);  
    DrawRatio("PbPb_0010",1);  


    legThermal->SetNColumns(6);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #gamma_{s}, #gamma_{q}, #chi^{2}/NDF", "0");    

    particlesToExcludeFromChi2="[313][1000010020][1000020030][1010010030]"; // do not consider K* and nuclei 
    shiftRatioDataModel = -0.3;
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_PbPb_0010_without_nuclei.txt" , kBlack      , kSolid      , "SHARE 3               , 156 #pm 4  , 4364 #pm 848  , 1 (fix)        , 1 (fix)        , 12.4/6");
    shiftRatioDataModel = -0.1;
    PlotThermusYields("NEW_fit_gamma_q_fixed_PbPb_0010_without_nuclei.txt"   , kRed+1      , kDashDotted , "SHARE 3               , 155 #pm 3  , 4406 #pm 766  , 1.07 #pm 0.05  , 1 (fix)        , 9.6/5");
    shiftRatioDataModel = 0.1;
    PlotThermusYields("NEW_fit_gamma_q_s_free_PbPb_0010_without_nuclei.txt"  , kCyan+2     , kDashed     , "SHARE 3               ,  138 #pm 6 , 3064 #pm 1319 ,  1.98 #pm 0.68 ,  1.63 #pm 0.38 , 3.1/4");
    shiftRatioDataModel = 0.3;
    particlesToExcludeFromChi2="[313]"; // do not consider K*  
    PlotThermusYields("NEW_fit_gamma_q_s_free_PbPb_0010_with_nuclei.txt"     , kOrange - 1 , kDotted     , "SHARE 3 (with nuclei) , 152 #pm 8  , 4445 #pm 743  , 1.16 #pm 0.20  ,  1.06 #pm 0.12 ,  9.0/7");


    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);
    DrawMarkerKStarNoFit() ;
    DrawMarkerNucleiNoFit();
    DrawExtrapolatedSymbolsYieldsPbPb0010();


    legThermal->SetX1(0.12249  );
    legThermal->SetY1(0.0454769);
    legThermal->SetX2(0.821285 );
    legThermal->SetY2(0.383481 );

    SaveCanvas("Fit_PbPb0010_SHARE_All");
    array =0;
  } else if( what == "fitGSI_PbPb6080") {
    DrawRatio("frame",1);  
    DrawRatio("PbPb_6080",1);  
    legThermal->SetNColumns(4);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #chi^{2}/NDF", "0");
    particlesToExcludeFromChi2 = "[313][1000020030][1010010030]";
    PlotGSIYields    ("data+therm_fit2_s2760_60-80qm14.dat", kBlack, kSolid, "GSI, 157 #pm 2, 210 #pm 20, 8.2/7");

    legThermal->SetX1(0.143574  );
    legThermal->SetY1(0.0731318 );
    legThermal->SetX2(0.659639  );
    legThermal->SetY2(0.245206  );
    
    myPadHisto->cd();
    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);

    DrawMarkerKStarNoFit(1) ;
    
    //    DrawExtrapolatedSymbolsYieldsPbPb0010(0.143574, 0.251352, 0.351406, 0.343535,0); 
    //DrawExtrapolatedSymbolsYieldsPbPb0010();

    SaveCanvas("Fit_PbPb6080_GSI");
    array =0;
  }else if( what == "fitGSI_PbPb2040") {
    DrawRatio("frame",1);  
    DrawRatio("PbPb_2040",1);  
    legThermal->SetNColumns(4);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #chi^{2}/NDF", "0");
    particlesToExcludeFromChi2 = "[313][1000020030][1010010030]";
    PlotGSIYields    ("data+therm_fit2_s2760_20-40qm14.dat", kBlack, kSolid, "GSI, 161 #pm 3, 1725 #pm 220, 24.6/7");

    legThermal->SetX1(0.143574  );
    legThermal->SetY1(0.0731318 );
    legThermal->SetX2(0.659639  );
    legThermal->SetY2(0.245206  );
    
    myPadHisto->cd();
    NewLegendQM(0.651606, 0.765993, 0.909639, 0.865951, 1);

    DrawMarkerKStarNoFit(1) ;
    //    DrawMarkerKStarNoFit() ;
    
    //    DrawExtrapolatedSymbolsYieldsPbPb0010(0.143574, 0.251352, 0.351406, 0.343535,0); 
    //DrawExtrapolatedSymbolsYieldsPbPb0010();

    SaveCanvas("Fit_PbPb2040_GSI");
    array =0;
  } else if( what == "fitGSI_pPb0005") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitGSI_pPb2040") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitThermus_GammaSFree_pPb0005") {

    DrawRatio("frame",1);  
    DrawRatio("pPb0005",1, 0.00001);  
    particlesToExcludeFromChi2="[313][1000020030][1010010030] "; // do not consider K* 
    legThermal->SetNColumns(6);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), #gamma_{s}, r_{C} (fm), r (fm), #chi^{2}/NDF", "0");    

    shiftRatioDataModel=-0.2;
    PlotThermusYields("lhc5020_final_0005_single_gc_output_gs1_wevc_nkst.txt" , kBlack    , kSolid      , "THERMUS 2.3 GC , 158 #pm 2 , 1 (fix)       ,  N/A          , 3.40 #pm 0.11 , 32.0/7");
    shiftRatioDataModel =0;
    PlotThermusYields("lhc5020_final_0005_single_gc_output_gsf_wevc_nkst.txt" , kOrange-1 , kDashed     , "THERMUS 2.3 GC , 159 #pm 2 , 0.98 #pm 0.03 , N/A           , 3.40 #pm 0.11 , 31.3/6");
    shiftRatioDataModel = 0.2;
    PlotThermusYields("lhc5020_final_0005_single_sc_output_gs1_woevc_nkst.txt", kCyan+2   , kDashDotted , "THERMUS 2.3 SC , 158 #pm 3 , 1 (fix)       , 4.61 #pm 3.77 , 3.07 #pm 0.13  , 29.6/6");
   //lhc5020_final_0005_single_gc_output_gsf.txt
    //lhc5020_final_0005_single_sc_output_gs1.txt
    
    NewLegendQM(0.650602, 0.694971, 0.909639, 0.865951, 1);
    DrawExtrapNotInFitpPb0005() ;

    legThermal->SetX1(0.121486 ); 
    legThermal->SetY1(0.0741793); 
    legThermal->SetX2(0.759036 ); 
    legThermal->SetY2(0.384575 ); 
 
     SaveCanvas("Fit_pPb0005_THERMUS");

    //DrawExtrapolatedSymbolsAndLegendPbPb0010();

    array =0;
  } else if( what == "fitThermus_GammaSFree_pPb2040") {
    DrawRatio("frame",1);  
    DrawRatio("pPb2040",1, 0.00001);  
    particlesToExcludeFromChi2="[313][1000020030][1010010030] "; // do not consider K* 
    legThermal->SetNColumns(6);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), #gamma_{s}, r_{C} (fm), r (fm), #chi^{2}/NDF", "0");    

    shiftRatioDataModel=-0.2;
    PlotThermusYields("lhc5020_final_2040_single_gc_output_gs1_wevc_nkst.txt" , kBlack    , kSolid      , "THERMUS 2.3 GC , 155 #pm 2 , 1 (fix)       ,  N/A          , 2.83 #pm 0.08 , 40.1/7");
    shiftRatioDataModel =0;
    PlotThermusYields("lhc5020_final_2040_single_gc_output_gsf_wevc_nkst.txt" , kOrange-1 , kDashed     , "THERMUS 2.3 GC , 156 #pm 2 , 0.93 #pm 0.03 , N/A           , 2.83 #pm 0.08 , 34.6/6");
    shiftRatioDataModel = 0.2;
    PlotThermusYields("lhc5020_final_2040_single_sc_output_gs1_woevc_nkst.txt", kCyan+2   , kDashDotted , "THERMUS 2.3 SC , 156 #pm 2 , 1 (fix)       , 3.45 #pm 0.77 , 2.54 #pm 0.09 , 35.5/6");
   //lhc5020_final_0005_single_gc_output_gsf.txt
    //lhc5020_final_0005_single_sc_output_gs1.txt
       
    NewLegendQM(0.650602, 0.694971, 0.909639, 0.865951, 1);
    DrawExtrapNotInFitpPb0005(0) ;

    legThermal->SetX1(0.121486 ); 
    legThermal->SetY1(0.0741793); 
    legThermal->SetX2(0.759036 ); 
    legThermal->SetY2(0.384575 ); 



    SaveCanvas("Fit_pPb2040_THERMUS");

    array =0;
  } else if( what == "fitThermus_GammaSFree_pPb6080") {

    DrawRatio("frame",1);  
    DrawRatio("pPb6080",1, 0.00001);  
    particlesToExcludeFromChi2="[313][1000020030][1010010030][1000010020]"; // do not consider K* 
    legThermal->SetNColumns(6);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), #gamma_{s}, r_{C} (fm), r (fm), #chi^{2}/NDF", "0");    

    shiftRatioDataModel=-0.2;

    PlotThermusYields("lhc5020_final_6080_single_gc_output_gs1_wevc_nkst.txt" , kBlack    , kSolid      , "THERMUS 2.3 GC , 152 #pm 2 , 1 (fix)    , N/A , 2.18 #pm 0.06 , 48.2/7");
    shiftRatioDataModel =0;
    PlotThermusYields("lhc5020_final_6080_single_gc_output_gsf_wevc_nkst.txt" , kOrange-1 , kDashed     , "THERMUS 2.3 GC , 154 #pm 2 , 0.88 #pm 3 , N/A , 2.21 #pm 0.07 , 28.8/6");
    shiftRatioDataModel = 0.2;
    PlotThermusYields("lhc5020_final_6080_single_sc_output_gs1_woevc_nkst.txt", kCyan+2   , kDashDotted , "THERMUS 2.3 SC , 154 #pm 2 ,1 (fix) ,3.18 #pm 0.57 ,1.96 #pm 0.07 ,40.5/6");


   //lhc5020_final_0005_single_gc_output_gsf.txt
    //lhc5020_final_0005_single_sc_output_gs1.txt
       
    NewLegendQM(0.650602, 0.694971, 0.909639, 0.865951, 1);
    DrawExtrapNotInFitpPb0005(0) ;
 
    legThermal->SetX1(0.121486  ); 
    legThermal->SetY1(0.0268308 ); 
    legThermal->SetX2(0.758032  ); 
    legThermal->SetY2(0.337226  ); 


    SaveCanvas("Fit_pPb6080_THERMUS");



    array =0;
  }  else if( what == "fitThermus_RC_pPb0005") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitThermus_RC_pPb2040") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitThermus_RC_pPb6080") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitShare_pPb0005") {


    DrawRatio("frame",1);  
    DrawRatio("pPb0005",1, 0.00001);  
    legThermal->SetNColumns(6);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #gamma_{s}, #gamma_{q}, #chi^{2}/NDF", "0");    
    particlesToExcludeFromChi2="[313][1000010020][1000020030][1010010030]"; // do not consider K* and nuclei 
 
    shiftRatioDataModel = -0.2;
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_pPb_0005_without_nuclei.txt" , kBlack    , kSolid      , "SHARE 3 , 158 #pm 3 , 121 #pm 18  , 1 (fix)       , 1 (fix) , 24.3/6");
    shiftRatioDataModel = 0;
    PlotThermusYields("NEW_fit_gamma_q_fixed_pPb_0005_without_nuclei.txt"   , kBlue+1   , kDashDotted , "SHARE 3 , 161 #pm 3 , 115 #pm 17  , 0.93 #pm 0.04 , 1 (fix) , 20.3/5");
    shiftRatioDataModel = 0.2;
    PlotThermusYields("NEW_fit_gamma_q_s_free_pPb_0005_without_nuclei.txt"  , kOrange+2 , kDashDotted , "SHARE 3 , 144 #pm 1 , 81 #pm 25   ,  1.599 #pm 0 , 1.71 #pm 0.06 , 11.4/4");

    NewLegendQM(0.650602, 0.694971, 0.909639, 0.865951, 1);
    DrawExtrapNotInFitpPb0005() ;
    DrawMarkerNucleiNoFit();
    legThermal->SetX1(0.124498 ); 
    legThermal->SetY1(0.0715488); 
    legThermal->SetX2(0.672691 ); 
    legThermal->SetY2(0.384575 ); 

    SaveCanvas("Fit_pPb0005_SHARE");


    array =0;

  } else if( what == "fitShare_pPb0005_NoOmega") {


    DrawRatio("frame",1);  
    DrawRatio("pPb0005",1, 0.00001);  
    legThermal->SetNColumns(6);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #gamma_{s}, #gamma_{q}, #chi^{2}/NDF", "0");    
    particlesToExcludeFromChi2="[313][1000010020][1000020030][1010010030]"; // do not consider K* and nuclei 
 
    shiftRatioDataModel = -0.3;
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_pPb_0005_without_nuclei.txt"                 , kBlack    , kSolid      , "SHARE 3                 , 158 #pm 3  , 121 #pm 18  , 1 (fix)       , 1 (fix), 24.3/6");
    shiftRatioDataModel = -.1;                                                                                                                                                                       
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_pPb_0005_without_nuclei_excluding_Omega.txt" , kBlue+1   , kDashDotted , "SHARE 3 (No #Omega)     , 163 #pm 4  ,  99 #pm 17  , 1 (fix)       , 1 (fix), 15.6/5");
    shiftRatioDataModel = 0.1;                                                                                                                                                                       
    PlotThermusYields("NEW_fit_gamma_q_fixed_pPb_0005_without_nuclei_excluding_Omega.txt"   , kOrange+2 , kDashed     , "SHARE 3 (No #Omega)     , 163 #pm 3  , 100 #pm 16  , 0.96 #pm 0.04 , 1 (fix), 14.7/4");
    shiftRatioDataModel = 0.3;                                                                                                                                                                       
    PlotThermusYields("NEW_fit_gamma_q_s_fixed_pPb_0005_with_nuclei_excluding_Omega.txt"    , kCyan+2   , kDotted     , "SHARE 3 (No #Omega + d) , 160 #pm 3  , 114 #pm 17  , 1 (fix)       , 1 (fix), 18.6/6");

    NewLegendQM(0.650602, 0.694971, 0.909639, 0.865951, 1);
    DrawExtrapNotInFitpPb0005() ;
    DrawMarkerNucleiNoFit();
    legThermal->SetX1(0.124498 ); 
    legThermal->SetY1(0.0715488); 
    legThermal->SetX2(0.672691 ); 
    legThermal->SetY2(0.384575 ); 

    SaveCanvas("Fit_pPb0005_SHARE_NoOmega");


    array =0;

  }  else if( what == "fitShare_pPb2040") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitShare_pPb6080") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitGSI_pp") {
    DrawRatio("frame",1);
    DrawRatio("pp7",1, 0.00001);
    legThermal->SetNColumns(5);
    AddLineToThermalLegend(legThermal, "Model,T (MeV), V (fm^{3}), #gamma_{s}, #chi^{2}/NDF", "0");
    particlesToExcludeFromChi2 = "[313][1000020030][1010010030]";
    shiftRatioDataModel=-0.1;
    PlotGSIYields    ("data+therm_fit2_s7000.dat", kBlack, kSolid, "GSI, 146 #pm 2, 25 #pm 2, 1 (fix),  78.2/7", 0);
    shiftRatioDataModel=0.1;
    PlotGSIYields    ("data+therm_fit2_s7000gs.dat", kRed+1, kDashed, "GSI, 150 #pm 2, 23 #pm 2, 0.88, 45.6/7", 0);
    legThermal->SetX1(0.143574 );
    legThermal->SetY1(0.0731318);
    legThermal->SetX2(0.659639 );
    legThermal->SetY2(0.374263 );
    myPadHisto->cd();
    
    NewLegendQM(0.725904, 0.761431, 0.874498, 0.841323, 1);
  
    DrawMarkerKStarNoFit();
    TLegend * leg  = new TLegend(0.143574, 0.380408, 0.351406, 0.4603, NULL, "brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(14);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    TLegendEntry * entry=leg->AddEntry("TMarker","Not in fit","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(markerNoFit);
    entry->SetMarkerSize(1.2);
    entry->SetTextFont(43);
    leg->Draw();

    SaveCanvas("Fit_pp7000_GSI");
    //DrawExtrapNotInFitpPb0005(0);
    array =0;
  }else if( what == "fitFlorence_pp") {
    DrawRatio("frame",1);
    DrawRatio("pp7",1, 0.00001);
    legThermal->SetNColumns(1);
    AddLineToThermalLegend(legThermal, "Model", "0");
    shiftRatioDataModel = -0.1;
    particlesToExcludeFromChi2 = "[1000020030][1010010030][1000010020]";
    PlotFlorenceYields("pp7000_Florence.txt", kBlack, kSolid, "Becattini et al.");
    shiftRatioDataModel = 0.1;
    particlesToExcludeFromChi2 = "[313][1000020030][1010010030]";
    PlotGSIYields    ("data+therm_fit2_s7000gs.dat", kRed+1, kDashed, "GSI GC",0);
    myPadHisto->cd();
    NewLegendQM(0.725904, 0.761431, 0.874498, 0.841323, 1);

    SaveCanvas("Fit_pp7000_Florence");
    //DrawExtrapNotInFitpPb0005(0);
    array =0;
  } else if( what == "fitGSI_fullCanonical") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitThermus_RC_pp") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitShare_pp") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  }



  else {
    std::cout << "Unknown ratio " << what.Data() << std::endl;
    return;
  }

  if(!correlatedUnc) {
    std::cout << "correlatedUnc not set!" << std::endl;
    
  }
  std::cout << "CORR: " << correlatedUnc[1] << std::endl;

  if(array) {
    if(isYield) {
      TGraphErrors * hstat = 
      AliPWGHistoTools::GetGraphFromHisto(GetHistoYields(array,  system,  energy, centrality, label, color, marker, kStatError, shiftloc)
                                          ,0);
      hstat->SetMarkerSize(1.2);
      hstat->Draw("PZ");
      AliPWGHistoTools::GetGraphFromHisto(GetHistoYields(array,  system,  energy, centrality, label+"NoLegend", color, marker, kSystError, shiftloc)
                                          ,0)->Draw("[]");
    } else {
      AliPWGHistoTools::GetGraphFromHisto(GetHistoRatios(array,  system,  energy, centrality, label, color, marker, kStatError, shiftloc)
                                          ,0)->Draw("PZ");
      AliPWGHistoTools::GetGraphFromHisto(GetHistoRatios(array,  system,  energy, centrality, label+"NoLegend", color, marker, kSystError, shiftloc)
                                          ,0)->Draw("[]");
    }
  }
  

}

void DrawFrame(Bool_t isYield) {

  myCan = new TCanvas("myCan","ThermalFits",50,10,1000,isYield ? 950 : 650);
  myCan->Draw();
  myCan->cd();
  // Set the Pads
  Double_t boundaryLabels = isYield ? 0.88 : 0.85;//0.92;
  myPadLabel = new TPad("myPadLabel","myPadLabel",0.0, boundaryLabels,1.0,1.0);
  myPadSetUp(myPadLabel);
  myPadLabel->Draw();

  myPadHisto = new TPad("myPadHisto","myPadHisto",0.0,isYield ? 0.4 : 0.05 ,1.0, boundaryLabels,0);
  myPadSetUp(myPadHisto);
  myPadHisto->Draw();

  myPadStdDev = new TPad("myPadStdDev","myPadStdDev",0.0,0.035,1.0,0.215,0);
  myPadSetUp(myPadStdDev);
  if(isYield)  myPadStdDev->Draw();
  myPadStdDev->SetGridx();
  myPadStdDev->SetGridy();

  // This pad is for the ratios data/model
  myPadRatio = new TPad("myPadRatio","myPadRatio",0.0,0.22,1.0,0.4,0);
  myPadSetUp(myPadRatio);
  if(isYield)  myPadRatio->Draw();
  myPadRatio->SetGridx();
  myPadRatio->SetGridy();



  myPadLabel->cd();

  double xLabelPosition[nratio] = {0.124498, 0.211847, 0.31, 0.38, 0.465, 0.575, 0.644, 0.72, 0.82, 0.905 }; 
  double xLabelYields[npart]    = {0.115, 0.185, 0.270, 0.330, 0.422, 0.485, 0.570, 0.625, 0.695, 0.785, 0.835, 0.915};
  double yLabelPosition     = 0.40;

  // labels
  if(isYield) {
    for(Int_t ipart = 0; ipart < npart; ipart++){
      TLatex *myPart = new TLatex(xLabelYields[ipart],yLabelPosition, yieldsLabels[ipart]);
      myLatexDraw(myPart,20);    
    }    
  }  else {
    for(Int_t iratio = 0; iratio < nratio; iratio++){
      TLatex *myRatio = new TLatex(xLabelPosition[iratio],yLabelPosition, ratiosLabels[iratio]);
      myLatexDraw(myRatio,20);    
    }
  }
  // Xi's and Omega's bar (there was no way to convince root to draw it properly)
  if(isYield) {
    TLine *line = new TLine(0.653,0.660,0.663,0.660);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.725,0.660,0.738,0.660);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.362,0.660,0.374,0.660);
    line->SetLineWidth(2);
    line->Draw();


  }
  else {
    TLine *line = new TLine(0.408,0.68,0.418,0.68);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.498,0.68,0.513,0.68);
    line->SetLineWidth(2);
    line->Draw();
   

  }

  
  if(isYield) {
    myPadStdDev->cd();
    myPadStdDev->SetTopMargin(0.0);
    myPadStdDev->SetTicks(1,1);

    Float_t devMax = 4.7;
    
    TH2F *myBlankStdDev = new TH2F("myBlankStdDev","myBlankStdDev",npart,1,npart+1,10,-devMax,+devMax);
    myBlankStdDev->GetXaxis()->SetLabelFont(43); // precision 3: size will be in pixels
    myBlankStdDev->GetYaxis()->SetLabelFont(43);
    myBlankStdDev->GetYaxis()->SetTitleFont(43);
    myBlankStdDev->SetLabelSize(23,"xy");
    myBlankStdDev->SetTitleSize(20,"y");
    myBlankStdDev->SetNdivisions(20,"x");
    myBlankStdDev->SetNdivisions(605,"y");
    myBlankStdDev->SetLabelOffset(0.012,"xy");
    myBlankStdDev->SetYTitle("(mod.-data)/#sigma_{data}");
    myBlankStdDev->SetTitleOffset(1.65,"y");
    myBlankStdDev->Draw();

    TH2F *myBlankRatio = new TH2F("myBlankRatio","myBlankRatio",npart,1,npart+1,10,-0.8,0.8);
    myBlankRatio->GetXaxis()->SetLabelFont(43); // precision 3: size will be in pixels
    myBlankRatio->GetYaxis()->SetLabelFont(43);
    myBlankRatio->GetYaxis()->SetTitleFont(43);
    myBlankRatio->SetLabelSize(23,"xy");
    myBlankRatio->SetTitleSize(20,"y");
    myBlankRatio->SetNdivisions(20,"x");
    myBlankRatio->SetNdivisions(605,"y");
    myBlankRatio->SetLabelOffset(0.012,"xy");
    myBlankRatio->SetYTitle("(mod.-data)/mod.");
    myBlankRatio->SetTitleOffset(1.65,"y");
    myPadRatio->cd();
    //    myPadRatio->SetLogy();
    myBlankRatio->Draw();

  }

  myPadHisto->cd();
  myPadHisto->SetBottomMargin(isYield ? 0.01 : 0.04);
  //  myPadHisto->SetLogy();
  myPadHisto->SetTicks(1,1);
  TH1 *myBlankHisto = 0;
  if(isYield) {
    myBlankHisto =  new TH2F("NoLegend","NoLegend",npart,1,npart+1,10,  0.00002, maxy+0.01  );
    gPad->SetLogy();

  } else {
   myBlankHisto =  new TH2F("NoLegend","NoLegend",nratio,1,nratio+1,10,  0, maxy+0.01  );
  }
  myBlankHisto->GetXaxis()->SetLabelFont(43); // precision 3: size will be in pixels
  myBlankHisto->GetYaxis()->SetLabelFont(43);
  myBlankHisto->GetYaxis()->SetTitleFont(43);
  myBlankHisto->SetLabelSize(23,"y");
  myBlankHisto->SetLabelSize(0,"x");
  myBlankHisto->SetTitleSize(26,"y");
  myBlankHisto->SetMaximum(10);
  myBlankHisto->SetMinimum(0);
  myBlankHisto->SetNdivisions(isYield? 20 :10,"x");
  myBlankHisto->SetNdivisions(505,"y");
  if(isYield) myBlankHisto->SetYTitle("d#it{N}/d#it{y}");
  myBlankHisto->SetLabelOffset(0.012,"xy");
  myBlankHisto->SetTitleOffset(isYield ? 1.3 : 1,"y");
  myBlankHisto->Draw();

  if(isYield) {
    legThermal = new TLegend(0.144578, 0.0702247, 0.659639, 0.383226);
    legThermal->SetBorderSize(1);
    legThermal->SetTextFont(43);
    legThermal->SetTextSize(18);
    //    legThermal->SetNColumns(6);
    legThermal->SetFillColor(0);
    legThermal->SetLineWidth(1);
    legThermal->Draw();
  }
}

void myLatexDraw(TLatex *currentLatex, Float_t currentSize, Int_t currentColor){
  currentLatex->SetTextFont(43);
  currentLatex->SetTextSizePixels(Int_t(currentSize));
  //  currentLatex->SetTextAngle(0);
  currentLatex->Draw();
  return;
}

void myPaveSetup(float rRatio, float rRange3, float rRange5,
		 int rFillColor){

  float cHiRange = 0, cLoRange = 0;
  if(rRange3<rRange5) {cHiRange=rRange5;cLoRange=rRange3;}
  else {cHiRange=rRange3;cLoRange=rRange5;}

  TPave *cPave= new TPave(rRatio-0.25,cLoRange,rRatio+0.25,cHiRange,0,"br");
  cPave->SetFillColor(rFillColor);
  cPave->SetLineColor(1);
  cPave->Draw();
}

void myPadSetUp(TPad *currentPad){
  currentPad->SetLeftMargin(0.10);
  currentPad->SetRightMargin(0.02);
  currentPad->SetTopMargin(0.02);  
  currentPad->SetBottomMargin(0.02);
  return;
}
TGraphErrors*  PlotThermusYields(const char * filename, Int_t color, Int_t lineStyle,
                                 const char * tag) {

  Int_t lw = lineStyle == kSolid ? 2 : 3; // Set line width


  std::map<Int_t,Double_t> mapYields;
  std::map<Int_t,Double_t> mapStdDev;

  Int_t pdg;
  Double_t yield, stddev;
  ifstream thermusFile(filename);
  TString line;
  Double_t chi2 = 0;
  std::cout << "---"<<tag<<"---" << std::endl;

  //  std::istream is(thermusFile);
  // Read the std dev and the ratio in 2 maps, then plot them in a graph.
  while(line.ReadLine(thermusFile, kTRUE)) {
    if(line.BeginsWith("#")) continue;
    TObjArray * tokens = line.Tokenize(" \t");
    if(tokens->GetEntries() != 3) continue;// not a line with data
    //    thermusFile >> pdg >> yield >> stddev;
    pdg    = ((TObjString*)tokens->At(0))->String().Atof();
    yield  = ((TObjString*)tokens->At(1))->String().Atof();
    stddev = ((TObjString*)tokens->At(2))->String().Atof();

    if( thermusFile.eof() ) break;
    std::cout << "PDG " << pdg << " " << yield << " " << stddev << std::endl;
    
    mapYields[TMath::Abs(pdg)] += yield;
    mapStdDev[TMath::Abs(pdg)] += stddev;
    if(pdg < 0) {
      // found the antiparticle: now compute the mean
      mapYields[TMath::Abs(pdg)] /=2;
      mapStdDev[TMath::Abs(pdg)] /=2;
    }
    delete tokens;
  }

  // Now plot
  TGraphErrors * gThermus = new TGraphErrors;
  TGraphErrors * gThermusStdDev = new TGraphErrors;
  for(Int_t ipart = 0; ipart < npart; ipart++){
    gThermus->SetPoint(ipart, ipart+1.5, mapYields[particleYields[ipart]]);
    gThermus->SetPointError(ipart, 0.3, 0);
    gThermusStdDev->SetPoint(ipart, ipart+1.5, mapStdDev[particleYields[ipart]]);
    gThermusStdDev->SetPointError(ipart, 0.3, 0);
  }

  myPadHisto->cd();
  gThermus->Draw("PZ");
  gThermus->SetLineWidth(lw);

  gThermus->SetLineColor(color);
  gThermus->SetLineStyle(lineStyle);
  gThermus->SetTitle("NoLegend");


  //  myPadStdDev->cd();
  // gThermusStdDev->Draw("PZ");
  // gThermusStdDev->SetLineWidth(lw);
  // gThermusStdDev->SetLineColor(color);    
  // gThermusStdDev->SetLineStyle(lineStyle);
  TGraphErrors* gStdDev2 = 0;
  TGraphErrors* gRatio   = 0;
  std::cout << "CHI2: " 
	    <<  GetGraphRatioAndStdDev(gThermus, gRatio, gStdDev2)
	    << std::endl;
  myPadRatio->cd(); 
  gRatio->Draw("PZ");
  myPadStdDev->cd();
  gStdDev2->Draw("PZ");
  myPadHisto->cd();
  AddLineToThermalLegend(gThermus, tag, "l");
  return gThermus;
}

TGraphErrors*  PlotFlorenceYields(const char * filename, Int_t color, Int_t lineStyle,
                                 const char * tag) {

  Int_t lw = lineStyle == kSolid ? 2 : 3; // Set line width


  std::map<Int_t,Double_t> mapYields;
  std::map<Int_t,Double_t> mapStdDev;

  Int_t pdg;
  Double_t yield, yieldbar, stddev;
  ifstream thermusFile(filename);
  TString line;
  std::cout << "---"<<tag<<"---" << std::endl;

  //  std::istream is(thermusFile);
  // Read the std dev and the ratio in 2 maps, then plot them in a graph.
  while(line.ReadLine(thermusFile, kTRUE)) {
    if(line.BeginsWith("#")) continue;
    TObjArray * tokens = line.Tokenize(" \t");
    if(tokens->GetEntries() != 4) continue;// not a line with data
    //    thermusFile >> pdg >> yield >> stddev;
    pdg      = ((TObjString*)tokens->At(0))->String().Atof();
    yield    = ((TObjString*)tokens->At(2))->String().Atof();
    yieldbar = ((TObjString*)tokens->At(3))->String().Atof(); // Antiparticle yield
    if(pdg == 0) {
      // not a line with data
      delete tokens;
      continue ;
    }
    if( thermusFile.eof() ) break;
    std::cout << "PDG " << pdg << " " << yield << " " << yieldbar << std::endl;
    
    mapYields[TMath::Abs(pdg)] += yieldbar ? ((yield+yieldbar)/2) : yield; // If the antiparticle exists, use the average

    delete tokens;
  }

  // Now plot
  TGraphErrors * gFlorence = new TGraphErrors;
  for(Int_t ipart = 0; ipart < npart; ipart++){
    gFlorence->SetPoint(ipart, ipart+1.5, mapYields[particleYields[ipart]]);
    gFlorence->SetPointError(ipart, 0.3, 0);
  }

  myPadHisto->cd();
  gFlorence->Draw("PZ");
  gFlorence->SetLineWidth(lw);

  gFlorence->SetLineColor(color);
  gFlorence->SetLineStyle(lineStyle);
  gFlorence->SetTitle("NoLegend");



  TGraphErrors* gStdDev2 = 0;
  TGraphErrors* gRatio   = 0;
  std::cout << "CHI2: " 
	    <<  GetGraphRatioAndStdDev(gFlorence, gRatio, gStdDev2)
	    << std::endl;
  myPadRatio->cd(); 
  gRatio->Draw("PZ");
  myPadStdDev->cd();
  gStdDev2->Draw("PZ");
  myPadHisto->cd();
  AddLineToThermalLegend(gFlorence, tag, "l");
  return gFlorence;
}


TGraphErrors*  PlotGSIYields(const char * filename, Int_t color, Int_t lineStyle,
                             const char * tag, Bool_t isPbPb) {

  // tag is a comma separated list of elements to be added to the legend as diferent columns

  Int_t lw = lineStyle == kSolid ? 2 : 3; // Set line width

 
  const Int_t pdgPbPb0010[] = {211, -211, 321, -321, 310, 313, 333, 2212, -2212, 3122, 3312, -3312, 3334, -3334, 1000010020, 1000020030, 1010010030, -1010010030};
  const Int_t pdgpp[] = {211, 321, 310, 313, 333, 2212, 3122, 3312, 3334, 1000010020};

  const Int_t *pdgOrder = isPbPb ? pdgPbPb0010 : pdgpp;

  std::map<Int_t,Double_t> mapYields;
  std::map<Int_t,Double_t> mapStdDev;
  std::map<Int_t,Double_t> mapUncert;
  std::map<Int_t,Double_t> mapData;

  Double_t data, uncert, model;
  
  ifstream gsiFile(filename);
  //  std::istream is(thermusFile);
  std::cout << "----- "<<tag<<" -----" << std::endl;
  
  // Read the std dev and the ratio in 2 maps, then plot them in a graph.
  Int_t ipart = 0;
  while(gsiFile) {
    gsiFile >> data >> uncert >> model;
    if( gsiFile.eof() ) break;
    Int_t pdg = pdgOrder[ipart];
    std::cout << "PDG " << pdg << " " << data << std::endl;
    mapYields[TMath::Abs(pdg)] += model;
    mapUncert[TMath::Abs(pdg)] += uncert;
    mapData[TMath::Abs(pdg)]   += data;
    if(pdg < 0) {
      // found the antiparticle: now compute the mean
      mapYields[TMath::Abs(pdg)] /=2;
      mapData[TMath::Abs(pdg)] /=2;
      mapUncert[TMath::Abs(pdg)] /=2;
      
    }
    ipart++;
  }

  // Now plot
  TGraphErrors * gGsi = new TGraphErrors;
  TGraphErrors * gGsiStdDev = new TGraphErrors;
  std::cout << "PDG    \tmodel\tdata\tuncert\tstddev" << std::endl;  // header
  for(Int_t ipart = 0; ipart < npart; ipart++){ // Here we use npart, as this is what we wnat to plot!
    Int_t pdg = particleYields[ipart];
    mapStdDev[TMath::Abs(pdg)]  = ( mapYields[TMath::Abs(pdg)] - mapData[TMath::Abs(pdg)]) /  mapUncert[TMath::Abs(pdg)] ;
    std::cout << "PDG " << pdg <<"\t" 
	      << mapYields[TMath::Abs(pdg)] << "\t" << mapData[TMath::Abs(pdg)] <<"\t" 
	      << mapUncert[TMath::Abs(pdg)] << "\t" << mapStdDev[TMath::Abs(pdg)]  
	      << std::endl;

    if(!mapYields[particleYields[ipart]]) mapYields[particleYields[ipart]] = -10;
    if(!mapStdDev[particleYields[ipart]]) mapStdDev[particleYields[ipart]] = -10;
    

    gGsi->SetPoint(ipart, ipart+1.5, mapYields[particleYields[ipart]]);
    gGsi->SetPointError(ipart, 0.3, 0);
    gGsiStdDev->SetPoint(ipart, ipart+1.5, mapStdDev[particleYields[ipart]]);
    gGsiStdDev->SetPointError(ipart, 0.3, 0);
  }

  

  myPadHisto->cd();
  gGsi->Draw("PZ");
  gGsi->SetLineWidth(lw);

  gGsi->SetLineColor(color);    
  gGsi->SetLineStyle(lineStyle);
  gGsi->SetTitle("NoLegend");


  // myPadStdDev->cd();
  // gGsiStdDev->Draw("PZ");
  // gGsiStdDev->SetLineWidth(lw);
  // gGsiStdDev->SetLineColor(color);    
  // gGsiStdDev->SetLineStyle(lineStyle);
  // myPadHisto->cd();

  TGraphErrors* gStdDev2 = 0;
  TGraphErrors* gRatio   = 0;
  std::cout << "CHI2: " 
	    <<  GetGraphRatioAndStdDev(gGsi, gRatio, gStdDev2)
	    << std::endl;
  myPadRatio->cd(); 
  gRatio->Draw("PZ");
  myPadStdDev->cd();
  gStdDev2->Draw("PZ");

  AddLineToThermalLegend(gGsi, tag, "L");

  return gGsi;
}


void DrawExtrapolatedSymbolsAndLegendPbPb0010() {
    myPadHisto->cd();


    TLegend *leg = new TLegend(    0.149598, 0.782203, 0.415663, 0.858447,NULL,"pNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(18);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(2);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetMargin(0.1);

    TLegendEntry *entry=leg->AddEntry("TMarker","Extrapolated (Pb-Pb 0-10%)","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(27);
    entry->SetMarkerSize(1.2);
    entry->SetTextFont(43);
    leg->Draw();
    myPadLabel->cd();
    // Markers for extrapolated points

    TMarker *marker = new TMarker(0.666667,0.111825,27);
    marker->SetMarkerStyle(27);
    marker->SetMarkerSize(1.2);
    marker->Draw();
    marker = new TMarker(0.920683,0.111825,27);
    marker->SetMarkerStyle(27);
    marker->SetMarkerSize(1.2);
    marker->Draw();

    // BR for 3He

    myPadLabel->cd();
    TLatex *   tex = new TLatex(0.73,0.05,"BR = 25%");
    tex->SetNDC();
    tex->SetTextFont(43);
    tex->SetTextSize(14);
    tex->SetLineWidth(2);
    tex->Draw();



}
void DrawExtrapolatedSymbolsAndLegendpPb0005() {
    myPadHisto->cd();


    TLegend *leg = new TLegend(    0.149598, 0.709972, 0.415663, 0.786216,NULL,"pNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(18);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(2);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->SetMargin(0.1);
    TLegendEntry *entry=leg->AddEntry("TMarker","Extrapolated (p-Pb 0-5%)","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(28);
    entry->SetMarkerSize(1.2);
    entry->SetTextFont(43);
    leg->Draw();
    myPadLabel->cd();

   TMarker *marker = new TMarker(0.590361,0.111825,28);
   marker->SetMarkerStyle(28);
   marker->SetMarkerSize(1.2);
   marker->Draw();
   marker = new TMarker(0.938755,0.111825,28);
   marker->SetMarkerStyle(28);
   marker->SetMarkerSize(1.2);
   marker->Draw();
}

void AddLineToThermalLegend(TObject * obj, TString line, const char * optFirst) {

  // This should be a comma-separated list of text to be added to the
  // columns.  If the number of entries does not match the numer of
  // columns, it returns an error
  TObjArray * tokens = line.Tokenize(",");
  if(tokens->GetEntries() != legThermal->GetNColumns()) {
    std::cout << "Wrong number of columns (" << tokens->GetEntries() << ","<<legThermal->GetNColumns()<<") for the thermal legend, not adding " << line.Data() << std::endl;
    return;
  }
  
  TIter iter(tokens);
  TObjString * col;
  Bool_t first = 1;
  while((col = (TObjString*)iter.Next())) {
    // Add entry removing whitespaces
    legThermal->AddEntry(obj, col->String().Strip(TString::kBoth, ' ').Data(), first ? optFirst : "0");
    if (first) first = 0;
  }
}

void DrawExtrapolatedSymbolsYieldsPbPb0010(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Bool_t plotExtraploatedLegend){
    // Markers for extrapolated points
  myPadLabel->cd();
  TMarker * marker = new TMarker(0.36245,0.111825,markerExtrap);
  marker->SetMarkerStyle(markerExtrap);
  marker->SetMarkerSize(1.2);
  marker->Draw();
  marker = new TMarker(0.945783,0.111825,markerExtrap);
  marker->SetMarkerStyle(markerExtrap);
  marker->SetMarkerSize(1.2);
  marker->Draw();
  // BR for 3He
  myPadHisto->cd();
  TLatex *   tex = new TLatex(11.15, 2.5e-5,"BR = 25%");
  tex->SetTextFont(43);
  tex->SetTextSize(14);
  tex->SetLineWidth(2);
  tex->Draw();

  myPadHisto->cd();



    TLegend * leg = new TLegend(x1,y1,x2,y2,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
  leg->SetTextSize(14);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  TLegendEntry * entry=leg->AddEntry("TMarker","Not in fit","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(markerNoFit);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(43);
  if(plotExtraploatedLegend) {
    entry=leg->AddEntry("TMarker","Extrapolated","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(markerExtrap);
    entry->SetMarkerSize(1.2);
    entry->SetTextFont(43);
  }
  leg->Draw();
  leg->SetMargin(0.1);

  
}

void DrawMarkerKStarNoFit(Bool_t plotLegend) {
  myPadLabel->cd();
  TMarker *marker = new TMarker(0.344378,0.111825,markerNoFit);
  marker->SetMarkerStyle(markerNoFit);
  marker->SetMarkerSize(1.2);
  marker->Draw();
  myPadHisto->cd();
  
  if(plotLegend) {
    

    TLegend * leg = new TLegend(0.126506, 0.253051, 0.335341, 0.345118,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(14);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    TLegendEntry * entry=leg->AddEntry("TMarker","Not in fit","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(markerNoFit);
    entry->SetMarkerSize(1.2);
    entry->SetTextFont(43);
    leg->Draw();
  }
}

void DrawMarkerNucleiNoFit() {
  myPadLabel->cd();
  TMarker *marker = new TMarker(0.928715,0.111825,markerNoFit);
  marker->SetMarkerStyle(markerNoFit);
  marker->SetMarkerSize(1.2);
  marker->Draw();
  
  marker = new TMarker(0.791751,0.111825,markerNoFit);
  marker->SetMarkerStyle(markerNoFit);
  marker->SetMarkerSize(1.2);
  marker->Draw();

  marker = new TMarker(0.866466,0.111825,markerNoFit);
  marker->SetMarkerStyle(markerNoFit);
  marker->SetMarkerSize(1.2);
  marker->Draw();
  myPadHisto->cd();

}

Double_t GetGraphRatioAndStdDev(TGraphErrors * gModel, TGraphErrors * &gRatio, TGraphErrors *&gStdDev) {
  
  // Plots the ratios data/model

  // I recomputed the stddev and the chi2 here, because some values
  // changed slightly due to typos in some of the input files with
  // respect to the ones used for those fits. 


  if(!gModel) {
    std::cout << "EMPTY MODEL" << std::endl;
    return 0;
  }
  // 0. set up the graphs which we will need
  TGraphErrors * gTemp = 0;
  TGraphErrors * gStat = 0;
  TGraphErrors * gSyst = 0;
  // the cloning below takes care of the style
  gRatio = (TGraphErrors*)gModel->Clone();
  gRatio->Clear();
  gStdDev = (TGraphErrors*)gModel->Clone();
  gStdDev->Clear();

  // 1. Find the data graphs. We need both the stat and the syst, since we want the ratio with the total uncertainty (or not?)
  TVirtualPad * currentPad = gPad;
  myPadHisto->cd();
  TList * list = myPadHisto->GetListOfPrimitives(); 
  TIterator * iter = list->MakeIterator();

  while ((gTemp = (TGraphErrors*) iter->Next())){
    if(gTemp->InheritsFrom("TGraphErrors")) {
      // Found a graph, it is the correct one?
      TString name = gTemp->GetName();
      TString title = gTemp->GetTitle();
      std::cout << "name " << name.Data() << std::endl;
      
      if (name.Contains("hPart")) {
	// ok, it's the data
	if (title.Contains("NoLegend")) gSyst = gTemp; // it's the syst error FIXME: it would be better to add the error type to the data
	else gStat = gTemp;
      }
      if(gStat && gSyst) {
	std::cout << "Cannot find gStat or gSyst (" << gStat << "," << gSyst <<")" << std::endl;
	break; // found both stat and syst uncertainties
      }
    }
  }

  

  // Compute the ratio, the stddev and the chi2
  Int_t npoint = gModel->GetN(); // We are sure that data and model have the same number of points in the same order, because they were created using the particleYields array (FIXME: is this also true for GSI?)
  Double_t chi2 = 0;
  for(Int_t ipoint = 0; ipoint < npoint; ipoint++){

    Double_t yield = gStat->GetY()[ipoint];
    Double_t stat  = gStat->GetEY()[ipoint];
    Double_t syst  = gSyst->GetEY()[ipoint];
    Double_t error = TMath::Sqrt(stat*stat+syst*syst);    
    Double_t model = gModel->GetY()[ipoint];
    Double_t width = gModel->GetEX()[ipoint];
    Double_t ratio  = (model-yield)/model;
    Double_t stddev = (model-yield)/error;
    if(!particlesToExcludeFromChi2.Contains(Form("[%d]", particleYields[ipoint]))) chi2 += TMath::Power(stddev,2);
    else std::cout << "Ecluding PDG "<< particleYields[ipoint] <<" from chi2 calculation" << std::endl;
    Double_t errorRatio = error/model;
    gRatio->SetPoint(ipoint, gModel->GetX()[ipoint]+shiftRatioDataModel, ratio);
    gRatio->SetPointError(ipoint, 0, errorRatio);
    if(model) {
      gStdDev->SetPoint(ipoint, gModel->GetX()[ipoint], stddev);
      gStdDev->SetPointError(ipoint, width, 0);
    }
    // The commented block down here was used to compare to the estimate provided directly by Boris and Benjamin
    // gStdDev->SetPoint(ipoint, gModel->GetX()[ipoint]-0.4, stddev);
    // gStdDev->SetPointError(ipoint, 0.2, 0);
    // std::cout << "PDG " << particleYields[ipoint] <<"\t" 
    // 	      << model << "\t" << yield <<"\t" 
    // 	      << error << "\t" << stddev  
    // 	      << std::endl;
  }
  
  gRatio->SetLineStyle(kSolid);
  gRatio->SetMarkerStyle(kOpenSquare);
  gRatio->SetMarkerColor(gRatio->GetLineColor());
  gRatio->SetLineWidth(2);

  currentPad->cd();
  return chi2;
}


void SaveCanvas(const char * name) {

  if(!saveCanvas) return;

  std::cout << "Saving " << name << ".{eps,pdf,root,C}" << std::endl;
  
  myCan->Update();
  gSystem->ProcessEvents();
  myCan->Print(Form("%s.eps",name));
  myCan->Print(Form("%s.root",name));
  myCan->Print(Form("%s.C",name));
  gSystem->Exec(Form("epstopdf %s.eps", name));
  gSystem->Exec(Form("if [ \"$USER\" = \"mfloris\" ]; then cp %s.{eps,pdf,root,C} /Users/mfloris/Documents/PapersNTalks/ALICE/ThermalFits/img/; fi ",name));

}

void DrawExtrapNotInFitpPb0005(Bool_t drawExtrap) {
  myPadLabel->cd();
  TMarker *marker = new TMarker(0.344378,0.111825,28);
  marker->SetMarkerStyle(28);
  marker->SetMarkerSize(1.2);
  marker->Draw();
  // marker = new TMarker(0.7851406,0.111825,28);
  // marker->SetMarkerStyle(28);
  // marker->SetMarkerSize(1.2);
  // marker->Draw();
  if(drawExtrap) {
    marker = new TMarker(0.364458,0.111825,27);
    marker->SetMarkerStyle(27);
    marker->SetMarkerSize(1.2);
    marker->Draw();
    marker = new TMarker(0.792,0.111825,27);
    marker->SetMarkerStyle(27);
    marker->SetMarkerSize(1.2);
    marker->Draw();
  }
  myPadHisto->cd();
  
  

  TLegend * leg = 0;
  
  if (drawExtrap) leg = new TLegend(0.123494, 0.400358, 0.331325, 0.534512,NULL,"brNDC");
  else leg = new TLegend(0.123494, 0.395097, 0.332329, 0.474011, NULL, "brNDC");
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
  leg->SetTextSize(14);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  TLegendEntry * entry=leg->AddEntry("TMarker","Not in fit","p");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(markerNoFit);
  entry->SetMarkerSize(1.2);
  entry->SetTextFont(43);

  if(drawExtrap) {
    entry=leg->AddEntry("TMarker","Extrapolated","p");
    entry->SetLineColor(1);
    entry->SetLineStyle(1);
    entry->SetLineWidth(1);
    entry->SetMarkerColor(1);
    entry->SetMarkerStyle(markerExtrap);
    entry->SetMarkerSize(1.2);
    entry->SetTextFont(43);
  }
  leg->Draw();
  leg->SetMargin(0.1);
}
