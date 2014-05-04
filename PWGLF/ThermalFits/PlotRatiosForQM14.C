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

void LoadArrays() ;
//void AddLabel(Float_t x, Float_t y, TString text);
void myLatexDraw(TLatex *currentLatex, Float_t currentSize=0.5, Int_t currentColor=1);
void myPaveSetup(float rRatio=0, float rRange3=0, float rRange5=0,
		 int rFillColor=0);
void myPadSetUp(TPad *currentPad);
TGraphErrors*  PlotThermusYields(const char * filename, Int_t color, Int_t lineStyle,
                                 const char * tag,
                                 Double_t t = -1, Double_t terr=-1, Double_t v=-1, Double_t verr=-1, 
                                 Double_t gs=-1, Double_t gser=-1, Double_t gq=-1, Double_t gerr=-1,
                                 Double_t chi2=-1, Double_t ndf=-1) ;

TGraphErrors * PlotGSIYields(const char * fileName, Int_t color=kBlack, Int_t lineStyle = kSolid,
                             const char * tag = "GSI",
                             Double_t t = -1, Double_t terr=-1, Double_t v=-1, Double_t verr=-1, 
                             Double_t gs=-1, Double_t gser=-1, Double_t gq=-1, Double_t gerr=-1,
                             Double_t chi2=-1, Double_t ndf=-1) ;

// Ratios to be draw. Remember to change the labels in DrawFrame if you change this
const Int_t nratio = 10;
//  Int_t denum[nratio]    = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGPi       , kPDGDeuteron , kPDGPi          , kPDGK   , -kPDGK};

Int_t num  [nratio]            = {kPDGK  , kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron , kPDGHE3      , kPDGHyperTriton , kPDGPhi , kPDGKStar};
Int_t denum[nratio]            = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGProton   , kPDGDeuteron , kPDGPi          , kPDGK   , kPDGK};
Int_t isSum[nratio]            = {1      , 1          ,  1         ,   1     , 1         , 1            , 0            , 1               , 1       , 1      };
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
static const Double_t scale[]  = {1      , 3          ,  0.5       ,  30     ,  250      , 50           , 100          , 4e5             , 2       , 1      };
//static const Double_t scale[]  = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,};
const Int_t npart = 12;
Int_t particleYields  [npart] = {kPDGPi ,kPDGK   ,kPDGKS0, kPDGKStar, kPDGPhi, kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron, kPDGHyperTriton, kPDGHE3    };
Int_t isSumYields[npart]      = {1      ,1       ,0      , 1        , 0      , 1           ,1            ,1        ,1          ,0           , 1              , 0          };
//Int_t isSumInputFiles[npart]  = {1      ,1       ,0      , 1        , 0      , 1           ,1            ,1        ,1          ,0           , 1              , 0          };
const char * yieldsLabels[]          = {"#frac{#pi^{+}+#pi^{-}}{2}",
                                        "#frac{K^{+}+K^{-}}{2}",
                                        "K_{S}^{0}",
                                        "#frac{K*+#bar{K}*}{2}",
                                        "#phi",
                                        "#frac{p+#bar{p}}{2}",
                                        "#Lambda", 
                                        "#frac{#Xi^{-}+#Xi^{+}}{2}",
                                        "#frac{#Omega^{-}+#Omega^{+}}{2}",
                                        "d",
                                        "#frac{{}^{3}_{#Lambda}H+{}^{3}_{#Lambda}#bar{H}}{2}",
                                        "{}^{3}He",
};

// Preferred colors and markers
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar,0};
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
Double_t correlatedUncLocalPbPb[14] = {0.0424 , 0.0424     ,  0.041     ,  0      , 0         , 0.0424       , 0.0424       , 0               , 0.05    , 0.05   };
Double_t correlatedUncLocalPP  [14] = {0.0424 , 0.0424     ,  0.041     ,  0      , 0         , 0.0424       , 0.0424       , 0               , 0.0424    , 0.0424   };
Double_t correlatedUncZero[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

TCanvas *myCan = 0;
TPad    *myPadStdDev =0;
TPad    *myPadHisto  =0;
TPad    *myPadLabel  =0;
TLegend * legThermal = 0;

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
  PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A0005", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A2040", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A6080", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M6080", /* separateCharges*/1);
  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M2030", /*separateCharges*/1);

  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/0);
  PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/0);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A0005", /*separateCharges*/0);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A2040", /*separateCharges*/0);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A6080", /*separateCharges*/0);
  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M6080", /*separateCharges*/0);  
  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M2030", /*separateCharges*/0);

  return 0;
#endif

  SetStyle();

  //DrawRatio("allpp");  
  //  DrawRatio("PbPbWithPP7TeV");
  // DrawRatio("allsyst");
  //  DrawRatio("PbPb6080andpPb0005");
  //DrawRatio("pp_vsRHIC");
  DrawRatio("PbPb_vsRHIC");
  //  DrawRatio("aliceall");


  // Yields and FITS
  // maxy=20000;
  // DrawRatio("fit_ReferenceFit_PbPb0010", 1);
  //  DrawRatio("fitShare_pPb0005", 1);
  //DrawRatio("fitShare_All_PbPb0010", 1);


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

  //  Double_t isSum = -1; // if this is -1, then the sum criterion is ignored
  for(Int_t ipart = 1; ipart <= npart; ipart++){
    AliParticleYield * part = AliParticleYield::FindParticle(arr, particleYields[ipart-1], system, energy, centrality,isSumYields[ipart-1]);
    if(!part && isSumYields[ipart-1]) {
      //Could not find the particle, but the sum was requested: build the sum!
      part = AliParticleYield::FindParticle(arr, particleYields[ipart-1], system, energy, centrality,  0);
      AliParticleYield * part2 = AliParticleYield::FindParticle(arr, -particleYields[ipart-1], system, energy, centrality,  0);
        if(part2 && part) part = AliParticleYield::Add(part, part2);        
        else if(part) part->Scale(2.); // If we only found a particle, we can scale it by a factor 2.
        else part = 0;
    }    
    if(!part){
      std::cout << "Cannot find " << particleYields[ipart-1] << std::endl;
      continue;
    }
    if(isSumYields[ipart-1]) part->Scale(0.5); // take average
    
    std::cout << " Plotting " << particleYields[ipart-1] << std::endl;

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
  h->SetMarkerSize(1.2);
  
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
  gStyle->SetLineWidth(2);
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
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

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
	  l->AddEntry(obj, obj->GetTitle(), style);
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
    l->SetTextSize(16);
    TMarker * m = new TMarker(0.1, 0.1,markerNoFit);
    m->SetMarkerSize(1.2);
    l->AddEntry( m, "Not in fit", "p");
    m = new TMarker(0.1,0.1,markerExtrap);
    m->SetMarkerSize(1.2);
    l->AddEntry( m, "Extrapolated (Pb-Pb 0-10%)", "p");
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
    
    TPaveText *pt = new TPaveText(    0.169679, 0.842881, 0.378514, 0.929595,"brNDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(0);
    pt->SetLineWidth(1);
    pt->SetTextFont(43);
    pt->SetTextSize(23);
    pt->AddText("ALICE Preliminary");
    pt->Draw();
    // pt = new TPaveText(    0.167671, 0.815206, 0.375502, 0.865021,"brNDC");
    // pt->SetBorderSize(0);
    // pt->SetFillColor(0);
    // pt->SetLineWidth(1);
    // pt->SetTextFont(43);
    // pt->SetTextSize(20);
    // pt->AddText("particle+antiparticle");
    // pt->Draw();

    gPad->SetGridx();

  }
  else if (what == "PbPb_0010") {
    array = arrPbPb;
    system = 2; energy = 2760.; centrality = "V0M0010";
    label = "Pb-Pb, #sqrt{s}_{NN} = 2.76 TeV, 0-10%";
    color  = kRed+1;
    marker = kFullCircle;
    if(!shift)    shiftloc =  0;
    correlatedUnc = correlatedUncLocalPbPb;

  }

  else if (what == "PbPb_6080") {
    array = arrPbPb;
    system = 2; energy = 2760.; centrality = "V0M6080";
    label = "Pb-Pb, #sqrt{s}_{NN} = 2.76 TeV, 60-80%";
    color  = kBlue+1;
    marker = kFullCircle;
    if(!shift)    shiftloc =  0.0;
    correlatedUnc = correlatedUncLocalPbPb;
  }
  else if (what == "pp7") {
    array = arrpp7;
    system = 0; energy = 7000.; centrality = "";
    label = "pp #sqrt{s} = 7 TeV";
    color  = kMagenta+1;
    marker = kFullCircle;
    if(!shift)    shiftloc =  0.2;
    correlatedUnc = correlatedUncLocalPP;
  }
  else if (what == "pp900") {
    array = arrpp900;
    system = 0; energy = 900.; centrality = "";
    label = "pp #sqrt{s} = 0.9 TeV";
    color  = kCyan+2;
    marker = kFullCircle;
    if(!shift)    shiftloc =  -0.2;
    correlatedUnc = correlatedUncLocalPP;  
  }
  else if (what == "pp276") {
    array = arrpp276;
    system = 0; energy = 2760.; centrality = "";
    label = "pp #sqrt{s} = 2.76 TeV";
    color  = kYellow+2;
    marker = kFullCircle;
    if(!shift)    shiftloc = 0;
    correlatedUnc = correlatedUncLocalPP;
  }
  else if (what == "pPb0005") {
    array = arrpPb;
    system = 1; energy = 5020.; centrality = "V0A0005";
    label = "p-Pb, #sqrt{s}_{NN} = 5.02 TeV, V0A 0-5%";
    color  = kBlack;
    marker = kFullCircle;
    if(!shift)    shiftloc = -0.2;
    correlatedUnc = correlatedUncLocalPP;
  } 
  else if (what == "PbPbSTAR") {
    array = arrSTARPbPb;
    system = 2; energy = 200.; centrality = "0005";
    label = "STAR, Au-Au, #sqrt{s}_{NN} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenStar;
    if(!shift)    shiftloc = +0.2;
    correlatedUnc = correlatedUncZero;
  }
  else if (what == "PbPbPHENIX") {
    array = arrPHENIXPbPb;
    system = 2; energy = 200.; centrality = "0005";
    label = "PHENIX, Au-Au, #sqrt{s}_{NN} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenSquare;
    if(!shift)    shiftloc = -0.15;
    correlatedUnc = correlatedUncZero;
  }
  else if (what == "PbPbBRAHMS") {
    array = arrBRAHMSPbPb;
    system = 2; energy = 200.; centrality = "0010";
    label = "BRAHMS, Au-Au, #sqrt{s}_{NN} = 0.2 TeV, 0-10%";
    color  = kBlack;
    marker = kOpenCross;
    if(!shift)    shiftloc = -0.3;
    correlatedUnc = correlatedUncZero;
  }
  else if (what == "ppSTAR") {
    array = arrSTARpp;
    system = 0; energy = 200.; centrality = "";
    label = "STAR, pp, #sqrt{s} = 0.2 TeV";
    color  = kBlack;
    marker = kOpenStar;
    if(!shift)    shiftloc = 0.;
    correlatedUnc = correlatedUncZero;
  }
  else if (what == "ppPHENIX") {
    array = arrPHENIXpp;
    system = 0; energy = 200.; centrality = "";
    label = "PHENIX, pp, #sqrt{s} = 0.2 TeV";
    color  = kBlack;
    marker = kOpenSquare;
    if(!shift)    shiftloc = -0.2;
    correlatedUnc = correlatedUncZero;
  } 
  // From here on, it's meta names, to draw several series of ratios
  else if (what == "allpp"){
    DrawRatio("frame");
    DrawRatio("pp7");
    DrawRatio("pp276");
    DrawRatio("pp900");
    array =0;
    NewLegendQM(0.62249, 0.635734, 0.910643, 0.94673);

    myCan->Update();
    gSystem->ProcessEvents();
    myCan->Print("Ratios_pponly.eps");
    gSystem->Exec("epstopdf Ratios_pponly.eps");
    gSystem->Exec("if [ \"$USER\" = \"mfloris\" ]; then cp Ratios_pponly.{eps,pdf} /Users/mfloris/Documents/PapersNTalks/ALICE/ThermalFits/img/; fi ");
  }
  else if (what == "PbPbWithPP7TeV"){
    
    DrawRatio("frame");
    DrawRatio("PbPb_0010");
    DrawRatio("pp7");
    array =0;

    NewLegendQM(0.413655, 0.748094, 0.910643, 0.948736);

    myCan->Update();
    gSystem->ProcessEvents();
    myCan->Print("Ratios_withpp7tev.eps");
    gSystem->Exec("epstopdf Ratios_withpp7tev.eps");
    gSystem->Exec("if [ \"$USER\" = \"mfloris\" ]; then cp Ratios_withpp7tev.{eps,pdf} /Users/mfloris/Documents/PapersNTalks/ALICE/ThermalFits/img/; fi ");
  } else if(what == "allsyst") {

    DrawRatio("frame");
    DrawRatio("PbPb_0010", 0, 0.2);
    DrawRatio("pp7", 0, -0.2);
    DrawRatio("pPb0005", 0, 0.00001);
    array =0;


    NewLegendQM(0.462851, 0.631722, 0.89257, 0.936697);
    DrawExtrapolatedSymbolsAndLegendPbPb0010();
    DrawExtrapolatedSymbolsAndLegendpPb0005();
    myCan->Update();
    gSystem->ProcessEvents();
    myCan->Print("Ratios_allsystems.eps");
    gSystem->Exec("epstopdf Ratios_allsystems.eps");
    gSystem->Exec("if [ \"$USER\" = \"mfloris\" ]; then cp Ratios_allsystems.{eps,pdf} /Users/mfloris/Documents/PapersNTalks/ALICE/ThermalFits/img/; fi ");

  }else if(what =="PbPb6080andpPb0005") {
    DrawRatio("frame");
    DrawRatio("PbPb_6080");
    DrawRatio("pPb0005");
    array=0;


    NewLegendQM(0.413655, 0.72803, 0.910643, 0.948736);
    DrawExtrapolatedSymbolsAndLegendpPb0005();
    myCan->Update();
    gSystem->ProcessEvents();
    myCan->Print("Ratios_6080vspPb.eps");
    gSystem->Exec("epstopdf Ratios_6080vspPb.eps");
    gSystem->Exec("if [ \"$USER\" = \"mfloris\" ]; then cp Ratios_6080vspPb.{eps,pdf} /Users/mfloris/Documents/PapersNTalks/ALICE/ThermalFits/img/; fi ");
    
    
  }else if(what =="pp_vsRHIC") {
    DrawRatio("frame");
    DrawRatio("pp7");
    DrawRatio("ppSTAR");
    DrawRatio("ppPHENIX");
    array=0;

    NewLegendQM(    0.554217, 0.677869, 0.910643, 0.948736);

    myCan->Update();
    gSystem->ProcessEvents();
    myCan->Print("Ratios_vsRHIC_pp.eps");
    gSystem->Exec("epstopdf Ratios_vsRHIC_pp.eps");
    gSystem->Exec("if [ \"$USER\" = \"mfloris\" ]; then cp Ratios_vsRHIC_pp.{eps,pdf} /Users/mfloris/Documents/PapersNTalks/ALICE/ThermalFits/img/; fi ");
        
  } else if (what =="PbPb_vsRHIC") {
    DrawRatio("frame");
    DrawRatio("PbPb_0010");
    DrawRatio("PbPbSTAR");
    DrawRatio("PbPbPHENIX");
    DrawRatio("PbPbBRAHMS");
    array = 0;



    NewLegendQM(    0.434739, 0.591593, 0.939759, 0.936697);
    DrawExtrapolatedSymbolsAndLegendPbPb0010();

    myCan->Update();
    gSystem->ProcessEvents();
    myCan->Print("Ratios_vsRHIC_PbPb.eps");
    gSystem->Exec("epstopdf Ratios_vsRHIC_PbPb.eps");
    gSystem->Exec("if [ \"$USER\" = \"mfloris\" ]; then cp Ratios_vsRHIC_PbPb.{eps,pdf} /Users/mfloris/Documents/PapersNTalks/ALICE/ThermalFits/img/; fi ");
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
    // FIXME: sistemare valori rapporti
    PlotThermusYields("test_outputfile.txt"               , kBlack , kSolid      , "Thermus 2.3" , 155  ,2  ,4220  ,500  ,-1  ,0 ,-1 ,0 ,39.5  ,10);
    PlotGSIYields("data+therm_s2760_0-10qm14.dat"         , kRed+1 , kDashed     , "GSI"         , 156.5  ,2  , 5380, 560  ,-1  ,0 ,-1 ,0 ,38.4  ,14 );
    PlotThermusYields("fit_gamma_q_s_fixed_PbPb_0010.txt" , kBlue+1, kDashDotted , "SHARE 3"     , 156 , 3 , 4387 , 767 , -1 ,0 ,-2 ,0 ,17.2 ,10);

    NewLegendQM(0.613454, 0.701578, 0.940763, 0.918272, 1);


    array =0;
  } else if( what == "fitGSI_NoPions_PbPb0010") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
    
  }  else if( what == "fitGSI_NoProtons_PbPb0010") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
    
  } else if( what == "fitThermus_GammaSFree_PbPb0010") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitShare_GammaSGammaQFree_PbPb0010") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitShare_All_PbPb0010") {
    DrawRatio("frame",1);  
    DrawRatio("PbPb_0010",1);  
    PlotThermusYields("fit_gamma_q_s_fixed_PbPb_0010.txt" , kBlack , kSolid , "SHARE 3" , 100 , 100 , 100  , 100 , -1   , -1  , -2 , -1 , -1  , -1);
    PlotThermusYields("fit_gamma_q_fixed_PbPb_0010.txt"   , kBlue  , kSolid , "SHARE 3" , 100 , 100 , 100  , 100 , -1   , -1  , -2 , -1 , -1  , -1);
    PlotThermusYields("fit_gamma_q_s_free_PbPb_0010.txt"  , kRed   , kSolid , "SHARE 3" , 100 , 100 , 100  , 100 , -1   , -1  , -2 , -1 , -1  , -1);
    array =0;
  } else if( what == "fitGSI_PbPb6080") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitGSI_pPb0005") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitGSI_pPb2040") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitThermus_GammaSFree_pPb0005") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitThermus_GammaSFree_pPb2040") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitThermus_GammaSFree_pPb6080") {
    std::cout << "MISSING DATA" << std::endl;
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
    DrawRatio("pPb0005",1);  
    PlotThermusYields("fit_gamma_q_fixed_pPb_0005.txt"    , kBlue  , kDashDotted , "SHARE 3"     , 161 , 3 , 114 , 17 , 0.93 ,0.04 ,-2 ,0 ,26.0 ,6);
    PlotThermusYields("fit_gamma_q_s_free_pPb_0005.txt"   , kRed   , kDashDotted , "SHARE 3"     , 162 , 3 , 111 , 16 , 0.93 ,0.03 ,-2 ,0 ,61.4 ,6);
    PlotThermusYields("fit_gamma_q_s_fixed_pPb_0005.txt"   , kGreen   , kDashDotted , "SHARE 3"    , 158 , 2 , 125 , 16 , -1 ,-1 ,-2 ,0 ,30.8 ,7);
    PlotThermusYields("fit_gamma_q_s_free_with_d_pPb_0005.txt", kBlack, kSolid, "SHARE 3", 163, 5, 124, 15, 0.85, 0.08, 0.92, 0.06, 29.3, 7);
    NewLegendQM(0.613454, 0.701578, 0.940763, 0.918272, 1);
    array =0;

  } else if( what == "fitShare_pPb2040") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitShare_pPb6080") {
    std::cout << "MISSING DATA" << std::endl;
    array =0;
  } else if( what == "fitGSI_pp") {
    std::cout << "MISSING DATA" << std::endl;
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
      AliPWGHistoTools::GetGraphFromHisto(GetHistoYields(array,  system,  energy, centrality, label, color, marker, kStatError, shiftloc)
                                          ,0)->Draw("PZ");
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

  myCan = new TCanvas("myCan","BD april.2014",50,10,1000,650);
  myCan->Draw();
  myCan->cd();
  // Set the Pads
  Double_t boundaryLabels = 0.85;//0.92;
  myPadLabel = new TPad("myPadLabel","myPadLabel",0.0, boundaryLabels,1.0,1.0);
  myPadSetUp(myPadLabel);
  myPadLabel->Draw();

  myPadHisto = new TPad("myPadHisto","myPadHisto",0.0,isYield ? 0.25 : 0.05 ,1.0, boundaryLabels,0);
  myPadSetUp(myPadHisto);
  myPadHisto->Draw();

  myPadStdDev = new TPad("myPadStdDev","myPadStdDev",0.0,0.0,1.0,0.25,0);
  myPadSetUp(myPadStdDev);
  if(isYield)  myPadStdDev->Draw();
  myPadStdDev->SetGridx();
  myPadStdDev->SetGridy();
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
    TLine *line = new TLine(0.653,0.75,0.663,0.75);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.728,0.75,0.741,0.75);
    line->SetLineWidth(2);
    line->Draw();

    // Markers for extrapolated points
    TMarker *marker = new TMarker(0.339357,0.111825,markerNoFit);
    marker->SetMarkerStyle(28);
    marker->SetMarkerSize(1.2);
    marker->Draw();
    marker = new TMarker(0.369478,0.111825,markerExtrap);
    marker->SetMarkerStyle(27);
    marker->SetMarkerSize(1.2);
    marker->Draw();
    marker = new TMarker(0.938755,0.111825,markerExtrap);
    marker->SetMarkerStyle(markerExtrap);
    marker->SetMarkerSize(1.2);
    marker->Draw();

  }
  else {
    TLine *line = new TLine(0.408,0.75,0.418,0.75);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(0.498,0.75,0.513,0.75);
    line->SetLineWidth(2);
    line->Draw();
   

  }

  
  if(isYield) {
    myPadStdDev->cd();
    myPadStdDev->SetTopMargin(0.0);
    myPadStdDev->SetTicks(1,1);

    Float_t devMax = 3.5;
    
    TH2F *myBlankStdDev = new TH2F("myBlankStdDev","myBlankStdDev",npart,1,npart+1,10,-devMax,+devMax);
    myBlankStdDev->GetXaxis()->SetLabelFont(43); // precision 3: size will be in pixels
    myBlankStdDev->GetYaxis()->SetLabelFont(43);
    myBlankStdDev->GetYaxis()->SetTitleFont(43);
    myBlankStdDev->SetLabelSize(23,"xy");
    myBlankStdDev->SetTitleSize(20,"y");
    myBlankStdDev->SetNdivisions(20,"x");
    myBlankStdDev->SetNdivisions(605,"y");
    myBlankStdDev->SetLabelOffset(0.012,"xy");
    myBlankStdDev->SetYTitle("(model-data)/#sigma_{data}");
    myBlankStdDev->SetTitleOffset(1.3,"y");
    myBlankStdDev->Draw();
  }

  myPadHisto->cd();
  myPadHisto->SetBottomMargin(0.01);
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
  myBlankHisto->SetLabelSize(23,"xy");
  myBlankHisto->SetTitleSize(26,"y");
  myBlankHisto->SetMaximum(10);
  myBlankHisto->SetMinimum(0);
  myBlankHisto->SetNdivisions(isYield? 20 :10,"x");
  myBlankHisto->SetNdivisions(505,"y");
  if(isYield) myBlankHisto->SetYTitle("d#it{N}/d#it{y}");
  myBlankHisto->SetLabelOffset(0.012,"xy");
  myBlankHisto->SetTitleOffset(1,"y");
  myBlankHisto->Draw();

  if(isYield) {
    legThermal= new TLegend(0.144578, 0.0702247, 0.659639, 0.383226);
    legThermal->SetBorderSize(1);
    legThermal->SetTextFont(43);
    legThermal->SetTextSize(14);
    legThermal->SetNColumns(6);
    legThermal->SetFillColor(0);
    legThermal->SetLineWidth(1);
    legThermal->Draw();
    // FIXME: make a method "Set header"?
    legThermal->AddEntry(myBlankHisto, "Model", "0");
    legThermal->AddEntry(myBlankHisto, "T (MeV)", "0");
    legThermal->AddEntry(myBlankHisto, "V (fm^{3})", "0");
    legThermal->AddEntry(myBlankHisto, "#gamma_{s}", "0");
    legThermal->AddEntry(myBlankHisto, "#gamma_{q}", "0");
    legThermal->AddEntry(myBlankHisto, "#chi^{2}/NDF", "0");


  } else {
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
                                 const char * tag,
                                 Double_t t, Double_t terr, Double_t v, Double_t verr, 
                                 Double_t gs, Double_t gserr, Double_t gq, Double_t gqerr,
                                 Double_t chi2, Double_t ndf) {

  Int_t lw = lineStyle == kSolid ? 2 : 3; // Set line width


  std::map<Int_t,Double_t> mapYields;
  std::map<Int_t,Double_t> mapStdDev;

  Int_t pdg;
  Double_t yield, stddev;
  ifstream thermusFile(filename);
  TString line;
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


  myPadStdDev->cd();
  gThermusStdDev->Draw("PZ");
  gThermusStdDev->SetLineWidth(lw);
  gThermusStdDev->SetLineColor(color);    
  gThermusStdDev->SetLineStyle(lineStyle);
  myPadHisto->cd();

  legThermal->AddEntry(gThermus, tag, "l");
  if (t>0) {
    legThermal->AddEntry(gThermus, Form("%0.0f #pm %0.0f",t,terr), "0");
  } else {
    legThermal->AddEntry(gThermus, " ", "0");
  }
  if(v>0) {
    legThermal->AddEntry(gThermus, Form("%0.0f #pm %0.0f",v,verr), "0");
  } else {
    legThermal->AddEntry(gThermus, " ", "0");
  }

  if(gs>0) {
    legThermal->AddEntry(gThermus, Form("%0.2f #pm %0.2f",gs,gserr), "0");
  } else {
    legThermal->AddEntry(gThermus, "fixed", "0");
  }
  if(gq>0) {
    legThermal->AddEntry(gThermus, Form("%0.2f #pm %0.2f",gq,gqerr), "0");
  } else if(gq==-1){ // -1 is thermus, -2 share
    legThermal->AddEntry(gThermus, "N/A", "0");
  } else {
    legThermal->AddEntry(gThermus, "fixed", "0");
  }
  if(chi2>0) {
    legThermal->AddEntry(gThermus, Form("%0.2f/%0.0f",chi2, ndf), "0");
  } else {
    legThermal->AddEntry(gThermus, " ", "0");
  }

  return gThermus;
}



TGraphErrors*  PlotGSIYields(const char * filename, Int_t color, Int_t lineStyle,
                             const char * tag,
                             Double_t t, Double_t terr, Double_t v, Double_t verr, 
                             Double_t gs, Double_t gserr, Double_t gq, Double_t gqerr,
                             Double_t chi2, Double_t ndf) {

  Int_t lw = lineStyle == kSolid ? 2 : 3; // Set line width

  const Int_t pdgPbPb0010[] = {211, -211, 321, -321, 310, 313, 333, 2212, -2212, 3122, 3312, -3312, 3334, -3334, 1000010020, 1000020030, 1010010030, -1010010030};
  const Int_t pdgPbPb6080[] = {211 , -211 , 321 , -321 , 310 , 313 , 333 , 2212 , -2212 , 3122 , 3312 , -3312 , 3334 , -3334 ,1000010020};
  const Int_t pdgpPb0005[]  = {211, 321, 310, 313, 333, 2212, 3122, 3312, 3334, 1000010020};

  std::map<Int_t,Double_t> mapYields;
  std::map<Int_t,Double_t> mapStdDev;
  std::map<Int_t,Double_t> mapUncert;
  std::map<Int_t,Double_t> mapData;

  Double_t data, uncert, model;
  
  ifstream gsiFile(filename);
  //  std::istream is(thermusFile);
  std::cout << "GSI" << std::endl;
  
  // Read the std dev and the ratio in 2 maps, then plot them in a graph.
  Int_t ipart = 0;
  while(gsiFile) {
    gsiFile >> data >> uncert >> model;
    if( gsiFile.eof() ) break;
    Int_t pdg = pdgPbPb0010[ipart];
    std::cout << "PDG " << pdg << " " << data << std::endl;;
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
    //    chi2 += 
    std::cout << "PDG " << pdg <<"\t" 
	      << mapYields[TMath::Abs(pdg)] << "\t" << mapData[TMath::Abs(pdg)] <<"\t" 
	      << mapUncert[TMath::Abs(pdg)] << "\t" << mapStdDev[TMath::Abs(pdg)]  
	      << std::endl;

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


  myPadStdDev->cd();
  gGsiStdDev->Draw("PZ");
  gGsiStdDev->SetLineWidth(lw);
  gGsiStdDev->SetLineColor(color);    
  gGsiStdDev->SetLineStyle(lineStyle);
  myPadHisto->cd();

  legThermal->AddEntry(gGsi, tag, "l");
  if (t>0) {
    legThermal->AddEntry(gGsi, Form("%0.0f #pm %0.0f",t,terr), "0");
  } else {
    legThermal->AddEntry(gGsi, " ", "0");
  }
  if(v>0) {
    legThermal->AddEntry(gGsi, Form("%0.0f #pm %0.0f",v,verr), "0");
  } else {
    legThermal->AddEntry(gGsi, " ", "0");
  }

  if(gs>0) {
    legThermal->AddEntry(gGsi, Form("%0.2f #pm %0.2f",gs,gserr), "0");
  } else {
    legThermal->AddEntry(gGsi, "fixed", "0");
  }
  if(gq>0) {
    legThermal->AddEntry(gGsi, Form("%0.2f #pm %0.2f",gq,gqerr), "0");
  } else {
    legThermal->AddEntry(gGsi, "N/A", "0");
  }
  if(chi2>0) {
    legThermal->AddEntry(gGsi, Form("%0.2f/%0.0f",chi2, ndf), "0");
  } else {
    legThermal->AddEntry(gGsi, "--", "0");
  }

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
    marker = new TMarker(0.920683,0.0904227,27);
    marker->SetMarkerStyle(27);
    marker->SetMarkerSize(1.2);
    marker->Draw();

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

   TMarker *marker = new TMarker(0.590361,0.0797218,28);
   marker->SetMarkerStyle(28);
   marker->SetMarkerSize(1.2);
   marker->Draw();
   marker = new TMarker(0.938755,0.0797218,28);
   marker->SetMarkerStyle(28);
   marker->SetMarkerSize(1.2);
   marker->Draw();
}
