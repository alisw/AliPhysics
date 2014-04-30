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



#endif

// Plots ratios for QM and saves input files for thermal models


enum MyParticles { kPDGPi = 211, kPDGK = 321, kPDGProton = 2212, kPDGKS0 = 310, kPDGLambda=3122, kPDGXi=3312,kPDGOmega=3334,kPDGPhi=333,kPDGKStar=313,kPDGDeuteron=1000010020,kPDGHE3 = 1000020030, kPDGHyperTriton = 1010010030, kPDGSigmaStarPlus=3224,kPDGSigmaStarMinus=3114,kPDGXiStar=3324};

typedef enum {kStatError, kSystError, kTotalError} myerror_t;

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle, Int_t icolor, Int_t imarker = kOpenSquare, Int_t errorsType = kTotalError, Float_t shift = 0) ;
TH1F * GetHistoYields(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle) ;
void   PrepareThermalModelsInputFiles(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, Bool_t separateCharges=0) ;
void SetStyle(Bool_t graypalette=0) ;
void NewLegendQM(Double_t x1, Double_t y1, Double_t x2, Double_t y2) ;
void DrawRatio(TString what);
void DrawYield(TString what);

void DrawFrame(Bool_t yields = 0) ;
void LoadArrays() ;
//void AddLabel(Float_t x, Float_t y, TString text);
void myLatexDraw(TLatex *currentLatex, Float_t currentSize=0.5, Int_t currentColor=1);
void myPaveSetup(float rRatio=0, float rRange3=0, float rRange5=0,
		 int rFillColor=0);
void myPadSetUp(TPad *currentPad);

// Ratios to be draw. Remember to change the labels in DrawFrame if you change this
const Int_t nratio = 10;
//  Int_t denum[nratio]    = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGPi       , kPDGDeuteron , kPDGPi          , kPDGK   , -kPDGK};

Int_t num  [nratio]            = {kPDGK  , kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron , kPDGHE3      , kPDGHyperTriton , kPDGPhi , kPDGKStar};
Int_t denum[nratio]            = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGProton   , kPDGDeuteron , kPDGPi          , kPDGK   , kPDGK};
Int_t isSum[nratio]            = {1      , 1          ,  1         ,   1     , 1         , 1            , 0            , 1               , 1       , 1      };
// const char * ratiosLabels[]          = {"K/^{}#pi", 
//                                         "p/^{}#pi",   
//                                         "#Lambda /^{}K_{s}^{0}", 
//                                         "#Xi/^{}#pi",
//                                         "#Omega/^{}#pi",
//                                         "d/^{}p",
//                                         "^{3}He^{}/d",
//                                         "{}^{3}_{#Lambda}H^{}/#pi",
//                                         "#phi /^{}K",
//                                         "K*/^{}K",};
const char * ratiosLabels[]          = {"#frac{K^{+}+K^{-}}{#pi^{+}+#pi^{-}}", 
                                        "#frac{p+#bar{p}}{#pi^{+}+#pi^{-}}", 
                                        "#frac{2#Lambda}{K_{s}^{0}}", 
                                        "#frac{#Xi^{-}+#Xi^{+}}{#pi^{+}+#pi^{-}}",
                                        "#frac{#Omega^{-}+#Omega^{+}}{#pi^{+}+#pi^{-}}",
                                        "#frac{d}{p+#bar{p}}",
                                        "#frac{{}^{3}He }{d}",
                                        "#frac{{}^{3}_{#Lambda}H+{}^{3}_{#Lambda}#bar{H} }{#pi^{+}+#pi^{-}}",
                                        "#frac{#phi}{K^{+}+K^{-}}",
                                        "#frac{K*+#bar{K}*}{K^{+}+K^{-}}",};
static const Double_t scale[]  = {1      , 3          ,  0.5       ,  30     ,  250      , 50           , 100          , 4e5             , 2       , 1      };
//static const Double_t scale[]  = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,};


// Preferred colors and markers
// const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7, kWhite}; // for syst bands
// const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2  , kWhite};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar,0};
//

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


TClonesArray * PlotRatiosForQM14() {
#if !(!defined (__CINT__) || (defined(__MAKECINT__)))
  LoadLibs();
#endif

  //
  LoadArrays();

  // Uncomment stuff in this section to save the inputs for thermal models
  //#define SAVE_INPUT_THERMAL_MODEL
#ifdef SAVE_INPUT_THERMAL_MODEL
  //  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/1);
  //PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A0005", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A2040", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A6080", /*separateCharges*/1);
  //  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M6080", /*separateCharges*/1);
  //PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M2030", /*separateCharges*/1);

  // PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/0);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A0005", /*separateCharges*/0);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A2040", /*separateCharges*/0);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A6080", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M6080", /*separateCharges*/0);  
  // PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M2030", /*separateCharges*/0);

  return 0;
#endif

  SetStyle();

  // TCanvas * c1 = new TCanvas("Ratios", "Ratios", 1400, 600);
  // c1->SetMargin( 0.0744986, 0.0329513, 0.225131, 0.83);
  //  c1->SetLogy();

  // CENTRAL
  //  DrawRatio("allpp");  
  //  DrawRatio("PbPbWithPP7TeV");
  //  DrawRatio("allsyst");
  //  DrawRatio("PbPb6080andpPb0005");
  //  DrawRatio("pp_vsRHIC");
  //  DrawRatio("PbPb_vsRHIC");
  //  DrawRatio("aliceall");
  return 0;
  // 
  TCanvas * c2 = new TCanvas("Yields", "Yields", 1400, 600);
  c2->SetMargin( 0.0744986, 0.0329513, 0.225131, 0.0593368);

  GetHistoYields(arrPbPb,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV, 0-10%")->Draw();
  GetHistoYields(arrThermus,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Thermus")->Draw("same");
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
  const Int_t npart = 12;
  Int_t particles  [npart] = {kPDGPi ,kPDGK   ,kPDGKS0, kPDGKStar, kPDGPhi, kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron, kPDGHyperTriton, kPDGHE3    };
  Int_t isSum[npart]       = {1      ,1       ,0      , 1        , 0      , 1           ,1            ,1        ,1          ,0           , 1              , 0          };

  Int_t ipartOut = 0; // Index for the array
  Int_t ipartOutGSI = 0; // Index for the array
  
  for(Int_t ipart = 0; ipart < npart; ipart++){
    if(!separateCharges) {
      AliParticleYield * part = AliParticleYield::FindParticle(arr, particles[ipart], system, energy, centrality,  isSum[ipart]);
      if(!part && isSum[ipart]) {
        //Could not find the particle, but the sum was requested: build the sum!
        part = AliParticleYield::FindParticle(arr, particles[ipart], system, energy, centrality,  0);
        AliParticleYield * part2 = AliParticleYield::FindParticle(arr, -particles[ipart], system, energy, centrality,  0);
        if(part2 && part) part = AliParticleYield::Add(part, part2);        
        else if(part) part->Scale(2.); // If we only found a particle, we can scale it by a factor 2.
        else part = 0;
      }
      // We want to save the average of particle and antiparticle in this case
      if(part) {
        if(isSum[ipart] && !part->IsTypeAverage()) part->Scale(0.5); // If it's not already an average, but just a sum, divide by 2
        new((*arrOut   )[ipartOut++]) AliParticleYield(*part);
        new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(*part);
      } else { // Add dummy particle to the GSI list
        new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(particles[ipart], system, energy, -10, -10, -10, -10, -10, -10, 5, 256, "DUMMY", 1, "ALICE");        
      }
    }
    else {
      // ignore isSum and try to find both particles
      Bool_t notFound = 0;
      AliParticleYield * part = AliParticleYield::FindParticle(arr, particles[ipart], system, energy, centrality,  0);
      if(part) new((*arrOut)[ipartOut++]) AliParticleYield(*part);
      else notFound=1;
      // Try to find antiparticle (-pdg code)
      part = AliParticleYield::FindParticle(arr, -particles[ipart], system, energy, centrality,  0);
      if(part) {
        new((*arrOut)[ipartOut++]) AliParticleYield(*part);
        new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(*part);
      }
      else if (notFound) {
        // If neither charge was found, check if we at least have the sum 
        part = AliParticleYield::FindParticle(arr, abs(particles[ipart]), system, energy, centrality,  1);
        if (part) {
          part->Scale(0.5);
          new((*arrOut)[ipartOut++]) AliParticleYield(*part);
          new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(*part);
        }
        else {
          new((*arrOutGSI)[ipartOutGSI++]) AliParticleYield(particles[ipart], system, energy, -10, -10, -10, -10, -10, -10, 5, 256, "DUMMY", 1, "ALICE");        
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


TH1F * GetHistoYields(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle) {

  const Int_t npart = 11;
  Int_t pdg  [npart]    = {kPDGPi, kPDGK  , kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron , kPDGHE3      , kPDGHyperTriton , kPDGPhi , kPDGKStar};
  Int_t isSum[npart]    = {1      ,1      ,1           ,0           ,1        ,1          ,0             ,0             ,1                ,0        ,1      };
  //Int_t isSum[npart]    = {0,0,0,0,0,0,0,0,0,0,1};
  //Double_t scale[npart] = {1      ,1      ,3           ,1           ,30       ,250         ,50            ,10            ,4e5              ,2       ,2      };
  Double_t scale[npart] = {1,5,30,30,200,1000,4000,2e6,2e6,20,20,};
  TH1F * h = new TH1F(Form("hPart_%d_%0.0f_%s",system,energy,centrality.Data()), histotitle, npart, 1, npart+1);

  //  Double_t isSum = -1; // if this is -1, then the sum criterion is ignored
  for(Int_t ipart = 1; ipart <= npart; ipart++){
    AliParticleYield * part = AliParticleYield::FindParticle(arr, pdg[ipart-1], system, energy, centrality,isSum[ipart-1]);
    if(!part) continue;
    part->Scale(scale[ipart-1]);
    h->SetBinContent(ipart, part->GetYield());
    h->SetBinError  (ipart, part->GetTotalError(0/* 0 = no normalization error */));
    h->GetXaxis()->SetBinLabel(ipart, Form("#splitline{%s}{%s}",Form("#times%2.0g",  scale[ipart-1]), part->GetLatexName()));

  }
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

void NewLegendQM(Double_t x1, Double_t y1, Double_t x2, Double_t y2) {

  const char * style = "lp";
  const char ** labels=0;
  //  Bool_t beautify=kFALSE;
  Bool_t useTitle=kTRUE;

  TLegend * l = new TLegend(x1, y1, x2, y2);
  l->SetFillColor(kWhite);
  l->SetTextFont(43);
  l->SetTextSize(25);
  l->SetBorderSize(1);
  l->SetLineWidth(1);
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

  l->Draw();

}


void DrawYield(TString what) {
  // This is used to simplify the code above
  // In order to draw syst error bars, we need to convert to graphs the syst errors histos
  if (what == "frame" ) {
    correlatedUnc = correlatedUncZero;
    DrawFrame(1);
  }

}

void DrawRatio(TString what) {
  // This is used to simplify the code above
  // In order to draw syst error bars, we need to convert to graphs the syst errors histos

  // Sample colors
  //  const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2  , kWhite};

  TClonesArray * array = 0;
  Int_t system,  color, marker;
  Float_t energy = 0, shift = 0;
  TString centrality, label;
  // FIXME: move this in the different sections below
  correlatedUnc = 0;
  std::cout << "Plotting " << what.Data() << std::endl;
  

  if (what == "frame" ) {
    correlatedUnc = correlatedUncZero;
    DrawFrame();
    // TH1 * h = GetHistoRatios(arrPbPb,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "NoLegend", kWhite);
    // h->Draw();
    // h->GetYaxis()->SetDecimals(1);
    // h->GetYaxis()->SetNdivisions(505);
    //    h->GetXaxis()->CenterLabels(1);
    //    Int_t nratio = h->GetNbinsX();
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
        Double_t shift = scale[iratio] < 50 ? 0.3 : 0.2;
        TLatex * scaleLabel = new TLatex(iratio+1+shift,  0.005, Form("#times %g", scale[iratio]));
        scaleLabel->SetTextFont(43);
        scaleLabel->SetTextSize(20);
        scaleLabel->Draw();
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
    shift =  0;
    correlatedUnc = correlatedUncLocalPbPb;

  }

  else if (what == "PbPb_6080") {
    array = arrPbPb;
    system = 2; energy = 2760.; centrality = "V0M6080";
    label = "Pb-Pb, #sqrt{s}_{NN} = 2.76 TeV, 60-80%";
    color  = kBlue+1;
    marker = kFullCircle;
    shift =  0.0;
    correlatedUnc = correlatedUncLocalPbPb;
  }
  else if (what == "pp7") {
    array = arrpp7;
    system = 0; energy = 7000.; centrality = "";
    label = "pp #sqrt{s} = 7 TeV";
    color  = kMagenta+1;
    marker = kFullCircle;
    shift =  0.2;
    correlatedUnc = correlatedUncLocalPP;
  }
  else if (what == "pp900") {
    array = arrpp900;
    system = 0; energy = 900.; centrality = "";
    label = "pp #sqrt{s} = 0.9 TeV";
    color  = kCyan+2;
    marker = kFullCircle;
    shift =  -0.2;
    correlatedUnc = correlatedUncLocalPP;  
  }
  else if (what == "pp276") {
    array = arrpp276;
    system = 0; energy = 2760.; centrality = "";
    label = "pp #sqrt{s} = 2.76 TeV";
    color  = kYellow+2;
    marker = kFullCircle;
    shift = 0;
    correlatedUnc = correlatedUncLocalPP;
  }
  else if (what == "pPb0005") {
    array = arrpPb;
    system = 1; energy = 5020.; centrality = "V0A0005";
    label = "p-Pb, #sqrt{s}_{NN} = 5.02 TeV, V0A 0-5%";
    color  = kBlack;
    marker = kFullCircle;
    shift = -0.2;
    correlatedUnc = correlatedUncLocalPP;
  } 
  else if (what == "PbPbSTAR") {
    array = arrSTARPbPb;
    system = 2; energy = 200.; centrality = "0005";
    label = "STAR, Au-Au, #sqrt{s}_{NN} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenStar;
    shift = +0.2;
    correlatedUnc = correlatedUncZero;
  }
  else if (what == "PbPbPHENIX") {
    array = arrPHENIXPbPb;
    system = 2; energy = 200.; centrality = "0005";
    label = "PHENIX, Au-Au, #sqrt{s}_{NN} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenSquare;
    shift = -0.15;
    correlatedUnc = correlatedUncZero;
  }
  else if (what == "PbPbBRAHMS") {
    array = arrBRAHMSPbPb;
    system = 2; energy = 200.; centrality = "0010";
    label = "BRAHMS, Au-Au, #sqrt{s}_{NN} = 0.2 TeV, 0-10%";
    color  = kBlack;
    marker = kOpenCross;
    shift = -0.3;
    correlatedUnc = correlatedUncZero;
  }
  else if (what == "ppSTAR") {
    array = arrSTARpp;
    system = 0; energy = 200.; centrality = "";
    label = "STAR, pp, #sqrt{s} = 0.2 TeV";
    color  = kBlack;
    marker = kOpenStar;
    shift = 0.;
    correlatedUnc = correlatedUncZero;
  }
  else if (what == "ppPHENIX") {
    array = arrPHENIXpp;
    system = 0; energy = 200.; centrality = "";
    label = "PHENIX, pp, #sqrt{s} = 0.2 TeV";
    color  = kBlack;
    marker = kOpenSquare;
    shift = -0.2;
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
    DrawRatio("PbPb_0010");
    DrawRatio("pp7");
    DrawRatio("pPb0005");
    array =0;

    NewLegendQM(0.413655, 0.641754, 0.910643, 0.948736);

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

    NewLegendQM(0.413655, 0.677869, 0.910643, 0.948736);

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


    NewLegendQM(0.380522, 0.589587, 0.975904, 0.936697);

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



  else {
    std::cout << "Unknown ratio " << what.Data() << std::endl;
    return;
  }

  if(!correlatedUnc) {
    std::cout << "correlatedUnc not set!" << std::endl;
    
  }
  std::cout << "CORR: " << correlatedUnc[1] << std::endl;

  if(array) {
    AliPWGHistoTools::GetGraphFromHisto(GetHistoRatios(array,  system,  energy, centrality, label, color, marker, kStatError, shift)
                                        ,0)->Draw("PZ");
    AliPWGHistoTools::GetGraphFromHisto(GetHistoRatios(array,  system,  energy, centrality, label+"NoLegend", color, marker, kSystError, shift)
                                        ,0)->Draw("[]");
  }
  

}
void DrawFrame(Bool_t isYields) {

  myCan = new TCanvas("myCan","BD april.2014",50,10,1000,650);
  myCan->Draw();
  myCan->cd();
  // Set the Pads
  Double_t boundaryLabels = 0.85;//0.92;
  TPad    *myPadLabel = new TPad("myPadLabel","myPadLabel",0.0, boundaryLabels,1.0,1.0);
  myPadSetUp(myPadLabel);
  myPadLabel->Draw();

  TPad    *myPadHisto = new TPad("myPadHisto","myPadHisto",0.0,isYields ? 0.25 : 0.05 ,1.0, boundaryLabels,0);
  myPadSetUp(myPadHisto);
  myPadHisto->Draw();

  TPad    *myPadStdDev = new TPad("myPadStdDev","myPadStdDev",0.0,0.0,1.0,0.25,0);
  myPadSetUp(myPadStdDev);
  if(isYields)  myPadStdDev->Draw();

  myPadLabel->cd();

  double xLabelPosition[nratio] = {0.124498, 0.211847, 0.31, 0.38, 0.465, 0.575, 0.644, 0.72, 0.82, 0.905 }; 
  double yLabelPosition     = 0.40;

  // labels
  for(Int_t iratio = 0; iratio < nratio; iratio++){
    TLatex *myRatio = new TLatex(xLabelPosition[iratio],yLabelPosition,ratiosLabels[iratio]);
    myLatexDraw(myRatio,20);    
  }
  // Xi's and Omega's bar (there was no way to convince root to draw it properly)
  TLine *line = new TLine(0.408,0.75,0.418,0.75);
  line->SetLineWidth(2);
  line->Draw();
  line = new TLine(0.498,0.75,0.513,0.75);
  line->SetLineWidth(2);
  line->Draw();

  
  if(isYields) {
    myPadStdDev->cd();
    myPadStdDev->SetTopMargin(0.0);
    myPadStdDev->SetTicks(1,1);

    Float_t devMax = 3.5;
    
    TH2F *myBlankStdDev = new TH2F("myBlankStdDev","myBlankStdDev",18,0,18,10,-devMax,+devMax);
    myBlankStdDev->GetXaxis()->SetLabelFont(43); // precision 3: size will be in pixels
    myBlankStdDev->GetYaxis()->SetLabelFont(43);
    myBlankStdDev->GetYaxis()->SetTitleFont(43);
    myBlankStdDev->SetLabelSize(23,"xy");
    myBlankStdDev->SetTitleSize(26,"y");
    myBlankStdDev->SetNdivisions(505,"x");
    myBlankStdDev->SetNdivisions(605,"y");
    myBlankStdDev->SetLabelOffset(0.012,"xy");
    myBlankStdDev->SetYTitle("std. dev. ^{  }");
    myBlankStdDev->SetTitleOffset(1,"y");
    myBlankStdDev->Draw();
  }

  myPadHisto->cd();
  myPadHisto->SetBottomMargin(0.01);
  //  myPadHisto->SetLogy();
  myPadHisto->SetTicks(1,1);
  
  TH2F *myBlankHisto = new TH2F("NoLegend","NoLegend",nratio,1,nratio+1,10, isYields ? -0.01 : 0, maxy+0.01  );
  myBlankHisto->GetXaxis()->SetLabelFont(43); // precision 3: size will be in pixels
  myBlankHisto->GetYaxis()->SetLabelFont(43);
  myBlankHisto->GetYaxis()->SetTitleFont(43);
  myBlankHisto->SetLabelSize(23,"xy");
  myBlankHisto->SetTitleSize(26,"y");
  myBlankHisto->SetMaximum(10);
  myBlankHisto->SetMinimum(0);
  myBlankHisto->SetNdivisions(10,"x");
  myBlankHisto->SetNdivisions(505,"y");
  if(isYields) myBlankHisto->SetYTitle("d#it{N}/d#it{y}");
  myBlankHisto->SetLabelOffset(0.012,"xy");
  myBlankHisto->SetTitleOffset(1,"y");
  myBlankHisto->Draw();
  
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

