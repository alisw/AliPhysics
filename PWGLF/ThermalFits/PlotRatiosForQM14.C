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



#endif

// Plots ratios for QM and saves input files for thermal models


enum MyParticles { kPDGPi = 211, kPDGK = 321, kPDGProton = 2212, kPDGKS0 = 310, kPDGLambda=3122, kPDGXi=3312,kPDGOmega=3334,kPDGPhi=333,kPDGKStar=313,kPDGDeuteron=1000010020,kPDGHE3 = 1000020030, kPDGHyperTriton = 1010010030, kPDGSigmaStarPlus=3224,kPDGSigmaStarMinus=3114,kPDGXiStar=3324};

typedef enum {kStatError, kSystError, kTotalError} myerror_t;

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle, Int_t icolor, Int_t imarker = kOpenSquare, Int_t errorsType = kTotalError, Float_t shift = 0) ;
TH1F * GetHistoYields(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle) ;
void   PrepareThermalModelsInputFiles(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, Bool_t separateCharges=0) ;
void SetStyle(Bool_t graypalette=0) ;
void NewLegendQM() ;
void DrawRatio(TString what);
void LoadArrays() ;

// Preferred colors and markers
// const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7, kWhite}; // for syst bands
// const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2  , kWhite};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar,0};

Double_t maxy = 0.4;

// Data arrays;
TClonesArray *arrPbPb=0, *arrpp7=0, *arrpPb=0, * arrpp276=0, * arrpp900=0, * arrThermus=0;
TClonesArray *arrSTARPbPb=0, *arrPHENIXPbPb=0, *arrBRAHMSPbPb=0;
TClonesArray *arrSTARpp  =0, *arrPHENIXpp=0;

const Double_t *scaleRatios = 0;
TClonesArray * PlotRatiosForQM14() {
#if !(!defined (__CINT__) || (defined(__MAKECINT__)))
  LoadLibs();
#endif

  //
  LoadArrays();

  // Uncomment stuff in this section to save the inputs for thermal models
  //#define SAVE_INPUT_THERMAL_MODEL
#ifdef SAVE_INPUT_THERMAL_MODEL
  //PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/0);
  //  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/1);
  // PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/1);


  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A0005", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A2040", /*separateCharges*/1);
  PrepareThermalModelsInputFiles(arrpPb, AliParticleYield::kCSpPb, 5020, "V0A6080", /*separateCharges*/1);
  //  PrepareThermalModelsInputFiles(arrPbPb, AliParticleYield::kCSPbPb, 2760, "V0M6080", /*separateCharges*/1);

  return 0;
#endif

  SetStyle();

  TCanvas * c1 = new TCanvas("Ratios", "Ratios", 1400, 600);
  c1->SetMargin( 0.0744986, 0.0329513, 0.225131, 0.83);
  //  c1->SetLogy();

  // CENTRAL
  DrawRatio("frame");
  
  DrawRatio("PbPb_0010");
  DrawRatio("PbPbSTAR");
  //  DrawRatio("PbPbPHENIX");
  //  DrawRatio("PbPbBRAHMS");
  // DrawRatio("PbPb_6080");
  
  //  DrawRatio("pp7");
  //DrawRatio("pPb0005");
  //DrawRatio("pp276");
  //DrawRatio("pp900");

  NewLegendQM();

  return 0;
  // 
  TCanvas * c2 = new TCanvas("Yields", "Yields", 1400, 600);
  c2->SetMargin( 0.0744986, 0.0329513, 0.225131, 0.0593368);

  GetHistoYields(arrPbPb,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV, 0-10%")->Draw();
  GetHistoYields(arrThermus,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Thermus")->Draw("same");
  NewLegendQM();
  return arrPbPb;
}

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle, Int_t icolor, Int_t imarker, Int_t errorType, Float_t shift) {
  // FIXME: THIS SHOULD BE REVIEWED TO MAKE SURE THE PLOTS ARE LABELLED CORRECTLY

  const Int_t nratio = 10;
  Int_t num  [nratio]    = {kPDGK  , kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron , kPDGHE3      , kPDGHyperTriton , kPDGPhi , kPDGKStar};
  //  Int_t denum[nratio]    = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGPi       , kPDGDeuteron , kPDGPi          , kPDGK   , -kPDGK};
  Int_t denum[nratio]    = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGProton   , kPDGDeuteron , kPDGPi          , kPDGK   , -kPDGK};
  Int_t isSum[nratio]    = {1      ,1           ,1           ,1        ,1          ,0             ,0             ,1                ,0        ,1      };
  static const Double_t scale[] = {1      ,3           ,1           ,30       ,250         ,50            ,100           ,4e5              ,2       ,2      };
  scaleRatios = scale;
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
      if(part1) part1 = new AliParticleYield(*part1); // We need to clone it to avoid a mess if we need to use this particle again later
      // Try with the !sum, if part 1 is not found
      if(!part1) {
        std::cout << "  Looking for " <<  num[iratio-1] << " ("<<!isSum[iratio-1]<<")"<<std::endl;
        part1 = AliParticleYield::FindParticle(arr, num[iratio-1], system, energy, centrality,!isSum[iratio-1]);
        if(part1) 
          {
            part1 = new AliParticleYield(*part1); // We need to clone it to avoid a mess if we need to use this particle again later
            // If the sum was requested, try to recover it!
            if(isSum[iratio-1]  && TDatabasePDG::Instance()->GetParticle(-denum[iratio])) { // Before looking for anti particle, check if it makes sense (antiparticle is different from particle)
              std::cout << "  Looking for " <<  -num[iratio-1] <<std::endl;
              AliParticleYield * part1_bar = AliParticleYield::FindParticle(arr, -num[iratio-1], system, energy, centrality,0);
              if(part1 && part1_bar) {
                std::cout << "Adding " << part1_bar->GetPartName() << " " << part1->GetPartName() << std::endl;            
                part1 = AliParticleYield::Add(part1, part1_bar);
                
              } else {
                std::cout << "WARNING: Sum requested but not found, scaling x2 " << part1->GetName() << std::endl;
                part1->Scale(2);
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
      if(part2) part2 = new AliParticleYield(*part2); // We need to clone it to avoid a mess if we need to use this particle again later
      if(!part2) {// Try with the !sum, if part 2 is not found
        std::cout << "  Looking for " <<  denum[iratio-1] << " ("<<!isSum[iratio-1]<<")"<<std::endl;
        part2 = AliParticleYield::FindParticle(arr, denum[iratio-1], system, energy, centrality,!isSum[iratio-1]);
        if(part2) 
          {
            part2 = new AliParticleYield(*part2); // We need to clone it to avoid a mess if we need to use this particle again later
            if(isSum[iratio-1] && TDatabasePDG::Instance()->GetParticle(-denum[iratio])) { // Before looking for anti particle, check if it makes sense 
              std::cout << "  Looking for " <<  -denum[iratio-1] << std::endl;
              AliParticleYield * part2_bar = AliParticleYield::FindParticle(arr, -denum[iratio-1], system, energy, centrality,0);
              if(part2 && part2_bar){
                std::cout << "Adding " << part2_bar->GetPartName() << " " << part2->GetPartName() << std::endl;            
                part2 = AliParticleYield::Add(part2, part2_bar);
              } else {
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
      ratio = AliParticleYield::Divide(part1, part2, 0, "YQ"); // Assume by that the systematics of part1 and part2 are uncorrelated.
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

  AliParticleYield::SaveAsASCIIFile(arrRatios, TString("ratios_")+h->GetName());
  return h;



}

void   PrepareThermalModelsInputFiles(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, Bool_t separateCharges)  {
  // If "Separate charges" is true, tries to dump both charges are dumped
  TClonesArray * arrOut = new TClonesArray("AliParticleYield");
  const Int_t npart = 12;
  Int_t particles  [npart] = {kPDGPi ,kPDGK   ,kPDGKS0, kPDGKStar, kPDGPhi, kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron, kPDGHyperTriton, kPDGHE3    };
  Int_t isSum[npart]       = {1      ,1       ,0      , 1        , 0      , 1           ,0           ,1        ,1          ,0           , 1              , 0          };

  Int_t ipartOut = 0; // Index for the array
  for(Int_t ipart = 0; ipart < npart; ipart++){
    if(!separateCharges) {
      AliParticleYield * part = AliParticleYield::FindParticle(arr, particles[ipart], system, energy, centrality,  isSum[ipart]);
      if(!part && isSum[ipart]) {
        //Could not find the particle, but the sum was requested: build the sum!
        part = AliParticleYield::FindParticle(arr, particles[ipart], system, energy, centrality,  0);
        AliParticleYield * part2 = AliParticleYield::FindParticle(arr, -particles[ipart], system, energy, centrality,  0);
        if(part2 && part) part = AliParticleYield::Add(part, part2);
        else part = 0;
      }
      if(part) new((*arrOut)[ipartOut++]) AliParticleYield(*part);
    }
    else {
      // ignore isSum and try to find both particles
      Bool_t notFound = 0;
      AliParticleYield * part = AliParticleYield::FindParticle(arr, particles[ipart], system, energy, centrality,  0);
      if(part) new((*arrOut)[ipartOut++]) AliParticleYield(*part);
      else notFound=1;
      // Try to find antiparticle (-pdg code)
      part = AliParticleYield::FindParticle(arr, -particles[ipart], system, energy, centrality,  0);
      if(part) new((*arrOut)[ipartOut++]) AliParticleYield(*part);
      else if (notFound) {
        // If neither charge was found, check if we at least have the sum 
        part = AliParticleYield::FindParticle(arr, abs(particles[ipart]), system, energy, centrality,  1);
        if (part) new((*arrOut)[ipartOut++]) AliParticleYield(*part);
      }
      
    }
  }
  std::cout << "Particles for thermal model fits:" << std::endl; 
  arrOut->Print("short");
  std::cout << "" << std::endl;
  // Write GSI input file
  TIter it(arrOut);
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
  //  Int_t isSum[npart]    = {1      ,1      ,1           ,0           ,1        ,1          ,0             ,0             ,1                ,0        ,1      };
  Int_t isSum[npart]    = {0,0,0,0,0,0,0,0,0,0,1};
  //  Double_t scale[npart] = {1      ,1      ,3           ,1           ,30       ,250         ,50            ,10            ,4e5              ,2       ,2      };
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
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
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

void NewLegendQM() {

  const char * style = "lp";
  const char ** labels=0;
  //  Bool_t beautify=kFALSE;
  Bool_t useTitle=kTRUE;
  
  TLegend * l = new TLegend(  0.0985145, 0.733119, 0.301016, 0.869775);
  l->SetFillColor(kWhite);

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


void DrawRatio(TString what) {
  // This is used to simplify the code above
  // In order to draw syst error bars, we need to convert to graphs the syst errors histos

  TClonesArray * array = 0;
  Int_t system,  color, marker;
  Float_t energy = 0, shift = 0;
  TString centrality, label;

  if (what == "frame" ) {
    // This is a bit of an hack: since the particle labels come
    // directly from AliPArticleYield, and since the PbPb sample is
    // the only one where we have all the ratios, we draw the PbPb
    // ratio here and then we set lines and markers to white. We also
    // add the "NoLegend" flag, so that it does not show up in the legend 
    TH1 * h = GetHistoRatios(arrPbPb,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "NoLegend", kWhite);
    h->Draw();
    Int_t nratio = h->GetNbinsX();
    for(Int_t iratio = 0; iratio < nratio; iratio++){
      Double_t exp = TMath::Floor(TMath::Log10(TMath::Abs(scaleRatios[iratio])));  
      Double_t man = scaleRatios[iratio] / TMath::Power(10, exp);
      if(exp > 2) {
        TLatex * scaleLabel = new TLatex(iratio+1+0.2,maxy*1.01, Form("#times %0.0f 10^{%0.0f}", man, exp));
        scaleLabel->Draw();
      } else {
        TLatex * scaleLabel = new TLatex(iratio+1+0.2,maxy*1.01, Form("#times %g", scaleRatios[iratio]));
        scaleLabel->Draw();
      }      
    }

    TLatex *   tex = new TLatex(8.8,0.037,"ALICE Preliminary");
    tex->SetTextFont(52);
    tex->SetLineWidth(2);
    tex->Draw();

    h->GetYaxis()->SetDecimals(1);
    h->GetYaxis()->SetNdivisions(505);


  }
  else if (what == "PbPb_0010") {
    array = arrPbPb;
    system = 2; energy = 2760.; centrality = "V0M0010";
    label = "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV, 0-10%";
    color  = kRed+1;
    marker = kFullCircle;
    shift =  0;
  }

  else if (what == "PbPb_6080") {
    array = arrPbPb;
    system = 2; energy = 2760.; centrality = "V0M6080";
    label = "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV, 60-80%";
    color  = kBlue+1;
    marker = kFullCircle;
    shift =  0.1;
  }
  else if (what == "pp7") {
    array = arrpp7;
    system = 0; energy = 7000.; centrality = "";
    label = "pp #sqrt{s} = 7 TeV";
    color  = kCyan-8;
    marker = kFullCircle;
    shift =  0.2;
  }
  else if (what == "pp900") {
    array = arrpp900;
    system = 0; energy = 900.; centrality = "";
    label = "pp #sqrt{s} = 0.9 TeV";
    color  = kMagenta-9;
    marker = kFullCircle;
    shift =  -0.2;
  }
  else if (what == "pp276") {
    array = arrpp276;
    system = 0; energy = 2760.; centrality = "";
    label = "pp #sqrt{s} = 2.76 TeV";
    color  = kYellow-7;
    marker = kFullCircle;
    shift = 0;
  }
  else if (what == "pPb0005") {
    array = arrpPb;
    system = 0; energy = 5020.; centrality = "V0A0005";
    label = "p-Pb, #sqrt{s_{NN}} = 5.02 TeV, V0A 0-5%";
    color  = kBlack;
    marker = kFullCircle;
    shift = -0.2;
  } 
  else if (what == "PbPbSTAR") {
    array = arrSTARPbPb;
    system = 2; energy = 200.; centrality = "00";
    label = "STAR, Pb-Pb, #sqrt{s_{NN}} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenStar;
    shift = +0.2;
  }
  else if (what == "PbPbPHENIX") {
    array = arrPHENIXPbPb;
    system = 2; energy = 200.; centrality = "00";
    label = "PHENIX, Pb-Pb, #sqrt{s_{NN}} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenSquare;
    shift = -0.2;
  }
  else if (what == "PbPbBRAHMS") {
    array = arrBRAHMSPbPb;
    system = 2; energy = 200.; centrality = "00";
    label = "BRAHMS, Pb-Pb, #sqrt{s_{NN}} = 0.2 TeV, 0-5%";
    color  = kBlack;
    marker = kOpenCross;
    shift = -0.4;
  }


  else {
    std::cout << "Unknown ratio " << what.Data() << std::endl;
  }

  if(array) {
    GetHistoRatios(array,  system,  energy, centrality, label, color, marker, kStatError, shift)->Draw("same");
    AliPWGHistoTools::GetGraphFromHisto(GetHistoRatios(array,  system,  energy, centrality, label+"NoLegend", color, marker, kSystError, shift)
                                        ,0)->Draw("[]");
  }
  

}
