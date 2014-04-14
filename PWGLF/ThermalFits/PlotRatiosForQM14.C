#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TClonesArray.h"
#include "AliParticleYield.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <fstream>
#endif

enum MyParticles { kPDGPi = 211, kPDGK = 321, kPDGProton = 2212, kPDGKS0 = 310, kPDGLambda=3122, kPDGXi=3312,kPDGOmega=3334,kPDGPhi=333,kPDGKStar=313,kPDGDeuteron=1000010020,kPDGHE3 = 1000020030, kPDGHyperTriton = 1010010030, kPDGSigmaStarPlus=3224,kPDGSigmaStarMinus=3114,kPDGXiStar=3324};

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle) ;
TH1F * GetHistoYields(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle) ;
void   PrepareThermalModelsInputFiles(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, Bool_t separateCharges=0) ;

// Plots ratios for QM and saves input files for thermal models

TClonesArray * PlotRatiosForQM14() {
#if !(!defined (__CINT__) || (defined(__MAKECINT__)))
   LoadLibs();
#endif
  TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("PbPb_2760_Cascades.txt");
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_DeuHelium3.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_Hypertriton.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_Kstar892.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_LambdaK0.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_PiKaPr.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_phi1020.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_AveragedNumbers.txt"));

  TClonesArray * arrpp7   = AliParticleYield::ReadFromASCIIFile("pp_7000.txt");

  TClonesArray * arrpp276 = AliParticleYield::ReadFromASCIIFile("pp_2760.txt");
  TClonesArray * arrpp900 = AliParticleYield::ReadFromASCIIFile("pp_900.txt");

  TClonesArray * arrpPb   = AliParticleYield::ReadFromASCIIFile("pPb_5020_MultiStrange.txt");
  arrpPb->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("pPb_5020_PiKaPrLamndaK0.txt"));
  arrpPb->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("pPb_5020_deuteron.txt"));

  TClonesArray * arrThermus = AliParticleYield::ReadFromASCIIFile("PbPb_2760_Thermus_Boris_20140407.txt");

  // PrepareThermalModelsInputFiles(arr, AliParticleYield::kCSPbPb, 2760, "V0M0010", separateCharges0);
  // PrepareThermalModelsInputFiles(arr, AliParticleYield::kCSPbPb, 2760, "V0M0010", separateCharges1);
  // PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arrpp7, AliParticleYield::kCSpp, 7000, "", /*separateCharges*/1);

  TCanvas * c1 = new TCanvas("Ratios", "Ratios", 1400, 600);
  c1->SetMargin( 0.0744986, 0.0329513, 0.225131, 0.0593368);

  //  c1->SetLogy();
  // CENTRAL
  TH1 * h =  GetHistoRatios(arr,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV, 0-10%");
  h->GetYaxis()->SetRangeUser(0, 0.4);
  h->Draw();
  // //GetHistoRatios(arrThermus,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Thermus")->Draw("same");
  // GetHistoRatios(arrpp7, AliParticleYield::kCSpp,   7000, ""       , "pp #sqrt{s} = 7 TeV"                   )->Draw("same");
  // GetHistoRatios(arrpPb,    AliParticleYield::kCSpPb,  5020, "V0A0005", "p-Pb, #sqrt{s_{NN}} = 5.02 TeV, V0A 0-5%")->Draw("same");
  // // GetHistoRatios(arrpp276, AliParticleYield::kCSpp,   2760, ""       , "pp #sqrt{s} = 2.76 TeV"                   )->Draw("same");
  // // GetHistoRatios(arrpp900, AliParticleYield::kCSpp,   900, ""       , "pp #sqrt{s} = 0.9 TeV"                   )->Draw("same");
  // NewLegend("", "lp",0,1,1);

  // Peripheral
  GetHistoRatios(arr,       AliParticleYield::kCSPbPb, 2760, "V0M6080", "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV, 60-80%")->Draw("same");
  //GetHistoRatios(arrThermus,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Thermus")->Draw("same");
  // GetHistoRatios(arrpp7, AliParticleYield::kCSpp,   7000, ""       , "pp #sqrt{s} = 7 TeV"                   )->Draw("same");
  GetHistoRatios(arrpPb,    AliParticleYield::kCSpPb,  5020, "V0A0005", "p-Pb, #sqrt{s_{NN}} = 5.02 TeV, V0A 0-5%")->Draw("same");
  // GetHistoRatios(arrpp276, AliParticleYield::kCSpp,   2760, ""       , "pp #sqrt{s} = 2.76 TeV"                   )->Draw("same");
  // GetHistoRatios(arrpp900, AliParticleYield::kCSpp,   900, ""       , "pp #sqrt{s} = 0.9 TeV"                   )->Draw("same");
  NewLegend("", "lp",0,1,1);

  //return;
  TCanvas * c2 = new TCanvas("Yields", "Yields", 1400, 600);
  c2->SetMargin( 0.0744986, 0.0329513, 0.225131, 0.0593368);

  GetHistoYields(arr,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV, 0-10%")->Draw();
  GetHistoYields(arrThermus,       AliParticleYield::kCSPbPb, 2760, "V0M0010", "Thermus")->Draw("same");
  NewLegend("", "lp",0,1,1);
  return arr;
}

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, const char * histotitle) {
  // FIXME: THIS SHOULD BE REVIEWED TO MAKE SURE THE PLOTS ARE LABELLED CORRECTLY

  const Int_t nratio = 10;
  Int_t num  [nratio]    = {kPDGK  , kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron , kPDGHE3      , kPDGHyperTriton , kPDGPhi , kPDGKStar};
  //  Int_t denum[nratio]    = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGPi       , kPDGDeuteron , kPDGPi          , kPDGK   , -kPDGK};
  Int_t denum[nratio]    = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGProton   , kPDGDeuteron , kPDGPi          , kPDGK   , -kPDGK};
  Int_t isSum[nratio]    = {1      ,1           ,0           ,1        ,1          ,0             ,0             ,1                ,0        ,1      };
  Double_t scale[nratio] = {1      ,3           ,1           ,30       ,250         ,50            ,100           ,4e5              ,2       ,2      };
  TH1F * h = new TH1F(Form("hRatio_%d_%0.0f_%s",system,energy,centrality.Data()), histotitle, nratio, 1, nratio+1);

  //  Double_t isSum = -1; // if this is -1, then the sum criterion is ignored
  for(Int_t iratio = 1; iratio <= nratio; iratio++){
    AliParticleYield * ratio = AliParticleYield::FindRatio(arr, num[iratio-1], denum[iratio-1], system, energy, centrality,isSum[iratio-1]);
    std::cout << num[iratio-1] << " " <<  denum[iratio-1]<< " " ;
    if(ratio)ratio->Print("short");


    if(!ratio) {
      // If the ratio is not found, try to build it!
      AliParticleYield * part1 = AliParticleYield::FindParticle(arr, num[iratio-1], system, energy, centrality,  isSum[iratio-1]);
      // Try with the !sum, if part 1 is not found
      if(!part1) {
        part1 = AliParticleYield::FindParticle(arr, num[iratio-1], system, energy, centrality,!isSum[iratio-1]);
        // If the sum was requested, try to recover it!
        if(isSum[iratio-1]) { 
          AliParticleYield * part1_bar = AliParticleYield::FindParticle(arr, -num[iratio-1], system, energy, centrality,0);
          if(part1 && part1_bar) {
            std::cout << "Adding " << part1_bar->GetPartName() << " " << part1->GetPartName() << std::endl;
            
            part1 = AliParticleYield::Add(part1, part1_bar);

          }
        } else if(part1) {
          // if the not sum was requested, but the sum is found, divide by 2 so that it is comparable
          part1->Scale(0.5);
        }
 
      }
      AliParticleYield * part2 = AliParticleYield::FindParticle(arr, denum[iratio-1], system, energy, centrality,isSum[iratio-1]);
      if(!part2) {// Try with the !sum, if part 2 is not found
        part2 = AliParticleYield::FindParticle(arr, denum[iratio-1], system, energy, centrality,!isSum[iratio-1]);
        if(isSum[iratio-1]) { 
          AliParticleYield * part2_bar = AliParticleYield::FindParticle(arr, -denum[iratio-1], system, energy, centrality,0);
          if(part2 && part2_bar) part2 = AliParticleYield::Add(part2, part2_bar);
        } else if(part2){
          // if the not sum was requested, but the sum is found, divide by 2 so that it is comparable
          part2->Scale(0.5);
        } 

      }
      ratio = AliParticleYield::Divide(part1, part2, 0, "YQ");
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
      h->SetBinError  (iratio, ratio->GetTotalError(0/* 0 = no normalization error */));
      h->GetXaxis()->SetBinLabel(iratio, Form("#splitline{%s}{%s}",Form("#times%2.2f",  scale[iratio-1]), ratio->GetLatexName()));
    }
    else {
      h->GetXaxis()->SetBinLabel(iratio, Form("#frac{%d}{%d}",num[iratio-1], denum[iratio-1]));
      
    }
  }
  
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
