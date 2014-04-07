#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TClonesArray.h"
#include "AliParticleYield.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <fstream>
#endif

enum MyParticles { kPDGPi = 211, kPDGK = 321, kPDGProton = 2212, kPDGKS0 = 310, kPDGLambda=3122, kPDGXi=3312,kPDGOmega=3334,kPDGPhi=333,kPDGKStar=313,kPDGDeuteron=1000010020,kPDGHE3 = 1000020030, kPDGHyperTriton = 1010010030, kPDGSigmaStarPlus=3224,kPDGSigmaStarMinus=3114,kPDGXiStar=3324};

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality) ;
void   PrepareThermalModelsInputFiles(TClonesArray * arr, Int_t system, Float_t energy, TString centrality, Bool_t separateCharges=0) ;

// Plots ratios for QM and saves input files for thermal models

TClonesArray * PlotRatiosForQM14() {
  LoadLibs();
  TClonesArray * arr = AliParticleYield::ReadFromASCIIFile("PbPb_2760_Cascades.txt");
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_DeuHelium3.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_Hypertriton.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_Kstar892.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_LambdaK0.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_PiKaPr.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_phi1020.txt"));
  arr->AbsorbObjects(  AliParticleYield::ReadFromASCIIFile("PbPb_2760_AveragedNumbers.txt"));

  TClonesArray * arrpp7 = AliParticleYield::ReadFromASCIIFile("pp_7000.txt");

  // PrepareThermalModelsInputFiles(arr, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/0);
  // PrepareThermalModelsInputFiles(arr, AliParticleYield::kCSPbPb, 2760, "V0M0010", /*separateCharges*/1);

  TCanvas * c1 = new TCanvas("Ratios", "Ratios", 1400, 600);
  c1->SetLogy();
  //  GetHistoRatios(arr, AliParticleYield::kCSPbPb, 2760, "V0M0010")->Draw();
  GetHistoRatios(arr, AliParticleYield::kCSpp,   7000, ".*")->Draw("");
  return arr;
}

TH1F * GetHistoRatios(TClonesArray * arr, Int_t system, Float_t energy, TString centrality) {

  const Int_t nratio = 10;
  Int_t num  [nratio] = {kPDGK  , kPDGProton , kPDGLambda , kPDGXi  , kPDGOmega , kPDGDeuteron , kPDGHE3      , kPDGHyperTriton , kPDGPhi , kPDGKStar};
  Int_t denum[nratio] = {kPDGPi , kPDGPi     , kPDGKS0    ,  kPDGPi , kPDGPi    , kPDGProton   , kPDGDeuteron , kPDGPi          , kPDGK   , -kPDGK};
  Int_t isSum[nratio] = {1      ,1           ,0           ,1        ,1          ,0             ,0             ,1                ,0        ,1};
  TH1F * h = new TH1F("hRatio", "hRatio", nratio, 1, nratio+1);

  //  Double_t isSum = -1; // if this is -1, then the sum criterion is ignored
  for(Int_t iratio = 1; iratio <= nratio; iratio++){
    AliParticleYield * ratio = AliParticleYield::FindRatio(arr, num[iratio-1], denum[iratio-1], system, energy, centrality,isSum[iratio-1]);
    std::cout << num[iratio-1] << " " <<  denum[iratio-1]<< " " ;
    if(ratio)ratio->Print("short");


    if(!ratio) {
      // If the ratio is not found, try to build it!
      AliParticleYield * part1 = AliParticleYield::FindParticle(arr, num[iratio-1], system, energy, centrality,  isSum[iratio-1]);
      if(!part1) {// Try with the !sum, if part 1 is not found
        part1 = AliParticleYield::FindParticle(arr, num[iratio-1], system, energy, centrality,!isSum[iratio-1]);
      }
      AliParticleYield * part2 = AliParticleYield::FindParticle(arr, denum[iratio-1], system, energy, centrality,isSum[iratio-1]);
      if(!part2) {// Try with the !sum, if part 2 is not found
        part2 = AliParticleYield::FindParticle(arr, denum[iratio-1], system, energy, centrality,!isSum[iratio-1]);
      }
      ratio = AliParticleYield::Divide(part1, part2);
      if(ratio) {
        std::cout << "WARNING: building ratio " << num[iratio-1] <<"/"<<denum[iratio-1]<<": Check uncertainties!!" << std::endl;
        ratio->Print("short");
      }
    }
    if(ratio){
      h->SetBinContent(iratio, ratio->GetYield());
      h->SetBinError  (iratio, ratio->GetTotalError(0/* 0 = no normalization error */));
      h->GetXaxis()->SetBinLabel(iratio, Form("#splitline{%s}{%s}",ratio->GetCentr().Data() , ratio->GetLatexName()));
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
