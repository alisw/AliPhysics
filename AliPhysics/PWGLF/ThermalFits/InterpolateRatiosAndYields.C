#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TClonesArray.h"
#include "AliParticleYield.h"

#endif
TClonesArray * arr =0;

void InterpolateRatios(Int_t pdg1, Int_t pdg2, TString centr1="V0M0005", TString centr2="V0M0510", TString centrfinal="V0M0010") ;
void Interpolate0010(Int_t pdg) ;
void Interpolate2040(Int_t pdg) ;
void Interpolate6080(Int_t pdg) ;

void ExtrapolateWithConstantRatioToPions(Int_t pdg, TString centrOrigin, TString centrDest);
Int_t collSystem = -1;
Float_t energy = -1;

void InterpolateRatiosAndYields() {
#if !(!defined (__CINT__) || (defined(__MAKECINT__)))
  LoadLibs();
#endif
  collSystem = 2; energy =2760;
  // *************** pi, K, pr *****************
  //  arr=   AliParticleYield::ReadFromASCIIFile("PbPb_2760_PiKaPr.txt");
  // Interpolate0010(211);
  // Interpolate0010(-211);
  // Interpolate0010(321);
  // Interpolate0010(-321);
  // Interpolate0010(2212);
  // Interpolate0010(-2212);
  // InterpolateRatios(2212,211, "V0M0005", "V0M0510", "V0M0010");  
  // InterpolateRatios(321,211 , "V0M0005", "V0M0510", "V0M0010");  
  // Interpolate6080(211);
  // Interpolate6080(-211);
  // Interpolate6080(2212);
  // Interpolate6080(-2212);
  // Interpolate6080(321);
  // Interpolate6080(-321);
  // InterpolateRatios(2212,211, "V0M6070", "V0M7080", "V0M6080");  
  //  InterpolateRatios(321,211, "V0M6070", "V0M7080", "V0M6080");    

  // Interpolate2040(211);
  // Interpolate2040(-211);
  // Interpolate2040(2212);
  // Interpolate2040(-2212);
  // Interpolate2040(321);
  // Interpolate2040(-321);

  // *************** Lambda and K0 *****************
  // arr=   AliParticleYield::ReadFromASCIIFile("PbPb_2760_LambdaK0.txt");
  // Interpolate0010(3122);
  // Interpolate0010(310);
  // *************** Helium 3 *****************
  arr = AliParticleYield::ReadFromASCIIFile("PbPb_2760_DeuHelium3.txt");
  arr->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("./PbPb_2760_AveragedNumbers.txt"));
  // --> 0-10

  //  ExtrapolateWithConstantRatioToPions(1000020030, "V0M0020", "V0M0010");
  // --> 10-20
  arr->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("./pbpb_2760_pikapr.txt"));
  ExtrapolateWithConstantRatioToPions(1000020030, "V0M0020", "V0M1020");
  // *************** Kstar *****************
  // arr = AliParticleYield::ReadFromASCIIFile("PbPb_2760_Kstar892.txt");
  // arr->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("./PbPb_2760_AveragedNumbers.txt"));
  // ExtrapolateWithConstantRatioToPions(313, "V0M0020", "V0M0010");

  // *************** pPb, deuteron *********************
  //  collSystem = 1; energy = 5020;
  // 1. Average pions
  //arr = AliParticleYield::ReadFromASCIIFile("pPb_5020_PiKaPrLamndaK0.txt");
  //  Interpolate0010(211);
  //  Interpolate0010(-211);
  // 2. Extrapolate the deuteron
  // arr = AliParticleYield::ReadFromASCIIFile("pPb_5020_deuteron.txt");
  // arr->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("pPb_5020_PiKaPrLamndaK0.txt"));
  //  ExtrapolateWithConstantRatioToPions(1000010020, "V0A0010", "V0A0005");
  // ExtrapolateWithConstantRatioToPions(-1000010020, "V0A0010", "V0A0005");
  // ExtrapolateWithConstantRatioToPions(1000010020, "V0A6000", "V0A6080");
  // ExtrapolateWithConstantRatioToPions(-1000010020, "V0A6000", "V0A6080");


}

void DUMP(){
  AliParticleYield * p0 = AliParticleYield::FindParticle(arr, 211,  2, 2760., "V0M0005");
  //  p0 = AliParticleYield::Add(p0, AliParticleYield::FindParticle(arr, -211,  2, 2760., "V0M0005"));
  AliParticleYield * p2 = AliParticleYield::FindParticle(arr, 2212, 2, 2760., "V0M0005");
  //  p2 = AliParticleYield::Add(p2,AliParticleYield::FindParticle(arr, -2212, 2, 2760., "V0M0005"));
  AliParticleYield *pratio = AliParticleYield::Divide(p2,p0);
  pratio->Print();
  AliParticleYield::FindRatio(arr, 2212, 211, 2, 2760, "V0M0005")->Print();

}

void Interpolate0010(Int_t pdg) {

  TString centrPrefix = collSystem == 2 ? "V0M" : "V0A";

  AliParticleYield * p0 = AliParticleYield::FindParticle(arr, pdg, collSystem, energy, centrPrefix+"0005");
  AliParticleYield * p1 = AliParticleYield::FindParticle(arr, pdg, collSystem, energy, centrPrefix+"0510");
  p0->Scale(0.5);
  p1->Scale(0.5);
  AliParticleYield * pav = AliParticleYield::Add(p0,p1);
  std::cout << pav->GetYield() << ", " << pav->GetStatError() << ", "<< pav->GetSystError() << std::endl;


} 
void Interpolate6080(Int_t pdg) {

  TString centrPrefix = collSystem == 2 ? "V0M" : "V0A";

  AliParticleYield * p0 = AliParticleYield::FindParticle(arr, pdg, collSystem, energy, centrPrefix+"6070");
  AliParticleYield * p1 = AliParticleYield::FindParticle(arr, pdg, collSystem, energy, centrPrefix+"7080");
  p0->Scale(0.5);
  p1->Scale(0.5);
  AliParticleYield * pav = AliParticleYield::Add(p0,p1);
  std::cout << pav->GetYield() << ", " << pav->GetStatError() << ", "<< pav->GetSystError() << std::endl;


} 

void Interpolate2040(Int_t pdg) {

  TString centrPrefix = collSystem == 2 ? "V0M" : "V0A";

  AliParticleYield * p0 = AliParticleYield::FindParticle(arr, pdg, collSystem, energy, centrPrefix+"2030");
  AliParticleYield * p1 = AliParticleYield::FindParticle(arr, pdg, collSystem, energy, centrPrefix+"3040");
  p0->Scale(0.5);
  p1->Scale(0.5);
  AliParticleYield * pav = AliParticleYield::Add(p0,p1);
  std::cout << pav->GetYield() << ", " << pav->GetStatError() << ", "<< pav->GetSystError() << std::endl;


} 

void InterpolateRatios(Int_t pdg1, Int_t pdg2, TString centr1, TString centr2, TString centrfinal) {

  // Get the ratios from the DB for the correct centralities
  AliParticleYield * ratio[2];
  ratio[0] = AliParticleYield::FindRatio(arr, pdg1, pdg2, 2, 2760., centr1, 1);
  ratio[1] = AliParticleYield::FindRatio(arr, pdg1, pdg2, 2, 2760., centr2, 1);

  AliParticleYield * average = AliParticleYield::Add(ratio[0], ratio[1]);
  average->Scale(0.5);

  AliParticleYield * pi[2];
  pi[0] = AliParticleYield::FindParticle(arr, pdg2, 2, 2760., centr1, 0);
  pi[0] = AliParticleYield::Add(pi[0],AliParticleYield::FindParticle(arr, -pdg2, 2, 2760., centr1, 0));
  pi[1] = AliParticleYield::FindParticle(arr, pdg2, 2, 2760., centr2, 0);
  pi[1] = AliParticleYield::Add(pi[1],AliParticleYield::FindParticle(arr, -pdg2, 2, 2760., centr2, 0));
  

  // Scale to get protons with errors corresponding to the ratio
  ratio[0]->Scale(pi[0]->GetYield()) ;
  ratio[1]->Scale(pi[1]->GetYield()) ;

  ratio[0] = AliParticleYield::Add(ratio[0], ratio[1]);
  pi[0]    = AliParticleYield::Add(pi[0],pi[1]);
  pi[0]->SetNormError(0);
  pi[0]->SetStatError(0);
  pi[0]->SetSystError(0);
  
  ratio[0]->Scale(1./pi[0]->GetYield());
  ratio[0]->SetCentr(centrfinal);
  cout << "*** "<< ratio[0]->GetPartName() << " " <<  centrfinal << " ***"<< std::endl;
  std::cout << "RATIO OF AVERAGE: " ;
  ratio[0]->Print("justvalue");
  //  average->Dump();
  std::cout << "AVERAGE OF RATIO: " ;
  average->Print("justvalue");

    
}

void ExtrapolateWithConstantRatioToPions(Int_t pdg, TString centrOrigin, TString centrDest) {

  AliParticleYield * part       =  AliParticleYield::FindParticle(arr, pdg,collSystem, energy, centrOrigin);
  AliParticleYield * pionOrigin =  AliParticleYield::FindParticle(arr, 211,collSystem, energy, centrOrigin);
  AliParticleYield * pionDest   =  AliParticleYield::FindParticle(arr, 211,collSystem, energy, centrDest);
  if(!part || !pionOrigin | !pionDest) {
    return;
  }
  // The pion yield is takes
  part->Scale (1./pionOrigin->GetYield());
  part->Scale (pionDest->GetYield());
  part->SetCentr(centrDest);
  part->SetMeasurementType(part->GetMeasurementType() | AliParticleYield::kTypeExtrPionRatio);
  part->Print();
  //  std::cout << "1" << std::endl;
  TClonesArray * arrOut = new TClonesArray("AliParticleYield");
  //  std::cout << "2" << std::endl;
  new((*arrOut)[0]) AliParticleYield(*part) ;

  //  std::cout << "3" << std::endl;
  //  std::cout << "4" << std::endl;
  AliParticleYield::SaveAsASCIIFile(arrOut, "temp.txt");
}
