#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "TClonesArray.h"
#include "AliParticleYield.h"

#endif
TClonesArray * arr =0;

void InterpolateRatios(Int_t pdg1, Int_t pdg2) ;
void Interpolate0010(Int_t pdg) ;
void ExtrapolateWithConstantRatioToPions(Int_t pdg, TString centrOrigin, TString centrDest);


void InterpolateRatiosAndYields() {
#if !(!defined (__CINT__) || (defined(__MAKECINT__)))
  LoadLibs();
#endif
  // *************** pi, K, pr *****************
  //  arr=   AliParticleYield::ReadFromASCIIFile("PbPb_2760_PiKaPr.txt");
  // Interpolate0010(211);
  // Interpolate0010(-211);
  // Interpolate0010(321);
  // Interpolate0010(-321);
  // Interpolate0010(2212);
  // Interpolate0010(-2212);
  //  InterpolateRatios(2212,211);  
  //  InterpolateRatios(321,211);  
  // *************** Lambda and K0 *****************
  // arr=   AliParticleYield::ReadFromASCIIFile("PbPb_2760_LambdaK0.txt");
  // Interpolate0010(3122);
  // Interpolate0010(310);
  // *************** Helium 3 *****************
  // arr = AliParticleYield::ReadFromASCIIFile("PbPb_2760_DeuHelium3.txt");
  // arr->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("./PbPb_2760_AveragedNumbers.txt"));
  // ExtrapolateWithConstantRatioToPions(1000020030, "V0M0020", "V0M0010");
  // *************** Kstar *****************
  arr = AliParticleYield::ReadFromASCIIFile("PbPb_2760_Kstar892.txt");
  arr->AbsorbObjects(AliParticleYield::ReadFromASCIIFile("./PbPb_2760_AveragedNumbers.txt"));
  ExtrapolateWithConstantRatioToPions(313, "V0M0020", "V0M0010");


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

  AliParticleYield * p0 = AliParticleYield::FindParticle(arr, pdg, 2, 2760., "V0M0005");
  AliParticleYield * p1 = AliParticleYield::FindParticle(arr, pdg, 2, 2760., "V0M0510");
  p0->Scale(0.5);
  p1->Scale(0.5);
  AliParticleYield * pav = AliParticleYield::Add(p0,p1);
  std::cout << pav->GetYield() << ", " << pav->GetStatError() << ", "<< pav->GetSystError() << std::endl;


} 

void InterpolateRatios(Int_t pdg1, Int_t pdg2) {

  AliParticleYield * ratio[2];
  ratio[0] = AliParticleYield::FindRatio(arr, pdg1, pdg2, 2, 2760., "V0M0005", 1);
  ratio[1] = AliParticleYield::FindRatio(arr, pdg1, pdg2, 2, 2760., "V0M0510", 1);
  AliParticleYield * average = AliParticleYield::Add(ratio[0], ratio[1]);
  average->Scale(0.5);
  AliParticleYield * pi[2];
  pi[0] = AliParticleYield::FindParticle(arr, pdg2, 2, 2760., "V0M0005", 0);
  pi[0] = AliParticleYield::Add(pi[0],AliParticleYield::FindParticle(arr, -pdg2, 2, 2760., "V0M0005", 0));
  pi[1] = AliParticleYield::FindParticle(arr, pdg2, 2, 2760., "V0M0510", 0);
  pi[1] = AliParticleYield::Add(pi[1],AliParticleYield::FindParticle(arr, -pdg2, 2, 2760., "V0M0510", 0));
  

  // Scale to get protons with errors corresponding to the ratio
  ratio[0]->Scale(pi[0]->GetYield()) ;
  ratio[1]->Scale(pi[1]->GetYield()) ;

  ratio[0]->Add(ratio[0], ratio[1]);
  pi[0]->Add(pi[0],pi[1]);
  pi[0]->SetNormError(0);
  pi[0]->SetStatError(0);
  pi[0]->SetSystError(0);
  
  ratio[0]->Scale(1./pi[0]->GetYield());
    ratio[0]->SetCentr("V0M0010");

  ratio[0]->Print();
  //  average->Dump();
  average->Print();

    
}

void ExtrapolateWithConstantRatioToPions(Int_t pdg, TString centrOrigin, TString centrDest) {

  AliParticleYield * part       =  AliParticleYield::FindParticle(arr, pdg, 2, 2760., centrOrigin);
  AliParticleYield * pionOrigin =  AliParticleYield::FindParticle(arr, 211, 2, 2760., centrOrigin);
  AliParticleYield * pionDest   =  AliParticleYield::FindParticle(arr, 211, 2, 2760., centrDest);
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
