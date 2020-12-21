/*************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//    Class to smear pt,phi,eta for HF cocktail (mother analysis task)   //
//                                        (description in .h file)       //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// c++ includes
#include <iostream>
//using namespace std;
// ROOT includes
#include <TMath.h>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TObjArray.h"
#include "AliVParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
// this task
#include "AliCocktailSmearing.h"

ClassImp(AliCocktailSmearing)

//_______________________________________________________________________________________________
AliCocktailSmearing::AliCocktailSmearing():
fUseRelPResolution(kFALSE),
  fPResArr(0x0),
  fOpeningAngleResArr(0x0),
  fArrResoPt(0x0),
  fArrResoEta(0x0),
  fArrResoPhi_Pos(0x0),
  fArrResoPhi_Neg(0x0),
  fOton(kFALSE)
{
  // Constructor
}

//_______________________________________________________________________________________________
AliCocktailSmearing::~AliCocktailSmearing()
{
  // Destructor
}


//_______________________________________________________________________________________________
void AliCocktailSmearing::ReadResoFile(TFile *fRes)
{
  //
  // Set resolution arrays
  //
  
  TObjArray *ArrResoPt = 0x0;
  TObjArray *ArrResoEta = 0x0;
  TObjArray *ArrResoPhi_Pos = 0x0;
  TObjArray *ArrResoPhi_Neg = 0x0;
  
  // Set resolutions: default Run 2
  if(fRes && fRes->IsOpen()){
    // Oton
    ArrResoPt = (TObjArray *)fRes->Get("RelPtResArrCocktail");
    ArrResoEta = (TObjArray *)fRes->Get("EtaResArrVsPt");
    if (ArrResoEta == 0x0) {
      ArrResoEta = (TObjArray *)fRes->Get("EtaResArr");
    }
    ArrResoPhi_Pos = (TObjArray *)fRes->Get("PhiPosResArrVsPt");
    if (ArrResoPhi_Pos == 0x0) {
      ArrResoPhi_Pos = (TObjArray *)fRes->Get("PhiPosResArr");
      std::cout << "Use PhiPosResArr" << std::endl;
    }
    ArrResoPhi_Neg = (TObjArray *)fRes->Get("PhiEleResArrVsPt");
    if (ArrResoPhi_Neg == 0x0) {
      ArrResoPhi_Neg = (TObjArray *)fRes->Get("PhiEleResArr");
      std::cout << "Use PhiEleResArr" << std::endl;
    }
  }
  else Printf("no resolution file given!");

  fArrResoPt = ArrResoPt;
  fArrResoEta = ArrResoEta;
  fArrResoPhi_Pos = ArrResoPhi_Pos;
  fArrResoPhi_Neg = ArrResoPhi_Neg;

  if(fArrResoPt && fArrResoEta && fArrResoPhi_Pos && fArrResoPhi_Neg) fOton = kTRUE;

  // if not found resolution Run 2, set resolution histos Run 1 old if there
  if(!fOton) {
    // Set resolutions
    if(fRes && fRes->IsOpen()){
      SetResolutionP           ((TObjArray*) fRes->Get("RelPResArr"),kTRUE);
      if (!fPResArr) {
	SetResolutionP         ((TObjArray*) fRes->Get("DeltaPResArr"),kFALSE);
      }
      if (fRes->Get("OpeningAngleResArr"))
	SetResolutionOpeningAngle((TObjArray*) fRes->Get("OpeningAngleResArr"));
      else
	SetResolutionOpeningAngle((TObjArray*) fRes->Get("OpAngleResArr_US"));
    }
  }

  
}


//_______________________________________________________________________________________________
void AliCocktailSmearing::SetSeed(UInt_t rndmseed)
{
  //
  // Set seed
  //
  TRandom3 *rand = new TRandom3(0);
  if (rndmseed>0) rand->SetSeed(rndmseed);
  gRandom = rand;
}


//_______________________________________________________________________________________________
void AliCocktailSmearing::Print() {
  //
  // Print info
  //

  if(!fPResArr)
    printf("\n\n========================================\n  No array for smearing! \n========================================\n\n");
  else
    printf("\n\n========================================\n  Smearing arrays: \n========================================\n\n");
  
  std::cout << "  pResArr:             " << fPResArr ;
  if(fUseRelPResolution) std::cout << "   relative momentum resolution (Prec/Pgen) will be used" << std::endl;
  else std::cout << "   momentum difference (or old input DeltaP/Pgen = (Pgen-Prec)/Pgen) will be used" << std::endl;
  std::cout << "  fOpeningAngleResArr: " << fOpeningAngleResArr << std::endl;
  std::cout << std::endl;
  std::cout << "  Smearing from Oton:         " << fOton << std::endl;
  std::cout << "  fArrResoPt:         " << fArrResoPt << std::endl;
  std::cout << "  fArrResoEta:           " << fArrResoEta << std::endl;
  std::cout << "  fArrResoPhi_Pos:        " << fArrResoPhi_Pos << std::endl;
  std::cout << "  fArrResoPhi_Neg:        " << fArrResoPhi_Neg << std::endl;
  std::cout << std::endl;
  std::cout << "  random number seed:  " << gRandom->GetSeed() << std::endl;
  std::cout << std::endl;

  
}
//_____________________________________________________________________________________________
TLorentzVector  AliCocktailSmearing::ApplySmearingOton(const TLorentzVector& vec, short ch)
{
  //
  // Smearing in pt, eta, phi: Run 2 method
  //
  
  //Double_t theta, phi, pt, p, mass, eta;
  Double_t phi, pt, mass, eta;
  TLorentzVector resvec;
  
  mass = 0.51099906e-3;
  pt = vec.Pt();
  //p = vec.P();
  //theta = vec.Theta();
  phi = vec.Phi();
  eta = vec.Eta();
  
  // smear pt
  Int_t ptbin     = ((TH2D *)(fArrResoPt->At(0)))->GetXaxis()->FindBin(pt);
  Int_t ptbin_max = ((TH2D *)(fArrResoPt->At(0)))->GetXaxis()->GetNbins();
  // make sure that no underflow or overflow bins are used
  if (ptbin < 1)
    ptbin = 1;
  else if (ptbin > ptbin_max)
    ptbin = ptbin_max;
  Double_t smearing = ((TH1D *)(fArrResoPt->At(ptbin)))->GetRandom() * pt;
  const Double_t sPt = pt - smearing;
  
  // smear eta
  ptbin     = ((TH2D *)(fArrResoEta->At(0)))->GetXaxis()->FindBin(pt);
  ptbin_max = ((TH2D *)(fArrResoEta->At(0)))->GetXaxis()->GetNbins();
  if (ptbin < 1)
    ptbin = 1;
  else if (ptbin > ptbin_max)
    ptbin = ptbin_max;
  smearing = ((TH1D *)(fArrResoEta->At(ptbin)))->GetRandom();
  const Double_t sEta = eta - smearing;
  
  // smear phi
  ptbin     = ((TH2D *)(fArrResoPhi_Pos->At(0)))->GetXaxis()->FindBin(pt);
  ptbin_max = ((TH2D *)(fArrResoPhi_Pos->At(0)))->GetXaxis()->GetNbins();
  if (ptbin < 1)
    ptbin = 1;
  if (ptbin > ptbin_max)
    ptbin = ptbin_max;
  if (ch > 0) {
    smearing = ((TH1D *)(fArrResoPhi_Pos->At(ptbin)))->GetRandom();
  } else if (ch < 0) {
    smearing = ((TH1D *)(fArrResoPhi_Neg->At(ptbin)))->GetRandom();
  }
  const Double_t sPhi = phi - smearing;
  
  // printf(" Original Pt = %f Phi %f Eta %f -> final pt = %f Phi %f Eta %f
  // \n",pt,phi,eta,sPt,sPhi,sEta);
  
  const Double_t sPx = sPt * cos(sPhi);
  const Double_t sPy = sPt * sin(sPhi);
  const Double_t sPz = sPt * sinh(sEta);
  const Double_t sP = sPt * cosh(sEta);
  const Double_t sE = sqrt(sP * sP + mass * mass);
  
  resvec.SetPxPyPzE(sPx, sPy, sPz, sE);
  return resvec;
}
//_______________________________________________________________________________________________
TLorentzVector AliCocktailSmearing::Smear(AliVParticle * part)
{
  //
  // Smear of |p|: Run 1 method
  //
  //double theta, phi, p, px, py, pz, pt, E, mass, sigma;
  double theta, phi, p, px, py, pz, pt, E, mass;
  TLorentzVector result;
  
  mass = 0.51099906e-3;   
  p  = part->P();
  pt = part->Pt();
  theta = part->Theta();
  phi   = part->Phi();
  
  if(!fPResArr) {
    //printf("No smearing\n");
    px   = p*sin(theta)*cos(phi);
    py   = p*sin(theta)*sin(phi);
    pz   = p*cos(theta);
    E    = sqrt(p*p + mass*mass);
    result.SetPxPyPzE(px,py,pz,E);
    return (result);
  }
  
  //printf("Smearing\n");
  TH2D * hDeltaPtvsPt = (TH2D*)fPResArr->At(0);
  // Get the momentum slice histogram for sampling of the smearing
  // (in some input file versions this also contains a langau fit )
  // since histogram bins start at 1, we can use the bin number directly for the array index (first slice stored in position 1).
  Int_t histIndex = TMath::Min( hDeltaPtvsPt->GetXaxis()->FindBin(pt), fPResArr->GetLast() );
  if (histIndex<1) histIndex=1; // just in case some track is below the first p-bin (which currently starts at 100 MeV).
  TH1D* hisSlice = (TH1D*)fPResArr->At(histIndex);
  
  if (fUseRelPResolution==kFALSE) { // old case
    Double_t SmearingPPercent = hisSlice->GetRandom();
    Double_t SmearingP = p * SmearingPPercent;
    p -= SmearingP;
  }
  else {
    if (hisSlice->GetEntries()>0)
      p *= hisSlice->GetRandom();
  }
  
  px   = p*sin(theta)*cos(phi);
  py   = p*sin(theta)*sin(phi);
  pz   = p*cos(theta);
  E    = sqrt(p*p + mass*mass);
  
  //cout << p << "  " << px << "   " << py << "  " << pz << "   " << E << endl;
  
  result.SetPxPyPzE(px,py,pz,E);
  return (result);
}


// from LightFlavorGenerator
//_______________________________________________________________________________________________
Double_t AliCocktailSmearing::GetSmearing(TObjArray *arr, Double_t x)
{
  //
  // Smearing of opening angle: Run 1 method
  //

  TH1D *hisSlice(0x0);
  TH2D *hDeltaXvsXgen = static_cast<TH2D*> (arr->At(0));
  // Get the slice histogram for sampling of the smearing.
  // Since histogram bins start at 1, we can use the bin number directly for the array index (first slice stored in position 1).
  Int_t histIndex = TMath::Min( hDeltaXvsXgen->GetXaxis()->FindBin(x), arr->GetLast() );
  if (histIndex<1) histIndex=1; // Just in case some values lie below the first bin.
  hisSlice = static_cast<TH1D*> (arr->At(histIndex));
  //printf("Smearing of %s. slice=%d, entries=%f\n", arr->GetName(), histIndex, hisSlice->GetEntries());
  // get smear parameter via random selection from the x slices retreived from the deltax plot
  Double_t smearing(0.);
  if(fUseRelPResolution) smearing = 1.;
  if(hisSlice) smearing = hisSlice->GetRandom();
  return smearing;
}
// from LightFlavorGenerator
//_______________________________________________________________________________________________
void AliCocktailSmearing::SmearOpeningAngle(TLorentzVector &lv1, TLorentzVector &lv2)
{
  //
  // Run 1 method
  //
  Double_t opAngle = lv1.Angle(lv2.Vect());
  Double_t smearVal(0.);
  if(fOpeningAngleResArr) smearVal = GetSmearing(fOpeningAngleResArr,opAngle);
  else return;
  TVector3 normalV = (lv1.Vect()).Cross(lv2.Vect());
  TVector3 rotAxis = normalV.Cross(fZaxis);
  Double_t rotAngle = normalV.Angle(fZaxis);
  lv1.Rotate(rotAngle,rotAxis);
  lv2.Rotate(rotAngle,rotAxis);
  Double_t phi1(lv1.Phi()),phi2(lv2.Phi());
  if(TMath::Abs(phi1-phi2)>TMath::Pi()){
    if(phi1>phi2){
      phi1 = phi1 - 0.5*smearVal;
      phi2 = phi2 + 0.5*smearVal;
    }
    else{
      phi1 = phi1 + 0.5*smearVal;
      phi2 = phi2 - 0.5*smearVal;
    }
  }
  else{
    if(phi1>phi2){
      phi1 = phi1 + 0.5*smearVal;
      phi2 = phi2 - 0.5*smearVal;
    }
    else{
      phi1 = phi1 - 0.5*smearVal;
      phi2 = phi2 + 0.5*smearVal;
    }
  }
  lv1.SetPtEtaPhiM(lv1.Pt(),lv1.Eta(),phi1,eMass());
  lv2.SetPtEtaPhiM(lv2.Pt(),lv2.Eta(),phi2,eMass());
  lv1.Rotate(-rotAngle,rotAxis);
  lv2.Rotate(-rotAngle,rotAxis);
}
