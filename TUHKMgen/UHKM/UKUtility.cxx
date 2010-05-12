//                                                                           
//                                                                            
//        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
//      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru
//                           November. 2, 2005                                
//
//

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"

#include "UKUtility.h" 

//_____________________________________________________________________
void IsotropicR3(Double_t r, Double_t *x, Double_t *y, Double_t *z) {
  //
  // return a random isotropic orientation
  //
  Double_t pZ  = 1. - 2.*(gRandom->Rndm());
  Double_t st  = TMath::Sqrt((1.-pZ)*(1.+pZ)) * r;
  Double_t phi = 2. * TMath::Pi() * (gRandom->Rndm());

  *x = st * cos(phi);
  *y = st * sin(phi);
  *z = pZ * r;
}

//_____________________________________________________________________
void IsotropicR3(Double_t r, TVector3 &pos) {
  //
  // return a random isotropic orientation
  //
  Double_t pZ  = 1. - 2.* (gRandom->Rndm());  
  Double_t st  = TMath::Sqrt((1.-pZ)*(1.+pZ)) * r;
  Double_t phi = 2. * TMath::Pi() * (gRandom->Rndm());

  pos.SetX(st * TMath::Cos(phi));
  pos.SetY(st * TMath::Sin(phi));
  pos.SetZ(pZ * r);
}

//_____________________________________________________________________
void MomAntiMom(TLorentzVector &mom, Double_t mass, TLorentzVector &antiMom, 
		Double_t antiMass, Double_t initialMass) {
  //
  // perform a 2 - body decay and orientate randomly the product particles momentum vectors
  //
  Double_t r = initialMass * initialMass - mass * mass - antiMass * antiMass;
  if (r * r - 4 * mass * mass * antiMass * antiMass < 0.) throw "MomAntiMom";
      
  Double_t pAbs = .5 * TMath::Sqrt(r * r - 4 * mass * mass * antiMass * antiMass) / initialMass;
  TVector3 mom3;
  IsotropicR3(pAbs, mom3);
  mom.SetVectM(mom3, mass);
  antiMom.SetVectM(- mom3, antiMass);
}
