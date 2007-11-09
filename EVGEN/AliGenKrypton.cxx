/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/*
$Log$
Revision 1.2  2007/10/23 09:27:16  hristov
Adding dependence on the dip angle (Marek)

Revision 1.1  2007/06/24 20:53:11  hristov
New generator for the krypton runs of TPC (Marek)

*/ 

//
// generates single Kr decay, in order to generate the calibration data
// one should use it together with the AliGenCocktail class
//
#include "AliGenKrypton.h"
#include "TPDGCode.h"
//
ClassImp(AliGenKrypton)

//________________________________________________________________________
AliGenKrypton::AliGenKrypton(){
  //
  // Default constructor
  //
}
//________________________________________________________________________
void AliGenKrypton::Generate(){
  Double_t eelectron[6];
  Double_t egamma[2];
  Float_t polar[3]={0.,0.,0.};
  Int_t nelectron, ngamma;
  Int_t nt;
  //
  // generate decay vertex within the gas TPC volume
  //
  Float_t rmin,rmax,zmax;
  zmax=249.7;
  rmin=78.8;
  rmax=258.;
  Float_t me=511.e-6;
  Double_t rnd;
  //
  rnd=gRandom->Rndm();
  Float_t r = (rmax-rmin)*rnd+rmin;
  rnd=gRandom->Rndm();
  Float_t phi=TMath::TwoPi()*rnd;
  //
  Float_t origin[3];
  //
  rnd=gRandom->Rndm();
  origin[2]=zmax*(2.*rnd-1.);
  origin[0]=r*TMath::Cos(phi);
  origin[1]=r*TMath::Sin(phi);
  //
  Float_t ptot,p[3];
  //
  // generate decay
  //
  KrDecay(nelectron,ngamma,eelectron,egamma);
  //
  // electrons
  //
  for(Int_t i=0;i<nelectron;i++){
    rnd=gRandom->Rndm();
    phi=TMath::TwoPi()*rnd; 
    rnd=gRandom->Rndm();
    Double_t theta = TMath::Pi()*rnd;   
    ptot=TMath::Sqrt(eelectron[i]*(eelectron[i]+2.*me));
    p[0]=ptot*TMath::Cos(phi)*TMath::Sin(theta);
    p[1]=ptot*TMath::Sin(phi)*TMath::Sin(theta);    
    p[2]=ptot*TMath::Cos(theta);
    //
    // her push particle
    //
    PushTrack(fTrackIt,-1,kElectron,p,origin,polar,0,kPPrimary,nt);
  }
  //
  // gammas
  //
  for(Int_t i=0;i<ngamma;i++){
    rnd=gRandom->Rndm();
    Double_t theta = TMath::Pi()*rnd;
    phi=TMath::TwoPi()*rnd;    
    ptot=egamma[i];
    p[0]=ptot*TMath::Cos(phi)*TMath::Sin(theta);
    p[1]=ptot*TMath::Sin(phi)*TMath::Sin(theta);
    rnd=gRandom->Rndm();
    p[2]=ptot*TMath::Cos(theta);
    //
    // her push particle
    //
    PushTrack(fTrackIt,-1,kGamma,p,origin,polar,0,kPPrimary,nt);
  }  
}
//________________________________________________________________________
void AliGenKrypton::KrDecay(Int_t &nelectron, Int_t &ngamma, Double_t *eelectron, Double_t *egamma)
{
  Double_t prob1[2]={0.76,0.88}; // 0.76, 0.12, 0.12
  Double_t prob2=0.95;           // 0.95, 0.05
  nelectron=0;
  ngamma=0;
  
  Double_t rnd;
  rnd = gRandom->Rndm();
  //
  //
  // first decay - 32 keV
  //
   if(rnd < prob1[0]) {
    // 2 electrons
    nelectron = 2;
    eelectron[0]=30.e-6;
    eelectron[1]=1.8e-6;
  }
  else if (rnd > prob1[1]){
    // 4 electrons
    nelectron=4;
    eelectron[0]=18.e-6;
    eelectron[1]=10.e-6;
    eelectron[2]=1.8e-6;
    eelectron[3]=1.8e-6;
  }
  else {
    // 2 electrons + 1 gamma
    nelectron = 2;
    ngamma = 1;
    eelectron[0]=18.e-6;
    eelectron[1]=1.8e-6;
    egamma[0]=12.e-6;
  }
  //
  //  second decay - 9 keV
  // 
  rnd=gRandom->Rndm();
  //
  if(rnd < prob2){
    // 2 electrons
    nelectron+=2;
    eelectron[nelectron-2]=7.6e-6;
    eelectron[nelectron-1]=1.8e-6;
  }
  else {
    ngamma++;
    egamma[ngamma-1]=9.e-6;
  }
  

}
//________________________________________________________________________
