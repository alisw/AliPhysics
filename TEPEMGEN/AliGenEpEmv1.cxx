/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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
 *                                                                        *
 *                                                                        *
 * Copyright(c) 1997, 1998, 2002, Adrian Alscher and Kai Hencken          *
 * See $ALICE_ROOT/EpEmGen/diffcross.f for full Copyright notice          *
 *                                                                        *
 *                                                                        *
 * Copyright(c) 2002 Kai Hencken, Yuri Kharlov, Serguei Sadovsky          *
 * See $ALICE_ROOT/EpEmGen/epemgen.f for full Copyright notice            *
 *                                                                        *
 **************************************************************************/

/* $Id$ */

// Event generator of single e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 TeV/nucleon.
// The generator is based on 5-dimentional differential cross section of the process.
//%
// References:
// [1] "Multiple electromagnetic electron positron pair production in
//      relativistic heavy ion collisions".
//      Adrian Alscher, Kai Hencken, Dirk Trautmann, and Gerhard Baur,
//      Phys. Rev. A55 (1997) 396.
// [2] K.Hencken, Yu.Kharlov, S.Sadovsky, Internal ALICE Note 2002-27.
//%
// Usage:
// Initialization:
//    AliGenEpEmv1 *gener = new AliGenEpEmv1();
//    gener->SetXXXRange(); // Set kinematics range
//    gener->Init();
// Event generation:
//    gener->Generate(); // Produce one e+e- pair with the event weight assigned 
//                       // to each track. The sum of event weights, divided by 
//                       // the total number of generated events, gives the 
//                       // integral cross section of the process of e+e- pair 
//                       // production in the above mentioned kinematics range.
//                       // Sum of the selected event weights, divided by the total 
//                       // number of generated events, gives the integral cross 
//                       // section corresponded to the set of selected events
//%
// The generator consists of several modules:
// 1) $ALICE_ROOT/EpEmGen/diffcross.f:
//    Exact calculation of the total differential e+ e- -pair production
//    in Relativistic Heavy Ion Collisions for a point particle in an
//    external field approach. See full comments in the mentioned file.
// 2) $ALICE_ROOT/EpEmGen/epemgen.f:
//    Generator of e+e- pairs produced in PbPb collisions at LHC
//    it generates events according to the parametrization of the
//    differential cross section. Produces events have weights calculated
//    by the exact differential cross section calculation (diffcross.f).
//    See full comments in the mentioned file.
// 3) Class TEpEmGen:
//    Interface from the fortran event generator to ALIROOT
// 4) Class AliGenEpEmv1:
//    The event generator to call within ALIROOT
//%
// Author of this module: Yuri.Kharlov@cern.ch
// 9 October 2002

#include "AliGenEpEmv1.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TEpEmGen.h>

ClassImp(AliGenEpEmv1)

//------------------------------------------------------------

AliGenEpEmv1::AliGenEpEmv1()
{
  // Default constructor
  // Avoid zero pt
  if (fPtMin == 0) fPtMin = 1.E-04;
}

//____________________________________________________________
AliGenEpEmv1::~AliGenEpEmv1()
{
  // Destructor
}

//____________________________________________________________
void AliGenEpEmv1::Init()
{
  // Initialisation:
  // 1) define a generator
  // 2) initialize the generator of e+e- pair production

  fMass = TDatabasePDG::Instance()->GetParticle(11)->Mass();

  SetMC(new TEpEmGen());
  fEpEmGen = (TEpEmGen*) fMCEvGen;
  fEpEmGen ->Initialize(fYMin,fYMax,fPtMin,fPtMax);
  fEvent = 0;
}

//____________________________________________________________
void AliGenEpEmv1::Generate()
{
  //
  // Generate one e+e- pair
  // Gaussian smearing on the vertex is done if selected. 
  //%
  // Each produced e+e- pair is defined by the following variables:
  // rapidities of e-, e+ (yElectron,yPositron)
  // log10(pt in MeV/c) of e-, e+ (xElectron,xPositron)
  // azymuth angles between e- and e+ (phi12)
  //%
  // On output an event weight is given (weight) which is assigned to each track.
  // The sum of event weights, divided by the total number of generated events, 
  // gives the integral cross section of the e+e- pair production in the   
  // selected kinematics range.	  
  //

  Float_t polar[3]= {0,0,0};
  Float_t origin[3];
  Float_t p[3];

  Double_t ptElectron,ptPositron, phiElectron,phiPositron, mt;
  Double_t phi12=0,xElectron=0,xPositron=0,yElectron=0,yPositron=0,weight=0;
  Int_t   j, nt, id;
  Float_t random[6];

  fEpEmGen->GenerateEvent(fYMin,fYMax,fPtMin,fPtMax,
	   yElectron,yPositron,xElectron,xPositron,phi12,weight);
  if (fDebug == 1)
    printf("AliGenEpEmv1::Generate(): y=(%f,%f), x=(%f,%f), phi=%f\n",
	   yElectron,yPositron,xElectron,xPositron,phi12);

  for (j=0;j<3;j++) origin[j]=fOrigin[j];
  if(fVertexSmear==kPerEvent) {
    Rndm(random,6);
    for (j=0;j<3;j++) {
      origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
	TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
    }
  }

  Rndm(random,1);
  ptElectron  = TMath::Power(10,xElectron) * 1.e-03;;
  ptPositron  = TMath::Power(10,xPositron) * 1.e-03;;
  phiElectron = fPhiMin + random[0] * (fPhiMax-fPhiMin);
  phiPositron = phiElectron + phi12;

  // Produce electron
  mt = TMath::Sqrt(ptElectron*ptElectron + fMass*fMass);
  p[0] = ptElectron*TMath::Cos(phiElectron);
  p[1] = ptElectron*TMath::Sin(phiElectron);
  p[2] = mt*TMath::SinH(yElectron);
  id =  11;
  if (fDebug == 2)
    printf("id=%+3d, p = (%+11.4e,%+11.4e,%+11.4e) GeV\n",id,p[0],p[1],p[2]);
  PushTrack(fTrackIt,-1, id,p,origin,polar,0,kPPrimary,nt,weight);

  // Produce positron
  mt = TMath::Sqrt(ptPositron*ptPositron + fMass*fMass);
  p[0] = ptPositron*TMath::Cos(phiPositron);
  p[1] = ptPositron*TMath::Sin(phiPositron);
  p[2] = mt*TMath::SinH(yPositron);
  id = -11;
  if (fDebug == 2)
    printf("id=%+3d, p = (%+11.4e,%+11.4e,%+11.4e) GeV\n",id,p[0],p[1],p[2]);
  PushTrack(fTrackIt,-1, id,p,origin,polar,0,kPPrimary,nt,weight);
  
  fEvent++;
  if (fEvent%1000 == 0) {
    printf("=====> AliGenEpEmv1::Generate(): \n   Event %d, sigma=%f +- %f kb\n",
  	   fEvent,fEpEmGen->GetXsection(),fEpEmGen->GetDsection());
  }
}

