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

// Event generator for background from e+e- pair production in ultraperipheral PbPb collisions
// at 5.5 TeV/nucleon, integrated over specific readout cycle
// Derived from AliGenEpEmv1
// Author: ruben.shahoyan@cern.ch
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
//    AliGenQEDBg *gener = new AliGenQEDBg();
//    gener->SetXXXRange(); // Set kinematics range
//    gener->SetLumiIntTime(double lumi, double sec); // luminosity and intergration time in seconds
//    gener->Init();
// Event generation:
//    gener->Generate(); // Produce poissonian number of e+e- pair with average number
//    corresponding to requested integration time with given beam luminosity
//    Each pair has its own vertex, and time randomly distributed at originT 
//    and originT+ integration time                    
//
//    For details of pair generation see AliGenEpEmv1.cxx by Yuri.Kharlov@cern.ch

/*
  The most useful way of using it is in the overlay with other generators (using the cocktail):
  In the Config.C 
  // load the library:
  gSystem->Load("libTEPEMGEN");
  //
  // add to AliGenCocktail as:
  AliGenCocktail *cocktail = new AliGenCocktail();
  // ... setup other stuff
  // 
  // QED background
  AliGenQEDBg*  genBg = new AliGenQEDBg();
  genBg->SetEnergyCMS(5500);
  genBg->SetProjectile("A", 208, 82);
  genBg->SetTarget    ("A", 208, 82);
  genBg->SetYRange(-6.,3);
  genBg->SetPtRange(1.e-3,1.0);      // Set pt limits (GeV) for e+-: 1MeV corresponds to max R=13.3mm at 5kGaus
  genBg->SetLumiIntTime(6.e27,20e-6); // luminosity and integration time
  //
  cocktail->AddGenerator(genBg,"QEDep",1);
  // We need to generate independent vertex for every event, this must be set after adding to cocktail
  // Note that the IO origin and sigma's is transferred from AliGenCocktail to daughter generators
  genBg->SetVertexSource(kInternal);
*/

#include "AliLog.h"
#include "AliGenQEDBg.h"
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>
#include <TEpEmGen.h>

ClassImp(AliGenQEDBg)

//------------------------------------------------------------

AliGenQEDBg::AliGenQEDBg()
:  fLumi(0)
  ,fXSection(0)
  ,fXSectionEps(1e-2)
  ,fIntTime(0)
  ,fPairsInt(-1)
  ,fMinXSTest(1e3)
  ,fMaxXSTest(1e7)
{
}

//____________________________________________________________
AliGenQEDBg::~AliGenQEDBg()
{
  // Destructor
}

//____________________________________________________________
void AliGenQEDBg::Init()
{
  // Initialisation:
  AliInfo(Form("Will estimate QED bg. for L=%e cm^-2*s^-1 and Integration Time of %e s.",fLumi,fIntTime));
  if (fLumi<=0 || fIntTime<=0) {
    AliWarning("One of parameters is not set properly, no pairs will be generated");
    return;
  }
  //
  // initialize the generator of e+e- pair production
  AliGenEpEmv1::Init();
  //
  fPairsInt = 0;
  int ngen = 0;
  AliInfo(Form("Estimating x-section with min.relative precision of %f and min/max test: %d/%d",
	       fXSectionEps,int(fMinXSTest),int(fMaxXSTest)));
  //
  double yElectron,yPositron,xElectron,xPositron,phi12,weight,err=0;
  fXSection = -1;
  do {
    fEpEmGen->GenerateEvent(fYMin,fYMax,fPtMin,fPtMax,yElectron,yPositron,xElectron,xPositron,phi12,weight);
    if (++ngen>fMinXSTest) { // ensure min number of tests
      fXSection = fEpEmGen->GetXsection();
      err = fEpEmGen->GetDsection();
    }
  } while(!((fXSection>0 && err/fXSection<fXSectionEps) || ngen>fMaxXSTest));
  //
  if (fXSection<=0) {
    AliError(Form("X-section = %e after %d trials, cannot generate",fXSection,ngen));
    return;
  }
  fPairsInt = fXSection*1e-21*fLumi*fIntTime; // xsestion is in kbarn!
  AliInfo(Form("Estimated x-secion: %e+-%ekb in %d tests, <Npairs>=%e per %e time interval",
	       fXSection,err,ngen,fPairsInt,fIntTime));
  //
}

//____________________________________________________________
void AliGenQEDBg::Generate()
{
  //
  // Generate poissian <fPairsInt> e+e- pairs, each one with its vertex
  //
  Float_t polar[3]= {0,0,0};
  Float_t origin[3];
  Float_t time = 0.;
  Float_t p[3];
  //
  int npairs=0,nt,id;;
  if (fPairsInt>0) npairs=gRandom->Poisson(fPairsInt);
  AliInfoF("<nQED>=%e -> %d pairs will be generated",fPairsInt,npairs);
  if (npairs<1) return;
  //
  Double_t ptElectron,ptPositron, phiElectron,phiPositron, mt, ms2 = fMass*fMass;
  Double_t phi12=0,xElectron=0,xPositron=0,yElectron=0,yPositron=0,weight=0;
  //
  for (int i=0;i<npairs;i++) {
    // each pair has its own vertex and time
    Vertex();
    for (int j=0;j<3;j++) origin[j] = fVertex[j];
    time = fTimeOrigin+gRandom->Rndm()*fIntTime;
    //
    fEpEmGen->GenerateEvent(fYMin,fYMax,fPtMin,fPtMax,yElectron,yPositron,xElectron,xPositron,phi12,weight);
    ptElectron  = TMath::Power(10,xElectron) * 1.e-03;;
    ptPositron  = TMath::Power(10,xPositron) * 1.e-03;;
    phiElectron = fPhiMin + gRandom->Rndm() * (fPhiMax-fPhiMin);
    phiPositron = phiElectron + phi12;
    // Produce electron
    mt = TMath::Sqrt(ptElectron*ptElectron + ms2);
    p[0] = ptElectron*TMath::Cos(phiElectron);
    p[1] = ptElectron*TMath::Sin(phiElectron);
    p[2] = mt*TMath::SinH(yElectron);
    id =  11;
    PushTrack(fTrackIt,-1, id,p,origin,polar,time,kPPrimary,nt,1);    
    //
    // Produce positron
    mt = TMath::Sqrt(ptPositron*ptPositron + ms2);
    p[0] = ptPositron*TMath::Cos(phiPositron);
    p[1] = ptPositron*TMath::Sin(phiPositron);
    p[2] = mt*TMath::SinH(yPositron);
    id = -11;
    PushTrack(fTrackIt,-1, id,p,origin,polar,time,kPPrimary,nt,1);
    //
  }
  fEvent++;
  //
  fHeader.SetNProduced(2*npairs);
  fHeader.SetEventWeight(1);
  fHeader.SetInteractionTime(fTimeOrigin);
  AddHeader(&fHeader);
}

//__________________________________________________________
void AliGenQEDBg::SetLumiIntTime(double lumi, double intTime)
{
  // assign luminosity and integration time
  if (lumi<=0)    AliFatal(Form("Luminosity must be positive, in cm^-2*s^-1, %e asked",lumi));
  if (intTime<=0) AliFatal(Form("Integration time must be positive, in seconnds, %e asked",intTime));
  fLumi = lumi;
  fIntTime = intTime;
  //
}

//__________________________________________________________
void AliGenQEDBg::SetMinMaxXSTest(double mn,double mx)
{
  // set min,max number of generator calls for xsection estimates
  fMinXSTest = mn>100 ? mn : 100.;
  fMaxXSTest = mx>mx ? mx : mx+100.;
}
