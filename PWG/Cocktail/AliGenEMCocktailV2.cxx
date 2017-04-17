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

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Implementation of AliGenEMCocktailV2 for electron, di-electron,         //
// and photon cocktail calculations.                                       //
// It is based on AliGenEMCocktail                                         //
//                                                                         //
// Responsible: Friederike Bock (friederike.bock@cern.ch)                  //
//              Lucas Altenkaemper (lucas.altenkamper@cern.ch)             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

// Class to create the cocktail for physics with electrons, di-electrons,
// and photons from the decay of the following sources:
// pizero, eta, rho, omega, etaprime, phi
// Kinematic distributions of the sources are taken from AliGenEMlibV2.
// Decay channels can be selected via the method SetDecayMode.
// Particles can be generated flat in pT with weights according to the
// chosen pT distributions from AliGenEMlibV2 (weighting mode: kNonAnalog),
// or they are generated according to the pT distributions themselves
// (weighting mode: kAnalog)


#include <TObjArray.h>
#include <TParticle.h>
#include <TF1.h>
#include <TVirtualMC.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include "AliGenCocktailEventHeader.h"

#include "AliGenCocktailEntry.h"
#include "AliGenEMCocktailV2.h"
#include "AliGenEMlibV2.h"
#include "AliGenBox.h"
#include "AliGenParam.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliGenCorrHF.h"

ClassImp(AliGenEMCocktailV2)

//________________________________________________________________________
AliGenEMCocktailV2::AliGenEMCocktailV2():AliGenCocktail(),
  fDecayer(0),
  fDecayMode(kAll),
  fWeightingMode(kNonAnalog),
  fNPart(1000),
  fParametrizationFile(""),
  fParametrizationDir(""),
  fYieldArray(),
  fCollisionSystem(AliGenEMlibV2::kpp7TeV),
  fCentrality(AliGenEMlibV2::kpp),
  fV2Systematic(AliGenEMlibV2::kNoV2Sys),
  fUseYWeighting(kFALSE),
  fDynPtRange(kFALSE),
  fForceConv(kFALSE),
  fSelectedParticles(kGenHadrons)
{
  // Constructor
}

// initialize static member
TF1*  AliGenEMCocktailV2::fPtParametrization[]    = {0x0};
TF1*  AliGenEMCocktailV2::fParametrizationProton  = NULL;
TH1D* AliGenEMCocktailV2::fMtScalingFactorHisto   = NULL;
TH2F* AliGenEMCocktailV2::fPtYDistribution[]      = {0x0};

//_________________________________________________________________________
AliGenEMCocktailV2::~AliGenEMCocktailV2()
{
  // Destructor
}

//_________________________________________________________________________
void AliGenEMCocktailV2::SetHeaviestHadron(ParticleGenerator_t part)
{
  Int_t val=kGenPizero;
  while(val<part) val|=val<<1;
  
  fSelectedParticles=val;
  return;
}

//_________________________________________________________________________
TF1* AliGenEMCocktailV2::GetPtParametrization(Int_t np) {

  if (np<18)
    return fPtParametrization[np];
  else if (np==18)
    return fParametrizationProton;
  else
    return NULL;
}

//_________________________________________________________________________
TH1D* AliGenEMCocktailV2::GetMtScalingFactors() {
  return fMtScalingFactorHisto;
}

//_________________________________________________________________________
TH2F* AliGenEMCocktailV2::GetPtYDistribution(Int_t np) {

  if (np<18)
    return fPtYDistribution[np];
  else
    return NULL;
}

//_________________________________________________________________________
void AliGenEMCocktailV2::GetPtRange(Double_t &ptMin, Double_t &ptMax) {
  ptMin = fPtMin;
  ptMax = fPtMax;
}

//-------------------------------------------------------------------
Double_t AliGenEMCocktailV2::GetMaxPtStretchFactor(Int_t pdgCode) {

  Double_t massParticle = TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass();
  Double_t massPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();

  Double_t factor = massParticle/massPi0;
  if (factor*fPtMax > 200) factor = 200./fPtMax; // so far the input pt parametrizations are defined up to pt = 200 GeV/c

  return factor;
}

//-------------------------------------------------------------------
Double_t AliGenEMCocktailV2::GetYWeight(Int_t np, TParticle* part) {

  if (!fUseYWeighting) return 1.;

  Double_t weight = 0.;
  if (fPtYDistribution[np]) {
    if (part->Pt() > fPtYDistribution[np]->GetXaxis()->GetXmin() && part->Pt() < fPtYDistribution[np]->GetXaxis()->GetXmin()) {
      if (part->Y() > fPtYDistribution[np]->GetYaxis()->GetXmin() && part->Y() < fPtYDistribution[np]->GetYaxis()->GetXmin()) {
        weight = fPtYDistribution[np]->GetBinContent(fPtYDistribution[np]->GetXaxis()->FindBin(part->Pt()), fPtYDistribution[np]->GetYaxis()->FindBin(part->Y()));
        if (weight)
          return weight;
        else
          return 1.;
      } else
        return 1.;
    } else
      return 1.;
  } else
    return 1.;
}

//_________________________________________________________________________
void AliGenEMCocktailV2::CreateCocktail()
{
  // create and add sources to the cocktail
  fDecayer->SetForceDecay(fDecayMode);
  fDecayer->ForceDecay();
  
  // Set kinematic limits
  Double_t ptMin  = fPtMin;
  Double_t ptMax  = fPtMax;
  Double_t yMin   = fYMin;;
  Double_t yMax   = fYMax;;
  Double_t phiMin = fPhiMin*180./TMath::Pi();
  Double_t phiMax = fPhiMax*180./TMath::Pi();
  AliInfo(Form("Ranges pT:%4.1f : %4.1f GeV/c, y:%4.2f : %4.2f, Phi:%5.1f : %5.1f degres",ptMin,ptMax,yMin,yMax,phiMin,phiMax));
  AliInfo(Form("the parametrised sources uses the decay mode %d",fDecayMode));
  AliInfo(Form("generating %d particles per source",fNPart));
  AliInfo(Form("Selected Params:collision system - %d , centrality - %d",fCollisionSystem, fCentrality));
  //Initialize user selection for Pt Parameterization and centrality:
  AliGenEMlibV2::SelectParams(fCollisionSystem, fCentrality,fV2Systematic);
  AliGenEMlibV2::SetMtScalingFactors(fParametrizationFile, fParametrizationDir);
  SetMtScalingFactors();
  AliGenEMlibV2::SetPtParametrizations(fParametrizationFile, fParametrizationDir);
  SetPtParametrizations();

  if (fDynPtRange)
    AliInfo("Dynamical adaption of pT range was chosen, the number of generated particles will also be adapted");

  if (fUseYWeighting) {
    AliInfo("Rapidity weighting will be used");
    AliGenEMlibV2::SetPtYDistributions(fParametrizationFile, fParametrizationDir);
    SetPtYDistributions();
  }

  // Create and add electron sources to the generator
  // pizero
  if(fSelectedParticles&kGenPizero){
    AliGenParam *genpizero=0;
    Char_t namePizero[10];
    snprintf(namePizero,10,"Pizero");
    //fNPart/0.925: increase number of particles so that we have the chosen number of particles in the chosen eta range
    // 	genpizero = new AliGenParam(fNPart/0.925, new AliGenEMlibV2(), AliGenEMlibV2::kPizero, "DUMMY");
    //fYMin/0.925: increase eta range, so that the electron yield is constant (<5% change) over the chosen eta range
    // genpizero->SetYRange(fYMin/0.925, fYMax/0.925);
    
    // NOTE Theo: fNPart/0.925: increase number of particles so that we have the chosen number of particles in the chosen eta range
    // NOTE Theo: fYMin/0.925: increase eta range, so that the electron yield is constant (<5% change) over the chosen eta range
    // NOTE Friederike: the additional factors here cannot be fixed numbers, if you need them
    // 					generate a setting which puts them for you but never do it hardcoded - electrons are not the only ones
    //					using the cocktail
    genpizero = new AliGenParam(fNPart, new AliGenEMlibV2(), AliGenEMlibV2::kPizero, "DUMMY");
    genpizero->SetYRange(fYMin, fYMax);

    AddSource2Generator(namePizero,genpizero);
    TF1 *fPtPizero = genpizero->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kPizero] = fPtPizero->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kPizero] = fPtPizero->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }
  
  // eta
  if(fSelectedParticles&kGenEta){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(221);

    AliGenParam *geneta=0;
    Char_t nameEta[10];
    snprintf(nameEta,10,"Eta");
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    geneta = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kEta, "DUMMY");
    geneta->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameEta,geneta,maxPtStretchFactor);
    TF1 *fPtEta = geneta->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kEta] = fPtEta->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kEta] = fPtEta->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t *)0,1.e-6);
#endif
  }
  
  // rho
  if(fSelectedParticles&kGenRho0){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(113);

    AliGenParam *genrho=0;
    Char_t nameRho[10];
    snprintf(nameRho,10,"Rho");
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genrho = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kRho0, "DUMMY");
    genrho->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameRho,genrho,maxPtStretchFactor);
    TF1 *fPtRho = genrho->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kRho0] = fPtRho->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kRho0] = fPtRho->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t *)0,1.e-6);
#endif
  }
  
  // omega
  if(fSelectedParticles&kGenOmega){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(223);

    AliGenParam *genomega=0;
    Char_t nameOmega[10];
    snprintf(nameOmega,10,"Omega");
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genomega = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kOmega, "DUMMY");
    genomega->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameOmega,genomega,maxPtStretchFactor);
    TF1 *fPtOmega = genomega->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kOmega] = fPtOmega->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kOmega] = fPtOmega->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t *)0,1.e-6);
#endif
  }
  
  // etaprime
  if(fSelectedParticles&kGenEtaprime){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(331);

    AliGenParam *genetaprime=0;
    Char_t nameEtaprime[10];
    snprintf(nameEtaprime,10,"Etaprime");
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genetaprime = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kEtaprime, "DUMMY");
    genetaprime->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameEtaprime,genetaprime,maxPtStretchFactor);
    TF1 *fPtEtaprime = genetaprime->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kEtaprime] = fPtEtaprime->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kEtaprime] = fPtEtaprime->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t *)0,1.e-6);
#endif
  }
  
  // phi
  if(fSelectedParticles&kGenPhi){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(333);

    AliGenParam *genphi=0;
    Char_t namePhi[10];
    snprintf(namePhi,10,"Phi");
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genphi = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kPhi, "DUMMY");
    genphi->SetYRange(fYMin, fYMax);

    AddSource2Generator(namePhi,genphi,maxPtStretchFactor);
    TF1 *fPtPhi = genphi->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kPhi] = fPtPhi->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kPhi] = fPtPhi->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t *)0,1.e-6);
#endif
  }
  
  // jpsi
  if(fSelectedParticles&kGenJpsi){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(443);

    AliGenParam *genjpsi=0;
    Char_t nameJpsi[10];
    snprintf(nameJpsi,10,"Jpsi");
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genjpsi = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kJpsi, "DUMMY");
    genjpsi->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameJpsi,genjpsi,maxPtStretchFactor);
    TF1 *fPtJpsi = genjpsi->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kJpsi] = fPtJpsi->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kJpsi] = fPtJpsi->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t *)0,1.e-6);
#endif
  }
  
  // sigma
  if(fSelectedParticles&kGenSigma0){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(3212);

    AliGenParam * gensigma=0;
    Char_t nameSigma[10];
    snprintf(nameSigma,10, "Sigma0");
    gensigma = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kSigma0, "DUMMY");
    gensigma->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameSigma,gensigma,maxPtStretchFactor);
    TF1 *fPtSigma = gensigma->GetPt();
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
    fYieldArray[kSigma0] = fPtSigma->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#else
    fYieldArray[kSigma0] = fPtSigma->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#endif
  }
  
  // k0short
  if(fSelectedParticles&kGenK0s){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(310);

    AliGenParam * genkzeroshort=0;
    Char_t nameK0short[10];
    snprintf(nameK0short, 10, "K0short");
    genkzeroshort = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kK0s, "DUMMY");
    genkzeroshort->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameK0short,genkzeroshort,maxPtStretchFactor);
    TF1 *fPtK0short = genkzeroshort->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kK0s] = fPtK0short->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kK0s] = fPtK0short->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }

  // k0long
  if(fSelectedParticles&kGenK0l){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(130);

    AliGenParam * genkzerolong=0;
    Char_t nameK0long[10];
    snprintf(nameK0long, 10, "K0long");
    genkzerolong = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kK0l, "DUMMY");
    genkzerolong->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameK0long,genkzerolong,maxPtStretchFactor);
    TF1 *fPtK0long = genkzerolong->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kK0l] = fPtK0long->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kK0l] = fPtK0long->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }
  
  // Lambda
  if(fSelectedParticles&kGenLambda){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(3122);

    AliGenParam * genLambda=0;
    Char_t nameLambda[10];
    snprintf(nameLambda, 10, "Lambda");
    genLambda = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kLambda, "DUMMY");
    genLambda->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameLambda,genLambda,maxPtStretchFactor);
    TF1 *fPtLambda = genLambda->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kLambda] = fPtLambda->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kLambda] = fPtLambda->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }

  // Delta++
  if(fSelectedParticles&kGenDeltaPlPl){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(2224);

    AliGenParam * genkdeltaPlPl=0;
    Char_t nameDeltaPlPl[10];
    snprintf(nameDeltaPlPl, 10, "DeltaPlPl");
    genkdeltaPlPl = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kDeltaPlPl, "DUMMY");
    genkdeltaPlPl->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameDeltaPlPl,genkdeltaPlPl,maxPtStretchFactor);
    TF1 *fPtDeltaPlPl = genkdeltaPlPl->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDeltaPlPl] = fPtDeltaPlPl->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kDeltaPlPl] = fPtDeltaPlPl->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }
  
  // Delta+
  if(fSelectedParticles&kGenDeltaPl){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(2214);

    AliGenParam * genkdeltaPl=0;
    Char_t nameDeltaPl[10];
    snprintf(nameDeltaPl, 10, "DeltaPl");
    genkdeltaPl = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kDeltaPl, "DUMMY");
    genkdeltaPl->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameDeltaPl,genkdeltaPl,maxPtStretchFactor);
    TF1 *fPtDeltaPl = genkdeltaPl->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDeltaPl] = fPtDeltaPl->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kDeltaPl] = fPtDeltaPl->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }
  
  // Delta-
  if(fSelectedParticles&kGenDeltaMi){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(1114);

    AliGenParam * genkdeltaMi=0;
    Char_t nameDeltaMi[10];
    snprintf(nameDeltaMi, 10, "DeltaMi");
    genkdeltaMi = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kDeltaMi, "DUMMY");
    genkdeltaMi->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameDeltaMi,genkdeltaMi,maxPtStretchFactor);
    TF1 *fPtDeltaMi = genkdeltaMi->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDeltaMi] = fPtDeltaMi->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kDeltaMi] = fPtDeltaMi->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }
  
  // Delta0
  if(fSelectedParticles&kGenDeltaZero){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(2114);

    AliGenParam * genkdeltaZero=0;
    Char_t nameDeltaZero[10];
    snprintf(nameDeltaZero, 10, "DeltaZero");
    genkdeltaZero = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kDeltaZero, "DUMMY");
    genkdeltaZero->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameDeltaZero,genkdeltaZero,maxPtStretchFactor);
    TF1 *fPtDeltaZero = genkdeltaZero->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDeltaZero] = fPtDeltaZero->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kDeltaZero] = fPtDeltaZero->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }
  
  // rho+
  if(fSelectedParticles&kGenRhoPl){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(213);

    AliGenParam * genkrhoPl=0;
    Char_t nameRhoPl[10];
    snprintf(nameRhoPl, 10, "RhoPl");
    genkrhoPl = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kRhoPl, "DUMMY");
    genkrhoPl->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameRhoPl,genkrhoPl,maxPtStretchFactor);
    TF1 *fPtRhoPl = genkrhoPl->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kRhoPl] = fPtRhoPl->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kRhoPl] = fPtRhoPl->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }
  
  // rho-
  if(fSelectedParticles&kGenRhoMi){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(213);

    AliGenParam * genkrhoMi=0;
    Char_t nameRhoMi[10];
    snprintf(nameRhoMi, 10, "RhoMi");
    genkrhoMi = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kRhoMi, "DUMMY");
    genkrhoMi->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameRhoMi,genkrhoMi,maxPtStretchFactor);
    TF1 *fPtRhoMi = genkrhoMi->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kRhoMi] = fPtRhoMi->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kRhoMi] = fPtRhoMi->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }
  
  // K0*
  if(fSelectedParticles&kGenK0star){
    Double_t maxPtStretchFactor = 1.;
    if (fDynPtRange) maxPtStretchFactor = GetMaxPtStretchFactor(313);

    AliGenParam * genkK0star=0;
    Char_t nameK0star[10];
    snprintf(nameK0star, 10, "K0star");
    genkK0star = new AliGenParam((Int_t)(maxPtStretchFactor*fNPart), new AliGenEMlibV2(), AliGenEMlibV2::kK0star, "DUMMY");
    genkK0star->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameK0star,genkK0star,maxPtStretchFactor);
    TF1 *fPtK0star = genkK0star->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kK0star] = fPtK0star->Integral(fPtMin,maxPtStretchFactor*fPtMax,1.e-6);
#else
    fYieldArray[kK0star] = fPtK0star->Integral(fPtMin,maxPtStretchFactor*fPtMax,(Double_t*)0,1.e-6);
#endif
  }
  
  TParticlePDG *elPDG=TDatabasePDG::Instance()->GetParticle(11);
  TDatabasePDG::Instance()->AddParticle("ForcedConversionElecton-","ForcedConversionElecton-",elPDG->Mass(),true,0,elPDG->Charge(),elPDG->ParticleClass(),220011,0);
  TDatabasePDG::Instance()->AddParticle("ForcedConversionElecton+","ForcedConversionElecton+",elPDG->Mass(),true,0,-elPDG->Charge(),elPDG->ParticleClass(),-220011,0);
  
  // direct gamma
  if(fDecayMode!=kGammaEM) return;
  
  TParticlePDG *gammaPDG=TDatabasePDG::Instance()->GetParticle(22);
  
  if(fSelectedParticles&kGenDirectRealGamma){
    TDatabasePDG::Instance()->AddParticle("DirectRealGamma","DirectRealGamma",0,true,0,0,gammaPDG->ParticleClass(),220000);
    AliGenParam *genDirectRealG=0;
    Char_t nameDirectRealG[10];
    snprintf(nameDirectRealG,10,"DirectRealGamma");
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genDirectRealG = new AliGenParam(fNPart, new AliGenEMlibV2(), AliGenEMlibV2::kDirectRealGamma, "DUMMY");
    genDirectRealG->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameDirectRealG,genDirectRealG);
    TF1 *fPtDirectRealG = genDirectRealG->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDirectRealGamma] = fPtDirectRealG->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kDirectRealGamma] = fPtDirectRealG->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }
  
  if(fSelectedParticles&kGenDirectVirtGamma){
    TDatabasePDG::Instance()->AddParticle("DirectVirtGamma","DirectVirtGamma",0,true,0,0,gammaPDG->ParticleClass(),220001);
    AliGenParam *genDirectVirtG=0;
    Char_t nameDirectVirtG[10];
    snprintf(nameDirectVirtG,10,"DirectVirtGamma");
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genDirectVirtG = new AliGenParam(fNPart, new AliGenEMlibV2(), AliGenEMlibV2::kDirectVirtGamma, "DUMMY");
    genDirectVirtG->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameDirectVirtG,genDirectVirtG);
    TF1 *fPtDirectVirtG = genDirectVirtG->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDirectVirtGamma] = fPtDirectVirtG->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kDirectVirtGamma] = fPtDirectVirtG->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }
}

//-------------------------------------------------------------------
void AliGenEMCocktailV2::AddSource2Generator(Char_t* nameSource,
                                             AliGenParam* const genSource,
                                             Double_t maxPtStretchFactor)
{
  // add sources to the cocktail
  Double_t phiMin = fPhiMin*180./TMath::Pi();
  Double_t phiMax = fPhiMax*180./TMath::Pi();
  
  genSource->SetPtRange(fPtMin, maxPtStretchFactor*fPtMax);
  genSource->SetPhiRange(phiMin, phiMax);
  genSource->SetWeighting(fWeightingMode);
  genSource->SetForceGammaConversion(fForceConv);
  if (!TVirtualMC::GetMC()) genSource->SetDecayer(fDecayer);
  genSource->Init();
		
  AddGenerator(genSource,nameSource,1.); // Adding Generator
}

//-------------------------------------------------------------------
void AliGenEMCocktailV2::Init()
{
  // Initialisation
  TIter next(fEntries);
  AliGenCocktailEntry *entry;
  if (fStack) {
    while((entry = (AliGenCocktailEntry*)next())) {
      entry->Generator()->SetStack(fStack);
      ((AliGenParam*)entry->Generator())->SetDecayer(fDecayer);
      ((AliGenParam*)entry->Generator())->SetParamsExplicitly(new AliGenEMlibV2(), ((AliGenParam*)entry->Generator())->GetParam(), "DUMMY");
    }
  }
}

//_________________________________________________________________________
Bool_t AliGenEMCocktailV2::SetPtParametrizations() {

  TF1* tempFct = NULL;
  for(Int_t i=0; i<19; i++) {
    tempFct = AliGenEMlibV2::GetPtParametrization(i);
    if (!tempFct) return kFALSE;
    if (i<18)
      fPtParametrization[i] = new TF1(*tempFct);
    else
      fParametrizationProton = new TF1(*tempFct);
  }
  return kTRUE;
}

//_________________________________________________________________________
void AliGenEMCocktailV2::SetMtScalingFactors() {
  
  TH1D* tempMtFactorHisto = AliGenEMlibV2::GetMtScalingFactors();
  fMtScalingFactorHisto = new TH1D(*tempMtFactorHisto);
}

//_________________________________________________________________________
Bool_t AliGenEMCocktailV2::SetPtYDistributions() {

  TH2F* tempPtY = NULL;
  for(Int_t i=0; i<18; i++) {
    tempPtY = AliGenEMlibV2::GetPtYDistribution(i);
    if (tempPtY)
      fPtYDistribution[i] = new TH2F(*tempPtY);
    else
      fPtYDistribution[i] = NULL;
  }

  return kTRUE;
}

//_________________________________________________________________________
void AliGenEMCocktailV2::Generate()
{
  // Generate event
  TIter next(fEntries);
  AliGenCocktailEntry *entry = 0;
  AliGenerator* gen = 0;
  
  if (fHeader) delete fHeader;
  fHeader = new AliGenCocktailEventHeader("Electromagnetic Cocktail Header");
  
  const TObjArray *partArray = gAlice->GetMCApp()->Particles();
		
  // Generate the vertex position used by all generators
  if(fVertexSmear == kPerEvent) Vertex();
  
  //Reseting stack
  AliRunLoader * runloader = AliRunLoader::Instance();
  if (runloader)
    if (runloader->Stack())
      runloader->Stack()->Clean();
  
  // Loop over generators and generate events
  Int_t igen = 0;
  Float_t evPlane;
  Rndm(&evPlane,1);
  evPlane*=TMath::Pi()*2;
  while((entry = (AliGenCocktailEntry*)next())) {
    gen = entry->Generator();
    gen->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2));
    
    if (fNPart > 0) {
      igen++;
      if (igen == 1) entry->SetFirst(0);
      else  entry->SetFirst((partArray->GetEntriesFast())+1);
      gen->SetEventPlane(evPlane);
      gen->Generate();
      entry->SetLast(partArray->GetEntriesFast());
    }
  }
  next.Reset();
  
  // Setting weights for proper absolute normalization
  Int_t iPart, iMother, iGrandMother;
  Int_t pdgMother = 0;
  Int_t pdgGrandMother = 0;
  Double_t weight = 0.;
  Double_t dNdy = 0.;
  Double_t yWeight = 0.;
  Int_t maxPart = partArray->GetEntriesFast();
  for(iPart=0; iPart<maxPart; iPart++){
    TParticle *part = gAlice->GetMCApp()->Particle(iPart);
    iMother = part->GetFirstMother();
    TParticle *mother = 0;
    TParticle *grandmother = 0;
    if (iMother>=0){
      mother = gAlice->GetMCApp()->Particle(iMother);
      pdgMother = mother->GetPdgCode();
      iGrandMother = mother->GetFirstMother();
      if (iGrandMother>=0){
        grandmother = gAlice->GetMCApp()->Particle(iGrandMother);
        pdgGrandMother = grandmother->GetPdgCode();
      }
      if(abs(part->GetPdgCode())==220011){
        // handle electrons from forced conversion
        part->SetPdgCode(TMath::Sign(abs(part->GetPdgCode())-220000,part->GetPdgCode()));
        if(pdgMother!=220000){
          int iGrandMother = mother->GetFirstMother();
          TParticle *grandmother = gAlice->GetMCApp()->Particle(iGrandMother);
          pdgMother = grandmother->GetPdgCode();
        }
      } else if (part->GetPdgCode()==22 && iMother>=0 && iGrandMother>=0){
          pdgMother = grandmother->GetPdgCode();
      }
    } else pdgMother = part->GetPdgCode();

    switch (pdgMother){
      case 111:
        dNdy = fYieldArray[kPizero];
        yWeight = GetYWeight(kPizero, part);
        break;
      case 221:
        dNdy = fYieldArray[kEta];
        yWeight = GetYWeight(kEta, part);
        break;
      case 113:
        dNdy = fYieldArray[kRho0];
        yWeight = GetYWeight(kRho0, part);
        break;
      case 223:
        dNdy = fYieldArray[kOmega];
        yWeight = GetYWeight(kOmega, part);
        break;
      case 331:
        dNdy = fYieldArray[kEtaprime];
        yWeight = GetYWeight(kEtaprime, part);
        break;
      case 333:
        dNdy = fYieldArray[kPhi];
        yWeight = GetYWeight(kPhi, part);
        break;
      case 443:
        dNdy = fYieldArray[kJpsi];
        yWeight = GetYWeight(kJpsi, part);
        break;
      case 220000:
        dNdy = fYieldArray[kDirectRealGamma];
        yWeight = 0.;
        break;
      case 220001:
        dNdy = fYieldArray[kDirectVirtGamma];
        yWeight = 0.;
        break;
      case 3212:
        dNdy = fYieldArray[kSigma0];
        yWeight = GetYWeight(kSigma0, part);
        break;
      case 310:
        dNdy = fYieldArray[kK0s];
        yWeight = GetYWeight(kK0s, part);
        break;
      case 130:
        dNdy = fYieldArray[kK0l];
        yWeight = GetYWeight(kK0l, part);
        break;
      case 3122:
        dNdy = fYieldArray[kLambda];
        yWeight = GetYWeight(kLambda, part);
        break;
      case 2224:
        dNdy = fYieldArray[kDeltaPlPl];
        yWeight = GetYWeight(kDeltaPlPl, part);
        break;
      case 2214:
        dNdy = fYieldArray[kDeltaPl];
        yWeight = GetYWeight(kDeltaPl, part);
        break;
      case 1114:
        dNdy = fYieldArray[kDeltaMi];
        yWeight = GetYWeight(kDeltaMi, part);
        break;
      case 2114:
        dNdy = fYieldArray[kDeltaZero];
        yWeight = GetYWeight(kDeltaZero, part);
        break;
      case 213:
        dNdy = fYieldArray[kRhoPl];
        yWeight = GetYWeight(kRhoPl, part);
        break;
      case -213:
        dNdy = fYieldArray[kRhoMi];
        yWeight = GetYWeight(kRhoMi, part);
        break;
      case 313:
        dNdy = fYieldArray[kK0star];
        yWeight = GetYWeight(kK0star, part);
        break;
      default:
        dNdy = 0.;
        yWeight = 0.;
    }
    
    if (fUseYWeighting && yWeight)
      weight = yWeight*dNdy*part->GetWeight();
    else
      weight = dNdy*part->GetWeight();
    part->SetWeight(weight);
  }	
  
  fHeader->SetNProduced(maxPart);

  TArrayF eventVertex;
  eventVertex.Set(3);
  for (Int_t j=0; j < 3; j++) eventVertex[j] = fVertex[j];
  
  fHeader->SetPrimaryVertex(eventVertex);
  
  gAlice->SetGenEventHeader(fHeader);
}
