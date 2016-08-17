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

/* $Id: AliGenEMCocktail.cxx 40702 2010-04-26 13:09:52Z morsch $ */

// Class to create the cocktail for physics with electrons, di-electrons,
// and photons from the decay of the following sources:
// pizero, eta, rho, omega, etaprime, phi
// Kinematic distributions of the sources are taken from AliGenEMlib.
// Decay channels can be selected via the method SetDecayMode.
// Particles can be generated flat in pT with weights according to the
// chosen pT distributions from AliGenEMlib (weighting mode: kNonAnalog),
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
#include "AliGenEMCocktail.h"
#include "AliGenEMlib.h"
#include "AliGenBox.h"
#include "AliGenParam.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliGenCorrHF.h"

ClassImp(AliGenEMCocktail)  
  
//________________________________________________________________________
AliGenEMCocktail::AliGenEMCocktail()
:AliGenCocktail(),
  fDecayer(0),
  fDecayMode(kAll),
  fWeightingMode(kNonAnalog),
  fNPart(1000),
  fYieldArray(),
  fCollisionSystem(AliGenEMlib::kpp7TeV),
  fPtSelectPi0(AliGenEMlib::kPizeroParam),
  fPtSelectEta(AliGenEMlib::kEtaParampp),
  fPtSelectOmega(AliGenEMlib::kOmegaParampp),
  fPtSelectPhi(AliGenEMlib::kPhiParampp),
  fCentrality(AliGenEMlib::kpp),
  fV2Systematic(AliGenEMlib::kNoV2Sys),
  fForceConv(kFALSE),
  fSelectedParticles(kGenHadrons)
{
  // Constructor
}

//_________________________________________________________________________
AliGenEMCocktail::~AliGenEMCocktail()
{
  // Destructor
}

//_________________________________________________________________________
void AliGenEMCocktail::SetHeaviestHadron(ParticleGenerator_t part)
{
  Int_t val=kGenPizero;
  while(val<part) val|=val<<1;

  fSelectedParticles=val;
  return;
}

//_________________________________________________________________________
void AliGenEMCocktail::CreateCocktail()
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
  AliInfo(Form("Selected Params:collision system - %d , centrality - %d, pi0 param - %d, eta param - %d, omega param - %d, phi param - %d",fCollisionSystem, fCentrality, fPtSelectPi0, fPtSelectEta, fPtSelectOmega, fPtSelectPhi));
  //Initialize user selection for Pt Parameterization and centrality:
  AliGenEMlib::SelectParams(fCollisionSystem, fPtSelectPi0, fPtSelectEta, fPtSelectOmega, fPtSelectPhi, fCentrality,fV2Systematic);

  // Create and add electron sources to the generator
  // pizero
  if(fSelectedParticles&kGenPizero){
    AliGenParam *genpizero=0;
    Char_t namePizero[10];    
    snprintf(namePizero,10,"Pizero");    
    //fNPart/0.925: increase number of particles so that we have the chosen number of particles in the chosen eta range
    // 	genpizero = new AliGenParam(fNPart/0.925, new AliGenEMlib(), AliGenEMlib::kPizero, "DUMMY");
    //fYMin/0.925: increase eta range, so that the electron yield is constant (<5% change) over the chosen eta range
    // genpizero->SetYRange(fYMin/0.925, fYMax/0.925);
		
    // NOTE Theo: fNPart/0.925: increase number of particles so that we have the chosen number of particles in the chosen eta range
    // NOTE Theo: fYMin/0.925: increase eta range, so that the electron yield is constant (<5% change) over the chosen eta range
    // NOTE Friederike: the additional factors here cannot be fixed numbers, if you need them 
    // 					generate a setting which puts them for you but never do it hardcoded - electrons are not the only ones 
    //					using the cocktail
    genpizero = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kPizero, "DUMMY");
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
    AliGenParam *geneta=0;
    Char_t nameEta[10];    
    snprintf(nameEta,10,"Eta");    
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    geneta = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kEta, "DUMMY");
    geneta->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameEta,geneta);
    TF1 *fPtEta = geneta->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kEta] = fPtEta->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kEta] = fPtEta->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }

  // rho  
  if(fSelectedParticles&kGenRho0){
    AliGenParam *genrho=0;
    Char_t nameRho[10];    
    snprintf(nameRho,10,"Rho");    
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genrho = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kRho0, "DUMMY");
    genrho->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameRho,genrho);
    TF1 *fPtRho = genrho->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kRho0] = fPtRho->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kRho0] = fPtRho->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }
	
  // omega
  if(fSelectedParticles&kGenOmega){
    AliGenParam *genomega=0;
    Char_t nameOmega[10];    
    snprintf(nameOmega,10,"Omega");    
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genomega = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kOmega, "DUMMY");
    genomega->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameOmega,genomega);
    TF1 *fPtOmega = genomega->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kOmega] = fPtOmega->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kOmega] = fPtOmega->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }

  // etaprime
  if(fSelectedParticles&kGenEtaprime){
    AliGenParam *genetaprime=0;
    Char_t nameEtaprime[10];    
    snprintf(nameEtaprime,10,"Etaprime");    
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genetaprime = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kEtaprime, "DUMMY");
    genetaprime->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameEtaprime,genetaprime);
    TF1 *fPtEtaprime = genetaprime->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kEtaprime] = fPtEtaprime->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kEtaprime] = fPtEtaprime->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }

  // phi  
  if(fSelectedParticles&kGenPhi){
    AliGenParam *genphi=0;
    Char_t namePhi[10];    
    snprintf(namePhi,10,"Phi");    
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genphi = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kPhi, "DUMMY");
    genphi->SetYRange(fYMin, fYMax);
    AddSource2Generator(namePhi,genphi);
    TF1 *fPtPhi = genphi->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kPhi] = fPtPhi->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kPhi] = fPtPhi->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }

  // jpsi  
  if(fSelectedParticles&kGenJpsi){
    AliGenParam *genjpsi=0;
    Char_t nameJpsi[10];    
    snprintf(nameJpsi,10,"Jpsi");    
    // NOTE: the additional factors are set back to one as they are not the same for photons and electrons
    genjpsi = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kJpsi, "DUMMY");
    genjpsi->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameJpsi,genjpsi);
    TF1 *fPtJpsi = genjpsi->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kJpsi] = fPtJpsi->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kJpsi] = fPtJpsi->Integral(fPtMin,fPtMax,(Double_t *)0,1.e-6);
#endif
  }

  // sigma  
  if(fSelectedParticles&kGenSigma0){
    AliGenParam * gensigma=0;
    Char_t nameSigma[10];    
    snprintf(nameSigma,10, "Sigma0");    
    gensigma = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kSigma0, "DUMMY");
    gensigma->SetYRange(fYMin, fYMax);

    AddSource2Generator(nameSigma,gensigma);
    TF1 *fPtSigma = gensigma->GetPt();
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
    fYieldArray[kSigma0] = fPtSigma->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
#else
    fYieldArray[kSigma0] = fPtSigma->Integral(fPtMin,fPtMax,1.e-6);
#endif
  }
	
  // k0short
  if(fSelectedParticles&kGenK0s){
    AliGenParam * genkzeroshort=0;
    Char_t nameK0short[10];    
    snprintf(nameK0short, 10, "K0short");    
    genkzeroshort = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kK0s, "DUMMY");
    genkzeroshort->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameK0short,genkzeroshort);
    TF1 *fPtK0short = genkzeroshort->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kK0s] = fPtK0short->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kK0s] = fPtK0short->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
#endif
  }	
	
  // Delta++
  if(fSelectedParticles&kGenDeltaPlPl){
    AliGenParam * genkdeltaPlPl=0;
    Char_t nameDeltaPlPl[10];    
    snprintf(nameDeltaPlPl, 10, "DeltaPlPl");    
    genkdeltaPlPl = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kDeltaPlPl, "DUMMY");
    genkdeltaPlPl->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameDeltaPlPl,genkdeltaPlPl);
    TF1 *fPtDeltaPlPl = genkdeltaPlPl->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDeltaPlPl] = fPtDeltaPlPl->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kDeltaPlPl] = fPtDeltaPlPl->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
#endif
  }
	
  // Delta+
  if(fSelectedParticles&kGenDeltaPl){
    AliGenParam * genkdeltaPl=0;
    Char_t nameDeltaPl[10];    
    snprintf(nameDeltaPl, 10, "DeltaPl");    
    genkdeltaPl = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kDeltaPl, "DUMMY");
    genkdeltaPl->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameDeltaPl,genkdeltaPl);
    TF1 *fPtDeltaPl = genkdeltaPl->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDeltaPl] = fPtDeltaPl->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kDeltaPl] = fPtDeltaPl->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
#endif
  }
	
  // Delta-
  if(fSelectedParticles&kGenDeltaMi){
    AliGenParam * genkdeltaMi=0;
    Char_t nameDeltaMi[10];    
    snprintf(nameDeltaMi, 10, "DeltaMi");    
    genkdeltaMi = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kDeltaMi, "DUMMY");
    genkdeltaMi->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameDeltaMi,genkdeltaMi);
    TF1 *fPtDeltaMi = genkdeltaMi->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDeltaMi] = fPtDeltaMi->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kDeltaMi] = fPtDeltaMi->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
#endif
  }
	
  // Delta0
  if(fSelectedParticles&kGenDeltaZero){	
    AliGenParam * genkdeltaZero=0;
    Char_t nameDeltaZero[10];    
    snprintf(nameDeltaZero, 10, "DeltaZero");    
    genkdeltaZero = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kDeltaZero, "DUMMY");
    genkdeltaZero->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameDeltaZero,genkdeltaZero);
    TF1 *fPtDeltaZero = genkdeltaZero->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kDeltaZero] = fPtDeltaZero->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kDeltaZero] = fPtDeltaZero->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
#endif
  }
	
  // rho+
  if(fSelectedParticles&kGenRhoPl){	
    AliGenParam * genkrhoPl=0;
    Char_t nameRhoPl[10];    
    snprintf(nameRhoPl, 10, "RhoPl");    
    genkrhoPl = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kRhoPl, "DUMMY");
    genkrhoPl->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameRhoPl,genkrhoPl);
    TF1 *fPtRhoPl = genkrhoPl->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kRhoPl] = fPtRhoPl->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kRhoPl] = fPtRhoPl->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
#endif
  }

  // rho-
  if(fSelectedParticles&kGenRhoMi){	
    AliGenParam * genkrhoMi=0;
    Char_t nameRhoMi[10];    
    snprintf(nameRhoMi, 10, "RhoMi");    
    genkrhoMi = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kRhoMi, "DUMMY");
    genkrhoMi->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameRhoMi,genkrhoMi);
    TF1 *fPtRhoMi = genkrhoMi->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kRhoMi] = fPtRhoMi->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kRhoMi] = fPtRhoMi->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
#endif
  }

  // K0*
  if(fSelectedParticles&kGenK0star){	
    AliGenParam * genkK0star=0;
    Char_t nameK0star[10];    
    snprintf(nameK0star, 10, "K0star");    
    genkK0star = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kK0star, "DUMMY");
    genkK0star->SetYRange(fYMin, fYMax);
    AddSource2Generator(nameK0star,genkK0star);
    TF1 *fPtK0star = genkK0star->GetPt();
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,99,0)
    fYieldArray[kK0star] = fPtK0star->Integral(fPtMin,fPtMax,1.e-6);
#else
    fYieldArray[kK0star] = fPtK0star->Integral(fPtMin,fPtMax,(Double_t*)0,1.e-6);
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
    genDirectRealG = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kDirectRealGamma, "DUMMY");
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
    genDirectVirtG = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kDirectVirtGamma, "DUMMY");
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
void AliGenEMCocktail::AddSource2Generator(Char_t* nameSource, 
					   AliGenParam* const genSource)
{
  // add sources to the cocktail
  Double_t phiMin = fPhiMin*180./TMath::Pi();
  Double_t phiMax = fPhiMax*180./TMath::Pi();

  genSource->SetPtRange(fPtMin, fPtMax);  
  genSource->SetPhiRange(phiMin, phiMax);
  genSource->SetWeighting(fWeightingMode);
  genSource->SetForceGammaConversion(fForceConv);
  if (!TVirtualMC::GetMC()) genSource->SetDecayer(fDecayer);  
  genSource->Init();
		
  AddGenerator(genSource,nameSource,1.); // Adding Generator    
}

//-------------------------------------------------------------------
void AliGenEMCocktail::Init()
{
  // Initialisation
  TIter next(fEntries);
  AliGenCocktailEntry *entry;
  if (fStack) {
    while((entry = (AliGenCocktailEntry*)next())) {
      entry->Generator()->SetStack(fStack);
    }
  }
}

//_________________________________________________________________________
void AliGenEMCocktail::Generate()
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
  Int_t iPart, iMother;
  Int_t pdgMother = 0;
  Double_t weight = 0.;
  Double_t dNdy = 0.;
  Int_t maxPart = partArray->GetEntriesFast();
  for(iPart=0; iPart<maxPart; iPart++){      
    TParticle *part = gAlice->GetMCApp()->Particle(iPart);
    iMother = part->GetFirstMother();
    TParticle *mother = 0;
    if (iMother>=0){
      mother = gAlice->GetMCApp()->Particle(iMother);
      pdgMother = mother->GetPdgCode();
      if(abs(part->GetPdgCode())==220011){
	// handle electrons from forced conversion
	part->SetPdgCode(TMath::Sign(abs(part->GetPdgCode())-220000,part->GetPdgCode()));
	if(pdgMother!=220000){
	  int iGrandMother = mother->GetFirstMother();
	  TParticle *grandmother = gAlice->GetMCApp()->Particle(iGrandMother);
	  pdgMother = grandmother->GetPdgCode();
	}
      }
    } else pdgMother = part->GetPdgCode();

    switch (pdgMother){
    case 111:
      dNdy = fYieldArray[kPizero];
      break;
    case 221:
      dNdy = fYieldArray[kEta];
      break;
    case 113:
      dNdy = fYieldArray[kRho0];
      break;
    case 223:
      dNdy = fYieldArray[kOmega];
      break;
    case 331:
      dNdy = fYieldArray[kEtaprime];
      break;
    case 333:
      dNdy = fYieldArray[kPhi];
      break;
    case 443:
      dNdy = fYieldArray[kJpsi];
      break;
    case 220000:
      dNdy = fYieldArray[kDirectRealGamma];
      break;
    case 220001:
      dNdy = fYieldArray[kDirectVirtGamma];
      break;
    case 3212:
      dNdy = fYieldArray[kSigma0];
      break;				
    case 310:
      dNdy = fYieldArray[kK0s];
      break;
    case 2224:
      dNdy = fYieldArray[kDeltaPlPl];
      break;
    case 2214:
      dNdy = fYieldArray[kDeltaPl];
      break;
    case 1114:
      dNdy = fYieldArray[kDeltaMi];
      break;
    case 2114:
      dNdy = fYieldArray[kDeltaZero];
      break;
    case 213:	
      dNdy = fYieldArray[kRhoPl];
      break;
    case -213:	
      dNdy = fYieldArray[kRhoMi];
      break;
    case 313:	
      dNdy = fYieldArray[kK0star];
      break;
    default:
      dNdy = 0.;
    }

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
