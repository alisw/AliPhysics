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

// Class to create cocktails for physics with electrons, di-electrons,
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
#include "AliGenCocktailEventHeader.h"

#include "AliGenCocktailEntry.h"
#include "AliGenEMCocktail.h"
#include "AliGenEMlib.h"
#include "AliGenBox.h"
#include "AliGenParam.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliDecayer.h"
#include "AliDecayerPythia.h"
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
   fYieldArray()
{
// Constructor

}

//_________________________________________________________________________
AliGenEMCocktail::~AliGenEMCocktail()
{
// Destructor

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

// Create and add electron sources to the generator

// pizero
  AliGenParam * genpizero=0;
  Char_t namePizero[10];    
  snprintf(namePizero,10, "Pizero");    
  genpizero = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kPizero, "DUMMY");
  AddSource2Generator(namePizero,genpizero);
  TF1 *fPtPizero = genpizero->GetPt();
  fYieldArray[kGenPizero] = fPtPizero->Integral(fPtMin,fPtMax,(Double_t *) 0x0,1.e-6);
// eta  
  AliGenParam * geneta=0;
  Char_t nameEta[10];    
  snprintf(nameEta,10, "Eta");    
  geneta = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kEta, "DUMMY");
  AddSource2Generator(nameEta,geneta);
  TF1 *fPtEta = geneta->GetPt();
  fYieldArray[kGenEta] = fPtEta->Integral(fPtMin,fPtMax,(Double_t *) 0x0,1.e-6);
// rho  
  AliGenParam * genrho=0;
  Char_t nameRho[10];    
  snprintf(nameRho,10, "Rho");    
  genrho = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kRho, "DUMMY");
  AddSource2Generator(nameRho,genrho);
  TF1 *fPtRho = genrho->GetPt();
  fYieldArray[kGenRho] = fPtRho->Integral(fPtMin,fPtMax,(Double_t *) 0x0,1.e-6);
// omega
  AliGenParam * genomega=0;
  Char_t nameOmega[10];    
  snprintf(nameOmega,10, "Omega");    
  genomega = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kOmega, "DUMMY");
  AddSource2Generator(nameOmega,genomega);
  TF1 *fPtOmega = genomega->GetPt();
  fYieldArray[kGenOmega] = fPtOmega->Integral(fPtMin,fPtMax,(Double_t *) 0x0,1.e-6);
// etaprime
  AliGenParam * genetaprime=0;
  Char_t nameEtaprime[10];    
  snprintf(nameEtaprime,10, "Etaprime");    
  genetaprime = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kEtaprime, "DUMMY");
  AddSource2Generator(nameEtaprime,genetaprime);
  TF1 *fPtEtaprime = genetaprime->GetPt();
  fYieldArray[kGenEtaprime] = fPtEtaprime->Integral(fPtMin,fPtMax,(Double_t *) 0x0,1.e-6);
// phi  
  AliGenParam * genphi=0;
  Char_t namePhi[10];    
  snprintf(namePhi, 10, "Phi");    
  genphi = new AliGenParam(fNPart, new AliGenEMlib(), AliGenEMlib::kPhi, "DUMMY");
  AddSource2Generator(namePhi,genphi);
  TF1 *fPtPhi = genphi->GetPt();
  fYieldArray[kGenPhi] = fPtPhi->Integral(fPtMin,fPtMax,(Double_t *) 0x0,1.e-6);
}

//-------------------------------------------------------------------
void AliGenEMCocktail::AddSource2Generator(Char_t* nameSource, 
					 AliGenParam* const genSource)
{
// add sources to the cocktail
  Double_t phiMin = fPhiMin*180./TMath::Pi();
  Double_t phiMax = fPhiMax*180./TMath::Pi();

  genSource->SetPtRange(fPtMin, fPtMax);  
  genSource->SetYRange(fYMin, fYMax);
  genSource->SetPhiRange(phiMin, phiMax);
  genSource->SetWeighting(fWeightingMode);
  if (!gMC) genSource->SetDecayer(fDecayer);  
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
  AliGenCocktailEntry *preventry = 0;
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
  while((entry = (AliGenCocktailEntry*)next())) {
    gen = entry->Generator();
    gen->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2));
    gen->SetTime(fTime);
    
    if (fNPart > 0) {
      igen++;	
      if (igen == 1) entry->SetFirst(0);		
      else  entry->SetFirst((partArray->GetEntriesFast())+1);
      gen->SetNumberParticles(fNPart);		
      gen->Generate();
      entry->SetLast(partArray->GetEntriesFast());
      preventry = entry;
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
    }
    else
      pdgMother = part->GetPdgCode();
    switch (pdgMother){
    case 111:
      dNdy = fYieldArray[kGenPizero];
      break;
    case 221:
      dNdy = fYieldArray[kGenEta];
      break;
    case 113:
      dNdy = fYieldArray[kGenRho];
      break;
    case 223:
      dNdy = fYieldArray[kGenOmega];
      break;
    case 331:
      dNdy = fYieldArray[kGenEtaprime];
      break;
    case 333:
      dNdy = fYieldArray[kGenPhi];
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
  fHeader->SetInteractionTime(fTime);

  gAlice->SetGenEventHeader(fHeader);
}

    
