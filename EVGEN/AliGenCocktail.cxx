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

/* $Id$ */

// Container class for AliGenerator through recursion.
// Container is itself an AliGenerator.
// What is stored are not the pointers to the generators directly but to objects of type
// AliGenCocktail entry.   
// The class provides also iterator functionality.  
// Author: andreas.morsch@cern.ch 
//

#include <TList.h>
#include <TObjArray.h>
#include <TFormula.h>

#include "AliGenCocktail.h"
#include "AliGenCocktailEntry.h"
#include "AliCollisionGeometry.h"
#include "AliRun.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliGenCocktailEventHeader.h"

ClassImp(AliGenCocktail)

AliGenCocktail::AliGenCocktail()
    :AliGenerator(), 
     fNGenerators(0),
     fTotalRate(0.),
     fSRandom(kFALSE),
     fUsePerEventRate(kFALSE),
     fProb(0),
     fEntries(0),
     flnk1(0),
     flnk2(0), 
     fHeader(0),
     fSeed(0)
{
// Constructor
    fName = "Cocktail";
    fTitle= "Particle Generator using cocktail of generators";
}

AliGenCocktail::~AliGenCocktail()
{
// Destructor
    delete fEntries;
    fEntries = 0;
    //    delete fHeader; // It is removed in AliRunLoader
    fHeader = 0;
}

void AliGenCocktail::
AddGenerator(AliGenerator *Generator, const char* Name, Float_t RateExp, TFormula* formula, Int_t ntimes)
{
//
// Add a generator to the list 
// First check that list exists
    if (!fEntries) fEntries = new TList();
    fTotalRate += RateExp;
//
//  Forward parameters to the new generator
    if(TestBit(kPtRange) && !(Generator->TestBit(kPtRange)) && !(Generator->TestBit(kMomentumRange))) 
	Generator->SetPtRange(fPtMin,fPtMax);
    if(TestBit(kMomentumRange) && !(Generator->TestBit(kPtRange)) && !(Generator->TestBit(kMomentumRange)))
	Generator->SetMomentumRange(fPMin,fPMax);
    
    if (TestBit(kYRange) && !(Generator->TestBit(kYRange)))
	Generator->SetYRange(fYMin,fYMax);
    if (TestBit(kPhiRange) && !(Generator->TestBit(kPhiRange)))
	Generator->SetPhiRange(fPhiMin*180/TMath::Pi(),fPhiMax*180/TMath::Pi());
    if (TestBit(kThetaRange) && !(Generator->TestBit(kThetaRange)) && !(Generator->TestBit(kEtaRange)))
	Generator->SetThetaRange(fThetaMin*180/TMath::Pi(),fThetaMax*180/TMath::Pi());
    if (!(Generator->TestBit(kVertexRange))) {
	Generator->SetOrigin(fOrigin[0], fOrigin[1], fOrigin[2]);
	Generator->SetSigma(fOsigma[0], fOsigma[1], fOsigma[2]);
	Generator->SetVertexSmear(fVertexSmear);
	Generator->SetVertexSource(kContainer);
    }
    Generator->SetTrackingFlag(fTrackIt);
    Generator->SetContainer(this);

        
//
//  Add generator to list   
    char theName[256];
    snprintf(theName, 256, "%s_%d",Name, fNGenerators);
    Generator->SetName(theName);

    AliGenCocktailEntry *entry = 
	new AliGenCocktailEntry(Generator, Name, RateExp);
    if (formula) entry->SetFormula(formula);  
    entry->SetNTimes(ntimes);
     fEntries->Add(entry);
     fNGenerators++;
     flnk1 = 0;
     flnk2 = 0;
     fSRandom  = kFALSE;
     fHeader  = 0;
}

  void AliGenCocktail::Init()
{
// Initialisation
    TIter next(fEntries);
    AliGenCocktailEntry *entry;
    //
    // Loop over generators and initialize
    while((entry = (AliGenCocktailEntry*)next())) {
	if (fStack)  entry->Generator()->SetStack(fStack);
	if (fSeed)   entry->Generator()->SetSeed(fSeed);
	entry->Generator()->Init();
    }  

    next.Reset();

    if (fSRandom) {
	fProb.Set(fNGenerators);
	next.Reset();
	Float_t sum = 0.;
	while((entry = (AliGenCocktailEntry*)next())) {
	    sum += entry->Rate();
	} 

	next.Reset();
	Int_t i = 0;
	Float_t psum = 0.;
	while((entry = (AliGenCocktailEntry*)next())) {
	    psum +=  entry->Rate() / sum;
	    fProb[i++] = psum;
	}
    }
	next.Reset();
}

  void AliGenCocktail::FinishRun()
{
// Initialisation
    TIter next(fEntries);
    AliGenCocktailEntry *entry;
    //
    // Loop over generators and initialize
    while((entry = (AliGenCocktailEntry*)next())) {
	entry->Generator()->FinishRun();
    }  
}

 void AliGenCocktail::Generate()
{
//
// Generate event 
    TIter next(fEntries);
    AliGenCocktailEntry *entry = 0;
    AliGenCocktailEntry *preventry = 0;
    AliGenCocktailEntry *collentry = 0;
    AliGenerator* gen = 0;
    if (fHeader) delete fHeader;

    
    fHeader = new AliGenCocktailEventHeader("Cocktail Header");

    const TObjArray *partArray = gAlice->GetMCApp()->Particles();

//
//  Generate the vertex position used by all generators
//    
    if(fVertexSmear == kPerEvent) Vertex();

    TArrayF eventVertex;
    eventVertex.Set(3);
    for (Int_t j=0; j < 3; j++) eventVertex[j] = fVertex[j];

    if (!fSRandom) {
	//
	// Loop over generators and generate events
	Int_t igen   = 0;
	while((entry = (AliGenCocktailEntry*)next())) {
          Int_t ntimes = entry->NTimes();
	  if (fUsePerEventRate && (gRandom->Rndm() > entry->Rate())) continue;
	  
	  igen++;
	  if (igen ==1) {
	    entry->SetFirst(0);
	  } else {
	    entry->SetFirst((partArray->GetEntriesFast())+1);
	  }
	  gen = entry->Generator();
	  if (gen->ProvidesCollisionGeometry()) collentry = entry; 
	  //
	  //      Handle case in which current generator needs collision geometry from previous generator
          //
	  if (gen->NeedsCollisionGeometry() && (entry->Formula() == 0))
	    {
	      if (preventry && preventry->Generator()->ProvidesCollisionGeometry())
		{
		  gen->SetCollisionGeometry(preventry->Generator()->CollisionGeometry());
		} else {
		Fatal("Generate()", "No Collision Geometry Provided");
	      }
	    }
	  //
	  //      Number of signals is calculated from Collision Geometry
	  //      and entry with given centrality bin is selected
	  //
	  if (entry->Formula() != 0)
	    {
	      if (!collentry) {
		Fatal("Generate()", "No Collision Geometry Provided");
		return;
	      }
	      AliCollisionGeometry* coll = (collentry->Generator())->CollisionGeometry();
	      Float_t b  = coll->ImpactParameter();
	      Int_t nsig = Int_t(entry->Formula()->Eval(b));
	      Int_t bin = entry->Bin() - 100;
	      if (bin > 0) {
		if (bin != nsig) continue;
	      } else {
		if (nsig < 1) nsig = 1;
		AliInfo(Form("Signal Events %13.3f %5d %5d\n", b, coll->HardScatters(), nsig));
		ntimes = nsig;
	      }
	    }
	  gen->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2), fTime);
	  
	  gen->GenerateN(ntimes);
	  entry->SetLast(partArray->GetEntriesFast());
	  preventry = entry;
	}
    } else if (fSRandom) {
	//
	// Select a generator randomly
	//
	Int_t i;
	Float_t p0 =  gRandom->Rndm();

	for (i = 0; i < fNGenerators; i++) {
	    if (p0 < fProb[i]) break;
	}

	entry = (AliGenCocktailEntry*) fEntries->At(i);
	entry->SetFirst(0);
	gen = entry->Generator();
	gen->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2), fTime);
	gen->Generate();
	entry->SetLast(partArray->GetEntriesFast());
    } 
    
    next.Reset();

    // Event Vertex
    fHeader->SetPrimaryVertex(eventVertex);
    fHeader->CalcNProduced();
    if (fContainer) {
      fHeader->SetName(fName);
      fContainer->AddHeader(fHeader);
    } else {
      gAlice->SetGenEventHeader(fHeader);	
    }
}

void AliGenCocktail::SetVertexSmear(VertexSmear_t smear)
{
// Set vertex smearing and propagate it to the generators

  AliGenerator::SetVertexSmear(smear);
  TIter next(fEntries);
  while (AliGenCocktailEntry* entry = (AliGenCocktailEntry*)next()) {
    entry->Generator()->SetVertexSmear(smear);
  }
}

AliGenCocktailEntry *  AliGenCocktail::FirstGenerator()
{
// Iterator over generators: Initialisation
    flnk1 = fEntries->FirstLink();
    if (flnk1) {
	return (AliGenCocktailEntry*) (flnk1->GetObject());
    } else {
	return 0;
    }
}

AliGenCocktailEntry*  AliGenCocktail::NextGenerator()
{
// Iterator over generators: Increment
    flnk1 = flnk1->Next();
    if (flnk1) {
	return (AliGenCocktailEntry*) (flnk1->GetObject());
    } else {
	return 0;
    }
}

void AliGenCocktail::
FirstGeneratorPair(AliGenCocktailEntry*& e1, AliGenCocktailEntry*& e2)
{
// Iterator over generator pairs: Initialisation
    flnk2 = flnk1 = fEntries->FirstLink();
    if (flnk1) {
	e2 = e1 = (AliGenCocktailEntry*) (flnk1->GetObject());
    } else {
	e2= e1 = 0;
    }
}

void AliGenCocktail::
NextGeneratorPair(AliGenCocktailEntry*& e1, AliGenCocktailEntry*& e2)
{
// Iterator over generators: Increment
    flnk2 = flnk2->Next();
    if (flnk2) {
	e1 = (AliGenCocktailEntry*) (flnk1->GetObject());
	e2 = (AliGenCocktailEntry*) (flnk2->GetObject());	
    } else {
	flnk2 = flnk1 = flnk1->Next();
	if (flnk1) {
	    e1 = (AliGenCocktailEntry*) (flnk1->GetObject());
	    e2 = (AliGenCocktailEntry*) (flnk2->GetObject());
	} else {
	    e1=0;
	    e2=0;
	}
    }
}

void AliGenCocktail::AddHeader(AliGenEventHeader* header)
{
// Add a header to the list 
    if (fHeader) fHeader->AddHeader(header);
}




