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

#include "AliGenCocktail.h"
#include "AliGenCocktailEntry.h"
#include "AliCollisionGeometry.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliGenCocktail)

AliGenCocktail::AliGenCocktail()
                 :AliGenerator()
{
// Constructor
    fName = "Cocktail";
    fTitle= "Particle Generator using cocktail of generators";
    flnk1 = 0;
    flnk2 = 0;
    fNGenerators=0;
    fEntries = 0;
    
}

AliGenCocktail::AliGenCocktail(const AliGenCocktail & cocktail):
    AliGenerator(cocktail)
{
// Copy constructor
    cocktail.Copy(*this);
}

AliGenCocktail::~AliGenCocktail()
{
// Destructor
    delete fEntries;
}

void AliGenCocktail::
AddGenerator(AliGenerator *Generator, const char* Name, Float_t RateExp)
{
//
// Add a generator to the list 
// First check that list exists
    if (!fEntries) fEntries = new TList();

//
//  Forward parameters to the new generator
    if(TestBit(kPtRange) && !(Generator->TestBit(kPtRange)) && !(Generator->TestBit(kMomentumRange))) 
	Generator->SetPtRange(fPtMin,fPtMax);
    if(TestBit(kMomentumRange) && !(Generator->TestBit(kPtRange)) && !(Generator->TestBit(kMomentumRange)))
	Generator->SetMomentumRange(fPMin,fPMax);
    
    if (!(Generator->TestBit(kYRange)))    
	Generator->SetYRange(fYMin,fYMax);
    if (!(Generator->TestBit(kPhiRange)))   
	Generator->SetPhiRange(fPhiMin*180/TMath::Pi(),fPhiMax*180/TMath::Pi());
    if (!(Generator->TestBit(kThetaRange))) 
	Generator->SetThetaRange(fThetaMin*180/TMath::Pi(),fThetaMax*180/TMath::Pi());
    if (!(Generator->TestBit(kVertexRange))) 
	Generator->SetOrigin(fOrigin[0], fOrigin[1], fOrigin[2]);

    Generator->SetSigma(fOsigma[0], fOsigma[1], fOsigma[2]);
    Generator->SetVertexSmear(fVertexSmear);
    Generator->SetVertexSource(kContainer);
    Generator->SetTrackingFlag(fTrackIt);
        
//
//  Add generator to list   
    char theName[256];
    sprintf(theName, "%s_%d",Name, fNGenerators);
    Generator->SetName(theName);

    AliGenCocktailEntry *entry = 
	new AliGenCocktailEntry(Generator, Name, RateExp);
    
     fEntries->Add(entry);
     fNGenerators++;
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
	entry->Generator()->Init();
    }  
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
    AliGenerator* gen = 0;

    TObjArray *partArray = gAlice->GetMCApp()->Particles();

//
//  Generate the vertex position used by all generators
//    
    if(fVertexSmear == kPerEvent) Vertex();

    
  //
    // Loop over generators and generate events
    Int_t igen=0;
    
    while((entry = (AliGenCocktailEntry*)next())) {
	igen++;
	if (igen ==1) {
	    entry->SetFirst(0);
	} else {
	    entry->SetFirst((partArray->GetEntriesFast())+1);
	}
//
//      Handle case in which current generator needs collision geometry from previous generator
//
	gen = entry->Generator();
	if (gen->NeedsCollisionGeometry())
	{
	    if (preventry && preventry->Generator()->ProvidesCollisionGeometry())
	    {
		gen->SetCollisionGeometry(preventry->Generator()->CollisionGeometry());
	    } else {
		Fatal("Generate()", "No Collision Geometry Provided");
	    }
	}
	entry->Generator()->SetVertex(fVertex.At(0), fVertex.At(1), fVertex.At(2));
	entry->Generator()->Generate();
	entry->SetLast(partArray->GetEntriesFast());
	preventry = entry;
    }  
    next.Reset();
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

AliGenCocktail& AliGenCocktail::operator=(const  AliGenCocktail& rhs)
{
// Assignment operator
    rhs.Copy(*this); 
    return (*this);
}

void AliGenCocktail::Copy(TObject &) const
{
    Fatal("Copy","Not implemented!\n");
}



