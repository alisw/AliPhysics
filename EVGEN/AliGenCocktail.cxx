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
Revision 1.8  2000/10/27 13:53:29  morsch
AddGenerator: check testbit before setting the pT and momentum range
(D.Y. Peressounko)

Revision 1.7  2000/10/02 15:15:41  morsch
Use default streamer for AliGenCocktail

Revision 1.6  2000/09/07 17:04:31  morsch
In Streamer: TIter() after R__b << fEntries (Dmitri Yurevitch Peressounko)

Revision 1.5  2000/06/09 20:28:51  morsch
All coding rule violations except RS3 corrected (AM)

Revision 1.4  1999/09/29 09:24:12  fca
Introduction of the Copyright and cvs Log

*/

// Container class for AliGenerator through recursion.
// Container is itself an AliGenerator.
// What is stored are not the pointers to the generators directly but to objects of type
// AliGenCocktail entry.   
// The class provides also iterator functionality.  
// Author: andreas.morsch@cern.ch 
//

#include "AliGenCocktail.h"
#include "AliGenCocktailEntry.h"
#include "AliRun.h"
#include <TList.h>

ClassImp(AliGenCocktail)

AliGenCocktail::AliGenCocktail()
                 :AliGenerator()
{
// Constructor
    fEntries = new TList;
    fNGenerators=0;
}

AliGenCocktail::AliGenCocktail(const AliGenCocktail & cocktail)
{
// copy constructor
}

AliGenCocktail::~AliGenCocktail()
{
// Destructor
    delete fEntries;
}

void AliGenCocktail::
AddGenerator(AliGenerator *Generator, char* Name, Float_t RateExp)
{
//
//  Forward parameters to the new generator
    if(TestBit(kPtRange)) 
    Generator->SetPtRange(fPtMin,fPtMax);
    if(TestBit(kMomentumRange))
    Generator->SetMomentumRange(fPMin,fPMax);

    Generator->SetYRange(fYMin,fYMax);
    Generator->
	SetPhiRange(fPhiMin*180/TMath::Pi(),fPhiMax*180/TMath::Pi());
    Generator->
	SetThetaRange(fThetaMin*180/TMath::Pi(),fThetaMax*180/TMath::Pi());
    Generator->
	SetOrigin(fOrigin[0], fOrigin[1], fOrigin[2]);
    Generator->
	SetSigma(fOsigma[0], fOsigma[1], fOsigma[2]);
    Generator->SetVertexSmear(fVertexSmear);
    Generator->SetTrackingFlag(fTrackIt);    
//
//  Add generator to list   
    
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
	entry->Generator()->Init();
    }  
}

 void AliGenCocktail::Generate()
{
//
// Generate event 
    TIter next(fEntries);
    AliGenCocktailEntry *entry;
    AliGenCocktailEntry *e1;
    AliGenCocktailEntry *e2;
    TClonesArray *partArray = gAlice->Particles();
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
	entry->Generator()->Generate();
	entry->SetLast(partArray->GetEntriesFast());
    }  
    next.Reset();
    while((entry = (AliGenCocktailEntry*)next())) {
	entry->PrintInfo();
    }
    for (entry=FirstGenerator();
	 entry;
	 entry=NextGenerator()
	) {
	entry->PrintInfo();
    }
    for (FirstGeneratorPair(e1,e2);
	 (e1&&e2);
	 NextGeneratorPair(e1,e2)
	){
	printf("\n -----------------------------");
	e1->PrintInfo();
	e2->PrintInfo();
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

/*
void AliGenCocktail::Streamer(TBuffer &R__b)
{
    // Stream an object of class AliGenCocktail.

    AliGenCocktailEntry *entry;
    
    if (R__b.IsReading()) {
	Version_t R__v = R__b.ReadVersion(); if (R__v) { }
	AliGenerator::Streamer(R__b);
	R__b >> fNGenerators;
	R__b >> fEntries;
	TIter next(fEntries);
// Stream generation related information
	while((entry = (AliGenCocktailEntry*)next())) {
	    entry->Streamer(R__b);
	}  
    } else {
	R__b.WriteVersion(AliGenCocktail::IsA());
	AliGenerator::Streamer(R__b);
	R__b << fNGenerators;
	R__b << fEntries;
	TIter next(fEntries); 
// Stream generation related information
	while((entry = (AliGenCocktailEntry*)next())) {
	    entry->Streamer(R__b);
	}  
    }
}
*/

AliGenCocktail& AliGenCocktail::operator=(const  AliGenCocktail& rhs)
{
// Assignment operator
    return *this;
}


