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
Revision 1.16  2003/01/14 10:50:18  alibrary
Cleanup of STEER coding conventions

Revision 1.15  2003/01/07 14:13:22  morsch
Communication between generators provising and requesting collision
geometries.

Revision 1.14  2002/02/08 16:50:50  morsch
Add name and title in constructor.

Revision 1.13  2001/10/21 18:35:56  hristov
Several pointers were set to zero in the default constructors to avoid memory management problems

Revision 1.12  2001/06/18 13:07:30  morsch
Forward kinematic ranges to entries only if not set by user.

Revision 1.11  2001/01/30 09:23:12  hristov
Streamers removed (R.Brun)

Revision 1.10  2001/01/26 19:55:49  hristov
Major upgrade of AliRoot code

Revision 1.9  2000/12/21 16:24:06  morsch
Coding convention clean-up

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

#include <TList.h>
#include <TObjArray.h>

#include "AliGenCocktail.h"
#include "AliGenCocktailEntry.h"
#include "AliCollisionGeometry.h"
#include "AliRun.h"

ClassImp(AliGenCocktail)

AliGenCocktail::AliGenCocktail()
                 :AliGenerator()
{
// Constructor
    fName = "Cocktail";
    fTitle= "Particle Generator using cocktail of generators";
    fEntries = new TList;
    flnk1 = 0;
    flnk2 = 0;
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

    TObjArray *partArray = gAlice->Particles();
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
		       (preventry->Generator()->CollisionGeometry())->NN());
		
		gen->SetCollisionGeometry(preventry->Generator()->CollisionGeometry());
	    } else {
		Fatal("Generate()", "No Collision Geometry Provided");
	    }
	}
	
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
    return *this;
}


