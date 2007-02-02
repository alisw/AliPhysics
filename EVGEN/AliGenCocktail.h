#ifndef ALIGENCOCKTAIL_H
#define ALIGENCOCKTAIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Container class for AliGenerator through recursion.
// (Container is itself an AliGenerator)
// Author: andreas.morsch@cern.ch 
//

#include "AliGenerator.h"
#include <TArrayF.h>
#include <TList.h>

class AliGenCocktailEntry;
class AliGenCocktailEventHeader;
class TArrayF;


class AliGenCocktail : public AliGenerator
{
 public:
    AliGenCocktail();
     
    virtual ~AliGenCocktail();
    virtual void Init();
    virtual void FinishRun();
    virtual void Generate();
    virtual void SetVertexSmear(VertexSmear_t smear);
    virtual void SetRandomise(Bool_t flag) {fRandom = flag;}
    virtual void UsePerEventRates() {fUsePerEventRate  = kTRUE;}
	    
    //
    // Add a new generator to the list
    virtual void AddGenerator
	(AliGenerator *Generator, const char* Name, Float_t RateExp );
    virtual TList* Entries() {return fEntries;}
    // Iterators
    AliGenCocktailEntry*  FirstGenerator();
    AliGenCocktailEntry*  NextGenerator();
    void FirstGeneratorPair(AliGenCocktailEntry*&e1, AliGenCocktailEntry*&e2);
    void NextGeneratorPair (AliGenCocktailEntry*&e1, AliGenCocktailEntry*&e2);
    virtual void AddHeader(AliGenEventHeader* header);
	    
 protected:
    Int_t fNGenerators;                 // Number of generators booked
    Bool_t fRandom;                     // Flag to select random generator from list
    Bool_t fUsePerEventRate;            // Flag to generate the events according to the rate per event    
    TArrayF  fProb;                     // Probability of an event (if fRandom == kTRUE)
    TList  *fEntries;                   // List of Generators
    TObjLink *flnk1;                    // ! Iterator for first generator
    TObjLink *flnk2;                    // ! Iterator for second generator
    AliGenCocktailEventHeader* fHeader; // !Header container  
			   
//
 private:
    AliGenCocktail(const AliGenCocktail &cocktail);
    AliGenCocktail & operator=(const AliGenCocktail & rhs);

    ClassDef(AliGenCocktail,1) // Particle cocktail generator a la SHAKER
};

#endif





