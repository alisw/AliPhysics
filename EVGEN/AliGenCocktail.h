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

class AliGenCocktailEntry;
class TArrayF;


class AliGenCocktail : public AliGenerator
{
 public:
    AliGenCocktail();
    AliGenCocktail(const AliGenCocktail &cocktail);
     
    virtual ~AliGenCocktail();
    virtual void Init();
    virtual void FinishRun();
    virtual void Generate();
    virtual void SetVertexSmear(VertexSmear_t smear);
    virtual void SetRandomise(Bool_t flag) {fRandom = flag;}
	    
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
    AliGenCocktail & operator=(const AliGenCocktail & rhs);

 protected:
    Int_t fNGenerators;   // Number of generators booked
    Bool_t fRandom;       // Flag to select random generator from list
    TArrayF  fProb;       // Probability of an event (if fRandom == kTRUE)
    TList  *fEntries;     // List of Generators
    TObjLink *flnk1;      // ! Iterator for first generator
    TObjLink *flnk2;      // ! Iterator for second generator
    
//
 private:
    void Copy(TObject &arun) const;
    ClassDef(AliGenCocktail,1) // Particle cocktail generator a la SHAKER
};

#endif





