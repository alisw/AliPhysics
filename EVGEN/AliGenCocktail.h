#ifndef ALIGENCOCKTAIL_H
#define ALIGENCOCKTAIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenerator.h"

class AliGenCocktailEntry;


class AliGenCocktail : public AliGenerator
{
 public:
    AliGenCocktail();
    AliGenCocktail(const AliGenCocktail &cocktail);
     
    virtual ~AliGenCocktail();
    virtual void Init();
    virtual void Generate();
    //
    // Add a new generator to the list
    virtual void AddGenerator
	(AliGenerator *Generator, TString Name, Float_t RateExp );
    virtual TList* Entries() {return fEntries;}
    // Iterators
    AliGenCocktailEntry*  FirstGenerator();
    AliGenCocktailEntry*  NextGenerator();
    void FirstGeneratorPair(AliGenCocktailEntry*&e1, AliGenCocktailEntry*&e2);
    void NextGeneratorPair (AliGenCocktailEntry*&e1, AliGenCocktailEntry*&e2);
    AliGenCocktail & operator=(const AliGenCocktail & rhs);
    
 protected:
    Int_t fNGenerators;   // Number of generators booked
    TList  *fEntries;     // List of Generators
    TObjLink *flnk1;      // Iterator for first generator
    TObjLink *flnk2;      // Iterator for second generator
//
    ClassDef(AliGenCocktail,1) // Particle cocktail generator a la SHAKER
};

#endif





