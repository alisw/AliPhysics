#ifndef ALIGENCOCKTAIL_H
#define ALIGENCOCKTAIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////
//                                                       //
//  Class to generate the particles for the MC           //
//  The base class is empty                              //
//                                                       //
///////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TF1.h"
#include "TArrayF.h"
#include "TTree.h"
#include "AliGenerator.h"
#include <TClass.h>

class AliGenCocktailEntry : public TObject
{
 protected:
    AliGenerator *fGenerator;
    Int_t fNGenerated;
    Int_t fFirst;
    Int_t fLast;
    Float_t fRate;
    Float_t fKineBias;
    Float_t fBias;
    TString fName;
 public:
    AliGenCocktailEntry()
	{
	    fGenerator =0;
	    fNGenerated=0;
	    fFirst=-1;
	    fLast=-1;
	    fRate=0;
	    fKineBias=1;
	    fBias=1;
	    fName="unknown";
	}
    
    AliGenCocktailEntry
	(AliGenerator* Generator, TString Name, Float_t RateExp)
	{
	    fGenerator=Generator;
	    fNGenerated=0;
	    fFirst=-1;
	    fLast=-1;
	    fRate=RateExp;
	    fName=Name;
// 	    
	    fKineBias=1;
	    fBias=1;

	}
    ~AliGenCocktailEntry(){;}
    AliGenerator* Generator() {return fGenerator;}
    void SetGenerator(AliGenerator* generator){fGenerator=generator;}
    void SetFirst(Int_t first){fFirst=first;}
    void SetLast (Int_t last ){fLast =last;}
    Int_t GetFirst(){return fFirst;}
    Int_t GetLast (){return fLast;}
    Float_t Rate(){return fRate;}
    void  PrintInfo();
 private:
    ClassDef(AliGenCocktailEntry,1)
};


class AliGenCocktail : public AliGenerator
{
 protected:
    //
    // Number of generators booked

    Int_t fNGenerators;
    //
    // List of Generators
    TList  *fEntries;
    // Iterators
    TObjLink *flnk1;
    TObjLink *flnk2;
 public:
    AliGenCocktail();
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
    void FirstGeneratorPair(AliGenCocktailEntry*&, AliGenCocktailEntry*&);
    void NextGeneratorPair (AliGenCocktailEntry*&, AliGenCocktailEntry*&);
    ClassDef(AliGenCocktail,1)
//
};

#endif





