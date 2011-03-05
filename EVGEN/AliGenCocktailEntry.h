#ifndef ALIGENCOCKTAILENTRY_H
#define ALIGENCOCKTAILENTRY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Entry for AleGenCocktail container class. 
// See also comments there
// Author: andreas.morsch@cern.ch

#include <TNamed.h>

class AliGenerator;


class AliGenCocktailEntry : public TNamed
{
 public:
    AliGenCocktailEntry();
    AliGenCocktailEntry(const AliGenCocktailEntry &entry);
    AliGenCocktailEntry
	(AliGenerator* pGenerator, const char* Name, Float_t RateExp);
    ~AliGenCocktailEntry(){;}
    AliGenerator* Generator() {return fGenerator;}
    void SetGenerator(AliGenerator* generator){fGenerator=generator;}
    void SetFirst(Int_t first){fFirst=first;}
    void SetLast (Int_t last ){fLast =last;}
    Int_t GetFirst() const {return fFirst;}
    Int_t GetLast () const {return fLast;}
    Float_t Rate()   const {return fRate;}
    void  PrintInfo() const;
    AliGenCocktailEntry & operator =(const AliGenCocktailEntry & rhs);
 protected:
    AliGenerator *fGenerator;   // Pointer to generator
    Int_t fNGenerated;          // Number of primaries generated
    Int_t fFirst;               // First index in list of primaries
    Int_t fLast;                // Last index in list of primaries
    Float_t fRate;              // Rate per event
    Float_t fKineBias;          // Bias due to kinematic selecion
    Float_t fBias;              // Bias
    void Copy(TObject&) const;
 private:
    ClassDef(AliGenCocktailEntry,1) // Generator entry of AliGenCocktail
};
#endif





