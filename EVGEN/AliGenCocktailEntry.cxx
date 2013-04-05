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

// Entry for AliGenCocktail container class. 
// See also comments there.
// In addition to the pointer to the generator the class provides some 
// bookkeeping functionality (weights, links to particle list, ...)
// Author: andreas.morsch@cern.ch

#include <TString.h>
#include <TFormula.h>
#include "AliGenCocktailEntry.h"
#include "AliGenerator.h"


ClassImp(AliGenCocktailEntry)

AliGenCocktailEntry::AliGenCocktailEntry():
    fGenerator(0),
    fNGenerated(0),
    fFirst(-1),
    fLast(-1),
    fRate(0),
    fNTimes(1),
    fKineBias(1),
    fBias(1),
    fFormula(0)
{
// Default constructor

}

AliGenCocktailEntry:: AliGenCocktailEntry
(AliGenerator* pGenerator, const char* Name, Float_t RateExp) : 
  TNamed(Name, "Generator Cocktail Entry"),
  fGenerator(pGenerator),
  fNGenerated(0),
  fFirst(-1),
  fLast(-1),
  fRate(RateExp),
  fNTimes(1),
  fKineBias(1),
  fBias(1),
  fFormula(0)
{
    // Constructor
}

AliGenCocktailEntry::AliGenCocktailEntry(const AliGenCocktailEntry &entry):
    TNamed(entry),
    fGenerator(0),
    fNGenerated(0),
    fFirst(-1),
    fLast(-1),
    fRate(0),
    fNTimes(1),
    fKineBias(1),
    fBias(1),
    fFormula(0)
{
// Dummy copy constructor
    entry.Copy(*this);
}


void AliGenCocktailEntry::PrintInfo() const
{
// Print out information about generator entry
    printf("\n Generator: %s Generated Events: %d First: %d Last: %d",
	   (const char *) fName, fGenerator->NumberParticles(), fFirst, fLast);
}
 
AliGenCocktailEntry& AliGenCocktailEntry::operator
=(const  AliGenCocktailEntry& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return (*this);
}

void AliGenCocktailEntry::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}
