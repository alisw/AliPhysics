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
$Log $
*/

#include "AliGenCocktailEntry.h"
#include "AliGenerator.h"


ClassImp(AliGenCocktailEntry)


    AliGenCocktailEntry::AliGenCocktailEntry()
{
// Default constructor
    fGenerator =0;
    fNGenerated=0;
    fFirst=-1;
    fLast=-1;
    fRate=0;
    fKineBias=1;
    fBias=1;
}

AliGenCocktailEntry:: AliGenCocktailEntry
(AliGenerator* Generator, char * Name, Float_t RateExp):TNamed(Name, "Generator Cocktail Entry")
{
// Constructor using generator type, name and rate per event
    fGenerator=Generator;
    fNGenerated=0;
    fFirst=-1;
    fLast=-1;
    fRate=RateExp;
// 	    
    fKineBias=1;
    fBias=1;
}

AliGenCocktailEntry::AliGenCocktailEntry(const AliGenCocktailEntry &entry)
{
// Dummy copy constructor
}


void AliGenCocktailEntry::PrintInfo()
{
// Print out information about generator entry
    printf("\n Generator: %s Generated Events: %d First: %d Last: %d",
	   (const char *) fName, fGenerator->NumberParticles(), fFirst, fLast);
}
 
AliGenCocktailEntry& AliGenCocktailEntry::operator
=(const  AliGenCocktailEntry& rhs)
{
// Assignment operator
    return *this;
}
