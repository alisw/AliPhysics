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

#include "AliLhcProcessBT.h"
#include "AliLHC.h"
#include "AliLhcIRegion.h"
#include "AliLhcBeam.h"

ClassImp(AliLhcProcessBT)

AliLhcProcessBT::AliLhcProcessBT(AliLHC* lhc, const char* name, const char* title)
    :AliLhcProcess(lhc,name,title),
     fCrossSection(0.),
     fIRegions(0),
     fBetaMin(0.)
{
// Constructor
}



AliLhcProcessBT::AliLhcProcessBT(const AliLhcProcessBT& bt):
    AliLhcProcess(bt),
    fCrossSection(0.),
    fIRegions(0),
    fBetaMin(0.)
{
// Copy Constructor

}

AliLhcProcessBT::~AliLhcProcessBT()
{
// Destructor
}

void AliLhcProcessBT::Init()
{
  // Initialization
   printf("\n Initializing Process %s", GetName());
   printf("\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
   fIRegions = fAccelerator->IRegions();
}

void AliLhcProcessBT::Evolve(Float_t dt)
{
  printf("\n Here process %s %f:", GetName(), dt);
  TIter next(fIRegions);
  AliLhcIRegion *region;
  //
  // Loop over generators and initialize
  while((region = (AliLhcIRegion*)next())) {
    Float_t betanew = 
      region->Luminosity()/region->InitialLumi()*region->BetaStar();
    if (betanew > fBetaMin) region->SetBetaStar(betanew);
  }  
}


AliLhcProcessBT& AliLhcProcessBT::operator=(const  AliLhcProcessBT & /*rhs*/)
{
// Assignment operator
    return *this;
}


