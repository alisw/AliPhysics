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

#include "AliLhcProcessBB.h"
#include "AliLHC.h"
#include "AliLhcIRegion.h"
#include "AliLhcBeam.h"

ClassImp(AliLhcProcessBB)

AliLhcProcessBB::AliLhcProcessBB(AliLHC* lhc, const char* name, const char* title)
    :AliLhcProcess(lhc,name,title)
{
// Constructor
}


AliLhcProcessBB::~AliLhcProcessBB()
{
// Destructor

}

void AliLhcProcessBB::Init()
{
  // Initialization
   printf("\n Initializing Process %s", GetName());
   printf("\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^");
   fIRegions = fAccelerator->IRegions();
   fBeam1 = fAccelerator->Beam(0);
   fBeam2 = fAccelerator->Beam(1);
}

void AliLhcProcessBB::Evolve(Float_t dt)
{
  printf("\n Here process %s %f:", GetName(), dt);
  TIter next(fIRegions);
  AliLhcIRegion *region;
  //
  // Loop over generators and initialize
  while((region = (AliLhcIRegion*)next())) {
      Float_t rate = region->Luminosity()*fCrossSection;
      Float_t loss = rate*dt;
      fBeam1->RemoveParticles(loss);
      fBeam2->RemoveParticles(loss);
  }  
}


AliLhcProcessBB& AliLhcProcessBB::operator=(const  AliLhcProcessBB & /*rhs*/)
{
// Assignment operator
    return *this;
}


