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

#include "AliEMCALTriggerPatch.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliEMCALTriggerPatch) ;
/// \endcond

///
/// Default constructor
//____________
AliEMCALTriggerPatch::AliEMCALTriggerPatch() : TObject(),
fPosition(0x0),
fSum(0),
fTime(0),
fPeaks(0)
{ }

///
/// Constructor
//____________
AliEMCALTriggerPatch::AliEMCALTriggerPatch(Int_t i, Int_t j,  Int_t k, Int_t l) : TObject(),
fPosition(new TVector2(i, j)),
fSum(k),
fTime(l),
fPeaks(0)
{ }

///
/// Copy constructor
//____________________________________________________________________
AliEMCALTriggerPatch::AliEMCALTriggerPatch(const AliEMCALTriggerPatch& other) : TObject(other), 
fPosition(new TVector2(*other.fPosition)),
fSum(other.fSum),
fTime(other.fTime),
fPeaks(other.fPeaks)
{ }

///
/// Destructor
//____________
AliEMCALTriggerPatch::~AliEMCALTriggerPatch()
{	
  if (fPosition) delete fPosition;
}

///
/// Peak (add explanation)
//____________
void AliEMCALTriggerPatch::SetPeak(Int_t x, Int_t y, Int_t sizeX, Int_t sizeY)
{
  if (sizeX * sizeY > 31) AliError("32b limit exceeded!");
  
  fPeaks = (fPeaks | (1 << (y * sizeX + x)));
}

///
/// Print patch info
//____________
void AliEMCALTriggerPatch::Print(const Option_t*) const
{
  printf("]> Patch at (%2d , %2d) w/ sum %3d time %2d\n",
         (int)fPosition->X(), (int)fPosition->Y(), fSum, fTime); 
}
