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

#include <TObjArray.h>
#include <TMath.h>

#include "AliMUONHitMapA1.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONDigit.h"
#include "AliLog.h"

ClassImp(AliMUONHitMapA1)

AliMUONHitMapA1::AliMUONHitMapA1()
  : AliHitMap()
{
    // Default constructor
    fNpx          = 0;
    fNpy          = 0;
    fMaxIndex     = 0;
    
    fHitMap       = 0;
    fDigits       = 0;
}

//________________________________________________________________________________
AliMUONHitMapA1::AliMUONHitMapA1(Int_t idDE, AliMUONGeometrySegmentation* seg, TObjArray* dig)
  : AliHitMap()
{
// Constructor with new segmentation
    fNpx  = seg->Npx(idDE)+1;
    fNpy  = seg->Npy(idDE)+1;
    fMaxIndex=2*(fNpx+1)*2*(fNpy+1)+2*fNpy;
    
    fHitMap = new Int_t[fMaxIndex];
    fDigits =  dig;
    Clear();
}
//______________________________________________________________
AliMUONHitMapA1::AliMUONHitMapA1(const AliMUONHitMapA1 & hitMap)
  : AliHitMap(hitMap)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}
//_________________________________
AliMUONHitMapA1::~AliMUONHitMapA1()
{
// Destructor
    if (fHitMap) delete[] fHitMap;
}
//______________________________________
void AliMUONHitMapA1::Clear(const char *)
{
// Clear hitmap
    memset(fHitMap,0,sizeof(int)*fMaxIndex);
}
//___________________________________________________
Bool_t AliMUONHitMapA1::ValidateHit(Int_t ix, Int_t iy)
{
    //
    // Check if pad coordinates are within boundaries
    //
//    printf("\n Validate %d %d %d %d", ix, iy, fNpx, fNpy);
    
    return (TMath::Abs(ix) <= fNpx && TMath::Abs(iy) <= fNpy); 
}
//_________________________________________________________
Int_t AliMUONHitMapA1::CheckedIndex(Int_t ix, Int_t iy) const
{
// Return checked indices ix, iy
    Int_t index=2*fNpy*(ix+fNpx)+(iy+fNpy);
    if (index >= fMaxIndex) {
// 	printf("\n \n \n Try to read/write outside array !!!! \n \n %d %d %d %d %d %d",
// 	       ix,iy, fMaxIndex, index, fNpx, fNpy);
	return  fMaxIndex-1;
    } else {
	return index;
    }
}
//_____________________________
void  AliMUONHitMapA1::FillHits()
{
// Fill hits from digits list  
    Int_t ndigits = fDigits->GetEntriesFast();
    //printf("\n Filling hits into HitMap\n");
    //printf("FindRawClusters -- ndigits %d \n",ndigits);
    if (!ndigits) return;
    AliMUONDigit *dig;
    for (Int_t ndig=0; ndig<ndigits; ndig++) {
	dig = (AliMUONDigit*)fDigits->UncheckedAt(ndig);
	SetHit(dig->PadX(),dig->PadY(),ndig);
    }
}
//___________________________________________________________
void  AliMUONHitMapA1::SetHit(Int_t ix, Int_t iy, Int_t idigit)
{
// Assign digit to hit cell ix,iy
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=idigit+1;
    fHitMap[CheckedIndex(ix, iy)]=idigit+1;
}
//_______________________________________________
void AliMUONHitMapA1::DeleteHit(Int_t ix, Int_t iy)
{
// Delete hit at cell ix,iy
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=0;
    fHitMap[CheckedIndex(ix, iy)]=0;
}
//_____________________________________________
void AliMUONHitMapA1::FlagHit(Int_t ix, Int_t iy)
{
// Flag hit as used
    fHitMap[CheckedIndex(ix, iy)]=
	-TMath::Abs(fHitMap[CheckedIndex(ix, iy)]);
}
//________________________________________________________
Int_t AliMUONHitMapA1::GetHitIndex(Int_t ix, Int_t iy) const
{
// Get absolute value of contents of hit cell ix,iy
    return TMath::Abs(fHitMap[CheckedIndex(ix, iy)])-1;
}
//_______________________________________________________
TObject* AliMUONHitMapA1::GetHit(Int_t ix, Int_t iy) const
{
    // Get pointer to object at hit cell ix, iy
    // Force crash if index does not exist ! (Manu)
    Int_t index=GetHitIndex(ix,iy);
    return (index <0) ? 0 : fDigits->UncheckedAt(GetHitIndex(ix,iy));
}
//_________________________________________________
FlagType AliMUONHitMapA1::TestHit(Int_t ix, Int_t iy)
{
// Check if hit cell is empty, used or unused
//
    Int_t index = CheckedIndex(ix, iy);
    if (index<0 || index >= fMaxIndex) return kEmpty;

    Int_t inf=fHitMap[index];
    if (inf < 0) {
	return kUsed;
    } else if (inf == 0) {
	return kEmpty;
    } else {
	return kUnused;
    }
}
//________________________________________________________________________
AliMUONHitMapA1 & AliMUONHitMapA1::operator = (const AliMUONHitMapA1 & rhs) 
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal( "Not implemented.");
    
  return *this;  
}





