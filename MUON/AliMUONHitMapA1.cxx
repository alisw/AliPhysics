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
$Log$
Revision 1.1.2.2  2000/06/12 07:58:06  morsch
include TMath.h

Revision 1.1.2.1  2000/06/09 22:01:09  morsch
Code from AliMUONHitMap.h
Most coding rule violations corrected.

*/

#include "AliMUONHitMapA1.h"
#include "AliMUONSegmentation.h"
#include "AliMUONResponse.h"
#include "AliMUONDigit.h"

#include <TObjArray.h>
#include <TMath.h>

ClassImp(AliMUONHitMapA1)


AliMUONHitMapA1::AliMUONHitMapA1(AliMUONSegmentation *seg, TObjArray *dig)
{
// Constructor
    fSegmentation = seg;
    fNpx  = fSegmentation->Npx();
    fNpy  = fSegmentation->Npy();
    fMaxIndex=2*(fNpx+1)*2*(fNpy+1)+2*fNpy;
    
    fHitMap = new Int_t[fMaxIndex];
    fDigits =  dig;
    fNdigits = fDigits->GetEntriesFast();
    Clear();
}

AliMUONHitMapA1::AliMUONHitMapA1(const AliMUONHitMapA1 & hitMap)
{
// Dummy copy constructor
    ;
}

 
AliMUONHitMapA1::~AliMUONHitMapA1()
{
// Destructor
//    if (fDigits) delete   fDigits;
    if (fHitMap) delete[] fHitMap;
}

void AliMUONHitMapA1::Clear()
{
// Clear hitmap
    memset(fHitMap,0,sizeof(int)*fMaxIndex);
}

Int_t AliMUONHitMapA1::CheckedIndex(Int_t ix, Int_t iy)
{
// Return checked indices ix, iy
    Int_t index=2*fNpy*(ix+fNpx)+(iy+fNpy);
    if (index > fMaxIndex) {
	printf("\n \n \n Try to read/write outside array !!!! \n \n %d %d %d %d %d %d",ix,iy, fMaxIndex, index, fNpx, fNpy);
	return  fMaxIndex-1;
    } else {
	return index;
    }
}

	
void  AliMUONHitMapA1::FillHits()
{
// Fill hits from digits list  
    Int_t ndigits = fDigits->GetEntriesFast();
    //printf("\n Filling hits into HitMap\n");
    //printf("FindRawClusters -- ndigits %d \n",ndigits);
    if (!ndigits) return;
    AliMUONDigit *dig;
    for (Int_t ndig=0; ndig<fNdigits; ndig++) {
	dig = (AliMUONDigit*)fDigits->UncheckedAt(ndig);
	SetHit(dig->fPadX,dig->fPadY,ndig);
    }
}


void  AliMUONHitMapA1::SetHit(Int_t ix, Int_t iy, Int_t idigit)
{
// Assign digit to hit cell ix,iy
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=idigit+1;
    fHitMap[CheckedIndex(ix, iy)]=idigit+1;
}

void AliMUONHitMapA1::DeleteHit(Int_t ix, Int_t iy)
{
// Delete hit at cell ix,iy
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=0;
    fHitMap[CheckedIndex(ix, iy)]=0;
}

void AliMUONHitMapA1::FlagHit(Int_t ix, Int_t iy)
{
// Flag hit as used
    fHitMap[CheckedIndex(ix, iy)]=
	-TMath::Abs(fHitMap[CheckedIndex(ix, iy)]);
}

Int_t AliMUONHitMapA1::GetHitIndex(Int_t ix, Int_t iy)
{
// Get absolute value of contents of hit cell ix,iy
    return TMath::Abs(fHitMap[CheckedIndex(ix, iy)])-1;
}

TObject* AliMUONHitMapA1::GetHit(Int_t ix, Int_t iy)
{
    // Get pointer to object at hit cell ix, iy
    // Force crash if index does not exist ! (Manu)
    Int_t index=GetHitIndex(ix,iy);
    return (index <0) ? 0 : fDigits->UncheckedAt(GetHitIndex(ix,iy));
}

FlagType AliMUONHitMapA1::TestHit(Int_t ix, Int_t iy)
{
// Check if hit cell is empty, used or unused
//
    Int_t inf=fHitMap[CheckedIndex(ix, iy)];
    if (inf < 0) {
	return kUsed;
    } else if (inf == 0) {
	return kEmpty;
    } else {
	return kUnused;
    }
}

AliMUONHitMapA1 & AliMUONHitMapA1::operator = (const AliMUONHitMapA1 & rhs) 
{
// Dummy assignment operator
    return *this;
}


