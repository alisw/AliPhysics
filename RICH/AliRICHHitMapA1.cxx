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


#include "AliRICHHitMapA1.h"
#include "AliSegmentation.h"
#include "AliRICHDigit.h"

#include <TObjArray.h>
#include <TMath.h>

ClassImp(AliRICHHitMapA1)


AliRICHHitMapA1::AliRICHHitMapA1(AliSegmentation *seg, TObjArray *dig)
{

// Constructor for AliRICHMapA1

    fSegmentation = seg;
    fNpx  = fSegmentation->Npx();
    fNpy  = fSegmentation->Npy();
    fMaxIndex=2*(fNpx+1)*2*(fNpy+1)+2*fNpy;
    
    fHitMap = new Int_t[fMaxIndex];
    fDigits =  dig;
    fNdigits = fDigits->GetEntriesFast();
    Clear();
}


AliRICHHitMapA1::~AliRICHHitMapA1()
{
// Destructor
//    if (fDigits) delete   fDigits;
    if (fHitMap) delete[] fHitMap;
}

void AliRICHHitMapA1::Clear(const char *opt)
{

// Clear contents of hit map

    memset(fHitMap,0,sizeof(int)*fMaxIndex);
}

Int_t AliRICHHitMapA1::CheckedIndex(Int_t ix, Int_t iy) const
{

// Check if index is valid

    Int_t index=2*fNpy*(ix+fNpx)+(iy+fNpy);
    if (index > fMaxIndex) {
	printf("\n \n \n Try to read/write outside array !!!! \n \n %d %d %d %d %d %d",ix,iy, fMaxIndex, index, fNpx, fNpy);
	return  fMaxIndex-1;
    } else {
	return index;
    }
}

	
void  AliRICHHitMapA1::FillHits()
{

// Fill hits into HitMap

    Int_t ndigits = fDigits->GetEntriesFast();
    //printf("\n Filling hits into HitMap\n");
    //printf("FindRawClusters -- ndigits %d \n",ndigits);
    if (!ndigits) return;
    AliRICHDigit *dig;
    for (Int_t ndig=0; ndig<fNdigits; ndig++) {
	dig = (AliRICHDigit*)fDigits->UncheckedAt(ndig);
	SetHit(dig->fPadX,dig->fPadY,ndig);
    }
}


void  AliRICHHitMapA1::SetHit(Int_t ix, Int_t iy, Int_t idigit)
{
//
// Set current hit
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=idigit+1;
    fHitMap[CheckedIndex(ix, iy)]=idigit+1;
}

void AliRICHHitMapA1::DeleteHit(Int_t ix, Int_t iy)
{
// Delete hit
//
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=0;
    fHitMap[CheckedIndex(ix, iy)]=0;
}

void AliRICHHitMapA1::FlagHit(Int_t ix, Int_t iy)
{
//
// Flag hit

    fHitMap[CheckedIndex(ix, iy)]=
	-TMath::Abs(fHitMap[CheckedIndex(ix, iy)]);
}

Int_t AliRICHHitMapA1::GetHitIndex(Int_t ix, Int_t iy) const
{

// Return hit coordinates from index

  //printf("ix:%d, iy:%d, index:%d\n",ix,iy,CheckedIndex(ix, iy));
   
  return TMath::Abs(fHitMap[CheckedIndex(ix, iy)])-1;
    
}

TObject* AliRICHHitMapA1::GetHit(Int_t ix, Int_t iy) const
{

// Return index from coordinates

    Int_t index=GetHitIndex(ix,iy);
    // Force crash if index does not exist ! (Manu)
    return (index <0) ? 0 : fDigits->UncheckedAt(GetHitIndex(ix,iy));
}

FlagType AliRICHHitMapA1::TestHit(Int_t ix, Int_t iy)
{

// Is there a hit?

    Int_t inf=fHitMap[CheckedIndex(ix, iy)];
    if (inf < 0) {
	return kUsed;
    } else if (inf == 0) {
	return kEmpty;
    } else {
	return kUnused;
    }
}


