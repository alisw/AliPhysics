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


#include "AliRICHTresholdMap.h"
#include "AliRICHSegmentation.h"
#include "AliRICHDigit.h"

#include <TObjArray.h>
#include <TMath.h>

ClassImp(AliRICHTresholdMap)


AliRICHTresholdMap::AliRICHTresholdMap(AliRICHSegmentation *seg)
{

// Constructor for AliRICHTresholdMap

    fSegmentation = seg;
    fNpx  = fSegmentation->Npx();
    fNpy  = fSegmentation->Npy();
    fMaxIndex=2*(fNpx+1)*2*(fNpy+1)+2*fNpy;
    
    fHitMap = new Int_t[fMaxIndex];
    Clear();
}


AliRICHTresholdMap::~AliRICHTresholdMap()
{
// Destructor
//    if (fDigits) delete   fDigits;
    if (fHitMap) delete[] fHitMap;
}

void AliRICHTresholdMap::Clear()
{

// Clear contents of hit map

    memset(fHitMap,0,sizeof(int)*fMaxIndex);
}

Int_t AliRICHTresholdMap::CheckedIndex(Int_t ix, Int_t iy)
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

	
void  AliRICHTresholdMap::FillHits()
{

// Dummy
}


void  AliRICHTresholdMap::SetHit(Int_t ix, Int_t iy, Int_t charge)
{
//
// Set current hit
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=idigit+1;

     fHitMap[CheckedIndex(ix, iy)]=charge;
}

void AliRICHTresholdMap::DeleteHit(Int_t ix, Int_t iy)
{
// Delete hit
//
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=0;
    fHitMap[CheckedIndex(ix, iy)]=0;
}

void AliRICHTresholdMap::FlagHit(Int_t ix, Int_t iy)
{
//
// Flag hit

    fHitMap[CheckedIndex(ix, iy)]=
	-TMath::Abs(fHitMap[CheckedIndex(ix, iy)]);
}

Int_t AliRICHTresholdMap::GetHitIndex(Int_t ix, Int_t iy)
{

// Return hit coordinates from index

  //printf("ix:%d, iy:%d, index:%d\n",ix,iy,CheckedIndex(ix, iy));
   
  return TMath::Abs(fHitMap[CheckedIndex(ix, iy)]);
    
}

TObject* AliRICHTresholdMap::GetHit(Int_t ix, Int_t iy)
{

// Return index from coordinates

    Int_t index=GetHitIndex(ix,iy);
    if (index>0)
      return((TObject*) index);
    else
      return((TObject*) 0);
}

FlagType AliRICHTresholdMap::TestHit(Int_t ix, Int_t iy)
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

