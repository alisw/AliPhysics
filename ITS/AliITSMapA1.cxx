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


#include "AliITSMapA1.h"
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"
#include "AliITSdigit.h"

#include <TObjArray.h>
#include <TMath.h>


ClassImp(AliITSMapA1)

AliITSMapA1::AliITSMapA1()
{
  // default constructor
  fSegmentation = 0;
  fNpz=0;
  fNpx=0;
  fMaxIndex=0;         
  
  fHitMap = 0;
  fObjects = 0;
  fNobjects = 0;
  fMapThreshold=0;
}

AliITSMapA1::AliITSMapA1(AliITSsegmentation *seg)
{
  //constructor
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Int_t[fMaxIndex];
  fObjects = 0;
  fNobjects = 0;
  fMapThreshold=0;
  ClearMap();
}

AliITSMapA1::AliITSMapA1(AliITSsegmentation *seg, TObjArray *obj)
{
  //constructor
  fNobjects = 0;
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Int_t[fMaxIndex];
  fObjects =  obj;
  if (fObjects) fNobjects = fObjects->GetEntriesFast();
  fMapThreshold=0;
  ClearMap();
}

AliITSMapA1::AliITSMapA1(AliITSsegmentation *seg, TObjArray *obj, Int_t thr)
{
  //constructor
  fNobjects = 0;
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Int_t[fMaxIndex];
  fObjects =  obj;
  if (fObjects) fNobjects = fObjects->GetEntriesFast();
  fMapThreshold=thr;
  ClearMap();
}


AliITSMapA1::~AliITSMapA1()
{
  //destructor
  if (fHitMap) delete[] fHitMap;
}

//__________________________________________________________________________
AliITSMapA1::AliITSMapA1(const AliITSMapA1 &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fNpx = source.fNpx;
  this->fNpz = source.fNpz;
  this->fObjects = source.fObjects;
  this->fNobjects = source.fNobjects;
  this->fMaxIndex = source.fMaxIndex;
  this->fHitMap = source.fHitMap;
  this->fMapThreshold = source.fMapThreshold;
  return;
}

//_________________________________________________________________________
AliITSMapA1& 
  AliITSMapA1::operator=(const AliITSMapA1 &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fNpx = source.fNpx;
  this->fNpz = source.fNpz;
  this->fObjects = source.fObjects;
  this->fNobjects = source.fNobjects;
  this->fMaxIndex = source.fMaxIndex;
  this->fHitMap = source.fHitMap;
  this->fMapThreshold = source.fMapThreshold;
  return *this;
}

void AliITSMapA1::ClearMap()
{
  //clear array
  memset(fHitMap,0,sizeof(int)*fMaxIndex);
}

void AliITSMapA1::SetArray(TObjArray *obj)
{
  // set array of objects
  fObjects =  obj;
  if (fObjects) fNobjects = fObjects->GetEntriesFast();
}


Int_t AliITSMapA1::CheckedIndex(Int_t iz, Int_t ix)
{
  //check boundaries and return an index in array
  Int_t index=fNpx*iz+ix;

  //if (index > fMaxIndex) {
  if (index > fMaxIndex || index < 0) {
    printf("\n \n \n Try to read/write outside array !!!! \n \n %d %d %d %d %d %d \n",iz,ix, fMaxIndex, index, fNpz, fNpx);
    // force crash
    return  -1;
  } else {
    return index;
  }
}


void  AliITSMapA1::FillMap()
{
  // fill array with digits indices

  //printf("fMapThreshold %d\n",fMapThreshold);
  Int_t ndigits = fObjects->GetEntriesFast();
  if (!ndigits) return;
  
  AliITSdigit *dig;
  for (Int_t ndig=0; ndig<ndigits; ndig++) {
    dig = (AliITSdigit*)fObjects->UncheckedAt(ndig);
    if(dig->fSignal > fMapThreshold) {
      SetHit(dig->fCoord1,dig->fCoord2,ndig);
    }
  }
  
}

void  AliITSMapA1::SetHit(Int_t iz, Int_t ix, Int_t idigit)
{
  // set the digit index at a certain position in array
  fHitMap[CheckedIndex(iz, ix)]=idigit+1;
}

void AliITSMapA1::DeleteHit(Int_t iz, Int_t ix)
{
  // delete an entry in array
  fHitMap[CheckedIndex(iz, ix)]=0;
}

void AliITSMapA1::FlagHit(Int_t iz, Int_t ix)
{
  // flag an entry in array
  fHitMap[CheckedIndex(iz, ix)]=
      -TMath::Abs(fHitMap[CheckedIndex(iz, ix)]);

}

Int_t AliITSMapA1::GetHitIndex(Int_t iz, Int_t ix)
{
  // return the digit index from a specific entry in array
  return TMath::Abs(fHitMap[CheckedIndex(iz, ix)])-1;
}

TObject* AliITSMapA1::GetHit(Int_t iz, Int_t ix)
{
  // return the pointer to the digit 
  Int_t index=GetHitIndex(iz,ix);
  // Force crash if index does not exist ! 
  return (index <0) ? 0 : fObjects->UncheckedAt(GetHitIndex(iz,ix));
}

Double_t AliITSMapA1::GetSignal(Int_t iz, Int_t ix)
{
  // get a pad signal 
  Double_t signal;
  AliITSdigit *dig = (AliITSdigit*)GetHit(iz,ix);
  if(dig) signal=(Double_t)dig->fSignal;
  else signal=0.;
  return signal;

}

FlagType AliITSMapA1::TestHit(Int_t iz, Int_t ix)
{
  // check whether the digit has already been flagged

  if (CheckedIndex(iz, ix) < 0) return kEmpty;
  Int_t inf=fHitMap[CheckedIndex(iz, ix)]; 
  if (inf < 0) {
      return kUsed;
  } else if (inf == 0) {
      return kEmpty;
  } else {
      return kUnused;
  }
}
