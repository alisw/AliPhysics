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

#include <TH1.h>
#include "AliITSMap.h"

ClassImp(AliITSMap)

ClassImp(AliITSMapA1)

AliITSMapA1::AliITSMapA1(AliITSsegmentation *seg)
{
  //constructor
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Int_t[fMaxIndex];
  fObjects = 0;
  ClearMap();
}

AliITSMapA1::AliITSMapA1(AliITSsegmentation *seg, TObjArray *obj)
{
  //constructor
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Int_t[fMaxIndex];
  fObjects =  obj;
  if (fObjects) fNobjects = fObjects->GetEntriesFast();
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
  if (index > fMaxIndex) {
    printf("\n \n \n Try to read/write outside array !!!! \n \n %d %d %d %d %d %d",iz,ix, fMaxIndex, index, fNpz, fNpx);
    // force crash
    return  -1;
  } else {
    return index;
  }
}


void  AliITSMapA1::FillMap()
{
  // fill array with digits indices
  Int_t ndigits = fObjects->GetEntriesFast();
  //printf("MapA1: ndigits fNobjects %d %d \n",ndigits,fNobjects);
  if (!ndigits) return;
  
  AliITSdigit *dig;
  Int_t ndig;
  for(ndig=0; ndig<ndigits; ndig++) {
    dig = (AliITSdigit*)fObjects->UncheckedAt(ndig);
    //printf("MapA1: ndig fCoord1 fCoord2 %d %d %d \n",dig->fCoord1,dig->fCoord2,ndig);
    SetHit(dig->fCoord1,dig->fCoord2,ndig);
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


Flag_t AliITSMapA1::TestHit(Int_t iz, Int_t ix)
{
  // check whether the digit has already been flagged
  Int_t inf=fHitMap[CheckedIndex(iz, ix)]; 
  if (inf < 0) {
    return kUsed;
  } else if (inf == 0) {
	return kEmpty;
  } else {
    return kUnused;
    }
}
//_______________________________________________________________________
void AliITSMapA1::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliITSMapA1.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliITSMap::Streamer(R__b);
      R__b >> fSegmentation;
      R__b >> fNpx;
      R__b >> fNpz;
      R__b >> fObjects;
      R__b >> fNobjects;
      R__b >> fMaxIndex;
      R__b.ReadArray(fHitMap);
   } else {
      R__b.WriteVersion(AliITSMapA1::IsA());
      AliITSMap::Streamer(R__b);
      R__b << fSegmentation;
      R__b << fNpx;
      R__b << fNpz;
      R__b << fObjects;
      R__b << fNobjects;
      R__b << fMaxIndex;
      R__b.WriteArray(fHitMap,fMaxIndex);
   }
}

//========================================================================
ClassImp(AliITSMapA2)

  AliITSMapA2::AliITSMapA2(AliITSsegmentation *seg)
{
  //constructor
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Double_t[fMaxIndex];
  fMapThreshold=0.;
  ClearMap();
}

//--------------------------------------
AliITSMapA2::AliITSMapA2(AliITSsegmentation *seg, TObjArray *hist, Double_t thresh)
{
  //constructor
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Double_t[fMaxIndex];
    fObjects =  hist;
    if (fObjects) fNobjects = fObjects->GetEntriesFast();
    fMapThreshold = thresh;
    ClearMap();
}
//--------------------------------------


AliITSMapA2::~AliITSMapA2()
{
  //destructor
  if (fHitMap) delete[] fHitMap;
}
//--------------------------------------

//__________________________________________________________________________
AliITSMapA2::AliITSMapA2(const AliITSMapA2 &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fMapThreshold = source.fMapThreshold;
  this->fHitMap = source.fHitMap;
  return;
}

//_________________________________________________________________________
AliITSMapA2& 
  AliITSMapA2::operator=(const AliITSMapA2 &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fMapThreshold = source.fMapThreshold;
  this->fHitMap = source.fHitMap;
  return *this;
}

void AliITSMapA2::ClearMap()
{
  //clear array
  memset(fHitMap,0,sizeof(Double_t)*fMaxIndex);
}

//--------------------------------------
void  AliITSMapA2::FillMap()
{
  
  // fills signal map from digits - apply a threshold for signal
   
  if (!fObjects) return; 
  
  Int_t ndigits = fObjects->GetEntriesFast();
  printf("MapA2: ndigits fNobjects %d %d \n",ndigits,fNobjects);
  if (!ndigits) return;
  
  AliITSdigit *dig;
  Int_t ndig;
  for(ndig=0; ndig<ndigits; ndig++) {
    dig = (AliITSdigit*)fObjects->UncheckedAt(ndig);
    Double_t signal = (Double_t)(dig->fSignal);
    if (signal >= fMapThreshold) SetHit(dig->fCoord1,dig->fCoord2,signal);
  }
}

//--------------------------------------
void  AliITSMapA2::SetHit(Int_t iz, Int_t ix, Double_t signal)
{
  // set signal at a certain position in array
  fHitMap[CheckedIndex(iz, ix)]=signal;
  
}

//--------------------------------------
void AliITSMapA2::DeleteHit(Int_t iz, Int_t ix)
{
  //set the entry value to zero
  fHitMap[CheckedIndex(iz, ix)]=0;
}

//--------------------------------------
void AliITSMapA2::FlagHit(Int_t iz, Int_t ix)
{
  //flag an entry
  fHitMap[CheckedIndex(iz, ix)]=
    -1000.*TMath::Abs((Int_t)(fHitMap[CheckedIndex(iz, ix)])+1.);
  
}

//--------------------------------------
Int_t AliITSMapA2::GetHitIndex(Int_t iz, Int_t ix)
{
  //return the index of an entry in array 
  return CheckedIndex(iz, ix);
}

//--------------------------------------
TObject* AliITSMapA2::GetHit(Int_t i, Int_t dummy)
{
  
  //return a pointer to the 1D histogram
  if (fObjects) {
    
    return fObjects->UncheckedAt(i); 
    
  } else return NULL;
  
}

//--------------------------------------
Double_t AliITSMapA2::GetSignal(Int_t iz, Int_t ix)
{
  //get signal in a cell 
  Int_t index=GetHitIndex(iz,ix);
  return (index <0) ? 0. : fHitMap[CheckedIndex(iz, ix)];
}

//--------------------------------------
Double_t AliITSMapA2::GetSignal(Int_t index)
{
  //get signal in a cell 
  if (index<fMaxIndex) return (index <0) ? 0. : fHitMap[index];
  else return 0.;
}
//--------------------------------------
Flag_t AliITSMapA2::TestHit(Int_t iz, Int_t ix)
{
  // check if the entry has already been flagged
    Int_t inf=(Int_t)fHitMap[CheckedIndex(iz, ix)];
    
    if (inf <= -1000) {
      return kUsed;
    } else if (inf == 0) {
      return kEmpty;
    } else {
      return kUnused;
    }
}

//--------------------------------------
void  AliITSMapA2::FillMapFromHist()
{
  
  // fills map from 1D histograms
  
  if (!fObjects) return; 
  
  // an example
  Int_t i,j;
  for(i=0; i<fNobjects; i++) {
    TH1F *hist =(TH1F *)fObjects->UncheckedAt(i);
    Int_t nsamples = hist->GetNbinsX();
    for(j=0; j<nsamples; j++) {
      Double_t signal = (Double_t)(hist->GetBinContent(j+1));
      if (signal >= fMapThreshold) SetHit(i,j,signal);
    }
  }
  
}
//--------------------------------------
void  AliITSMapA2::FillHist()
{
  
  // fill 1D histograms from map
  if (!fObjects) return; 
  
  // an example
  Int_t i,j;
  for(i=0; i<fNobjects; i++) {
    TH1F *hist =(TH1F *)fObjects->UncheckedAt(i);
    for(j=0; j<fNpx; j++) {
      Double_t signal=GetSignal(i,j);
      if (signal >= fMapThreshold) hist->Fill((Float_t)j,signal);
    }
  }
  
}
//--------------------------------------
void  AliITSMapA2::ResetHist()
{
  //
  // Reset histograms 
  //
  
  if (!fObjects) return; 
  
  Int_t i;
  for(i=0; i<fNobjects; i++) {
    if ((*fObjects)[i])    ((TH1F*)(*fObjects)[i])->Reset();
  }
  
}
//______________________________________________________________________________
void AliITSMapA2::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliITSMapA2.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliITSMapA1::Streamer(R__b);
      R__b.ReadArray(fHitMap);
      R__b >> fMapThreshold;
   } else {
      R__b.WriteVersion(AliITSMapA2::IsA());
      AliITSMapA1::Streamer(R__b);
      R__b.WriteArray(fHitMap, fMaxIndex); // fMaxIndex is from AliITSMapA1.
      R__b << fMapThreshold;
   }
}
