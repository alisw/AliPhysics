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
#include <TObjArray.h>
#include <TMath.h>

#include "AliITSMapA2.h"
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"
#include "AliITSdigit.h"


ClassImp(AliITSMapA2)

  AliITSMapA2::AliITSMapA2(AliITSsegmentation *seg)
{
  //constructor
  fScaleSizeZ=1;
  fScaleSizeX=1;
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Double_t[fMaxIndex];
  fMapThreshold=0.;
  ClearMap();
}
//--------------------------------------
  AliITSMapA2::AliITSMapA2(AliITSsegmentation *seg, Int_t scalesizeX, Int_t scalesizeZ)
{
  //constructor
  fSegmentation = seg;
  fScaleSizeX=scalesizeX;
  fScaleSizeZ=scalesizeZ;
  fNpz=fScaleSizeZ*fSegmentation->Npz();
  fNpx=fScaleSizeX*fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Double_t[fMaxIndex];
  fMapThreshold=0.;
  ClearMap();
}

//--------------------------------------
AliITSMapA2::AliITSMapA2(AliITSsegmentation *seg, TObjArray *obj, Double_t thresh)
{
  //constructor
  fScaleSizeZ=1;
  fScaleSizeX=1;
  fSegmentation = seg;
  fNpz=fSegmentation->Npz();
  fNpx=fSegmentation->Npx();
  fMaxIndex=fNpz*fNpx+fNpx;             // 2 halves of detector
  
  fHitMap = new Double_t[fMaxIndex];
  fObjects =  obj;
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
  this->fScaleSizeX = source.fScaleSizeX;
  this->fScaleSizeZ = source.fScaleSizeZ;
  this->fHitMap = source.fHitMap;
  return;
}

//_________________________________________________________________________
AliITSMapA2& 
  AliITSMapA2::operator=(const AliITSMapA2 &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fMapThreshold = source.fMapThreshold;
  this->fScaleSizeX = source.fScaleSizeX;
  this->fScaleSizeZ = source.fScaleSizeZ;
  this->fHitMap = source.fHitMap;
  return *this;
}

//_________________________________________________________________________
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
  if (!ndigits) return;
  
  AliITSdigit *dig;
  for (Int_t ndig=0; ndig<ndigits; ndig++) {
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
FlagType AliITSMapA2::TestHit(Int_t iz, Int_t ix)
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
  for( Int_t i=0; i<fNobjects; i++) {
    TH1F *hist =(TH1F *)fObjects->UncheckedAt(i);
    Int_t nsamples = hist->GetNbinsX();
    for( Int_t j=0; j<nsamples; j++) {
      Double_t signal = (Double_t)(hist->GetBinContent(j+1));
      if (signal >= fMapThreshold) SetHit(i,j,signal);
    }
  }
  
}
//--------------------------------------
void  AliITSMapA2::FillHist()
{
  
  // fill 1D histograms from map
  if (!fObjects || fScaleSizeX != 1) return; 
  
  // an example
  for( Int_t i=0; i<fNobjects; i++) {
    TH1F *hist =(TH1F *)fObjects->UncheckedAt(i);
    for( Int_t j=0; j<fNpx; j++) {
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
  
  for( Int_t i=0; i<fNobjects; i++) {
    if ((*fObjects)[i])    ((TH1F*)(*fObjects)[i])->Reset();
  }
  
}

