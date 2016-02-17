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
///////////////////////////////////////////////////////////////////////////////
// This is the base class for ITSMFT detector signal simulations. Data members //
///////////////////////////////////////////////////////////////////////////////
#include <TRandom.h>
#include <TArrayI.h>
#include "TSeqCollection.h"
#include "AliLog.h"
#include "AliITSMFTSimulation.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliITSMFTSDigit.h"
#include "AliITSMFTChip.h"
#include "AliITSMFTParamList.h"
using namespace TMath;

ClassImp(AliITSMFTSimulation)

//______________________________________________________________________
AliITSMFTSimulation::AliITSMFTSimulation()
:  fSeg(0)
  ,fCalibDead(0)
  ,fCalibNoisy(0)
  ,fSensMap(0)
  ,fSimuParam(0)
  ,fResponseParam(0)
  ,fChip(0)
  ,fReadOutCycleOffset(0)
  ,fReadOutCycleLength(25e-6)  
  ,fEvent(0)
  ,fDebug(0)
{
    // Default constructor
}
//______________________________________________________________________
AliITSMFTSimulation::AliITSMFTSimulation(AliITSMFTSimuParam* sim,AliITSMFTSensMap* map)
  :fSeg(0)
  ,fCalibDead(0)
  ,fCalibNoisy(0)
  ,fSensMap(map)
  ,fSimuParam(sim)
  ,fResponseParam(0)
  ,fChip(0)
  ,fReadOutCycleOffset(0)
  ,fReadOutCycleLength(25e-6)
  ,fEvent(0)
  ,fDebug(0)
{
    // Default constructor
}

//__________________________________________________________________________
AliITSMFTSimulation::AliITSMFTSimulation(const AliITSMFTSimulation &s) 
  :TObject(s)
  ,fSeg(s.fSeg)
  ,fCalibDead(s.fCalibDead)
  ,fCalibNoisy(s.fCalibNoisy)
  ,fSensMap(s.fSensMap)
  ,fSimuParam(s.fSimuParam)   
  ,fResponseParam(s.fResponseParam)
  ,fChip(s.fChip)
  ,fReadOutCycleOffset(s.fReadOutCycleOffset)
  ,fReadOutCycleLength(s.fReadOutCycleLength)
  ,fEvent(s.fEvent)
  ,fDebug(s.fDebug)
{
  //     Copy Constructor 
}

//_________________________________________________________________________
AliITSMFTSimulation&  AliITSMFTSimulation::operator=(const AliITSMFTSimulation &s)
{
  //    Assignment operator
  if(&s == this) return *this;
  fSeg       = s.fSeg;
  fCalibDead = s.fCalibDead;
  fCalibNoisy= s.fCalibNoisy;
  fSensMap   = s.fSensMap;
  fSimuParam = s.fSimuParam;
  fResponseParam = s.fResponseParam;
  fChip    = s.fChip;
  fReadOutCycleOffset = s.fReadOutCycleOffset;
  fReadOutCycleLength = s.fReadOutCycleLength;
  fEvent     = s.fEvent;
  return *this;
}

//______________________________________________________________________
void AliITSMFTSimulation::InitSimulationChip(AliITSMFTChip* mod, Int_t event, AliITSMFTSegmentationPix* seg, AliITSMFTParamList* resp)
{
  //  This function creates maps to build the list of tracks for each
  //  summable digit. Inputs defined by base class.
  //
  SetChip(mod);
  SetSegmentation(seg);
  SetResponseParam(resp);
  ClearMap();
  memset(fCyclesID,0,(1+2*kMaxROCycleAccept)*sizeof(Bool_t));
  //
  SetEvent(event);
  
}

//______________________________________________________________________
Bool_t AliITSMFTSimulation::AddSDigitsToChip(TSeqCollection *pItemArr,Int_t mask )
{
  // Add Summable digits to chip maps.
  // Inputs:
  //    pItemArr  Array of AliITSpListItems (SDigits).
  //    mask    Track number off set value 
  //
  Int_t nItems = pItemArr->GetEntries();
  Bool_t sig = kFALSE;
  // 
  for( Int_t i=0; i<nItems; i++ ) {
    AliITSMFTSDigit * pItem = (AliITSMFTSDigit *)(pItemArr->At( i ));
    if(pItem->GetChip() != int(fChip->GetIndex()) ) AliFatal(Form("SDigits chip %d != current chip %d: exit", pItem->GetChip(),fChip->GetIndex()));
    if(pItem->GetSumSignal()>0.0 ) sig = kTRUE;
    AliITSMFTSDigit* oldItem = (AliITSMFTSDigit*)fSensMap->GetItem(pItem);
    if (!oldItem) {
      oldItem = (AliITSMFTSDigit*)fSensMap->RegisterItem( new(fSensMap->GetFree()) AliITSMFTSDigit(*pItem) );
      if (mask) oldItem->ShiftIndices(mask);
    }
    else oldItem->AddTo(mask, pItem);
  }
  return sig;
}

//______________________________________________________________________
void AliITSMFTSimulation::UpdateMapSignal(UInt_t col,UInt_t row,Int_t trk,Int_t ht,Double_t signal, Int_t roCycle) 
{
  // update map with new hit
  // Note: roCycle can be anything between -kMaxROCycleAccept : kMaxROCycleAccept
  if (Abs(roCycle)>kMaxROCycleAccept) {
    AliError(Form("CycleID %d is outside of allowed +-%d range",roCycle,kMaxROCycleAccept));
    return;
  }
  UInt_t ind = fSensMap->GetIndex(col,row,roCycle);
  AliITSMFTSDigit* oldItem = (AliITSMFTSDigit*)fSensMap->GetItem(ind);  
  if (!oldItem) {    
    fSensMap->RegisterItem( new(fSensMap->GetFree()) AliITSMFTSDigit(trk,ht,fChip->GetIndex(),ind,signal,roCycle) );
    fCyclesID[roCycle+kMaxROCycleAccept] = kTRUE;
  }
  else oldItem->AddSignal(trk,ht,signal);
  //
}

//______________________________________________________________________
void AliITSMFTSimulation::UpdateMapNoise(UInt_t col,UInt_t row,Double_t noise, Int_t roCycle) 
{
  // update map with new hit
  if (Abs(roCycle)>kMaxROCycleAccept) {
    AliError(Form("CycleID %d is outside of allowed +-%d range",roCycle,kMaxROCycleAccept));
    return;
  }
  UInt_t ind = fSensMap->GetIndex(col,row,roCycle);
  AliITSMFTSDigit* oldItem = (AliITSMFTSDigit*)fSensMap->GetItem(ind);
  if (!oldItem) {
    fSensMap->RegisterItem( new(fSensMap->GetFree()) AliITSMFTSDigit(fChip->GetIndex(),ind,noise,roCycle) );
    fCyclesID[roCycle+kMaxROCycleAccept] = kTRUE;
  }
  else oldItem->AddNoise(noise);
}

//______________________________________________________________________
Int_t AliITSMFTSimulation::GenOrderedSample(UInt_t nmax,UInt_t ngen,TArrayI &vals,TArrayI &indx)
{
  // generate random sample [0:nmax] of ngen variables, and fill orreder indices 
  // return actual number of generated values 
  if (vals.GetSize()<(int)ngen) vals.Set(ngen);
  if (indx.GetSize()<(int)ngen) indx.Set(ngen);
  int* valA = vals.GetArray();
  int* indA = indx.GetArray();
  if (ngen>=nmax) {
    ngen = nmax-1;
    for (int i=(int)ngen;i--;) {valA[i]=indA[i]=i;}
    return ngen;
  }
  Bool_t rep;
  for (int i=0;i<(int)ngen;i++) {
    do { // exclude repetitions
      rep = kFALSE;
      valA[i] = gRandom->Rndm()*nmax;
      for (int j=i;j--;) if (valA[j]==valA[i]) {rep=kTRUE;break;}
    } while(rep);
  }
  Sort((int)ngen,valA,indA,kFALSE);
  return ngen;
}

//______________________________________________________________________
Double_t AliITSMFTSimulation::GenerateReadOutCycleOffset()
{
  // Generate randomly the strobe
  // phase w.r.t to the LHC clock
  return fReadOutCycleOffset = fReadOutCycleLength*gRandom->Rndm();
  // fReadOutCycleOffset = 25e-9*gRandom->Rndm(); // clm: I think this way we shift too much 10-30 us! The global shift should be between the BCs?!
  // RS: 25 ns is too small number, the staggering will not work. Let's at the moment keep fully random shift (still, no particle from correct
  // collision will be lost) untill real number is specified
 //
}
