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
//////////////////////////////////////////////////////////////////////////////
// This is the base class for ITS detector signal simulations. Data members //
// include are a pointer to the AliITSDetTypeSim clas in order to access    //
// segmentation and response objects                                        // 
// classes. See the detector specific implementations for the propper code. //
//////////////////////////////////////////////////////////////////////////////
#include <TRandom.h>
#include "TSeqCollection.h"
#include "AliITSUSimulation.h"
#include "AliITSUSDigit.h"
#include "AliITSUModule.h"
#include "AliITSUParamList.h"
using namespace TMath;

ClassImp(AliITSUSimulation)

//______________________________________________________________________
AliITSUSimulation::AliITSUSimulation()
:  fSeg(0)
  ,fCalibDead(0)
  ,fCalibNoisy(0)
  ,fSensMap(0)
  ,fSimuParam(0)
  ,fResponseParam(0)
  ,fModule(0)
  ,fEvent(0)
  ,fDebug(0)
{
    // Default constructor
}
//______________________________________________________________________
AliITSUSimulation::AliITSUSimulation(AliITSUSimuParam* sim,AliITSUSensMap* map)
  :fSeg(0)
  ,fCalibDead(0)
  ,fCalibNoisy(0)
  ,fSensMap(map)
  ,fSimuParam(sim)
  ,fResponseParam(0)
  ,fModule(0)
  ,fEvent(0)
  ,fDebug(0)
{
    // Default constructor
}

//__________________________________________________________________________
AliITSUSimulation::AliITSUSimulation(const AliITSUSimulation &s) 
  :TObject(s)
  ,fSeg(s.fSeg)
  ,fCalibDead(s.fCalibDead)
  ,fCalibNoisy(s.fCalibNoisy)
  ,fSensMap(s.fSensMap)
  ,fSimuParam(s.fSimuParam)   
  ,fResponseParam(s.fResponseParam)
  ,fModule(s.fModule)
  ,fEvent(s.fEvent)
  ,fDebug(s.fDebug)
{
  //     Copy Constructor 
}

//_________________________________________________________________________
AliITSUSimulation&  AliITSUSimulation::operator=(const AliITSUSimulation &s)
{
  //    Assignment operator
  if(&s == this) return *this;
  fSeg       = s.fSeg;
  fCalibDead = s.fCalibDead;
  fCalibNoisy= s.fCalibNoisy;
  fSensMap   = s.fSensMap;
  fSimuParam = s.fSimuParam;
  fResponseParam = s.fResponseParam;
  fModule    = s.fModule;
  fEvent     = s.fEvent;
  return *this;
}

//______________________________________________________________________
void AliITSUSimulation::InitSimulationModule(AliITSUModule* mod, Int_t event, AliITSsegmentation* seg, AliITSUParamList* resp)
{
  //  This function creates maps to build the list of tracks for each
  //  summable digit. Inputs defined by base class.
  //
  SetModule(mod);
  SetSegmentation(seg);
  SetResponseParam(resp);
  ClearMap();
  memset(fCyclesID,0,(1+2*kMaxROCycleAccept)*sizeof(Bool_t));
  //
  if (event != fEvent) GenerateReadOutCycleOffset(); 
  SetEvent(event);
  
}

//______________________________________________________________________
Bool_t AliITSUSimulation::AddSDigitsToModule(TSeqCollection *pItemArr,Int_t mask )
{
  // Add Summable digits to module maps.
  // Inputs:
  //    pItemArr  Array of AliITSpListItems (SDigits).
  //    mask    Track number off set value 
  //
  Int_t nItems = pItemArr->GetEntries();
  Bool_t sig = kFALSE;
  // 
  for( Int_t i=0; i<nItems; i++ ) {
    AliITSUSDigit * pItem = (AliITSUSDigit *)(pItemArr->At( i ));
    if(pItem->GetModule() != int(fModule->GetIndex()) ) AliFatal(Form("SDigits module %d != current module %d: exit", pItem->GetModule(),fModule->GetIndex()));
    if(pItem->GetSumSignal()>0.0 ) sig = kTRUE;
    AliITSUSDigit* oldItem = (AliITSUSDigit*)fSensMap->GetItem(pItem);
    if (!oldItem) {
      oldItem = (AliITSUSDigit*)fSensMap->RegisterItem( new(fSensMap->GetFree()) AliITSUSDigit(*pItem) );
      if (mask) oldItem->ShiftIndices(mask);
    }
    else oldItem->AddTo(mask, pItem);
  }
  return sig;
}

//______________________________________________________________________
void AliITSUSimulation::UpdateMapSignal(UInt_t col,UInt_t row,Int_t trk,Int_t ht,Double_t signal, Int_t roCycle) 
{
  // update map with new hit
  // Note: roCycle can be anything between -kMaxROCycleAccept : kMaxROCycleAccept
  if (Abs(roCycle)>kMaxROCycleAccept) {
    AliError(Form("CycleID %d is outside of allowed +-%d range",roCycle,kMaxROCycleAccept));
    return;
  }
  UInt_t ind = fSensMap->GetIndex(col,row,roCycle);
  AliITSUSDigit* oldItem = (AliITSUSDigit*)fSensMap->GetItem(ind);  
  if (!oldItem) {    
    fSensMap->RegisterItem( new(fSensMap->GetFree()) AliITSUSDigit(trk,ht,fModule->GetIndex(),ind,signal,roCycle) );
    fCyclesID[roCycle+kMaxROCycleAccept] = kTRUE;
  }
  else oldItem->AddSignal(trk,ht,signal);
  //
}

//______________________________________________________________________
void AliITSUSimulation::UpdateMapNoise(UInt_t col,UInt_t row,Double_t noise, Int_t roCycle) 
{
  // update map with new hit
  if (Abs(roCycle)>kMaxROCycleAccept) {
    AliError(Form("CycleID %d is outside of allowed +-%d range",roCycle,kMaxROCycleAccept));
    return;
  }
  UInt_t ind = fSensMap->GetIndex(col,row,roCycle);
  AliITSUSDigit* oldItem = (AliITSUSDigit*)fSensMap->GetItem(ind);
  if (!oldItem) {
    fSensMap->RegisterItem( new(fSensMap->GetFree()) AliITSUSDigit(fModule->GetIndex(),ind,noise,roCycle) );
    fCyclesID[roCycle+kMaxROCycleAccept] = kTRUE;
  }
  else oldItem->AddNoise(noise);
}

//______________________________________________________________________
Int_t AliITSUSimulation::GenOrderedSample(UInt_t nmax,UInt_t ngen,TArrayI &vals,TArrayI &indx)
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
