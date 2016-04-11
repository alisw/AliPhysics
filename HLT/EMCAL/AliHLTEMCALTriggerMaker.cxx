/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Markus Fasel                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include <vector>

#include "AliEMCALTriggerTypes.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerBitConfig.h"
#include "AliEMCALTriggerMapping.h"
#include "AliEMCALTriggerTypes.h"
#include "AliEMCALTriggerAlgorithm.h"
#include "AliEMCALTriggerDataGrid.h"
#include "AliEMCALTriggerPatchFinder.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerRawPatch.h"
#include "AliEMCALTriggerConstants.h"

#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloTriggerPatchContainerStruct.h"
#include "AliHLTCaloTriggerPatchDataStruct.h"
#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALTriggerMaker.h"

ClassImp(AliHLTEMCALTriggerMaker)

AliHLTEMCALTriggerMaker::AliHLTEMCALTriggerMaker() :
  TObject(),
  AliHLTLogging(),
  fTriggerPatchDataPtr(NULL),
  fkGeometryPtr(NULL),
  fPatchFinder(NULL),
  fL0PatchFinder(NULL),
  fADCValues(NULL),
  fADCOfflineValues(NULL),
  fL0Amplitudes(NULL),
  fTriggerBitMasks(NULL),
  fTriggerBitConfig(NULL),
  fJetPatchSize(8),
  fJetSubregionSize(4),
  fGammaPatchSize(2),
  fGammaSubregionSize(1),
  fBkgPatchSize(8),
  fBkgSubregionSize(4),
  fBufferSize(0),
  fBkgThresholdOnline(-0.1),
  fBkgThresholdOffline(-0.1),
  fRunBkgAlgorithm(kTRUE),
  fLevel0ThresholdOnline(0),
  fLevel0ThresholdOffline(0)
{
  fTriggerBitConfig = new AliEMCALTriggerBitConfigNew;
  memset(fJetThresholdOnline, 0, sizeof(Float_t) * kNthresholds);
  memset(fJetThresholdOffline, 0, sizeof(Float_t) * kNthresholds);
  memset(fGammaThresholdOnline, 0, sizeof(Float_t) * kNthresholds);
  memset(fGammaThresholdOffline, 0, sizeof(Float_t) * kNthresholds);
}

AliHLTEMCALTriggerMaker::~AliHLTEMCALTriggerMaker() {
  if(fTriggerBitConfig) delete fTriggerBitConfig;
  if(fPatchFinder) delete fPatchFinder;
  if(fL0PatchFinder) delete fL0PatchFinder;
  if(fADCValues) delete fADCValues;
  if(fADCOfflineValues) delete fADCOfflineValues;
  if(fL0Amplitudes) delete fL0Amplitudes;
  if(fTriggerBitMasks) delete fTriggerBitMasks;
}

void AliHLTEMCALTriggerMaker::ResetADC(){
  fADCValues->Reset();
  fADCOfflineValues->Reset();
  fL0Amplitudes->Reset();
  fTriggerBitMasks->Reset();
}

void AliHLTEMCALTriggerMaker::AddDigit(const AliHLTCaloDigitDataStruct *digit){
  /*
   * TODO Crosscheck
   */
  Int_t fastorIndex;
  fkGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetFastORIndexFromCellIndex(fkGeometryPtr->GetGeometryPtr()->GetAbsCellIdFromCellIndexes(digit->fModule, digit->fX, digit->fZ), fastorIndex);
  int globCol, globRow;
  fkGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetPositionInEMCALFromAbsFastORIndex(fastorIndex, globCol, globRow);
  (*fADCOfflineValues)(globCol, globRow) += digit->fEnergy;
}

void AliHLTEMCALTriggerMaker::SetADC(Int_t col, Int_t row, Float_t adc){
  (*fADCValues)(col, row) = adc;
}

void AliHLTEMCALTriggerMaker::SetL0Amplitude(Int_t col, Int_t row, Float_t amp){
  (*fL0Amplitudes)(col, row) = amp;
}

void AliHLTEMCALTriggerMaker::SetBitMask(Int_t col, Int_t row, Int_t bitMask){
  (*fTriggerBitMasks)(col, row) = bitMask;
}

Int_t AliHLTEMCALTriggerMaker::FindPatches(){
  /*
  if(availableSize < sizeof(AliHLTCaloTriggerPatchDataStruct)){
    HLTError("Not enough space to write new trigger patches");
    return -1;
  }
  */

  //AliHLTUInt32_t mysize = availableSize;
  std::vector<AliEMCALTriggerRawPatch> foundpatches = fPatchFinder->FindPatches(*fADCValues, *fADCOfflineValues);
  Int_t patchcount = 0;
AliHLTCaloTriggerPatchDataStruct *next = NULL;
  for(std::vector<AliEMCALTriggerRawPatch>::iterator patchiter = foundpatches.begin(); patchiter != foundpatches.end(); ++patchiter){
	if(fBufferSize < sizeof(AliHLTCaloTriggerPatchDataStruct)){
		HLTWarning("Buffer exceeded after %d trigger patches", patchcount);
		break;
	}
    next = fTriggerPatchDataPtr + 1;
    // Set offline bits
    // Note: trigger patch can contain more than one patch type
    Int_t offlinebits = 0;
    if(patchiter->GetPatchSize() == fGammaPatchSize){
      if(patchiter->GetADC() > fGammaThresholdOnline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchiter->GetOfflineADC() > fGammaThresholdOffline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchiter->GetADC() > fGammaThresholdOnline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetGammaLowBit());
      if(patchiter->GetOfflineADC() > fGammaThresholdOffline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetGammaLowBit());
    }
    if (patchiter->GetPatchSize() == fJetPatchSize){
      if(patchiter->GetADC() > fJetThresholdOnline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchiter->GetOfflineADC() > fJetThresholdOffline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchiter->GetADC() > fJetThresholdOnline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetJetLowBit());
      if(patchiter->GetOfflineADC() > fJetThresholdOffline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetJetLowBit());
    }
    if (patchiter->GetPatchSize() == fBkgPatchSize){
      if(patchiter->GetADC() > fBkgThresholdOnline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetBkgBit());
      if(patchiter->GetOfflineADC() > fBkgThresholdOffline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetBkgBit());
    }
    MakeHLTPatch(*patchiter, *fTriggerPatchDataPtr, offlinebits);
    fTriggerPatchDataPtr = next;
    patchcount++;
    fBufferSize -= sizeof(AliHLTCaloTriggerPatchDataStruct);
  }

  // Do Level0 patches as well
  std::vector<AliEMCALTriggerRawPatch> l0patches = fL0PatchFinder->FindPatches(*fL0Amplitudes, *fADCOfflineValues);
  for(std::vector<AliEMCALTriggerRawPatch>::iterator patchit = l0patches.begin(); patchit != l0patches.end(); ++patchit){
    if(!IsSameTRU(*patchit)) continue;
    Int_t offlinebits = 0;
    if(patchit->GetADC() > fLevel0ThresholdOnline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetLevel0Bit());
    if(patchit->GetOfflineADC() > fLevel0ThresholdOffline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetLevel0Bit());
    if(!offlinebits) continue;
    next = fTriggerPatchDataPtr + 1;
    MakeHLTPatch(*patchit, *fTriggerPatchDataPtr, offlinebits);
    fTriggerPatchDataPtr = next;
    patchcount++;
    fBufferSize -= sizeof(AliHLTCaloTriggerPatchDataStruct);

  }
  return patchcount;
}

void AliHLTEMCALTriggerMaker::Initialise(const AliHLTEMCALGeometry *geo){
  fkGeometryPtr = geo;

  // Allocate
  fADCValues = new AliEMCALTriggerDataGrid<float>;
  fADCValues->Allocate(48, fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2);

  fADCOfflineValues = new AliEMCALTriggerDataGrid<float>;
  fADCOfflineValues->Allocate(48, fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2);

  fL0Amplitudes = new AliEMCALTriggerDataGrid<float>;
  fL0Amplitudes->Allocate(48, fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2);

  fTriggerBitMasks = new AliEMCALTriggerDataGrid<int>;
  fTriggerBitMasks->Allocate(48, fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2);

  fPatchFinder = new AliEMCALTriggerPatchFinder<float>;
  InitializeLevel1PatchFinders(false);
  if(fkGeometryPtr->GetGeometryPtr()->GetNumberOfSuperModules() > 12){
    InitializeLevel1PatchFinders(true);
  }

  fL0PatchFinder = new AliEMCALTriggerPatchFinder<float>;
  InitializeLevel0PatchFinders(false);
  if(fkGeometryPtr->GetGeometryPtr()->GetNumberOfSuperModules() > 12){
	InitializeLevel0PatchFinders(true);
  }

}

void AliHLTEMCALTriggerMaker::InitializeLevel1PatchFinders(Bool_t isDCAL){
  Int_t rowMin = isDCAL ? 64 : 0, rowMax =  isDCAL ? 103 : 63;
  AliEMCALTriggerAlgorithm<float> *gammatrigger = new AliEMCALTriggerAlgorithm<float>(rowMin, rowMax, 0);
  gammatrigger->SetThresholds(0, 0);
  gammatrigger->SetPatchSize(fGammaPatchSize);
  gammatrigger->SetSubregionSize(fGammaSubregionSize);
  fPatchFinder->AddTriggerAlgorithm(gammatrigger);
  AliEMCALTriggerAlgorithm<float> *jettrigger = new AliEMCALTriggerAlgorithm<float>(rowMin, rowMax, 0);
  jettrigger->SetThresholds(0, 0);
  jettrigger->SetPatchSize(fJetPatchSize);
  jettrigger->SetSubregionSize(fJetSubregionSize);
  fPatchFinder->AddTriggerAlgorithm(jettrigger);
  if(fRunBkgAlgorithm){
    AliEMCALTriggerAlgorithm<float> *jetmedian = new AliEMCALTriggerAlgorithm<float>(rowMin, rowMax, 0);
    jetmedian->SetThresholds(-1, -1);
    jetmedian->SetPatchSize(fBkgPatchSize);
    jetmedian->SetSubregionSize(fBkgSubregionSize);
    fPatchFinder->AddTriggerAlgorithm(jetmedian);
  }
}

void AliHLTEMCALTriggerMaker::InitializeLevel0PatchFinders(Bool_t isDCAL){
  AliEMCALTriggerAlgorithm<float> *l0trigger =  new AliEMCALTriggerAlgorithm<float>(isDCAL ? 64 : 0, isDCAL ? 103 : 63, 0);
  l0trigger->SetThresholds(0, 0);
  l0trigger->SetPatchSize(2);
  l0trigger->SetSubregionSize(1);
  fL0PatchFinder->AddTriggerAlgorithm(l0trigger);
}

void AliHLTEMCALTriggerMaker::MakeHLTPatch(const AliEMCALTriggerRawPatch &input, AliHLTCaloTriggerPatchDataStruct &output, UInt_t offlinebits) const {
  output.fCol = input.GetColStart();
  output.fRow = input.GetRowStart();
  output.fSize = input.GetPatchSize();
  output.fADC = input.GetADC();
  output.fOfflineADC = input.GetOfflineADC();
  output.fBitMask = input.GetBitmask() | (*fTriggerBitMasks)(output.fCol, output.fRow) | offlinebits;
}

Bool_t AliHLTEMCALTriggerMaker::IsSameTRU(const AliEMCALTriggerRawPatch &patch) const {
  Bool_t isSameTRU = true;
  int myAbsID = -1, firstAbsID = -1;
  for(int irow = patch.GetRowStart(); irow < patch.GetRowStart() + patch.GetPatchSize(); ++irow){
    for(int icol = patch.GetColStart(); icol < patch.GetColStart() + patch.GetPatchSize(); ++icol){
      fkGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, myAbsID);
      if(firstAbsID < 0) firstAbsID = myAbsID;
      else {
        if(myAbsID != firstAbsID){
          // abs fastor IDs do not match, patch is not fully in same TRU
          isSameTRU = false;
          break;
        }
      }
    }
  }
  return isSameTRU;
}
