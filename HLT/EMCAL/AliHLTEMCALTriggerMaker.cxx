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
  fLevel0TimeMap(NULL),
  fTRUIndexMap(NULL),
  fTriggerBitConfig(NULL),
  fJetPatchSize(8),
  fJetSubregionSize(4),
  fGammaPatchSize(2),
  fGammaSubregionSize(1),
  fBkgPatchSize(8),
  fBkgSubregionSize(4),
  fL0MinTime(7),
  fL0MaxTime(10),
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
  if(fLevel0TimeMap) delete fLevel0TimeMap;
}

void AliHLTEMCALTriggerMaker::ResetADC(){
  fADCValues->Reset();
  fADCOfflineValues->Reset();
  fL0Amplitudes->Reset();
  fTriggerBitMasks->Reset();
  fLevel0TimeMap->Reset();
}

void AliHLTEMCALTriggerMaker::AddDigit(const AliHLTCaloDigitDataStruct *digit){
  /*
   * TODO Crosscheck
   */
  Int_t fastorIndex;
  fkGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetFastORIndexFromCellIndex(fkGeometryPtr->GetGeometryPtr()->GetAbsCellIdFromCellIndexes(digit->fModule, digit->fX, digit->fZ), fastorIndex);
  int globCol, globRow;
  fkGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetPositionInEMCALFromAbsFastORIndex(fastorIndex, globCol, globRow);
  (*fADCOfflineValues)(globCol, globRow) += digit->fEnergy/EMCALTrigger::kEMCL1ADCtoGeV;
}

void AliHLTEMCALTriggerMaker::SetADC(Int_t col, Int_t row, Float_t adc){
  (*fADCValues)(col, row) = adc;
}

void AliHLTEMCALTriggerMaker::SetL0Amplitude(Int_t col, Int_t row, Float_t amp){
  (*fL0Amplitudes)(col, row) = amp*4; // to compensate for the last two bits that are chopped away in the hardware chain
}

void AliHLTEMCALTriggerMaker::SetL0Time(Int_t col, Int_t row, UChar_t time){
  (*fLevel0TimeMap)(col, row) = time;
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
    UInt_t offlinebits = 0, onlinebitmask = 0;
    if(patchiter->GetPatchSize() == fGammaPatchSize){
      onlinebitmask = 1 << kL1GammaHigh | 1 << kL1GammaLow | 1 << (kL1GammaHigh + kTriggerTypeEnd) | 1 << (kL1GammaLow + kTriggerTypeEnd);
      if(patchiter->GetADC() > fGammaThresholdOnline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchiter->GetOfflineADC() > fGammaThresholdOffline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchiter->GetADC() > fGammaThresholdOnline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetGammaLowBit());
      if(patchiter->GetOfflineADC() > fGammaThresholdOffline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetGammaLowBit());
    }
    if (patchiter->GetPatchSize() == fJetPatchSize){
      onlinebitmask = 1 << kL1JetHigh | 1 << kL1JetLow | 1 << (kL1JetHigh + kTriggerTypeEnd) | 1 << (kL1JetLow + kTriggerTypeEnd);
      if(patchiter->GetADC() > fJetThresholdOnline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchiter->GetOfflineADC() > fJetThresholdOffline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchiter->GetADC() > fJetThresholdOnline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetJetLowBit());
      if(patchiter->GetOfflineADC() > fJetThresholdOffline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetJetLowBit());
    }
    if (patchiter->GetPatchSize() == fBkgPatchSize){
      if(patchiter->GetADC() > fBkgThresholdOnline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetBkgBit());
      if(patchiter->GetOfflineADC() > fBkgThresholdOffline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetBkgBit());
    }
    MakeHLTPatch(*patchiter, *fTriggerPatchDataPtr, offlinebits, onlinebitmask, 0);
    fTriggerPatchDataPtr = next;
    patchcount++;
    fBufferSize -= sizeof(AliHLTCaloTriggerPatchDataStruct);
  }

  // Do Level0 patches as well
  std::vector<AliEMCALTriggerRawPatch> l0patches = fL0PatchFinder->FindPatches(*fL0Amplitudes, *fADCOfflineValues);
  for(std::vector<AliEMCALTriggerRawPatch>::iterator patchit = l0patches.begin(); patchit != l0patches.end(); ++patchit){
    if(fBufferSize < sizeof(AliHLTCaloTriggerPatchDataStruct)){
      HLTWarning("Buffer exceeded after %d trigger patches", patchcount);
      break;
    }
    ELevel0TriggerStatus_t L0trigger = CheckForL0(patchit->GetColStart(), patchit->GetRowStart());
    // Check that it is a valid L0 patch candidate, i.e. all FastORs in the same TRU
    if (L0trigger == kNotLevel0) continue;
    Int_t onlinebits(0), offlinebits(0);
    if (L0trigger == kLevel0Fired) SETBIT(onlinebits, fTriggerBitConfig->GetLevel0Bit()+fTriggerBitConfig->GetTriggerTypesEnd());
    // No requirement on L0 time for offline and recalc patches, only they have to be in the same TRU
    if (patchit->GetADC() > fLevel0ThresholdOnline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetLevel0Bit());
    if (patchit->GetOfflineADC() > fLevel0ThresholdOffline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetLevel0Bit());
    if (!(offlinebits || onlinebits))  continue;
    next = fTriggerPatchDataPtr + 1;
    MakeHLTPatch(*patchit, *fTriggerPatchDataPtr, offlinebits, 0, onlinebits);
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

  fLevel0TimeMap = new AliEMCALTriggerDataGrid<unsigned char>;
  fLevel0TimeMap->Allocate(48, fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2);

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

  InitializeLookupTables();
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

void AliHLTEMCALTriggerMaker::MakeHLTPatch(const AliEMCALTriggerRawPatch &input, AliHLTCaloTriggerPatchDataStruct &output, UInt_t offlinebits, UInt_t onlinebitmask, UInt_t level0bits) const {
  output.fCol = input.GetColStart();
  output.fRow = input.GetRowStart();
  output.fSize = input.GetPatchSize();
  output.fADC = input.GetADC();
  output.fOfflineADC = input.GetOfflineADC();
  Int_t onlinebits = (*fTriggerBitMasks)(output.fCol, output.fRow) & onlinebitmask;
  output.fBitMask = input.GetBitmask() | onlinebits | offlinebits | level0bits;
}

void AliHLTEMCALTriggerMaker::InitializeLookupTables(){
  fTRUIndexMap = new AliEMCALTriggerDataGrid<int>;
  fTRUIndexMap->Allocate(48, fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2);
  int absFastor, truid, adc;
  for(int icol = 0; icol < 48; icol++){
    for(int irow = 0; irow < fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2; irow++){
      fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, absFastor);
      fkGeometryPtr->GetGeometryPtr()->GetTRUFromAbsFastORIndex(absFastor, truid, adc);
      (*fTRUIndexMap)(icol, irow) = truid;
    }
  }
}

AliHLTEMCALTriggerMaker::ELevel0TriggerStatus_t AliHLTEMCALTriggerMaker::CheckForL0(Int_t col, Int_t row) const {
  ELevel0TriggerStatus_t result = kLevel0Candidate;

  if(col < 0 || row < 0){
    AliError(Form("Patch outside range [col %d, row %d]", col, row));
    return kNotLevel0;
  }
  Int_t truref = (*fTRUIndexMap)(col, row), trumod(-1);
  const Int_t kColsEta = 48;
  const int kNRowsPhi = fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2;
  for(int ipos = 0; ipos < 2; ipos++){
    if(row + ipos >= kNRowsPhi) continue;    // boundary check
    for(int jpos = 0; jpos < 2; jpos++){
      if(col + jpos >= kColsEta) continue;  // boundary check
      // Check whether we are in the same TRU
      trumod = (*fTRUIndexMap)(col + jpos, row + ipos);
      if(trumod != truref) {
        result = kNotLevel0;
        return result;
      }
      if(col + jpos >= kColsEta) AliError(Form("Boundary error in col [%d, %d + %d]", col + jpos, col, jpos));
      if(row + ipos >= kNRowsPhi) AliError(Form("Boundary error in row [%d, %d + %d]", row + ipos, row, ipos));
    }
  }
  if (result == kLevel0Candidate) {
    Char_t l0times = (*fLevel0TimeMap)(col,row);
    if(l0times > fL0MinTime && l0times < fL0MaxTime) result = kLevel0Fired;
  }
  return result;
}
