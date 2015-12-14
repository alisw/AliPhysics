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
  fADCValues(NULL),
  fADCOfflineValues(NULL),
  fTriggerBitMasks(NULL),
  fTriggerBitConfig(NULL),
  fJetPatchSize(8),
  fJetSubregionSize(4),
  fGammaPatchSize(2),
  fGammaSubregionSize(1),
  fBkgPatchSize(8),
  fBkgSubregionSize(4),
  fBufferSize(0),
  fBkgThresholdOnline(0),
  fBkgThresholdOffline(0)
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
  if(fADCValues) delete fADCValues;
  if(fADCOfflineValues) delete fADCOfflineValues;
  if(fTriggerBitMasks) delete fTriggerBitMasks;
}

void AliHLTEMCALTriggerMaker::ResetADC(){
  fADCValues->Reset();
  fADCOfflineValues->Reset();
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
  for(std::vector<AliEMCALTriggerRawPatch>::iterator patchiter = foundpatches.begin(); patchiter != foundpatches.end(); ++patchiter){
	if(fBufferSize < sizeof(AliHLTCaloTriggerPatchDataStruct)){
		HLTWarning("Buffer exceeded after %d trigger patches", patchcount);
		break;
	}
    AliHLTCaloTriggerPatchDataStruct *next = fTriggerPatchDataPtr + 1;
    // Set offline bits
    Int_t offlinebits = 0;
    if(patchiter->GetPatchSize() == 2){
      if(patchiter->GetADC() > fGammaThresholdOnline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchiter->GetOfflineADC() > fGammaThresholdOffline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchiter->GetADC() > fGammaThresholdOnline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetGammaLowBit());
      if(patchiter->GetOfflineADC() > fGammaThresholdOffline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetGammaLowBit());
    } else if (patchiter->GetPatchSize() == 16){
      if(patchiter->GetADC() > fJetThresholdOnline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchiter->GetOfflineADC() > fJetThresholdOffline[kHighThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchiter->GetADC() > fJetThresholdOnline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetJetLowBit());
      if(patchiter->GetOfflineADC() > fJetThresholdOffline[kLowThreshold]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetJetLowBit());
    } else if (patchiter->GetPatchSize() == 8){
      if(patchiter->GetADC() > fBkgThresholdOnline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetBkgBit());
      if(patchiter->GetOfflineADC() > fBkgThresholdOffline) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetBkgBit());
    }
    MakeHLTPatch(*patchiter, *fTriggerPatchDataPtr, offlinebits);
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

  fTriggerBitMasks = new AliEMCALTriggerDataGrid<int>;
  fTriggerBitMasks->Allocate(48, fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2);

  fPatchFinder = new AliEMCALTriggerPatchFinder<float>;
  InitializeEMCALPatchFinders();
  if(fkGeometryPtr->GetGeometryPtr()->GetNumberOfSuperModules() > 12){
    InitializeDCALPatchFinders();
  }
}

void AliHLTEMCALTriggerMaker::InitializeEMCALPatchFinders(){
  AliEMCALTriggerAlgorithm<float> *gammatrigger = new AliEMCALTriggerAlgorithm<float>(0, 63, 0);
  gammatrigger->SetThresholds(0, 0);
  gammatrigger->SetPatchSize(fGammaPatchSize);
  gammatrigger->SetSubregionSize(fGammaSubregionSize);
  fPatchFinder->AddTriggerAlgorithm(gammatrigger);
  AliEMCALTriggerAlgorithm<float> *jettrigger = new AliEMCALTriggerAlgorithm<float>(0, 63, 0);
  jettrigger->SetThresholds(0, 0);
  jettrigger->SetPatchSize(fJetPatchSize);
  jettrigger->SetSubregionSize(fJetSubregionSize);
  fPatchFinder->AddTriggerAlgorithm(jettrigger);
  AliEMCALTriggerAlgorithm<float> *jetmedian = new AliEMCALTriggerAlgorithm<float>(0, 63, 0);
  jetmedian->SetThresholds(-1, -1);
  jetmedian->SetPatchSize(fBkgPatchSize);
  jetmedian->SetSubregionSize(fBkgSubregionSize);
  fPatchFinder->AddTriggerAlgorithm(jetmedian);
}

void AliHLTEMCALTriggerMaker::InitializeDCALPatchFinders(){
  AliEMCALTriggerAlgorithm<float> *gammatrigger = new AliEMCALTriggerAlgorithm<float>(64, 103, 0);
  gammatrigger->SetThresholds(0, 0);
  gammatrigger->SetPatchSize(fGammaPatchSize);
  gammatrigger->SetSubregionSize(fGammaSubregionSize);
  fPatchFinder->AddTriggerAlgorithm(gammatrigger);
  AliEMCALTriggerAlgorithm<float> *jettrigger = new AliEMCALTriggerAlgorithm<float>(64, 103, 0);
  jettrigger->SetThresholds(0, 0);
  jettrigger->SetPatchSize(fJetPatchSize);
  jettrigger->SetSubregionSize(fJetSubregionSize);
  fPatchFinder->AddTriggerAlgorithm(jettrigger);
  AliEMCALTriggerAlgorithm<float> *jetmedian = new AliEMCALTriggerAlgorithm<float>(64, 103, 0);
  jetmedian->SetThresholds(-1, -1);
  jetmedian->SetPatchSize(fBkgPatchSize);
  jetmedian->SetSubregionSize(fBkgSubregionSize);
  fPatchFinder->AddTriggerAlgorithm(jetmedian);
}

void AliHLTEMCALTriggerMaker::MakeHLTPatch(const AliEMCALTriggerRawPatch &input, AliHLTCaloTriggerPatchDataStruct &output, UInt_t offlinebits) const {
  output.fCol = input.GetColStart();
  output.fRow = input.GetRowStart();
  output.fSize = input.GetPatchSize();
  output.fADC = input.GetADC();
  output.fOfflineADC = input.GetOfflineADC();
  output.fBitMask = input.GetBitmask() | (*fTriggerBitMasks)(output.fCol, output.fRow) | offlinebits;
}
