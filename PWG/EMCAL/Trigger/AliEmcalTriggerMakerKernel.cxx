/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include <iostream>
#include <vector>
#include <cstring>

#include <TArrayI.h>
#include <TObjArray.h>

#include "AliAODCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliEmcalTriggerDataGridAP.h"
#include "AliEmcalTriggerPatchInfoAPV1.h"
#include "AliEmcalTriggerPatchFinderAP.h"
#include "AliEmcalTriggerAlgorithmAP.h"
#include "AliEmcalTriggerRawPatchAP.h"
#include "AliEmcalTriggerMakerKernel.h"
#include "AliEmcalTriggerSetupInfo.h"
#include "AliLog.h"
#include "AliVCaloCells.h"
#include "AliVCaloTrigger.h"
#include "AliVEvent.h"
#include "AliVVZERO.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTriggerMakerKernel)
/// \endcond

AliEmcalTriggerMakerKernel::AliEmcalTriggerMakerKernel():
  TObject(),
  fBadChannels(),
  fTriggerBitConfig(NULL),
  fGeometry(NULL),
  fPatchAmplitudes(NULL),
  fPatchADCSimple(NULL),
  fPatchADC(NULL),
  fLevel0TimeMap(NULL),
  fTriggerBitMap(NULL),
  fPatchFinder(NULL),
  fLevel0PatchFinder(NULL),
  fBkgThreshold(-1),
  fIsMC(kFALSE),
  fDebugLevel(0)
{
  memset(fThresholdConstants, 0, sizeof(Int_t) * 12);
  memset(fL1ThresholdsOffline, 0, sizeof(ULong64_t) * 4);
}

AliEmcalTriggerMakerKernel::~AliEmcalTriggerMakerKernel() {
  delete fPatchAmplitudes;
  delete fPatchADCSimple;
  delete fPatchADC;
  delete fLevel0TimeMap;
  delete fTriggerBitMap;
  delete fPatchFinder;
  delete fLevel0PatchFinder;
  if(fTriggerBitConfig) delete fTriggerBitConfig;
}

void AliEmcalTriggerMakerKernel::Init(){
  fPatchAmplitudes = new AliEmcalTriggerDataGridAP<double>;
  fPatchADCSimple = new AliEmcalTriggerDataGridAP<double>;
  fPatchADC = new AliEmcalTriggerDataGridAP<double>;
  fLevel0TimeMap = new AliEmcalTriggerDataGridAP<char>;
  fTriggerBitMap = new AliEmcalTriggerDataGridAP<int>;

  // Allocate containers for the ADC values
  int nrows = fGeometry->GetNTotalTRU() * 2;
  std::cout << "Allocating channel grid with 48 columns in eta and " << nrows << " rows in phi" << std::endl;
  fPatchAmplitudes->Allocate(48, nrows);
  fPatchADC->Allocate(48, nrows);
  fPatchADCSimple->Allocate(48, nrows);
  fLevel0TimeMap->Allocate(48, nrows);
  fTriggerBitMap->Allocate(48, nrows);

  // Initialize patch finder
  fPatchFinder = new AliEmcalTriggerPatchFinderAP<double>;
  fPatchFinder->AddTriggerAlgorithm(new AliEmcalJetTriggerAlgorithmAP<double>(0, 63, 0));
  fPatchFinder->AddTriggerAlgorithm(new AliEmcalGammaTriggerAlgorithmAP<double>(0, 63, 0));
  fPatchFinder->AddTriggerAlgorithm(new AliEmcalBkgTriggerAlgorithmAP<double>(0, 63, 0));
  if(nrows > 64){
    // Add trigger algorithms for DCAL
    fPatchFinder->AddTriggerAlgorithm(new AliEmcalJetTriggerAlgorithmAP<double>(64, nrows, 0));
    fPatchFinder->AddTriggerAlgorithm(new AliEmcalGammaTriggerAlgorithmAP<double>(64, nrows, 0));
    fPatchFinder->AddTriggerAlgorithm(new AliEmcalBkgTriggerAlgorithmAP<double>(64, nrows, 0));
  }

  fLevel0PatchFinder = new AliEmcalGammaTriggerAlgorithmAP<double>(0, nrows, 0);
}

void AliEmcalTriggerMakerKernel::Reset(){
  fPatchAmplitudes->Reset();
  fPatchADC->Reset();
  fPatchADCSimple->Reset();
  fLevel0TimeMap->Reset();
  fTriggerBitMap->Reset();
  memset(fL1ThresholdsOffline, 0, sizeof(ULong64_t) * 4);
}

void AliEmcalTriggerMakerKernel::ReadTriggerData(AliVCaloTrigger *trigger){
  trigger->Reset();
  Int_t globCol=-1, globRow=-1;
  Int_t adcAmp=-1, bitmap = 0;
  while(trigger->Next()){
    // get position in global 2x2 tower coordinates
    // A0 left bottom (0,0)
    trigger->GetPosition(globCol, globRow);
    // exclude channel completely if it is masked as hot channel
    if(fBadChannels.HasChannel(globCol, globRow)) continue;
    // for some strange reason some ADC amps are initialized in reconstruction
    // as -1, neglect those
    trigger->GetL1TimeSum(adcAmp);
    if (adcAmp>-1) (*fPatchADC)(globCol,globRow) = adcAmp;
    trigger->GetTriggerBits(bitmap);
    (*fTriggerBitMap)(globCol, globRow) = bitmap;

    // Handling for L0 triggers
    // For the ADC value we use fCaloTriggers->GetAmplitude()
    // In data, all patches which have 4 TRUs with proper level0 times are
    // valid trigger patches. Therefore we need to check all neighbors for
    // the level0 times, not only the bottom left. In order to obtain this
    // information, a lookup table with the L0 times for each TRU is created
    Float_t amplitude(0);
    trigger->GetAmplitude(amplitude);
    if(amplitude < 0) amplitude = 0;
    (*fPatchAmplitudes)(globCol,globRow) = amplitude*4; // values are shifted by 2 bits to fit in a 10 bit word (on the hardware side)
    Int_t nl0times(0);
    trigger->GetNL0Times(nl0times);
    if(nl0times){
      TArrayI l0times(nl0times);
      trigger->GetL0Times(l0times.GetArray());
      for(int itime = 0; itime < nl0times; itime++){
        if(l0times[itime] >7 && l0times[itime] < 10){
          (*fLevel0TimeMap)(globCol,globRow) = static_cast<Char_t>(l0times[itime]);
          break;
        }
      }
    }
  }
}

void AliEmcalTriggerMakerKernel::ReadCellData(AliVCaloCells *cells){
  // fill the patch ADCs from cells
  Int_t nCell = cells->GetNumberOfCells();
  for(Int_t iCell = 0; iCell < nCell; ++iCell) {
    // get the cell info, based in index in array
    Short_t cellId = cells->GetCellNumber(iCell);
    Double_t amp = cells->GetAmplitude(iCell);
    // get position
    Int_t absId=-1;
    fGeometry->GetFastORIndexFromCellIndex(cellId, absId);
    Int_t globCol=-1, globRow=-1;
    fGeometry->GetPositionInEMCALFromAbsFastORIndex(absId, globCol, globRow);
    // add
    (*fPatchADCSimple)(globCol,globRow) += amp/EmcalTriggerAP::kEMCL1ADCtoGeV;
  }
}

void AliEmcalTriggerMakerKernel::BuildL1ThresholdsOffline(const AliVVZERO *vzerodata){
  // get the V0 value and compute and set the offline thresholds
  // get V0, compute thresholds and save them as global parameters
  ULong64_t v0S = vzerodata->GetTriggerChargeA() + vzerodata->GetTriggerChargeC();
  for (Int_t i = 0; i < 4; ++i) {
    // A*V0^2/2^32+B*V0/2^16+C
    fL1ThresholdsOffline[i]= ( ((ULong64_t)fThresholdConstants[i][0]) * v0S * v0S ) >> 32;
    fL1ThresholdsOffline[i] += ( ((ULong64_t)fThresholdConstants[i][1]) * v0S ) >> 16;
    fL1ThresholdsOffline[i] += ((ULong64_t)fThresholdConstants[i][2]);
  }
}

TObjArray *AliEmcalTriggerMakerKernel::CreateTriggerPatches(const AliVEvent *inputevent){

  AliEmcalTriggerPatchInfoAPV1 *trigger, *triggerMainJet, *triggerMainGamma, *triggerMainLevel0;
  AliEmcalTriggerPatchInfoAPV1 *triggerMainJetSimple, *triggerMainGammaSimple;

  Double_t vertexpos[3];
  inputevent->GetPrimaryVertex()->GetXYZ(vertexpos);
  TVector3 vertexvec(vertexpos);

  Int_t isMC = fIsMC ? 1 : 0;
  Int_t offset = (1 - isMC) * fTriggerBitConfig->GetTriggerTypesEnd();

  // Create trigger bit masks. They are needed later to remove
  // trigger bits from the trigger bit mask for non-matching patch types
  Int_t jetPatchMask =  1 << fTriggerBitConfig->GetJetHighBit()
      | 1 << fTriggerBitConfig->GetJetLowBit()
      | 1 << (fTriggerBitConfig->GetJetHighBit() + fTriggerBitConfig->GetTriggerTypesEnd())
      | 1 << (fTriggerBitConfig->GetJetLowBit() + fTriggerBitConfig->GetTriggerTypesEnd()),
      gammaPatchMask = 1 << fTriggerBitConfig->GetGammaHighBit()
      | 1 << fTriggerBitConfig->GetGammaLowBit()
      | 1 << (fTriggerBitConfig->GetGammaHighBit() + fTriggerBitConfig->GetTriggerTypesEnd())
      | 1 << (fTriggerBitConfig->GetGammaLowBit() + fTriggerBitConfig->GetTriggerTypesEnd()),
      bkgPatchMask = 1 << fTriggerBitConfig->GetBkgBit();

  std::vector<AliEmcalTriggerRawPatchAP> patches = fPatchFinder->FindPatches(*fPatchADC, *fPatchADCSimple);
  TObjArray *result = new TObjArray(1000);
  for(std::vector<AliEmcalTriggerRawPatchAP>::iterator patchit = patches.begin(); patchit != patches.end(); ++patchit){
    // Apply offline and recalc selection
    Int_t offlinebits = 0;
    if(patchit->GetPatchSize() == 2){
      if(patchit->GetADC() > fL1ThresholdsOffline[1]) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kRecalcOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchit->GetOfflineADC() > fL1ThresholdsOffline[1]) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kOfflineOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchit->GetADC() > fL1ThresholdsOffline[3]) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kRecalcOffset + fTriggerBitConfig->GetGammaLowBit());
      if(patchit->GetOfflineADC() > fL1ThresholdsOffline[3]) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kOfflineOffset + fTriggerBitConfig->GetGammaLowBit());
    } else if (patchit->GetPatchSize() == 16){
      if(patchit->GetADC() > fL1ThresholdsOffline[0]) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kRecalcOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchit->GetOfflineADC() > fL1ThresholdsOffline[0]) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kOfflineOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchit->GetADC() > fL1ThresholdsOffline[2]) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kRecalcOffset + fTriggerBitConfig->GetJetLowBit());
      if(patchit->GetOfflineADC() > fL1ThresholdsOffline[2]) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kOfflineOffset + fTriggerBitConfig->GetJetLowBit());
    } else if (patchit->GetPatchSize() == 8){
      if(patchit->GetADC() > fBkgThreshold) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kRecalcOffset + fTriggerBitConfig->GetBkgBit());
      if(patchit->GetOfflineADC() > fBkgThreshold) SETBIT(offlinebits, AliEmcalTriggerPatchInfoAPV1::kOfflineOffset + fTriggerBitConfig->GetBkgBit());
    }
    // Remove unwanted bits from the online bits (gamma bits from jet patches and vice versa)
    Int_t onlinebits = (*fTriggerBitMap)(patchit->GetColStart(), patchit->GetRowStart());
    if(patchit->GetPatchSize() == 2) onlinebits &= gammaPatchMask;
    else if(patchit->GetPatchSize() == 16) onlinebits &= jetPatchMask;
    else if(patchit->GetPatchSize() == 8) onlinebits &= bkgPatchMask;
    else onlinebits = 0;
    if(!(onlinebits || offlinebits)) continue;
    // convert
    AliEmcalTriggerPatchInfoAPV1 *fullpatch = AliEmcalTriggerPatchInfoAPV1::CreateAndInitialize(patchit->GetColStart(), patchit->GetRowStart(),
        patchit->GetPatchSize(), patchit->GetADC(), patchit->GetOfflineADC(), patchit->GetOfflineADC() * EmcalTriggerAP::kEMCL1ADCtoGeV,
        patchit->GetBitmask() | onlinebits | offlinebits, vertexvec, fGeometry);
    fullpatch->SetTriggerBitConfig(fTriggerBitConfig);
    fullpatch->SetOffSet(offset);
    result->Add(fullpatch);
  }

  // Find Level0 patches
  std::vector<AliEmcalTriggerRawPatchAP> l0patches = fLevel0PatchFinder->FindPatches(*fPatchAmplitudes, *fPatchADCSimple);
  for(std::vector<AliEmcalTriggerRawPatchAP>::iterator patchit = patches.begin(); patchit != patches.end(); ++patchit){
    if(!CheckForL0(patchit->GetColStart(), patchit->GetRowStart())) continue;
    AliEmcalTriggerPatchInfoAPV1 *fullpatch = AliEmcalTriggerPatchInfoAPV1::CreateAndInitialize(patchit->GetColStart(), patchit->GetRowStart(),
        patchit->GetPatchSize(), patchit->GetADC(), patchit->GetOfflineADC(), patchit->GetOfflineADC() * EmcalTriggerAP::kEMCL1ADCtoGeV,
        1, vertexvec, fGeometry);
    fullpatch->SetTriggerBitConfig(fTriggerBitConfig);
    result->Add(fullpatch);

  }
  return result;
}

Bool_t AliEmcalTriggerMakerKernel::CheckForL0(Int_t col, Int_t row) const {
  if(col < 0 || row < 0){
    AliError(Form("Patch outside range [col %d, row %d]", col, row));
    return kFALSE;
  }
  Int_t truref(-1), trumod(-1), absFastor(-1), adc(-1);
  fGeometry->GetAbsFastORIndexFromPositionInEMCAL(col, row, absFastor);
  fGeometry->GetTRUFromAbsFastORIndex(absFastor, truref, adc);
  int nvalid(0);
  const int kNRowsPhi = fGeometry->GetNTotalTRU() * 2;
  for(int ipos = 0; ipos < 2; ipos++){
    if(row + ipos >= kNRowsPhi) continue;    // boundary check
    for(int jpos = 0; jpos < 2; jpos++){
      if(col + jpos >= kColsEta) continue;  // boundary check
      // Check whether we are in the same TRU
      trumod = -1;
      fGeometry->GetAbsFastORIndexFromPositionInEMCAL(col+jpos, row+ipos, absFastor);
      fGeometry->GetTRUFromAbsFastORIndex(absFastor, trumod, adc);
      if(trumod != truref) continue;
      if(col + jpos >= kColsEta) AliError(Form("Boundary error in col [%d, %d + %d]", col + jpos, col, jpos));
      if(row + ipos >= kNRowsPhi) AliError(Form("Boundary error in row [%d, %d + %d]", row + ipos, row, ipos));
      Char_t l0times = (*fLevel0TimeMap)(col + jpos,row + ipos);
      if(l0times > 7 && l0times < 10) nvalid++;
    }
  }
  if (nvalid != 4) return false;
  return true;
}
