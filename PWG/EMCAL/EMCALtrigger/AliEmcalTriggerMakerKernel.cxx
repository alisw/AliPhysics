/**************************************************************************************
 * Copyright (C) 2014, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>

#include <TArrayI.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TRandom.h>

#include "AliAODCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerConstants.h"
#include "AliEMCALTriggerDataGrid.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALTriggerPatchFinder.h"
#include "AliEMCALTriggerAlgorithm.h"
#include "AliEMCALTriggerRawPatch.h"
#include "AliEmcalTriggerMakerKernel.h"
#include "AliEmcalTriggerSetupInfo.h"
#include "AliLog.h"
#include "AliVCaloCells.h"
#include "AliVCaloTrigger.h"
#include "AliVEvent.h"
#include "AliVVZERO.h"

ClassImp(AliEmcalTriggerMakerKernel)

AliEmcalTriggerMakerKernel::AliEmcalTriggerMakerKernel():
  TObject(),
  fBadChannels(),
  fOfflineBadChannels(),
  fFastORPedestal(5000),
  fTriggerBitConfig(nullptr),
  fPatchFinder(nullptr),
  fLevel0PatchFinder(nullptr),
  fL0MinTime(7),
  fL0MaxTime(10),
  fApplyL0TimeCut(true),
  fMinCellAmp(0),
  fMinL0FastORAmp(0),
  fMinL1FastORAmp(0),
  fBkgThreshold(-1),
  fL0Threshold(0),
  fIsMC(kFALSE),
  fDebugLevel(0),
  fMinCellAmplitude(0.),
  fApplyOnlineBadChannelsToOffline(kFALSE),
  fApplyOnlineBadChannelsToSmeared(kFALSE),
  fConfigured(kFALSE),
  fSmearModelMean(nullptr),
  fSmearModelSigma(nullptr),
  fSmearThreshold(0.1),
  fScaleShift(0.),
  fScaleMult(1.),
  fConstNoiseFEESmear(0.),
  fMeanNoiseFEESmear(0.),
  fSigmaNoiseFEESmear(0.),
  fAddConstantNoiseFEESmear(false),
  fAddGaussianNoiseFEESmear(false),
  fUseNegPartGaussNoise(false),
  fDoBackgroundSubtraction(false),
  fGeometry(nullptr),
  fPatchAmplitudes(nullptr),
  fPatchADCSimple(nullptr),
  fPatchADC(nullptr),
  fPatchEnergySimpleSmeared(nullptr),
  fLevel0TimeMap(nullptr),
  fTriggerBitMap(nullptr),
  fADCtoGeV(1.)
{
  memset(fThresholdConstants, 0, sizeof(Int_t) * 12);
  memset(fL1ThresholdsOffline, 0, sizeof(ULong64_t) * 4);
  fCellTimeLimits[0] = -10000.;
  fCellTimeLimits[1] = 10000.;
  memset(fRhoValues, 0, sizeof(Double_t) * kNIndRho);
}

AliEmcalTriggerMakerKernel::~AliEmcalTriggerMakerKernel() {
  delete fPatchAmplitudes;
  delete fPatchADCSimple;
  delete fPatchADC;
  delete fPatchEnergySimpleSmeared;
  delete fLevel0TimeMap;
  delete fTriggerBitMap;
  delete fPatchFinder;
  delete fLevel0PatchFinder;
  if(fTriggerBitConfig) delete fTriggerBitConfig;
}

void AliEmcalTriggerMakerKernel::Init(){
  if (!fTriggerBitConfig) {
    AliWarning("Trigger bit configuration was not provided! Assuming new bit configuration (>= 2013).");
    AliEMCALTriggerBitConfig* triggerBitConfig = new AliEMCALTriggerBitConfigNew();
    SetTriggerBitConfig(triggerBitConfig);
  }

  fPatchAmplitudes = new AliEMCALTriggerDataGrid<double>;
  fPatchADCSimple = new AliEMCALTriggerDataGrid<double>;
  fPatchADC = new AliEMCALTriggerDataGrid<double>;
  fLevel0TimeMap = new AliEMCALTriggerDataGrid<char>;
  fTriggerBitMap = new AliEMCALTriggerDataGrid<int>;

  // Allocate containers for the ADC values
  int nrows = fGeometry->GetNTotalTRU() * 2;
  std::cout << "Allocating channel grid with 48 columns in eta and " << nrows << " rows in phi" << std::endl;
  fPatchAmplitudes->Allocate(48, nrows);
  fPatchADC->Allocate(48, nrows);
  fPatchADCSimple->Allocate(48, nrows);
  fLevel0TimeMap->Allocate(48, nrows);
  fTriggerBitMap->Allocate(48, nrows);

  if(fSmearModelMean && fSmearModelSigma){
    // Allocate container for energy smearing (if enabled)
    fPatchEnergySimpleSmeared = new AliEMCALTriggerDataGrid<double>;
    fPatchEnergySimpleSmeared->Allocate(48, nrows);
  }
}

void AliEmcalTriggerMakerKernel::AddL1TriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t bitmask, Int_t patchSize, Int_t subregionSize)
{
  if (!fPatchFinder) fPatchFinder = new AliEMCALTriggerPatchFinder<double>;
  AliEMCALTriggerAlgorithm<double> *trigger = new AliEMCALTriggerAlgorithm<double>(rowmin, rowmax, bitmask);
  trigger->SetPatchSize(patchSize);
  trigger->SetSubregionSize(subregionSize);
  fPatchFinder->AddTriggerAlgorithm(trigger);
}

void AliEmcalTriggerMakerKernel::SetL0TriggerAlgorithm(Int_t rowmin, Int_t rowmax, UInt_t bitmask, Int_t patchSize, Int_t subregionSize)
{
  if (!fLevel0PatchFinder) delete fLevel0PatchFinder;
  fLevel0PatchFinder = new AliEMCALTriggerAlgorithm<double>(rowmin, rowmax, bitmask);
  fLevel0PatchFinder->SetPatchSize(patchSize);
  fLevel0PatchFinder->SetSubregionSize(subregionSize);
}

void AliEmcalTriggerMakerKernel::ConfigureForPbPb2015()
{
  AliEMCALTriggerBitConfig* triggerBitConfig = new AliEMCALTriggerBitConfigNew();
  SetTriggerBitConfig(triggerBitConfig);

  // Initialize patch finder
  if (fPatchFinder) delete fPatchFinder;
  fPatchFinder = new AliEMCALTriggerPatchFinder<double>;

  SetL0TriggerAlgorithm(0, 103, 1<<fTriggerBitConfig->GetLevel0Bit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetGammaHighBit() | 1<<fTriggerBitConfig->GetGammaLowBit(), 2, 1);
  AddL1TriggerAlgorithm(64, 103, 1<<fTriggerBitConfig->GetGammaHighBit() | 1<<fTriggerBitConfig->GetGammaLowBit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetJetHighBit() | 1<<fTriggerBitConfig->GetJetLowBit() | 1<<fTriggerBitConfig->GetBkgBit(), 8, 4);
  AddL1TriggerAlgorithm(64, 103, 1<<fTriggerBitConfig->GetJetHighBit() | 1<<fTriggerBitConfig->GetJetLowBit() | 1<<fTriggerBitConfig->GetBkgBit(), 8, 4);
  fDoBackgroundSubtraction = true;
  fConfigured = true;
}

void AliEmcalTriggerMakerKernel::ConfigureForPP2015()
{
  AliEMCALTriggerBitConfig* triggerBitConfig = new AliEMCALTriggerBitConfigNew();
  SetTriggerBitConfig(triggerBitConfig);

  // Initialize patch finder
  if (fPatchFinder) delete fPatchFinder;
  fPatchFinder = new AliEMCALTriggerPatchFinder<double>;

  SetL0TriggerAlgorithm(0, 103, 1<<fTriggerBitConfig->GetLevel0Bit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetGammaHighBit() | 1<<fTriggerBitConfig->GetGammaLowBit(), 2, 1);
  AddL1TriggerAlgorithm(64, 103, 1<<fTriggerBitConfig->GetGammaHighBit() | 1<<fTriggerBitConfig->GetGammaLowBit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetJetHighBit() | 1<<fTriggerBitConfig->GetJetLowBit(), 16, 4);
  AddL1TriggerAlgorithm(64, 103, 1<<fTriggerBitConfig->GetJetHighBit() | 1<<fTriggerBitConfig->GetJetLowBit(), 8, 4);
  fConfigured = true;
}

void AliEmcalTriggerMakerKernel::ConfigureForPP20158x8()
{
  AliEMCALTriggerBitConfig* triggerBitConfig = new AliEMCALTriggerBitConfigNew();
  SetTriggerBitConfig(triggerBitConfig);

  // Initialize patch finder
  if (fPatchFinder) delete fPatchFinder;
  fPatchFinder = new AliEMCALTriggerPatchFinder<double>;

  SetL0TriggerAlgorithm(0, 103, 1<<fTriggerBitConfig->GetLevel0Bit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetGammaHighBit() | 1<<fTriggerBitConfig->GetGammaLowBit(), 2, 1);
  AddL1TriggerAlgorithm(64, 103, 1<<fTriggerBitConfig->GetGammaHighBit() | 1<<fTriggerBitConfig->GetGammaLowBit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetJetHighBit() | 1<<fTriggerBitConfig->GetJetLowBit(), 8, 4);
  AddL1TriggerAlgorithm(64, 103, 1<<fTriggerBitConfig->GetJetHighBit() | 1<<fTriggerBitConfig->GetJetLowBit(), 8, 4);
  fConfigured = true;
}

void AliEmcalTriggerMakerKernel::ConfigureForPPb2013()
{
  AliEMCALTriggerBitConfig* triggerBitConfig = new AliEMCALTriggerBitConfigNew();
  SetTriggerBitConfig(triggerBitConfig);

  // Initialize patch finder
  if (fPatchFinder) delete fPatchFinder;
  fPatchFinder = new AliEMCALTriggerPatchFinder<double>;

  SetL0TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetLevel0Bit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetGammaHighBit() | 1<<fTriggerBitConfig->GetGammaLowBit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetJetHighBit() | 1<<fTriggerBitConfig->GetJetLowBit(), 16, 4);
  fConfigured = true;
}

void AliEmcalTriggerMakerKernel::ConfigureForPP2012()
{
  AliEMCALTriggerBitConfig* triggerBitConfig = new AliEMCALTriggerBitConfigOld();
  SetTriggerBitConfig(triggerBitConfig);

  // Initialize patch finder
  if (fPatchFinder) delete fPatchFinder;
  fPatchFinder = new AliEMCALTriggerPatchFinder<double>;

  SetL0TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetLevel0Bit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetGammaHighBit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetJetHighBit(), 16, 4);
  fConfigured = true;
}

void AliEmcalTriggerMakerKernel::ConfigureForPbPb2011()
{
  AliEMCALTriggerBitConfig* triggerBitConfig = new AliEMCALTriggerBitConfigOld();
  SetTriggerBitConfig(triggerBitConfig);

  // Initialize patch finder
  if (fPatchFinder) delete fPatchFinder;
  fPatchFinder = new AliEMCALTriggerPatchFinder<double>;

  SetL0TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetLevel0Bit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetGammaHighBit(), 2, 1);
  AddL1TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetJetHighBit(), 16, 4);
  fConfigured = true;
}

void AliEmcalTriggerMakerKernel::ConfigureForPP2011()
{
  AliEMCALTriggerBitConfig* triggerBitConfig = new AliEMCALTriggerBitConfigOld();
  SetTriggerBitConfig(triggerBitConfig);

  // Initialize patch finder
  if (fPatchFinder) delete fPatchFinder;
  fPatchFinder = new AliEMCALTriggerPatchFinder<double>;

  SetL0TriggerAlgorithm(0, 63, 1<<fTriggerBitConfig->GetLevel0Bit(), 2, 1);
  fConfigured = true;
}

void AliEmcalTriggerMakerKernel::ReadOfflineBadChannelFromStream(std::istream& stream)
{
  Short_t absId = 0;

  while (stream.good()) {
    stream >> absId;
    AddOfflineBadChannel(absId);
  }
}

void AliEmcalTriggerMakerKernel::ReadOfflineBadChannelFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadOfflineBadChannelFromStream(file);
}

void AliEmcalTriggerMakerKernel::ReadFastORBadChannelFromStream(std::istream& stream)
{
  Short_t absId = -1;

  while (stream.good()) {
    stream >> absId;
    AddFastORBadChannel(absId);
  }
}

void AliEmcalTriggerMakerKernel::ReadFastORBadChannelFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadFastORBadChannelFromStream(file);
}

void AliEmcalTriggerMakerKernel::SetFastORPedestal(Short_t absId, Float_t ped)
{
  if (absId < 0 || absId >= fFastORPedestal.GetSize()) {
    AliWarning(Form("Abs. ID %d out of range (0,5000)", absId));
    return;
  }
  fFastORPedestal[absId] = ped;
}

void AliEmcalTriggerMakerKernel::ReadFastORPedestalFromStream(std::istream& stream)
{
  Short_t absId = 0;
  Float_t ped = 0;
  while (stream.good()) {
    stream >> ped;
    SetFastORPedestal(absId, ped);
    absId++;
  }
}

void AliEmcalTriggerMakerKernel::ReadFastORPedestalFromFile(const char* fname)
{
  std::ifstream file(fname);
  ReadFastORPedestalFromStream(file);
}

void AliEmcalTriggerMakerKernel::Reset(){
  fPatchAmplitudes->Reset();
  fPatchADC->Reset();
  fPatchADCSimple->Reset();
  fLevel0TimeMap->Reset();
  fTriggerBitMap->Reset();
  if(fPatchEnergySimpleSmeared) fPatchEnergySimpleSmeared->Reset();
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
    Int_t absId = -1;
    fGeometry->GetAbsFastORIndexFromPositionInEMCAL(globCol, globRow, absId);

    // trigger bits can also occur on online masked fastors. Therefore trigger
    // bits are handled before ADC values, and independently whether fastor is
    // masked or not
    bitmap = 0;
    trigger->GetTriggerBits(bitmap);
    if(bitmap){
      // protection against duplicate entries in the AliVCaloTriggers object:
      // the grid is anyhow initialized with 0, so in case of a 0 entry don't overwrite
      // existing entries
      try {
        (*fTriggerBitMap)(globCol, globRow) = bitmap;
      }
      catch (AliEMCALTriggerDataGrid<int>::OutOfBoundsException &e) {
        std::string dirstring = e.GetDirection() == AliEMCALTriggerDataGrid<int>::OutOfBoundsException::kColDir ? "Col" : "Row";
        AliErrorStream() << "Trigger maker task - filling trigger bit grid - index out-of-bounds in " << dirstring << ": " << e.GetIndex() << std::endl;
      }
    }

    // also Level0 times need to be handled without masking of the fastor ...
    // @TODO cross check
    Int_t nl0times(0);
    trigger->GetNL0Times(nl0times);
    if(nl0times){
      TArrayI l0times(nl0times);
      trigger->GetL0Times(l0times.GetArray());
      for(int itime = 0; itime < nl0times; itime++){
        try{
          (*fLevel0TimeMap)(globCol,globRow) = static_cast<Char_t>(l0times[itime]);
          break;
        }
        catch (AliEMCALTriggerDataGrid<char>::OutOfBoundsException &e) {
          std::string dirstring = e.GetDirection() == AliEMCALTriggerDataGrid<char>::OutOfBoundsException::kColDir ? "Col" : "Row";
          AliErrorStream() << "Trigger maker task - filling trigger bit grid - index out-of-bounds in " << dirstring << ": " << e.GetIndex() << std::endl;
        }
      }
    }

    // exclude channel completely if it is masked as hot channel
    if (fBadChannels.find(absId) != fBadChannels.end()){
      AliDebugStream(1) << "Found ADC for masked fastor " << absId << ", rejecting" << std::endl;
      continue;
    }
    // for some strange reason some ADC amps are initialized in reconstruction
    // as -1, neglect those
    trigger->GetL1TimeSum(adcAmp);
    if (adcAmp < 0) adcAmp = 0;

    if (adcAmp >= std::max(fMinL1FastORAmp, 1)) {
      // protection against duplicate entries: the grid is anyhow
      // initialized with 0, so in case of an entry with negative or
      // 0 ADC time sum don't overwrite the existing one
      try {
        (*fPatchADC)(globCol,globRow) = adcAmp;
      }
      catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e) {
        std::string dirstring = e.GetDirection() == AliEMCALTriggerDataGrid<double>::OutOfBoundsException::kColDir ? "Col" : "Row";
        AliErrorStream() << "Trigger maker task - filling trigger bit grid - index out-of-bounds in " << dirstring << ": " << e.GetIndex() << std::endl;
      }
    }

    // Handling for L0 triggers
    // For the ADC value we use fCaloTriggers->GetAmplitude()
    // In data, all patches which have 4 TRUs with proper level0 times are
    // valid trigger patches. Therefore we need to check all neighbors for
    // the level0 times, not only the bottom left. In order to obtain this
    // information, a lookup table with the L0 times for each TRU is created
    Float_t amplitude(0);
    trigger->GetAmplitude(amplitude);
    amplitude *= 4; // values are shifted by 2 bits to fit in a 10 bit word (on the hardware side)
    amplitude -= fFastORPedestal[absId];
    if(amplitude < 0) amplitude = 0;
    if (amplitude > std::max(0., double(fMinL0FastORAmp))) {
      // protection against duplicate entries: the grid is anyhow
      // initialized with 0, so in case of an entry with negative or
      // 0 ADC time sum don't overwrite the existing one
      try{
        (*fPatchAmplitudes)(globCol,globRow) = amplitude;
      }
      catch (AliEMCALTriggerDataGrid<int>::OutOfBoundsException &e) {
        std::string dirstring = e.GetDirection() == AliEMCALTriggerDataGrid<int>::OutOfBoundsException::kColDir ? "Col" : "Row";
        AliErrorStream() << "Trigger maker task - filling trigger bit grid - index out-of-bounds in " << dirstring << ": " << e.GetIndex() << std::endl;
      }
    }
  }

  // Reading of the rho values (only PbPb)
  if(fDoBackgroundSubtraction) {
    AliInfoStream() << "Reading median values: EMCAL " << fRhoValues[kIndRhoEMCAL] << ", DCAL " << fRhoValues[kIndRhoDCAL] << std::endl;
    fRhoValues[kIndRhoEMCAL] = trigger->GetMedian(0);     // EMCAL STU at position 0
    fRhoValues[kIndRhoDCAL] = trigger->GetMedian(1);      // DCAL STU at position 1
  }
}

void AliEmcalTriggerMakerKernel::ReadCellData(AliVCaloCells *cells){
  // fill the patch ADCs from cells
  Int_t nCell = cells->GetNumberOfCells();
  for(Int_t iCell = 0; iCell < nCell; ++iCell) {
    // get the cell info, based in index in array
    Short_t cellId = cells->GetCellNumber(iCell);

    // Check bad channel map
    if (fOfflineBadChannels.find(cellId) != fOfflineBadChannels.end()) {
      AliDebugStream(1) << "Cell " << cellId << " masked as bad channel, rejecting." << std::endl;
      continue;
    }

    Double_t amp = cells->GetAmplitude(iCell),
             celltime = cells->GetTime(iCell);
    if(celltime < fCellTimeLimits[0] || celltime > fCellTimeLimits[1]) continue;
    amp *= fScaleMult;
    if(amp < fMinCellAmplitude) continue;
    if(fScaleShift) amp += fScaleShift;
     
    amp = TMath::Max(amp, 0.);      // never go negative in energy
    // get position
    Int_t absId=-1;
    fGeometry->GetFastORIndexFromCellIndex(cellId, absId);
    if(fApplyOnlineBadChannelsToOffline){
      // Exclude FEE amplitudes from cells which are within a TRU which is masked at
      // online level. Using this the online acceptance can be applied to offline
      // patches as well.
      if(fBadChannels.find(absId) != fBadChannels.end()){
        AliDebugStream(1) << "Cell " << cellId << " corresponding to masked fastor " << absId << ", rejecting." << std::endl;
        continue;
      }
    }
    Int_t globCol=-1, globRow=-1;
    fGeometry->GetPositionInEMCALFromAbsFastORIndex(absId, globCol, globRow);
    // add
    amp /= fADCtoGeV;
    try {
      if (amp >= fMinCellAmp) (*fPatchADCSimple)(globCol,globRow) += amp;
    }
    catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e) {
    }
  }

  // Apply energy smearing (if enabled)
  if(fPatchEnergySimpleSmeared){
    AliDebugStream(1) << "Trigger Maker: Apply energy smearing" << std::endl;
    for(int icol = 0; icol < fPatchADCSimple->GetNumberOfCols(); icol++){
      for(int irow = 0; irow < fPatchADCSimple->GetNumberOfRows(); irow++){
        bool doChannel = true;
        if(irow >= 64 && irow < 100 && icol >= 16 && icol < 32) {
          // in PHOS Hole, continue
          continue;
        }
        if(fApplyOnlineBadChannelsToSmeared) {
          int absFastor = -1;
          fGeometry->GetAbsFastORIndexFromPositionInEMCAL(icol, irow, absFastor);
          if(absFastor > -1) {
            if(fBadChannels.find(absFastor) != fBadChannels.end()){
              AliDebugStream(1) << "In smearing, FastOR " << absFastor << " masked, rejecting." << std::endl;
              doChannel = false;
            }
          }
        }
        if(!doChannel) continue;
        double energyorig = (*fPatchADCSimple)(icol, irow) * fADCtoGeV;          // Apply smearing in GeV
        double energysmear = energyorig;
        if(energyorig > fSmearThreshold){
          double mean = fSmearModelMean->Eval(energyorig), sigma = fSmearModelSigma->Eval(energyorig);
          energysmear =  gRandom->Gaus(mean, sigma);
          if(energysmear < 0) energysmear = 0;      // only accept positive or 0 energy values
          AliDebugStream(1) << "Original energy " << energyorig << ", mean " << mean << ", sigma " << sigma << ", smeared " << energysmear << std::endl;
        }
        (*fPatchEnergySimpleSmeared)(icol, irow) += energysmear;
        // check whether to handle noise
        if(fAddConstantNoiseFEESmear) {
          (*fPatchEnergySimpleSmeared)(icol, irow) += fConstNoiseFEESmear;
        }
        if(fAddGaussianNoiseFEESmear) {
          // Accept also the negative part of the gaussian to simulate underfluctuations
          double noisevalue = gRandom->Gaus(fMeanNoiseFEESmear, fSigmaNoiseFEESmear);
          if(noisevalue < 0. && !fUseNegPartGaussNoise) 
            noisevalue = 0.;
          (*fPatchEnergySimpleSmeared)(icol, irow) += noisevalue;
        }
        // Truncate to 0
        (*fPatchEnergySimpleSmeared)(icol, irow) = TMath::Max((*fPatchEnergySimpleSmeared)(icol, irow), 0.);
      }
    }
    AliDebugStream(1) << "Smearing done" << std::endl;
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

void AliEmcalTriggerMakerKernel::CreateTriggerPatches(const AliVEvent *inputevent, std::vector<AliEMCALTriggerPatchInfo> &outputcont, Bool_t useL0amp){
  //std::cout << "Finding trigger patches" << std::endl;
  //AliEMCALTriggerPatchInfo *trigger, *triggerMainJet, *triggerMainGamma, *triggerMainLevel0;
  //AliEMCALTriggerPatchInfo *triggerMainJetSimple, *triggerMainGammaSimple;

  if (useL0amp) {
    fADCtoGeV = EMCALTrigger::kEMCL0ADCtoGeV_AP;
  }
  else {
    fADCtoGeV = EMCALTrigger::kEMCL1ADCtoGeV;
  }

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
      //l0PatchMask = 1 << fTriggerBitConfig->GetLevel0Bit();

  std::vector<AliEMCALTriggerRawPatch> patches;
  if (fPatchFinder) {
    if (useL0amp) {
      patches = fPatchFinder->FindPatches(*fPatchAmplitudes, *fPatchADCSimple);
    }
    else {
      patches = fPatchFinder->FindPatches(*fPatchADC, *fPatchADCSimple);
    }
  }
  outputcont.clear();
  for(std::vector<AliEMCALTriggerRawPatch>::iterator patchit = patches.begin(); patchit != patches.end(); ++patchit){
    // Apply offline and recalc selection
    // Remove unwanted bits from the online bits (gamma bits from jet patches and vice versa)
    Int_t offlinebits = 0, onlinebits = (*fTriggerBitMap)(patchit->GetColStart(), patchit->GetRowStart());
    if(HasPHOSOverlap(*patchit)) continue;
    if(IsGammaPatch(*patchit)){
      if(patchit->GetADC() > fL1ThresholdsOffline[1]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchit->GetOfflineADC() > fL1ThresholdsOffline[1]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetGammaHighBit());
      if(patchit->GetADC() > fL1ThresholdsOffline[3]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetGammaLowBit());
      if(patchit->GetOfflineADC() > fL1ThresholdsOffline[3]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetGammaLowBit());
      onlinebits &= gammaPatchMask;
    }
    if (IsJetPatch(*patchit)){
      if(patchit->GetADC() > fL1ThresholdsOffline[0]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchit->GetOfflineADC() > fL1ThresholdsOffline[0]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetJetHighBit());
      if(patchit->GetADC() > fL1ThresholdsOffline[2]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetJetLowBit());
      if(patchit->GetOfflineADC() > fL1ThresholdsOffline[2]) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetJetLowBit());
      onlinebits &= jetPatchMask;
    }
    if (IsBkgPatch(*patchit)){
      if(patchit->GetADC() > fBkgThreshold) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetBkgBit());
      if(patchit->GetOfflineADC() > fBkgThreshold) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetBkgBit());
      onlinebits &= bkgPatchMask;
    }
    if(fDoBackgroundSubtraction) {
      double area = TMath::Power(static_cast<double>(patchit->GetPatchSize())/8., 2);
      double rhoval = (patchit->GetRowStart() >= 64) ? fRhoValues[kIndRhoDCAL] : fRhoValues[kIndRhoEMCAL];  // Rho values are for a detector measured in the opposite arm
      AliDebugStream(1) << "Subtracting background for area " << area << ": " << rhoval  << " -> " << (area * rhoval) << std::endl;
      patchit->SetADC(patchit->GetADC() - area * rhoval);
    }
    // convert
    AliEMCALTriggerPatchInfo fullpatch;
    fullpatch.Initialize(patchit->GetColStart(), patchit->GetRowStart(),
        patchit->GetPatchSize(), patchit->GetADC(), patchit->GetOfflineADC(), patchit->GetOfflineADC() * fADCtoGeV,
        onlinebits | offlinebits, vertexvec, fGeometry);
    fullpatch.SetTriggerBitConfig(fTriggerBitConfig);
    fullpatch.SetOffSet(offset);
    if(fPatchEnergySimpleSmeared){
      // Add smeared energy
      double energysmear = 0;
      for(int icol = 0; icol < fullpatch.GetPatchSize(); icol++){
        for(int irow = 0; irow < fullpatch.GetPatchSize(); irow++){
          energysmear += (*fPatchEnergySimpleSmeared)(fullpatch.GetColStart() + icol, fullpatch.GetRowStart() + irow);
        }
      }
      AliDebugStream(1) << "Patch size(" << fullpatch.GetPatchSize() <<") energy " << fullpatch.GetPatchE() << " smeared " << energysmear << std::endl;
      fullpatch.SetSmearedEnergy(energysmear);
    }
    outputcont.push_back(fullpatch);
  }

  // Find Level0 patches
  std::vector<AliEMCALTriggerRawPatch> l0patches;
  if (fLevel0PatchFinder) l0patches = fLevel0PatchFinder->FindPatches(*fPatchAmplitudes, *fPatchADCSimple);
  for(std::vector<AliEMCALTriggerRawPatch>::iterator patchit = l0patches.begin(); patchit != l0patches.end(); ++patchit){
    Int_t offlinebits = 0, onlinebits = 0;
    if(HasPHOSOverlap(*patchit)) continue;
    ELevel0TriggerStatus_t L0status = CheckForL0(patchit->GetColStart(), patchit->GetRowStart());
    if (L0status == kNotLevel0) continue;
    if (L0status == kLevel0Fired) SETBIT(onlinebits, fTriggerBitConfig->GetLevel0Bit());
    if (patchit->GetADC() > fL0Threshold) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kRecalcOffset + fTriggerBitConfig->GetLevel0Bit());
    if (patchit->GetOfflineADC() > fL0Threshold) SETBIT(offlinebits, AliEMCALTriggerPatchInfo::kOfflineOffset + fTriggerBitConfig->GetLevel0Bit());

    AliEMCALTriggerPatchInfo fullpatch;
    fullpatch.Initialize(patchit->GetColStart(), patchit->GetRowStart(),
        patchit->GetPatchSize(), patchit->GetADC(), patchit->GetOfflineADC(), patchit->GetOfflineADC() * fADCtoGeV,
        onlinebits | offlinebits, vertexvec, fGeometry);
    fullpatch.SetTriggerBitConfig(fTriggerBitConfig);
    if(fPatchEnergySimpleSmeared){
      // Add smeared energy
      double energysmear = 0;
      for(int icol = 0; icol < fullpatch.GetPatchSize(); icol++){
        for(int irow = 0; irow < fullpatch.GetPatchSize(); irow++){
          energysmear += (*fPatchEnergySimpleSmeared)(fullpatch.GetColStart() + icol, fullpatch.GetRowStart() + irow);
        }
      }
      fullpatch.SetSmearedEnergy(energysmear);
    }
    outputcont.push_back(fullpatch);
  }
  // std::cout << "Finished finding trigger patches" << std::endl;
}

double AliEmcalTriggerMakerKernel::GetL0TriggerChannelAmplitude(Int_t col, Int_t row) const{
  double amp = 0;
  try {
    amp = (*fPatchAmplitudes)(col, row);
  } catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e) {

  }
  return amp;
}

double AliEmcalTriggerMakerKernel::GetTriggerChannelADC(Int_t col, Int_t row) const{
  double adc = 0;
  try {
    adc = (*fPatchADC)(col, row);
  } catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e) {

  }
  return adc;
}

double AliEmcalTriggerMakerKernel::GetTriggerChannelEnergyRough(Int_t col, Int_t row) const{
  double adc = 0;
  try {
    adc = (*fPatchADC)(col, row) * EMCALTrigger::kEMCL1ADCtoGeV;
  } catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e) {

  }
  return adc;
}

double AliEmcalTriggerMakerKernel::GetTriggerChannelADCSimple(Int_t col, Int_t row) const{
  double adc = 0;
  try {
    adc = (*fPatchADCSimple)(col, row);
  } catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e) {

  }
  return adc;
}

double AliEmcalTriggerMakerKernel::GetTriggerChannelEnergy(Int_t col, Int_t row) const {
  double adc = 0;
  try {
    adc = (*fPatchADCSimple)(col, row) * fADCtoGeV;
  } catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e) {

  }
  return adc;
}

double AliEmcalTriggerMakerKernel::GetTriggerChannelEnergySmeared(Int_t col, Int_t row) const {
  double adc = 0;
  if(fPatchEnergySimpleSmeared){
	  try {
		  adc = (*fPatchEnergySimpleSmeared)(col, row);
	  } catch (AliEMCALTriggerDataGrid<double>::OutOfBoundsException &e) {

	  }
  }
  return adc;
}

double AliEmcalTriggerMakerKernel::GetDataGridDimensionRows() const{
  return fPatchADC->GetNumberOfRows();
}

AliEmcalTriggerMakerKernel::ELevel0TriggerStatus_t AliEmcalTriggerMakerKernel::CheckForL0(Int_t col, Int_t row) const {
  ELevel0TriggerStatus_t result = kLevel0Candidate;

  if(col < 0 || row < 0){
    AliError(Form("Patch outside range [col %d, row %d]", col, row));
    return kNotLevel0;
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
      if(trumod != truref) {
        result = kNotLevel0;
        return result;
      }
      if(col + jpos >= kColsEta) AliError(Form("Boundary error in col [%d, %d + %d]", col + jpos, col, jpos));
      if(row + ipos >= kNRowsPhi) AliError(Form("Boundary error in row [%d, %d + %d]", row + ipos, row, ipos));
      if(fApplyL0TimeCut) {
        Char_t l0times = (*fLevel0TimeMap)(col + jpos,row + ipos);
        if(l0times > fL0MinTime && l0times < fL0MaxTime) nvalid++;
      } else {
        nvalid++;
      }
    }
  }
  if (nvalid == 4) result = kLevel0Fired;
  return result;
}

void AliEmcalTriggerMakerKernel::ClearFastORBadChannels(){
  fBadChannels.clear();
}

void AliEmcalTriggerMakerKernel::ClearOfflineBadChannels() {
  fOfflineBadChannels.clear();
}

Bool_t AliEmcalTriggerMakerKernel::IsGammaPatch(const AliEMCALTriggerRawPatch &patch) const {
  ULong_t bitmask = patch.GetBitmask(), testmask = 1 << fTriggerBitConfig->GetGammaHighBit() | 1 << fTriggerBitConfig->GetGammaLowBit();
  return bitmask & testmask;
}

Bool_t AliEmcalTriggerMakerKernel::IsJetPatch(const AliEMCALTriggerRawPatch &patch) const {
  ULong_t bitmask = patch.GetBitmask(), testmask = 1 << fTriggerBitConfig->GetJetHighBit() | 1 << fTriggerBitConfig->GetJetLowBit();
  return bitmask & testmask;
}

Bool_t AliEmcalTriggerMakerKernel::IsBkgPatch(const AliEMCALTriggerRawPatch &patch) const {
  ULong_t bitmask = patch.GetBitmask(), testmask = 1 << fTriggerBitConfig->GetBkgBit();
  return bitmask & testmask;
}

void AliEmcalTriggerMakerKernel::SetTriggerBitConfig(const AliEMCALTriggerBitConfig *const config) {
  if (config == fTriggerBitConfig) return;
  if (fTriggerBitConfig) delete fTriggerBitConfig;
  fTriggerBitConfig = config;
}

bool AliEmcalTriggerMakerKernel::HasPHOSOverlap(const AliEMCALTriggerRawPatch &patch) const {
  const int kEtaMinPhos = 16, kEtaMaxPhos = 31, kPhiMinPhos = 64, kPhiMaxPhos = 99;
  if(patch.GetRowStart() + patch.GetPatchSize() -1 < kPhiMinPhos) return false;   // EMCAL Patch
  if(patch.GetRowStart() > kPhiMaxPhos) return false;         // DCAL 1/3 supermodule
  if(patch.GetColStart() + patch.GetPatchSize() -1 < kEtaMinPhos) return false;
  if(patch.GetColStart() > kEtaMaxPhos) return false;
  return true;
}
