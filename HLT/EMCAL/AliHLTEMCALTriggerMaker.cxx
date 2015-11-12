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

#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerMapping.h"
#include "AliEMCALTriggerTypes.h"
#include "AliEMCALTriggerAlgorithm.h"
#include "AliEMCALTriggerDataGrid.h"
#include "AliEMCALTriggerPatchFinder.h"
#include "AliEMCALTriggerRawPatch.h"

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
  fOrigin(kOriginCELLS)
{
}

AliHLTEMCALTriggerMaker::~AliHLTEMCALTriggerMaker() {
  if(fPatchFinder) delete fPatchFinder;
  if(fADCValues) delete fADCValues;
}

void AliHLTEMCALTriggerMaker::ResetADC(){
  fADCValues->Reset();
}

void AliHLTEMCALTriggerMaker::AddDigit(const AliHLTCaloDigitDataStruct *digit){
  /*
   * TODO Crosscheck
   */
  Int_t fastorIndex;
  fkGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetFastORIndexFromCellIndex(fkGeometryPtr->GetGeometryPtr()->GetAbsCellIdFromCellIndexes(digit->fModule, digit->fX, digit->fZ), fastorIndex);
  int globCol, globRow;
  fkGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetPositionInEMCALFromAbsFastORIndex(fastorIndex, globCol, globRow);
  (*fADCValues)(globCol, globRow) += digit->fAmplitude;
}

void AliHLTEMCALTriggerMaker::SetADC(Int_t col, Int_t row, Float_t adc){
  (*fADCValues)(col, row) = adc;
}

Int_t AliHLTEMCALTriggerMaker::FindPatches(){
  /*
  if(availableSize < sizeof(AliHLTCaloTriggerPatchDataStruct)){
    HLTError("Not enough space to write new trigger patches");
    return -1;
  }
  */
  //AliHLTUInt32_t mysize = availableSize;
  std::vector<AliEMCALTriggerRawPatch> foundpatches = fPatchFinder->FindPatches(*fADCValues);
  Int_t patchcount = 0;
  for(std::vector<AliEMCALTriggerRawPatch>::iterator patchiter = foundpatches.begin(); patchiter != foundpatches.end(); ++patchiter){
    /*
    if(mysize < sizeof(AliHLTCaloTriggerPatchDataStruct)){
      HLTError("Not enough space - some patches will be lost");
      return patchcount;
    }
    */
    AliHLTCaloTriggerPatchDataStruct *next = fTriggerPatchDataPtr + 1;
    MakeHLTPatch(*patchiter, *fTriggerPatchDataPtr, fOrigin);
    fTriggerPatchDataPtr = next;
    patchcount++;
  }
  return patchcount;
}

void AliHLTEMCALTriggerMaker::Initialise(const AliHLTEMCALGeometry *geo){
  fkGeometryPtr = geo;

  // Allocate
  fADCValues = new AliEMCALTriggerDataGrid<float>;
  fADCValues->Allocate(48, fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU() * 2);

  fPatchFinder = new AliEMCALTriggerPatchFinder<float>;
  InitializeEMCALPatchFinders();
  if(fkGeometryPtr->GetGeometryPtr()->GetNumberOfSuperModules() > 12){
    InitializeDCALPatchFinders();
  }
}

void AliHLTEMCALTriggerMaker::InitializeEMCALPatchFinders(){
  AliEMCALTriggerAlgorithm<float> *gammatrigger = new AliEMCALGammaTriggerAlgorithm<float>(0, 63, 1<<kL1GammaHigh);
  fPatchFinder->AddTriggerAlgorithm(gammatrigger);
  AliEMCALTriggerAlgorithm<float> *jettrigger = new AliEMCALJetTriggerAlgorithm<float>(0, 63, 1<<kL1JetHigh);
  fPatchFinder->AddTriggerAlgorithm(jettrigger);
}

void AliHLTEMCALTriggerMaker::InitializeDCALPatchFinders(){
  AliEMCALTriggerAlgorithm<float> *gammatrigger = new AliEMCALGammaTriggerAlgorithm<float>(64, 103, 1<<kL1GammaHigh);
  fPatchFinder->AddTriggerAlgorithm(gammatrigger);
  AliEMCALTriggerAlgorithm<float> *jettrigger = new AliEMCALJetTriggerAlgorithm<float>(64, 103, 1<<kL1JetHigh);
  fPatchFinder->AddTriggerAlgorithm(jettrigger);
}

void AliHLTEMCALTriggerMaker::MakeHLTPatch(const AliEMCALTriggerRawPatch &input, AliHLTCaloTriggerPatchDataStruct &output, UChar_t origin) const {
  output.fCol = input.GetColStart();
  output.fRow = input.GetRowStart();
  output.fSize = input.GetPatchSize();
  output.fADC = input.GetADC();
  output.fBitMask = input.GetBitmask();
  output.fOrigin = origin;
}
