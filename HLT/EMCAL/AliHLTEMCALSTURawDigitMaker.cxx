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
#include "AliHLTEMCALSTURawDigitMaker.h"
#include "AliCaloConstants.h"
#include "AliCaloRawAnalyzerFactory.h"
#include "AliCaloRawAnalyzerFakeALTRO.h"
#include "AliCaloRawStreamV3.h"
#include "AliDAQ.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEMCALTriggerData.h"
#include "AliEMCALTriggerRawDigit.h"
#include "AliEMCALTriggerSTURawStream.h"
#include "AliHLTEMCALGeometry.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawEvent.h"
#include "AliRawReader.h"
#include "AliRawVEquipment.h"

ClassImp(AliHLTEMCALSTURawDigitMaker)

AliHLTEMCALSTURawDigitMaker::AliHLTEMCALSTURawDigitMaker() :
TObject(),
AliHLTLogging(),
fkGeometryPtr(NULL),
fDCSConfigSTU(NULL),
fTriggerData(NULL),
fNRawDigits(0)
{
  fDCSConfigSTU = new AliEMCALTriggerSTUDCSConfig;
  fTriggerData = new AliEMCALTriggerData;
  for (Int_t i=0; i<fgkNRawDigits; i++) fRawDigitIndex[i] = -1;
}

AliHLTEMCALSTURawDigitMaker::~AliHLTEMCALSTURawDigitMaker() {
  if(fDCSConfigSTU) delete fDCSConfigSTU;
  if(fTriggerData) delete fTriggerData;
}

void AliHLTEMCALSTURawDigitMaker::ProcessSTUStream(AliEMCALTriggerSTURawStream *stustream){

  //std::cout << "Start Processing STU data" << std::endl;
  HLTDebug("Start post processing the raw digit maker");
  Int_t idx;

  AliHLTCaloTriggerRawDigitDataStruct *hltdig(NULL);

  //std::cout << "Before check for payload" << std::endl;
  if (stustream && stustream->ReadPayLoad()) {
    //std::cout << "STU stream has payload" << std::endl;
    fTriggerData->SetL1DataDecoded(1);

    for (int i = 0; i < 2; i++) {
      fTriggerData->SetL1GammaThreshold(i, stustream->GetL1GammaThreshold(i));
      fTriggerData->SetL1JetThreshold(  i, stustream->GetL1JetThreshold(i)  );
    }

    Int_t v0[2] = { static_cast<Int_t>(stustream->GetV0A()),  static_cast<Int_t>(stustream->GetV0C())};

    // Modify DCS config from STU payload content
    for(Int_t i = 0; i < 3; i++){
      for(Int_t j = 0; j < 2; j++){
        fDCSConfigSTU->SetG(i, j, stustream->GetG(i, j));
        fDCSConfigSTU->SetJ(i, j, stustream->GetJ(i, j));
      }
    }
    fDCSConfigSTU->SetRawData(stustream->GetRawData());
    fDCSConfigSTU->SetRegion(stustream->GetRegionEnable());
    fDCSConfigSTU->SetFw(stustream->GetFwVersion());

    fTriggerData->SetL1FrameMask(stustream->GetFrameReceived());
    fTriggerData->SetL1V0(v0);
    Int_t type[15] =
    {
      static_cast<Int_t>(stustream->GetG(0, 0)),
      static_cast<Int_t>(stustream->GetG(1, 0)),
      static_cast<Int_t>(stustream->GetG(2, 0)),
      static_cast<Int_t>(stustream->GetJ(0, 0)),
      static_cast<Int_t>(stustream->GetJ(1, 0)),
      static_cast<Int_t>(stustream->GetJ(2, 0)),
      static_cast<Int_t>(stustream->GetG(0, 1)),
      static_cast<Int_t>(stustream->GetG(1, 1)),
      static_cast<Int_t>(stustream->GetG(2, 1)),
      static_cast<Int_t>(stustream->GetJ(0, 1)),
      static_cast<Int_t>(stustream->GetJ(1, 1)),
      static_cast<Int_t>(stustream->GetJ(2, 1)),
      static_cast<Int_t>(stustream->GetRawData()),
      static_cast<Int_t>(stustream->GetRegionEnable()),
      static_cast<Int_t>(stustream->GetFwVersion())
    };
    fTriggerData->SetL1TriggerType(type);

    fTriggerData->SetL1RawData(stustream->GetRawData());

    Int_t iTRU, x, y;

    TVector2 sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch;
    fDCSConfigSTU->GetSegmentation(sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch);

    //std::cout << "Before check for raw data" << std::endl;
    if (stustream->GetRawData()) {
      //std::cout << "STU has raw data" << std::endl;
      HLTDebug("| STU => TRU raw data are there!\n");

      Int_t nTRU = fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU();
      for (Int_t i = 0; i < nTRU; i++) {
        iTRU = fkGeometryPtr->GetGeometryPtr()->GetTRUIndexFromSTUIndex(i, 0);

        UInt_t adc[96]; for (Int_t j = 0; j < 96; j++) adc[j] = 0;

        stustream->GetADC(i, adc);

        for (Int_t j = 0; j < 96; j++) {
          if (adc[j] <= 0) continue;
          HLTDebug("| STU => TRU# %2d raw data: ADC# %2d: %d\n", iTRU, j, adc[j]);
          fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromTRU(iTRU, j, idx);
          SetL1TimeSum(GetRawDigit(idx), adc[j]);
        }
      }
    }

    // List of patches in EMCal coordinate system

    for (Int_t i = 0; i < stustream->GetNL0GammaPatch(); i++) {
      stustream->GetL0GammaPatch(i, iTRU, x);

      iTRU = fkGeometryPtr->GetGeometryPtr()->GetTRUIndexFromSTUIndex(iTRU, 0);

      const Int_t sizePatchL0 = 4;

      HLTDebug("| STU => Found L0 patch id: %2d in TRU# %2d\n", x, iTRU);

      Int_t idFastOR[4];
      for (Int_t j = 0; j < 4; j++) idFastOR[j] = -1;

      if (fkGeometryPtr->GetGeometryPtr()->GetFastORIndexFromL0Index(iTRU, x, idFastOR, sizePatchL0)) {
        idx = idFastOR[1];

        Int_t px, py;
        if (fkGeometryPtr->GetGeometryPtr()->GetPositionInEMCALFromAbsFastORIndex(idx, px, py)){
          HLTDebug("| STU => Add L0 patch at (%2d , %2d)\n", px, py);
          SetTriggerBit(GetRawDigit(idx), kL0, 1);
        }
      }
    }

    for (int ithr = 0; ithr < 2; ithr++) {
      for (Int_t i = 0; i < stustream->GetNL1GammaPatch(ithr); i++) {
        if (stustream->GetL1GammaPatch(i, ithr, iTRU, x, y)) { // col (0..23), row (0..3)
          iTRU = fkGeometryPtr->GetGeometryPtr()->GetTRUIndexFromSTUIndex(iTRU, 0);

          HLTDebug("| STU => Found L1 gamma patch at (%2d , %2d) in TRU# %2d\n", x, y, iTRU);

          Int_t vx = 23 - x, vy = y + 4 * int(iTRU / 2); // Position in EMCal frame

          if (iTRU % 2) vx += 24; // C side

          vx = vx - int(sizeL1gsubr.X()) * int(sizeL1gpatch.X()) + 1;

          if (vx >= 0 && vy < 63) {
            if (fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromPositionInEMCAL(vx, vy, idx)) {
              HLTDebug("| STU => Add L1 gamma [%d] patch at (%2d , %2d)\n", ithr, vx, vy);
              SetTriggerBit(GetRawDigit(idx), kL1GammaHigh + ithr, 1);
            }
          }
        }
      }

      for (Int_t i = 0; i < stustream->GetNL1JetPatch(ithr); i++) {
        if (stustream->GetL1JetPatch(i, ithr, x, y)) { // col (0,15), row (0,11)
          HLTDebug("| STU => Found L1 jet [%d] patch at (%2d , %2d)\n", ithr, x, y);

          Int_t ix = int(sizeL1jsubr.X()) * (11 - y - int(sizeL1jpatch.X()) + 1);
          Int_t iy = int(sizeL1jsubr.Y()) * (15 - x - int(sizeL1jpatch.Y()) + 1);

          if (ix >= 0 && iy >= 0) {
            if (fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromPositionInEMCAL(ix, iy, idx)) {
              HLTDebug("| STU => Add L1 jet patch at (%2d , %2d)\n", ix, iy);
              SetTriggerBit(GetRawDigit(idx), kL1JetHigh + ithr, 1);
            }
          }
        }
      }
    }
  }
  /*
  for(int idig = 0; idig < fNRawDigits; idig++){
    PrintRawDigit(fRawDigitBuffer[idig]);
  }
  */
}

Int_t AliHLTEMCALSTURawDigitMaker::WriteRawDigitsBuffer(AliHLTCaloTriggerRawDigitDataStruct *bufferptr) const {
  Int_t outputsize = 0;
  for(Int_t idig = 0; idig < fNRawDigits; idig++){
    *bufferptr = fRawDigitBuffer[idig];
    bufferptr++;
    outputsize += sizeof(AliHLTCaloTriggerRawDigitDataStruct);
  }
  return outputsize;
}

void AliHLTEMCALSTURawDigitMaker::Reset() {
  for (Int_t i = 0; i < fgkNRawDigits; i++) fRawDigitIndex[i] = -1;
  fNRawDigits = 0;
  fTriggerData->Reset();
}

AliHLTCaloTriggerRawDigitDataStruct &AliHLTEMCALSTURawDigitMaker::GetRawDigit(Int_t index){
  if(fRawDigitIndex[index] >= 0){
    return fRawDigitBuffer[fRawDigitIndex[index]];
  }

  fRawDigitIndex[index] = fNRawDigits;
  AliHLTCaloTriggerRawDigitDataStruct &dig = fRawDigitBuffer[fNRawDigits];
  InitializeRawDigit(dig);
  SetRawDigitID(dig, index);
  fNRawDigits++;
  return dig;
}
