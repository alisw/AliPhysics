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

void AliHLTEMCALSTURawDigitMaker::ProcessSTUStream(AliEMCALTriggerSTURawStream *stustream, Int_t detector){

  HLTDebug("Start post processing the raw digit maker");
  Int_t idx;

  AliHLTCaloTriggerRawDigitDataStruct *hltdig(NULL);


  if (stustream && stustream->ReadPayLoad()) {
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

    Int_t iTRU, jTRU, x, y;

    TVector2 sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch;
    fDCSConfigSTU->GetSegmentation(sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch);

    if (stustream->GetRawData()) {
      HLTDebug("| STU => TRU raw data are there!\n");

      //Int_t nTRU = fkGeometryPtr->GetGeometryPtr()->GetNTotalTRU();
      Int_t nTRU = detector == 0 ? 32 : 14;
      for (Int_t i = 0; i < nTRU; i++) {
        iTRU = fkGeometryPtr->GetGeometryPtr()->GetTRUIndexFromSTUIndex(i, detector);

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

      //if(iTRU >= 32) iTRU -= 32;
      iTRU = fkGeometryPtr->GetGeometryPtr()->GetTRUIndexFromSTUIndex(iTRU, detector);

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

    Int_t vx, vy, lphi;
    for (int ithr = 0; ithr < 2; ithr++) {
      for (Int_t i = 0; i < stustream->GetNL1GammaPatch(ithr); i++) {
        if (stustream->GetL1GammaPatch(i, ithr, iTRU, x, y)) { // col (0..23), row (0..3)

          if (fkGeometryPtr->GetGeometryPtr()->GetTriggerMappingVersion() == 1) {
            // Run 1
            iTRU = fkGeometryPtr->GetGeometryPtr()->GetTRUIndexFromSTUIndex(iTRU, detector);

            HLTDebug("| STU => Found L1 gamma patch at (%2d , %2d) in TRU# %2d\n", x, y, iTRU);

            vx = 23 - x;
            vy = y + 4 * int(iTRU / 2); // Position in EMCal frame

            if (iTRU % 2) vx += 24; // C side

            vx = vx - int(sizeL1gsubr.X()) * int(sizeL1gpatch.X()) + 1;
            lphi = 64;
            if (vx >= 0 && vy < lphi) {
              if (fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromPositionInEMCAL(vx, vy, idx)) {
            	HLTDebug("| STU => Add L1 gamma [%d] patch at (%2d , %2d)\n", ithr, vx, vy, index);
            	SetTriggerBit(GetRawDigit(idx), kL1GammaHigh + ithr, 1);
              }
            }
          } else {
            // Run 2
            fkGeometryPtr->GetGeometryPtr()->GetTRUFromSTU(iTRU, x, y, jTRU, vx, vy, detector);
            fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromPositionInTRU(jTRU, vx, vy, idx);
            SetTriggerBit(GetRawDigit(idx), kL1GammaHigh + ithr, 1);
          }
        }
      }

      for (Int_t i = 0; i < stustream->GetNL1JetPatch(ithr); i++) {
        if (stustream->GetL1JetPatch(i, ithr, x, y)) { // col (0,15), row (0,11)
          HLTDebug("| STU => Found L1 jet [%d] patch at (%2d , %2d)\n", ithr, x, y);

          if (fkGeometryPtr->GetGeometryPtr()->GetTriggerMappingVersion() == 1) {
            vx = 11 - y - int(sizeL1jpatch.X()) + 1;
            vy = 15 - x - int(sizeL1jpatch.Y()) + 1;
          }
          else {
            vx = y;
            vy = x;
          }
          vx *= int(sizeL1jsubr.X());
          vy *= int(sizeL1jsubr.Y());

          if (vx >= 0 && vy >= 0) {
            if (fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromPositionInEMCAL(vx, vy, idx)) {
              HLTDebug("| STU => Add L1 jet patch at (%2d , %2d)\n", ix, iy);
              SetTriggerBit(GetRawDigit(idx), kL1JetHigh + ithr, 1);
            }
          }
        }
      }
    }

    if (detector == 1) {
      UInt_t sregion[36] = {0};
      stustream->GetPHOSSubregion(sregion);
      for (int isr=0;isr<36;isr++) {
        if (fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromPHOSSubregion(isr, idx)) {
          SetL1SubRegion(GetRawDigit(idx), sregion[isr]);
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

Int_t AliHLTEMCALSTURawDigitMaker::WriteRawDigitsBuffer(AliHLTCaloTriggerRawDigitDataStruct *bufferptr, AliHLTUInt32_t &availableSize) const {
  Int_t outputsize = 0;
  for(Int_t idig = 0; idig < fNRawDigits; idig++){
	if(availableSize < sizeof(AliHLTCaloTriggerRawDigitDataStruct)){
		HLTWarning("Buffer exceeded after %d digits", idig);
		break;
	}
	if (fRawDigitBuffer[idig].fID < 0 || fRawDigitBuffer[idig].fID >= fgkNRawDigits)
	{
		HLTWarning("Invalid TRU Index in Output %d", fRawDigitBuffer[idig].fID);
		continue;
	}
    *bufferptr = fRawDigitBuffer[idig];
    bufferptr++;
    outputsize += sizeof(AliHLTCaloTriggerRawDigitDataStruct);
    availableSize -= sizeof(AliHLTCaloTriggerRawDigitDataStruct);
  }
  return outputsize;
}

void AliHLTEMCALSTURawDigitMaker::Reset() {
  for (Int_t i = 0; i < fgkNRawDigits; i++) fRawDigitIndex[i] = -1;
  fNRawDigits = 0;
  fTriggerData->Reset();
}

AliHLTCaloTriggerRawDigitDataStruct &AliHLTEMCALSTURawDigitMaker::GetRawDigit(Int_t index){
  if (index < 0){
    HLTWarning("Invalid STU Index %d", index);
    AliHLTCaloTriggerRawDigitDataStruct &dig = fRawDigitBuffer[fgkNRawDigits - 1];
    InitializeRawDigit(dig);
    return dig;
  }
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
