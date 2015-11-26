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
#include "AliHLTEMCALTriggerRawDigitMaker.h"
#include "AliRawEquipmentHeader.h"
#include "AliRawEvent.h"
#include "AliRawReader.h"
#include "AliRawVEquipment.h"

ClassImp(AliHLTEMCALTriggerRawDigitMaker)

const Int_t AliHLTEMCALTriggerRawDigitMaker::fgkSTUEqId = 4652;
const Int_t AliHLTEMCALTriggerRawDigitMaker::fgkNRawDigits = 5952;

AliHLTEMCALTriggerRawDigitMaker::AliHLTEMCALTriggerRawDigitMaker() :
TObject(),
AliHLTLogging(),
fkGeometryPtr(NULL),
fRawReader(NULL),
fCaloRawStream(NULL),
fSTURawStream(NULL),
fRawAnalyzer(NULL),
fDCSConfigSTU(NULL),
fTriggerData(NULL)
{
  fDCSConfigSTU = new AliEMCALTriggerSTUDCSConfig;
  fRawAnalyzer =  (AliCaloRawAnalyzerFakeALTRO*)AliCaloRawAnalyzerFactory::CreateAnalyzer(kFakeAltro);

  for (Int_t i=0; i<fgkNRawDigits; i++) fRawDigitIndex[i] = -1;
}

AliHLTEMCALTriggerRawDigitMaker::~AliHLTEMCALTriggerRawDigitMaker() {
  if(fRawAnalyzer) delete fRawAnalyzer;
  if(fDCSConfigSTU) delete fDCSConfigSTU;
}

void AliHLTEMCALTriggerRawDigitMaker::SetIO(AliRawReader* reader, AliCaloRawStreamV3& in, AliEMCALTriggerSTURawStream& inSTU, AliEMCALTriggerData* data) {
  fRawReader     = reader;
  fCaloRawStream = &in;
  fSTURawStream  = &inSTU;
  fTriggerData   = data;
}

void AliHLTEMCALTriggerRawDigitMaker::Add(const std::vector<AliCaloBunchInfo> &bunchlist) {
  Int_t    hwAdd   = fCaloRawStream->GetHWAddress();
  UShort_t iRCU    = fCaloRawStream->GetDDLNumber() % 2; // 0/1
  Int_t    iSM     = fCaloRawStream->GetModule();

  Int_t iTRU = fkGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetTRUIndexFromOnlineHwAdd(hwAdd,iRCU,iSM);

  if (GetLocalLoggingDefault() & kHLTLogDebug) {
    UShort_t iBranch = ( hwAdd >> 11 ) & 0x1; // 0/1
    HLTDebug("===\n");
    HLTDebug("| Hw Adress: 0x%x => SM# %2d / RCU# %d / Branch# %d / TRU# %2d / ADC# %2d\n",
         hwAdd, fCaloRawStream->GetModule(), iRCU, iBranch, iTRU, fCaloRawStream->GetColumn());
  }

  Int_t idx;

  Int_t timeSamples[256]; memset(timeSamples, 0, sizeof(Int_t) * 256);
  UChar_t nSamples = 0;

  UInt_t iBin   = bunchlist.at(0).GetStartBin();
  Int_t iBunch = 0;

  for (UInt_t i = 0; i < bunchlist.size(); i++) {
    AliCaloBunchInfo bunch = bunchlist.at(i);

    if (iBin > bunch.GetStartBin()) {
      iBin   = bunch.GetStartBin();
      iBunch = i;
    }

    if (fCaloRawStream->GetColumn() < 96) {
      const UShort_t* sig = bunch.GetData();
      Int_t startBin = bunch.GetStartBin();

      for (Int_t iS = 0; iS < bunch.GetLength(); iS++) {
        Int_t time = startBin--;
        Int_t amp  = sig[iS];

        if (amp) timeSamples[nSamples++] = ((time << 16) & 0xFF0000) | (amp & 0xFFFF);
        HLTDebug("ADC# %2d / time: %2d amplitude: %d\n", fCaloRawStream->GetColumn(), time, amp);
      }
    }
  }

  if (fCaloRawStream->GetColumn() > 95 && fCaloRawStream->GetColumn() < 106) {
    Int_t nBits = (fCaloRawStream->GetColumn() == 105) ? 6 : 10;
    const UShort_t* sig = bunchlist.at(iBunch).GetData();
    HLTDebug("| L0 id in F-ALTRO => bunch length is: %d\n", bunchlist.at(iBunch).GetLength());

    for (Int_t i = 0; i < bunchlist.at(iBunch).GetLength(); i++) {
      HLTDebug("| sig[%3d]: %x\n",i,sig[i]);

      for (Int_t j = 0; j < nBits; j++) {
        if (sig[i] & ( 1 << j )) {
          HLTDebug("| Add L0 patch index in TRU# %2d position %2d\n",iTRU,(fCaloRawStream->GetColumn() - 96) * 10 + j);

          if (fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromTRU(iTRU, (fCaloRawStream->GetColumn() - 96) * 10 + j, idx)) {
            SetL0Time(GetRawDigit(idx), iBin);
          }
        }
      }

      if (fCaloRawStream->GetColumn() == 105 && (sig[i] & (1 << 6))) {
        fTriggerData->SetL0Trigger(1, iTRU, 1);
        HLTDebug("=======TRU# %2d has issued a L0\n",iTRU);
      }
      iBin--;
    }
  } else {
    if (nSamples && fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromTRU(iTRU, fCaloRawStream->GetColumn(), idx)) {
      SetTimeSamples(GetRawDigit(idx), nSamples, timeSamples);

      if (GetLocalLoggingDefault() & kHLTLogDebug) {
        HLTDebug("| Add TRG digit of id# %4d from TRU# %2d ADC# %2d\n", idx, iTRU, fCaloRawStream->GetColumn());

        PrintRawDigit(GetRawDigit(idx));
        Int_t iSm, iTru, iEta, iPhi, iD[4], iFor;
        if (fkGeometryPtr->GetGeometryPtr()->GetPositionInTRUFromAbsFastORIndex(idx, iTru, iEta, iPhi))
          HLTDebug("| Position => TRU: %2d Eta: %2d Phi: %2d\n", iTru, iEta, iPhi);

        if (fkGeometryPtr->GetGeometryPtr()->GetPositionInSMFromAbsFastORIndex(idx, iSm, iEta, iPhi))
          HLTDebug("| Position =>  SM: %2d Eta: %2d Phi: %2d\n", iSm, iEta, iPhi);

        if (fkGeometryPtr->GetGeometryPtr()->GetCellIndexFromFastORIndex(idx, iD)) {
          HLTDebug("| tower iDs: ");
          for (Int_t i = 0; i < 4; i++) HLTDebug ("%5d ",iD[i]);
          for (Int_t i = 0; i < 4; i++) {
            if (fkGeometryPtr->GetGeometryPtr()->GetFastORIndexFromCellIndex(iD[i], iFor)) {
              HLTDebug("| tower %d to F-OR %d\n",iD[i],iFor);
            }
          }
        }
      }
    }
  }
}

void AliHLTEMCALTriggerRawDigitMaker::PostProcess(){

  HLTDebug("Start post processing the raw digit maker");
  Int_t idx;

  AliHLTEMCALTriggerRawDigitDataStruct *hltdig(NULL);

  fRawReader->Reset();
  fRawReader->Select("EMCAL",AliDAQ::GetFirstSTUDDL());

  AliRawVEvent *event = (AliRawEvent*)fRawReader->GetEvent();
  if (!event) {
    HLTError("STU DDL# %d not available!", AliDAQ::GetFirstSTUDDL());
    return;
  }

  Bool_t isSTUin = kFALSE;

  Int_t nSubEv = fRawReader->GetEvent()->GetNSubEvents();

  for ( Int_t iSubEv=0; iSubEv<nSubEv; iSubEv++) {
    AliRawVEvent *subEv = ((AliRawEvent*)fRawReader->GetEvent())->GetSubEvent(iSubEv);
    if ( !subEv ) continue;

    for (Int_t iEquip = 0; iEquip < subEv->GetNEquipments(); iEquip++) {
      Int_t eqId = subEv->GetEquipment(iEquip)->GetEquipmentHeader()->GetId();
      if (eqId == fgkSTUEqId) isSTUin = kTRUE;
    }
  }

  fRawReader->Reset();

  if (isSTUin && fSTURawStream && fSTURawStream->ReadPayLoad()) {
    fTriggerData->SetL1DataDecoded(1);

    for (int i = 0; i < 2; i++) {
      fTriggerData->SetL1GammaThreshold(i, fSTURawStream->GetL1GammaThreshold(i));
      fTriggerData->SetL1JetThreshold(  i, fSTURawStream->GetL1JetThreshold(i)  );
    }

    Int_t v0[2] = { static_cast<Int_t>(fSTURawStream->GetV0A()),  static_cast<Int_t>(fSTURawStream->GetV0C())};

    // Modify DCS config from STU payload content
    for(Int_t i = 0; i < 3; i++){
      for(Int_t j = 0; j < 2; j++){
        fDCSConfigSTU->SetG(i, j, fSTURawStream->GetG(i, j));
        fDCSConfigSTU->SetJ(i, j, fSTURawStream->GetJ(i, j));
      }
    }
    fDCSConfigSTU->SetRawData(fSTURawStream->GetRawData());
    fDCSConfigSTU->SetRegion(fSTURawStream->GetRegionEnable());
    fDCSConfigSTU->SetFw(fSTURawStream->GetFwVersion());

    fTriggerData->SetL1FrameMask(fSTURawStream->GetFrameReceived());
    fTriggerData->SetL1V0(v0);
    Int_t type[15] =
    {
      static_cast<Int_t>(fSTURawStream->GetG(0, 0)),
      static_cast<Int_t>(fSTURawStream->GetG(1, 0)),
      static_cast<Int_t>(fSTURawStream->GetG(2, 0)),
      static_cast<Int_t>(fSTURawStream->GetJ(0, 0)),
      static_cast<Int_t>(fSTURawStream->GetJ(1, 0)),
      static_cast<Int_t>(fSTURawStream->GetJ(2, 0)),
      static_cast<Int_t>(fSTURawStream->GetG(0, 1)),
      static_cast<Int_t>(fSTURawStream->GetG(1, 1)),
      static_cast<Int_t>(fSTURawStream->GetG(2, 1)),
      static_cast<Int_t>(fSTURawStream->GetJ(0, 1)),
      static_cast<Int_t>(fSTURawStream->GetJ(1, 1)),
      static_cast<Int_t>(fSTURawStream->GetJ(2, 1)),
      static_cast<Int_t>(fSTURawStream->GetRawData()),
      static_cast<Int_t>(fSTURawStream->GetRegionEnable()),
      static_cast<Int_t>(fSTURawStream->GetFwVersion())
    };
    fTriggerData->SetL1TriggerType(type);

    fTriggerData->SetL1RawData(fSTURawStream->GetRawData());

    Int_t iTRU, x, y;

    TVector2 sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch;
    fDCSConfigSTU->GetSegmentation(sizeL1gsubr, sizeL1gpatch, sizeL1jsubr, sizeL1jpatch);

    if (fSTURawStream->GetRawData()) {
      HLTDebug("| STU => TRU raw data are there!\n");

      Int_t nTRU = 32;//fGeometry->GetNTotalTRU();
      for (Int_t i = 0; i < nTRU; i++) {
        iTRU = fkGeometryPtr->GetGeometryPtr()->GetTRUIndexFromSTUIndex(i, 0);

        UInt_t adc[96]; for (Int_t j = 0; j < 96; j++) adc[j] = 0;

        fSTURawStream->GetADC(i, adc);

        for (Int_t j = 0; j < 96; j++) {
          if (adc[j] <= 0) continue;
          HLTDebug("| STU => TRU# %2d raw data: ADC# %2d: %d\n", iTRU, j, adc[j]);
          fkGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromTRU(iTRU, j, idx);
          SetL1TimeSum(GetRawDigit(idx), adc[j]);
        }
      }
    }

    // List of patches in EMCal coordinate system

    for (Int_t i = 0; i < fSTURawStream->GetNL0GammaPatch(); i++) {
      fSTURawStream->GetL0GammaPatch(i, iTRU, x);

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
      for (Int_t i = 0; i < fSTURawStream->GetNL1GammaPatch(ithr); i++) {
        if (fSTURawStream->GetL1GammaPatch(i, ithr, iTRU, x, y)) { // col (0..23), row (0..3)
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

      for (Int_t i = 0; i < fSTURawStream->GetNL1JetPatch(ithr); i++) {
        if (fSTURawStream->GetL1JetPatch(i, ithr, x, y)) { // col (0,15), row (0,11)
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
}

void AliHLTEMCALTriggerRawDigitMaker::Reset() {
  fRawDigitVector.clear();
  for (Int_t i = 0; i < fgkNRawDigits; i++) fRawDigitIndex[i] = -1;
}

AliHLTEMCALTriggerRawDigitDataStruct &AliHLTEMCALTriggerRawDigitMaker::GetRawDigit(Int_t index){
  if(fRawDigitIndex[index] > 0) return fRawDigitVector[fRawDigitIndex[index]];

  fRawDigitIndex[index] = fRawDigitVector.size();
  fRawDigitVector.push_back(AliHLTEMCALTriggerRawDigitDataStruct());
  AliHLTEMCALTriggerRawDigitDataStruct &dig = fRawDigitVector[fRawDigitIndex[index]];
  InitializeRawDigit(dig);
  SetRawDigitID(dig, index);
  return dig;
}
