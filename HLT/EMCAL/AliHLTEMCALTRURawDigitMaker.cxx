#include "AliCaloRawStreamV3.h"
#include "AliCaloBunchInfo.h"
#include "AliHLTCaloTriggerRawDigitDataStruct.h"
#include "AliHLTEMCALGeometry.h"
#include "AliHLTEMCALTRURawDigitMaker.h"

ClassImp(AliHLTEMCALTRURawDigitMaker)

AliHLTEMCALTRURawDigitMaker::AliHLTEMCALTRURawDigitMaker():
AliHLTLogging(),
fCaloRawStream(),
fGeometryPtr(NULL),
fNRawDigits(0)
{
  for(int ndx = 0; ndx < fgkNRawDigits; ndx++)
    fRawDigitIndex[ndx] = -1;
}

void AliHLTEMCALTRURawDigitMaker::SetRawReader(AliCaloRawStreamV3 *reader){
  fCaloRawStream = reader;
}

AliHLTEMCALTRURawDigitMaker::~AliHLTEMCALTRURawDigitMaker() {
  if(fGeometryPtr) delete fGeometryPtr;
}

void AliHLTEMCALTRURawDigitMaker::Initialize(Int_t runno){
  fGeometryPtr = new AliHLTEMCALGeometry(runno);
}

void AliHLTEMCALTRURawDigitMaker::Add(const std::vector<AliCaloBunchInfo> &bunchlist) {
  Int_t    hwAdd   = fCaloRawStream->GetHWAddress();
  UShort_t iRCU    = fCaloRawStream->GetDDLNumber() % 2; // 0/1
  Int_t    iSM     = fCaloRawStream->GetModule();

  Int_t iTRU = fGeometryPtr->GetGeometryPtr()->GetTriggerMapping()->GetTRUIndexFromOnlineHwAdd(hwAdd,iRCU,iSM);

  if (GetLocalLoggingDefault() & kHLTLogDebug) {
    UShort_t iBranch = ( hwAdd >> 11 ) & 0x1; // 0/1
    HLTDebug("===\n");
    HLTDebug("| Hw Adress: 0x%x => SM# %2d / RCU# %d / Branch# %d / TRU# %2d / ADC# %2d\n",
         hwAdd, fCaloRawStream->GetModule(), iRCU, iBranch, iTRU, fCaloRawStream->GetColumn());
  }

  Int_t idx;

  Int_t timeSamples[15]; memset(timeSamples, 0, sizeof(timeSamples));
  UChar_t nSamples = 0;

  UInt_t iBin   = bunchlist.at(0).GetStartBin();
  Int_t iBunch = 0;

  for (UInt_t i = 0; i < bunchlist.size(); i++) {
    AliCaloBunchInfo bunch = bunchlist.at(i);

    if (iBin > bunch.GetStartBin()) {
      iBin   = bunch.GetStartBin();
      iBunch = i;
    }
    Int_t column = fCaloRawStream->GetColumn();

    if (column < 96) {
      const UShort_t* sig = bunch.GetData();
      Int_t startBin = bunch.GetStartBin();

      for (Int_t iS = 0; iS < bunch.GetLength(); iS++) {
        Int_t time = startBin--;
        Int_t amp  = sig[iS];

        if (amp){
          if(nSamples >= 15){
            HLTError("Buffer for time samples exceeded, not possible to store more");
          } else {
            timeSamples[nSamples++] = ((time << 16) & 0xFF0000) | (amp & 0xFFFF);
            HLTDebug("ADC# %2d / time: %2d amplitude: %d\n", fCaloRawStream->GetColumn(), time, amp);
          }
        }
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

          if (fGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromTRU(iTRU, (fCaloRawStream->GetColumn() - 96) * 10 + j, idx)) {
            SetL0Time(GetRawDigit(idx), iBin);
          }
        }
      }
      iBin--;
    }
  } else {
    if (nSamples && fGeometryPtr->GetGeometryPtr()->GetAbsFastORIndexFromTRU(iTRU, fCaloRawStream->GetColumn(), idx)) {
      SetTimeSamples(GetRawDigit(idx), nSamples, timeSamples);

      if (GetLocalLoggingDefault() & kHLTLogDebug) {
        HLTDebug("| Add TRG digit of id# %4d from TRU# %2d ADC# %2d\n", idx, iTRU, fCaloRawStream->GetColumn());

        //PrintRawDigit(GetRawDigit(idx));
        Int_t iSm, iTru, iEta, iPhi, iD[4], iFor;
        if (fGeometryPtr->GetGeometryPtr()->GetPositionInTRUFromAbsFastORIndex(idx, iTru, iEta, iPhi))
          HLTDebug("| Position => TRU: %2d Eta: %2d Phi: %2d\n", iTru, iEta, iPhi);

        if (fGeometryPtr->GetGeometryPtr()->GetPositionInSMFromAbsFastORIndex(idx, iSm, iEta, iPhi))
          HLTDebug("| Position =>  SM: %2d Eta: %2d Phi: %2d\n", iSm, iEta, iPhi);

        if (fGeometryPtr->GetGeometryPtr()->GetCellIndexFromFastORIndex(idx, iD)) {
          HLTDebug("| tower iDs: ");
          for (Int_t i = 0; i < 4; i++) HLTDebug ("%5d ",iD[i]);
          for (Int_t i = 0; i < 4; i++) {
            if (fGeometryPtr->GetGeometryPtr()->GetFastORIndexFromCellIndex(iD[i], iFor)) {
              HLTDebug("| tower %d to F-OR %d\n",iD[i],iFor);
            }
          }
        }
      }
    }
  }
  /*
  std::cout << "Found TRU  raw digits: " << std::endl;
  for(Int_t idig = 0; idig < fNRawDigits; idig++){
    PrintRawDigit(fRawDigitBuffer[idig]);
  }
  */
}

void AliHLTEMCALTRURawDigitMaker::Reset() {
  for (Int_t i = 0; i < fgkNRawDigits; i++) fRawDigitIndex[i] = -1;
  fNRawDigits = 0;
}

Int_t AliHLTEMCALTRURawDigitMaker::WriteRawDigitsBuffer(AliHLTCaloTriggerRawDigitDataStruct *bufferptr, AliHLTUInt32_t &availableSize) const {
  Int_t outputsize = 0;
  if(availableSize < sizeof(AliHLTCaloTriggerRawDigitDataStruct)){
	  HLTWarning("Not enough space in buffer in order to write digit");
	  return 0;
  }
  for(Int_t idig = 0; idig < fNRawDigits; idig++){
	if(availableSize < sizeof(AliHLTCaloTriggerRawDigitDataStruct)){
		HLTWarning("Buffer exceeded after %d digits", idig);
		break;
	}
    *bufferptr = fRawDigitBuffer[idig];
    bufferptr++;
    outputsize += sizeof(AliHLTCaloTriggerRawDigitDataStruct);
    availableSize -= sizeof(AliHLTCaloTriggerRawDigitDataStruct);
  }
  return outputsize;
}

AliHLTCaloTriggerRawDigitDataStruct &AliHLTEMCALTRURawDigitMaker::GetRawDigit(Int_t index){
  if(fRawDigitIndex[index] >= 0){
    return fRawDigitBuffer[fRawDigitIndex[index]];
  }

  fRawDigitIndex[index] = fNRawDigits;
  new(fRawDigitBuffer + fNRawDigits) AliHLTCaloTriggerRawDigitDataStruct;
  AliHLTCaloTriggerRawDigitDataStruct &dig = fRawDigitBuffer[fNRawDigits];
  InitializeRawDigit(dig);
  SetRawDigitID(dig, index);
  fNRawDigits++;
  return dig;
}
