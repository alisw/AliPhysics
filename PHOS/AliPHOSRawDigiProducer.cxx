#include "TClonesArray.h"

#include "AliPHOSRawDigiProducer.h"
#include "AliPHOSRawDecoder.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"

ClassImp(AliPHOSRawDigiProducer)

void AliPHOSRawDigiProducer::MakeDigits(TClonesArray *digits, AliPHOSRawDecoder* decoder) 
{
  digits->Clear();
 
  AliPHOSGeometry* geo = AliPHOSGeometry::GetInstance();
  if(!geo) geo = AliPHOSGeometry::GetInstance("IHEP");

  Int_t    iDigit   = 0 ;
  Double_t time     = 0. ;
  Int_t    iOldDigit;
  Bool_t   seen,lowGainFlag;
  Int_t relId[4], absId =0;

  while (decoder->NextDigit()) {

    lowGainFlag = decoder->IsLowGain();
    time = decoder->GetTime();

    relId[0] = decoder->GetModule();
    relId[1] = 0;
    relId[2] = decoder->GetRow();
    relId[3] = decoder->GetColumn();
    geo->RelToAbsNumbering(relId, absId);

    // Add low gain digit only
    //if the high gain digit does not exist in the digits array

    seen = kFALSE;

    if(lowGainFlag) {
      for (iOldDigit=iDigit-1; iOldDigit>=0; iOldDigit--) {
        if ((dynamic_cast<AliPHOSDigit*>(digits->At(iOldDigit)))->GetId() == absId) {
          seen = kTRUE;
          break;
        }
      }
      if (!seen) {
        new((*digits)[iDigit]) AliPHOSDigit(-1,absId,(Float_t)decoder->GetEnergy(),time);
        iDigit++;
      }
    }

    // Add high gain digit only if it is not saturated;
    // replace low gain digit by a high gain one
    else {
      if (decoder->GetEnergy() >= 1023) continue;
      for (iOldDigit=iDigit-1; iOldDigit>=0; iOldDigit--) {
        if ((dynamic_cast<AliPHOSDigit*>(digits->At(iOldDigit)))->GetId() == absId) {
          digits->RemoveAt(iOldDigit);
          new((*digits)[iOldDigit]) AliPHOSDigit(-1,absId,(Float_t)decoder->GetEnergy(),time);
          seen = kTRUE;
          break;
        }
      }
      if (!seen) {
        new((*digits)[iDigit]) AliPHOSDigit(-1,absId,(Float_t)decoder->GetEnergy(),time);
        iDigit++;
      }
    }

  }

  digits->Compress();
  digits->Sort();
}
