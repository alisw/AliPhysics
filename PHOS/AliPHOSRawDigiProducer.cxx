/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

/* $Id$ */

//This class produces PHOS digits of one event
//using AliPHOSRawDecoder. 
//
//   For example:
//   TClonesArray *digits = new TClonesArray("AliPHOSDigit",100);
//   AliRawReader* rawReader = new AliRawReaderDate("2006run2211.raw");
//   AliPHOSRawDecoder dc(rawReader);
//   while (rawReader->NextEvent()) {
//     AliPHOSRawDigiProducer producer;
//     producer.MakeDigits(digits,&dc);
//   }

// Author: Boris Polichtchouk

// --- ROOT system ---
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOSRawDigiProducer.h"
#include "AliPHOSRawDecoder.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"

ClassImp(AliPHOSRawDigiProducer)

//--------------------------------------------------------------------------------------
void AliPHOSRawDigiProducer::MakeDigits(TClonesArray *digits, AliPHOSRawDecoder* decoder) 
{
  //Makes the job.
  //TClonesArray *digits and raw data decoder should be provided by calling function.

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
