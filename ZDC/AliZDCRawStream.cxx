/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to ZDC digits in raw data.
///
/// It loops over all ZDC digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Getters provide information about the current digit.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCRawStream.h"
#include "AliRawReader.h"

ClassImp(AliZDCRawStream)


//_____________________________________________________________________________
AliZDCRawStream::AliZDCRawStream(AliRawReader* rawReader) :
  fRawReader(rawReader),
  fADCValue(-1)
{
// create an object to read ZDC raw digits

  fSector[0] = 1;
  fSector[1] = -1;
  fRawReader->Select(kDDLOffset / 0x100);
}

//_____________________________________________________________________________
AliZDCRawStream::AliZDCRawStream(const AliZDCRawStream& stream) :
  TObject(stream),
  fADCValue(-1)
{
// copy constructor

  Fatal("AliZDCRawStream", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliZDCRawStream& AliZDCRawStream::operator = (const AliZDCRawStream& 
					      /* stream */)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliZDCRawStream::~AliZDCRawStream()
{
// destructor

}


//_____________________________________________________________________________
Bool_t AliZDCRawStream::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  if (!fRawReader->ReadNextInt((UInt_t&) fRawADC)) return kFALSE;
  fIsADCDataWord = kFALSE;
  
  //ADC Header
  if (fRawADC & 0x2000000) {
    printf("This is the ADC Header\n");
    printf("%d data words will follow \n",2*((fRawADC & 0x3f00) >> 8));
  } 
  //ADC EOB
  else if (fRawADC & 0x4000000) {
    printf("This is the ADC End Of Block\n");
    printf("This was event number %d\n",(fRawADC & 0xffffff));
  } else 
  //ADC Data Words
  {
    printf("This is an ADC Data Word\n");
    printf("Channel %d range %d\n", 
           (fRawADC & 0x1e0000) >> 17, (fRawADC & 0x10000) >> 16);
    if(fRawADC & 0x1000) printf("Data = overflow\n");
    fADCGain = (fRawADC & 0x10000) >> 16;
    fADCValue = (fRawADC & 0xfff);   
    fIsADCDataWord = kTRUE;

    Int_t ADCChannel = (fRawADC & 0x1e0000) >> 17;
    if (ADCChannel >= 0 && ADCChannel <= 4) { 
      fSector[0] = 1;
      fSector[1] = ADCChannel;
    } else if (ADCChannel >= 8 && ADCChannel <= 12) {
      fSector[0] = 2;
      fSector[1] = ADCChannel-8;
    } else if (ADCChannel == 5 || ADCChannel == 13){
      fSector[0] = 3;
      fSector[1] = (ADCChannel-5)/8;
    }
  }
  return kTRUE;
}
