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
  // Create an object to read ZDC raw digits

  fSector[0] = 0;
  fSector[1] = -1;
  fADCModule = 0;
  fRawReader->Select("ZDC");
}

//_____________________________________________________________________________
AliZDCRawStream::AliZDCRawStream(const AliZDCRawStream& stream) :
  TObject(stream),
  fADCValue(-1)
{
  // Copy constructor
  Fatal("AliZDCRawStream", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliZDCRawStream& AliZDCRawStream::operator = (const AliZDCRawStream& 
					      /* stream */)
{
  // Assignment operator
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliZDCRawStream::~AliZDCRawStream()
{
// Destructor

}


//_____________________________________________________________________________
Bool_t AliZDCRawStream::Next()
{
  // Read the next raw digit
  // Returns kFALSE if there is no digit left

  if(!fRawReader->ReadNextInt((UInt_t&) fRawADC)) return kFALSE;
  fIsADCDataWord = kFALSE;
  
  //ADC Header
  if(fRawADC & 0x2000000){
      if(((fRawADC & 0x3f00) >> 8) == 24)      fADCModule=1; //fRawADC=2001800 -> 24 words -> 1st ADC module
      else if(((fRawADC & 0x3f00) >> 8) == 20) fADCModule=2; //fRawADC=2001400 -> 20 words -> 2nd ADC module
      //
      //printf(" **** This is the ADC Header - %d data words will follow \n",((fRawADC & 0x3f00) >> 8));
      //printf("  fRawADC = %x, fADCModule = %d\n",fRawADC, fADCModule);
  }
  else if((fRawADC & 0x4000000) || (fRawADC & 0x3000000)){
    fSector[0] = 0;
    //ADC EOB
    /*if(fRawADC & 0x4000000){
      printf(" **** This is the ADC End Of Block - event number %d\n",(fRawADC & 0xffffff));
    }*/
  } 
  //ADC Data Words
  else{
    //printf("This is an ADC Data Word -> channel %d range %d\n",(fRawADC & 0x1e0000) >> 17, (fRawADC & 0x10000) >> 16);
    if(fRawADC & 0x1000) printf("ZDCRawStream -> ADC overflow\n");
    if(fRawADC & 0x2000) printf("ZDCRawStream -> ADC underflow\n");
    //
    fADCGain = (fRawADC & 0x10000) >> 16;
    fADCValue = (fRawADC & 0xfff);   
    fIsADCDataWord = kTRUE;
    //
    Int_t vADCChannel = (fRawADC & 0x1e0000) >> 17;
    if(fADCModule==1){  //1st ADC module
      if(vADCChannel >= 0 && vADCChannel <= 4){ 
        fSector[0] = 1;
        fSector[1] = vADCChannel;
      } 
      else if(vADCChannel >= 8 && vADCChannel <= 12){
        fSector[0] = 2;
        fSector[1] = vADCChannel-8;
      } 
      else if(vADCChannel == 5 || vADCChannel == 13){
        fSector[0] = 3;
        fSector[1] = ((vADCChannel-5)/8)+1;
      }
    }
    else if(fADCModule==2){  //2nd ADC module
      if(vADCChannel >= 0 && vADCChannel <= 4){ 
        fSector[0] = 4;
        fSector[1] = vADCChannel;
      } 
      else if(vADCChannel >= 8 && vADCChannel <= 12){
        fSector[0] = 5;
        fSector[1] = vADCChannel-8;
      } 
    }
    else printf("\t AliZDCRawStreamer -> ERROR! No valid ADC module!\n");
    
  }
  return kTRUE;
}
