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
//									     //
// This class provides access to ZDC digits in raw data.		     //
//									     //
// It loops over all ZDC digits in the raw data given by the AliRawReader.   //
// The Next method goes to the next digit. If there are no digits left	     //
// it returns kFALSE.							     //
// Getters provide information about the current digit. 		     //
//									     //
///////////////////////////////////////////////////////////////////////////////

#include "AliZDCRawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"

ClassImp(AliZDCRawStream)


//_____________________________________________________________________________
AliZDCRawStream::AliZDCRawStream(AliRawReader* rawReader) :
  fRawReader(rawReader),
  fRawADC(0),	 
  fADCModule(0),	 
  fADCValue(-1),	 
  fADCGain(0),
  fIsADCDataWord(kFALSE)
{
  // Create an object to read ZDC raw digits

  fRawReader->Select("ZDC");
}

//_____________________________________________________________________________
AliZDCRawStream::AliZDCRawStream(const AliZDCRawStream& stream) :
  TObject(stream)
{
  // Copy constructor
  fRawADC = stream.GetADCRaw();	 
  for(Int_t j=0; j<2; j++) fSector[j] = stream.GetSector(j);	 
  fADCModule = stream.GetADCModule();	 
  fADCValue = stream.GetADCValue();	 
  fADCGain = stream.GetADCGain();	 
  fIsADCDataWord = stream.IsADCDataWord();
  
    
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
  
  if(fRawADC & 0x2000000){//ADC Header
    fADCModule++;
    //printf(" \t AliZDCRawStream -> Header ADC %d -> ADC datum contains %d data words \n",fADCModule,((fRawADC & 0x3f00) >> 8));
  }
  else if((fRawADC & 0x4000000)){//ADC EOB
    //printf(" \t AliZDCRawStream -> ADC %d End Of Block - event number %d\n\n",fADCModule, (fRawADC & 0xffffff));
    fSector[0] = fSector[1] = 99;
  } 
  else if((fRawADC & 0x6000000)){//Not valid datum
    printf(" \t AliZDCRawStream -> Not valid datum in ADC module %d\n",fADCModule);
  }
  else{//ADC data word
    fIsADCDataWord = kTRUE;
    //printf(" \t \t ADC Data Word");
    if(!(fRawADC & 0x1000) && !(fRawADC & 0x2000)){ // Valid ADC data
      fADCGain = (fRawADC & 0x10000) >> 16;
      fADCValue = (fRawADC & 0xfff);   
      //
      Int_t vADCChannel = (fRawADC & 0x1e0000) >> 17;
      if(fADCModule==1 || fADCModule==3){  //1st & 3rd ADC modules
        if(vADCChannel >= 0 && vADCChannel <= 4){ 
          fSector[0] = 1; // ZN1
          fSector[1] = vADCChannel;
        } 
        else if(vADCChannel >= 8 && vADCChannel <= 12){
          fSector[0] = 2; // ZP1
          fSector[1] = vADCChannel-8;
        } 
        else if(vADCChannel == 5 || vADCChannel == 13){
          fSector[0] = 3; // ZEM 1,2
          fSector[1] = ((vADCChannel-5)/8)+1;
        }
      }
      else if(fADCModule==2 || fADCModule==4){  //2nd & 4rth ADC modules
        if(vADCChannel >= 0 && vADCChannel <= 4){ 
          fSector[0] = 4; // ZN2
          fSector[1] = vADCChannel;
        } 
        else if(vADCChannel >= 8 && vADCChannel <= 12){
          fSector[0] = 5; // ZP2
          fSector[1] = vADCChannel-8;
        } 
        else if(vADCChannel == 5 || vADCChannel == 13){
          fSector[0] = (vADCChannel-5)*3/8+1; // PM Ref 1,2
          fSector[1] = 5;
        }
      }
      else{
        AliWarning(" \t AliZDCRawStream -> No valid ADC module!");
        printf(" ADCmod = %d\n", fADCModule);
        fRawReader->AddMajorErrorLog(kInvalidADCModule);
      }
      /*printf(" \t \tADC %d ch %d range %d -> det %d quad %d value %x\n",
           fADCModule, vADCChannel, fADCGain, fSector[0], fSector[1], fADCValue);//Chiara debugging
      */
    }
    else if(fRawADC & 0x1000) 
      printf(" \t \tADC %d ch %d range %d  overflow\n",fADCModule, (fRawADC & 0x1e0000) >> 17, fADCGain); // Overflow
    else if(fRawADC & 0x2000) 
      printf(" \t \tADC %d ch %d range %d  underflow\n",fADCModule, (fRawADC & 0x1e0000) >> 17, fADCGain); // Underflow
  }

  return kTRUE;
}
