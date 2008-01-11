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
  fADCChannel(-1),	 
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
  fADCChannel = stream.GetADCChannel();	 
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
  //
  // --- DARC header
  if((fRawADC == 0xe52b6300) || (fRawADC == 0xe52c0300) ||  (fRawADC == 0xffffffff) ||
     (fRawADC == 0xdeadface) || (fRawADC == 0xdeadbeef)){
    //printf("    This is a DARC header!!!\n");
  }
  // --- End of data
  else if(fRawADC == 0xcafefade){
    //printf("    End of data!\n");
  } 
  // --- ADC buffer
  else{
    Bool_t firstADCHeader = kTRUE;
    //
    if((fRawADC & 0x6000000) == 0x6000000){ // Not valid datum BEFORE valid data!
      firstADCHeader = kFALSE;
      //printf("    AliZDCRawStream -> Not valid datum in ADC module!!! %d\n",fADCModule);
    }
    else if((fRawADC & 0x4000000) == 0x4000000){ // ADC EOB
      //printf("    AliZDCRawStream -> ADC %d End Of Block - Event no. %d\n\n",fADCModule, (fRawADC & 0xffffff));
      fSector[0] = fSector[1] = 99;
    } 
    else if((fRawADC & 0x2000000) == 0x2000000){ // ADC Header
      //printf("  fRawADC %x\n",fRawADC);
      // Reading the GEO address to determine ADC module
      fADCModule = (fRawADC & 0xf8000000) >> 27;
      fADCModule ++;
      // If the not valid datum isn't followed by the 1st ADC header
      // the event is corrupted (i.e., 2 gates arrived before trigger)
      if(fADCModule == 1) firstADCHeader = kTRUE;
      //if(fADCModule == 1) printf(" fADCModule %d, firstADCHeader %d\n",fADCModule,firstADCHeader);
      //printf("  AliZDCRawStream -> HEADER: ADC mod. %d contains %d data words \n",fADCModule,((fRawADC & 0x3f00) >> 8));
    }
    else{ // ADC data word
      fIsADCDataWord = kTRUE;
      //printf(" fRawADC = %x firstADCHeader %d\n",fRawADC,firstADCHeader);
      //
      if(!(fRawADC & 0x1000) && !(fRawADC & 0x2000) && firstADCHeader){ // Valid ADC data
        fADCGain = (fRawADC & 0x10000) >> 16;
        fADCValue = (fRawADC & 0xfff);   
        //
        fADCChannel = (fRawADC & 0x1e0000) >> 17;
	//
	// If raw data corresponds to an ADC ch. without physical signal
	// fsector[0] = fSector[1] = -1
	for(Int_t i=0; i<2; i++)fSector[i] = -1;
        //if(fSector[0]!=-1) printf(" \t \t ADC Data Word: ADC ch. %d",fADCChannel);
	//
        if(fADCModule==1 || fADCModule==3){  //1st & 3rd ADC modules
          if(fADCChannel >= 0 && fADCChannel <= 4){ 
            fSector[0] = 1; // ZN1
            fSector[1] = fADCChannel;
          } 
          else if(fADCChannel >= 8 && fADCChannel <= 12){
            fSector[0] = 2; // ZP1
            fSector[1] = fADCChannel-8;
          } 
          else if(fADCChannel == 5 || fADCChannel == 13){
            fSector[0] = 3; // ZEM 1,2
            fSector[1] = ((fADCChannel-5)/8)+1;
          }
        }
        else if(fADCModule==2 || fADCModule==4){  //2nd & 4rth ADC modules
          if(fADCChannel >= 0 && fADCChannel <= 4){ 
            fSector[0] = 4; // ZN2
          fSector[1] = fADCChannel;
          } 
          else if(fADCChannel >= 8 && fADCChannel <= 12){
            fSector[0] = 5; // ZP2
            fSector[1] = fADCChannel-8;
          } 
          else if(fADCChannel == 5 || fADCChannel == 13){
            fSector[0] = (fADCChannel-5)*3/8+1; // PM Ref 1,2
            fSector[1] = 5;
          }
        }
        else{
          AliWarning("    AliZDCRawStream -> No valid ADC module\n");
          fRawReader->AddMajorErrorLog(kInvalidADCModule);
        }
        //printf("  ADC %d ch %d range %d -> det %d quad %d value %x\n",
        //   fADCModule, fADCChannel, fADCGain, fSector[0], fSector[1], fADCValue);//Chiara debugging
      
      }// Valid ADC data
      else if(fRawADC & 0x1000){ 
        //printf("  ADC %d ch %d range %d  overflow\n",fADCModule, (fRawADC & 0x1e0000) >> 17, fADCGain); // Overflow
      }
      else if(fRawADC & 0x2000){ 
        //printf("  ADC %d ch %d range %d  underflow\n",fADCModule, (fRawADC & 0x1e0000) >> 17, fADCGain); // Underflow
      }
    }// ADC word
  }//Not DARC Header

  return kTRUE;
}
