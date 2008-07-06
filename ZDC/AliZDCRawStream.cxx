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
#include "AliRawDataHeader.h"
#include "AliRawEventHeaderBase.h"
#include "AliLog.h"

ClassImp(AliZDCRawStream)


//_____________________________________________________________________________
AliZDCRawStream::AliZDCRawStream(AliRawReader* rawReader) :
  fRawReader(rawReader),
  fBuffer(0),
  fEvType(0),
  fPosition(0),
  fIsCalib(kFALSE),
  fIsDARCHeader(kFALSE),
  fIsChMapping(kFALSE),
  fIsADCDataWord(kFALSE),
  fIsADCHeader(kFALSE),
  fIsADCEOB(kFALSE),
  fDARCEvBlockLenght(0),  
  fDARCBlockAttributes(0),
  fDeadfaceOffset(0),
  fDeadbeefOffset(0),
  fDataOffset(0),
  fModType(-1),
  fADCModule(-1),
  fADCNChannels(-1),	 
  fADCChannel(-1),	 
  fADCValue(-1),	 
  fADCGain(-1),
  fNConnCh(-1),
  fCabledSignal(-1)
{
  // Create an object to read ZDC raw digits
  fRawReader->Reset();
  fRawReader->Select("ZDC");
  //
  for(Int_t i=0; i<48; i++){
    for(Int_t j=0; j<3; j++) fMapADC[i][j]=-1;
  }

}

//_____________________________________________________________________________
AliZDCRawStream::AliZDCRawStream(const AliZDCRawStream& stream) :
  TObject(stream),
  fRawReader(stream.fRawReader),
  fBuffer(stream.GetRawBuffer()),
  fEvType(stream.fEvType),
  fPosition(stream.fPosition),
  fIsCalib(stream.fIsCalib),
  fIsDARCHeader(stream.fIsDARCHeader), 
  fIsChMapping(stream.fIsChMapping),
  fIsADCDataWord(stream.fIsADCDataWord), 
  fIsADCHeader(stream.fIsADCHeader), 
  fIsADCEOB(stream.fIsADCEOB), 
  fDARCEvBlockLenght(stream.fDARCEvBlockLenght),  
  fDARCBlockAttributes(stream.fDARCBlockAttributes),
  fDeadfaceOffset(stream.GetDeadfaceOffset()),
  fDeadbeefOffset(stream.GetDeadbeefOffset()),
  fDataOffset(stream.GetDataOffset()),
  fModType(stream.GetModType()),
  fADCModule(stream.GetADCModule()),	 
  fADCNChannels(stream.GetADCNChannels()),	 
  fADCChannel(stream.GetADCChannel()),	 
  fADCValue(stream.GetADCValue()),	 
  fADCGain(stream.GetADCGain()),
  fNConnCh(stream.fNConnCh),
  fCabledSignal(stream.GetCabledSignal())
{
  // Copy constructor
  for(Int_t j=0; j<2; j++) fSector[j] = stream.GetSector(j);	 
  for(Int_t i=0; i<48; i++){
    for(Int_t j=0; j<3; j++) fMapADC[i][j] = stream.fMapADC[i][j];
  }
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
void AliZDCRawStream::ReadCDHHeader()
{
  // Reading CDH 
  const AliRawDataHeader* header = fRawReader->GetDataHeader();
  if(!header) {
      AliError("\t No CDH in raw data streaming\n");
      fRawReader->AddMajorErrorLog(kCDHError);
      //
      // For the moment to debug the classe the event is read
      // also if the CDH is not present in the data buffer
      // ******* TO BE CHANGED!!! ***************************
      //return;
  }
  else{
    //printf("\t AliZDCRawStream::ReadCDHHeader -> Data Size = %d\n",fRawReader->GetDataSize());

    fDARCEvBlockLenght = header->fSize;
    //printf("\t AliZDCRawStream::ReadCDHHeader -> fDARCEvBlockLenght = %d\n",fDARCEvBlockLenght);
    
    //UChar_t message = header->GetAttributes();
    //printf("\t AliZDCRawStream::ReadCDHHeader -> Attributes %x\n",message);
    
    /*if(message & 0x10){ // COSMIC RUN
       printf("\t STANDALONE_COSMIC RUN raw data found\n");
    }
    else if(message & 0x20){ // PEDESTAL RUN
       printf("\t STANDALONE_PEDESTAL RUN raw data found\n");
    }
    else if(message & 0x30){ // LASER RUN
       printf("\t STANDALONE_LASER RUN raw data found\n");
    }*/
    
    if(header->GetL1TriggerMessage() & 0x1){ // Calibration bit set in CDH
      fIsCalib = kTRUE;
    }
    //printf("\t AliZDCRawStream::ReadCDHHeader -> L1TriggerMessage %x\n",header->GetL1TriggerMessage());
    printf("\t AliZDCRawStream::ReadCDHHeader -> fIsCalib = %d\n",fIsCalib);
    
    
    UInt_t status = header->GetStatus();
    //printf("\t AliZDCRawStream::ReadCDHHeader -> status = %d\n",status);
    if(status & 0x000f == 0x0001){
      AliWarning("CDH -> DARC trg0 overlap error\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if(status & 0x000f == 0x0002){
      AliWarning("CDH -> DARC trg0 missing error\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if(status & 0x000f == 0x0004){
      AliWarning("CDH -> DARC data parity error\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if(status & 0x000f == 0x0008){
      AliWarning("CDH -> DARC ctrl parity error\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    //
    if(status & 0x00f0 == 0x0010){
      AliWarning("CDH -> DARC trg unavailable\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if(status & 0x00f0 == 0x0020){
      AliWarning("CDH -> DARC FEE error\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    //
    if(status & 0x0f00 == 0x0200){
      AliWarning("CDH -> DARC L1 time violation\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if(status & 0x0f00 == 0x0400){
      AliWarning("CDH -> DARC L2 time-out\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if(status & 0x0f00 == 0x0800){
      AliWarning("CDH -> DARC prepulse time violation\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    //
    if(status & 0xf000 == 0x1000){
      AliWarning("CDH -> DARC other error\n");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
  }
  //
  fIsDARCHeader = kTRUE;
}


//_____________________________________________________________________________
Bool_t AliZDCRawStream::Next()
{
  // Read the next raw digit
  // Returns kFALSE if there is no digit left

  if(!fRawReader->ReadNextInt((UInt_t&) fBuffer)) return kFALSE;
  fIsChMapping = kFALSE; fIsADCHeader = kFALSE; 
  fIsADCDataWord = kFALSE; fIsADCEOB = kFALSE;
  
  fEvType = fRawReader->GetType();
  //printf("\n\t AliZDCRawStream::Next() -> ev. type %d\n",fEvType);
  
  if(fPosition==0){
    //if(fEvType==7 || fEvType ==8){ //Physics or calibration event
      //ReadEventHeader();
      ReadCDHHeader();
    //}
    fNConnCh=0;
  }
  //printf("\n  AliZDCRawStream::Next() - fBuffer[%d] = %x\n",fPosition, fBuffer);
  
  // -------------------------------------------
  // --- DARC header
  // -------------------------------------------
  if(fIsDARCHeader){
    printf("\t ---- DARC header ----\n");
    if(fIsCalib){
      fDeadfaceOffset = 9;
      fDeadbeefOffset = 25;
    }
    else{
      fDeadfaceOffset = 1;
      fDeadbeefOffset = 7;
    }
    fDataOffset = 1+fDeadbeefOffset;
    fIsDARCHeader = kFALSE;
  }
  
    
  // -------------------------------------------
  // --- Start of data event
  // --- decoding mapping of connected ADC ch.
  // -------------------------------------------
  if(fEvType==10){
    if(fPosition>fDataOffset){
      if((fBuffer&0xff000000) == 0xff000000){
        if(fPosition==(fDataOffset+1)){ 
	   printf("\n\n\t Reading ZDC mapping from StartOfData event\n");
	   fNConnCh=0;	
        }
	else{
	  printf("\n\t End of StartOfData event\n\n");
	  return kTRUE;
	}
      }
      else if((fBuffer&0x80000000)>>31 == 1){
        // Mapping identification
	fADCModule = ((fBuffer & 0x7f000000)>>24);
	fModType = ((fBuffer & 0xfff000)>>8);
	fADCNChannels = (fBuffer & 0xff);
	//
	//printf("\tGEO %d, mod. type %d, #ch. %d\n",fADCModule,fModType,fADCNChannels);
      }
      else if(fModType==0 && (fBuffer&0x80000000)>>31 == 0){
        // Channel signal
	if((fBuffer&0x40000000)>>30==0){ // high range chain ADC
	  fIsChMapping = kTRUE;
	  fADCChannel = ((fBuffer & 0x3fff0000)>>16);
	  fCabledSignal = (fBuffer&0xffff);
	  fMapADC[fNConnCh][0] = fADCModule;
	  fMapADC[fNConnCh][1] = fADCChannel;
	  fMapADC[fNConnCh][2] = fCabledSignal;
	  //printf("  RawData%d mod. %d ch. %d signal %d\n",fNConnCh,fADCModule,
	   // fADCChannel, fBuffer&0xffff);
	  //
	  fNConnCh++;
	  if(fNConnCh>48){
	    // Protection manually set since it returns:
	    // RawData48 mod. 3 ch. 2048 signal 515
	    // to be checked with Pietro!!!!!!!!!!!!!!!!!!!!!!!
	    AliWarning("\t AliZDCRawStream -> ERROR: no. of cabled channels > 48!!!\n");
	    return kFALSE;
	  }
	}
      }
    }
    fPosition++;
    return kTRUE;
  }
  
  // -------------------------------------------
  // --- DARC data
  // -------------------------------------------
  if(fPosition<fDeadfaceOffset){
    fPosition++;
    return kTRUE;
  }
  else if(fPosition==fDeadfaceOffset){
    if(fBuffer != 0xdeadface){
      AliWarning("AliZDCRawStream -> NO deadface after DARC data\n");
      fRawReader->AddMajorErrorLog(kDARCError);  
    }
    else{
      fPosition++;
      return kTRUE;
    }
  }
  
  // -------------------------------------------
  // --- DARC global data
  // -------------------------------------------
  else if(fPosition>fDeadfaceOffset && fPosition<fDeadbeefOffset){
    fPosition++;
    return kTRUE;
  }
  else if(fPosition==fDeadbeefOffset){
    if(fBuffer != 0xdeadbeef){
      AliWarning("AliZDCRawStream -> NO deadbeef after DARC global data\n");
      fRawReader->AddMajorErrorLog(kDARCError);  
    }
    else{
      fPosition++;
      return kTRUE;
    }
  }

  // -------------------------------------------
  // --- ZDC data
  // --- ADC buffer + scaler
  // -------------------------------------------
  else if(fPosition>=fDataOffset){
  
    // Not valid datum before the event 
    // there MUST be a NOT valid datum before the event!!!
    if(fPosition==fDataOffset){
      printf("\t **** ZDC data begin ****\n");
      if((fBuffer & 0x07000000) == 0x06000000){
        //printf("    AliZDCRawStream -> Not valid datum in ADC %d,"
        //       "position %d in word data buffer\n",fADCModule,fPosition);
      }
      else fRawReader->AddMajorErrorLog(kZDCDataError);
    }
    
    // If the not valid datum isn't followed by the 1st ADC header
    // the event is corrupted (i.e., 2 gates arrived before trigger)
    else if(fPosition==fDataOffset+1){
      if((fBuffer & 0x07000000) != 0x02000000){
        AliWarning("ZDC ADC -> The not valid datum is NOT followed by an ADC header!\n");
        fRawReader->AddMajorErrorLog(kZDCDataError);
      }
    }
     
    // Get geo address of current word to determine
    // if it is a scaler word (geo address == kScalerAddress)
    // if it is an ADC word (geo address != 8)
    Int_t kScalerAddress=8;
    fADCModule = ((fBuffer & 0xf8000000)>>27);
    if(fADCModule == kScalerAddress){
      DecodeScaler();
    }
    else{//ADC module
      // *** End of event
      if(fBuffer == 0xcafefade){
        printf("  AliZDCRawStream ->  End of ZDC event!\n");
      }
      // *** ADC header
      else if((fBuffer & 0x07000000) == 0x02000000){
        fIsADCHeader = kTRUE;
    	fADCNChannels = ((fBuffer & 0x00003f00)>>8);
    	printf("  AliZDCRawStream -> HEADER: ADC mod.%d has %d ch. \n",fADCModule,fADCNChannels);
      }
      // *** ADC data word
      else if((fBuffer & 0x07000000) == 0x00000000){
        fIsADCDataWord = kTRUE;
        fADCChannel = ((fBuffer & 0x1e0000) >> 17);
        fADCGain = ((fBuffer & 0x10000) >> 16);       
        fADCValue = (fBuffer & 0xfff);  
    	
	printf("  AliZDCRawStream -> DATA: ADC mod. %d ch. %d gain %d value %d\n",
	  fADCModule,fADCChannel,fADCGain,fADCValue);

	// Valid ADC data (not underflow nor overflow)
        if(!(fBuffer & 0x1000) && !(fBuffer & 0x2000)){ 
	  Int_t indSig = -1;
	  for(Int_t k=0; k<48; k++){
	     if(fMapADC[k][0]==fADCModule && fMapADC[k][1]==fADCChannel){
	       //
	       for(Int_t ci=0; ci<48; ci++)
        	printf("  %d mod. %d ch. %d signal %d\n",ci,fMapADC[ci][0],
          	fMapADC[ci][1], fMapADC[ci][2]);

	       indSig=k;
	       break;
	     } 
	  }
	  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  //    To understand the assignment of fSector[i]
	  //   have a look at the enum in AliZDCRawStream.h
	  //   fSector[0] = 1(ZNC+PMRefC), 2(ZPC), 3(ZEM), 
	  //   		    4(ZNA+PMRefA), 5(ZPA)
	  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  if((fMapADC[indSig][2]>=2 && fMapADC[indSig][2]<=6) || (fMapADC[indSig][2]>=26 && fMapADC[indSig][2]<=30) 
	    || fMapADC[indSig][2]==24 || fMapADC[indSig][2]==48){
	    fSector[0] = 4; 
	    //
	    if(fMapADC[indSig][2]==2 || fMapADC[indSig][2]==26) fSector[1]=0;
	    else if(fMapADC[indSig][2]==3 || fMapADC[indSig][2]==27) fSector[1]=1;
	    else if(fMapADC[indSig][2]==4 || fMapADC[indSig][2]==28) fSector[1]=2;
	    else if(fMapADC[indSig][2]==5 || fMapADC[indSig][2]==29) fSector[1]=3;
	    else if(fMapADC[indSig][2]==6 || fMapADC[indSig][2]==30) fSector[1]=4;
	    else if(fMapADC[indSig][2]==24 || fMapADC[indSig][2]==48) fSector[1]=5;
	  }
	  else if((fMapADC[indSig][2]>=7 && fMapADC[indSig][2]<=11) || (fMapADC[indSig][2]>=31 && fMapADC[indSig][2]<=35)){
	    fSector[0] = 5; 
	    //
	    if(fMapADC[indSig][2]==7 || fMapADC[indSig][2]==31) fSector[1]=0;
	    else if(fMapADC[indSig][2]==8 || fMapADC[indSig][2]==32) fSector[1]=1;
	    else if(fMapADC[indSig][2]==9 || fMapADC[indSig][2]==33) fSector[1]=2;
	    else if(fMapADC[indSig][2]==10 || fMapADC[indSig][2]==34) fSector[1]=3;
	    else if(fMapADC[indSig][2]==11 || fMapADC[indSig][2]==35) fSector[1]=4;
	  }
	  else if((fMapADC[indSig][2]>=12 && fMapADC[indSig][2]<=16) || (fMapADC[indSig][2]>=36 && fMapADC[indSig][2]<=40) 
	    || fMapADC[indSig][2]==25 || fMapADC[indSig][2]==49){
	    fSector[0] = 1; 
	    //
	    if(fMapADC[indSig][2]==12 || fMapADC[indSig][2]==36) fSector[1]=0;
	    else if(fMapADC[indSig][2]==13 || fMapADC[indSig][2]==37) fSector[1]=1;
	    else if(fMapADC[indSig][2]==14 || fMapADC[indSig][2]==38) fSector[1]=2;
	    else if(fMapADC[indSig][2]==15 || fMapADC[indSig][2]==39) fSector[1]=3;
	    else if(fMapADC[indSig][2]==16 || fMapADC[indSig][2]==40) fSector[1]=4;
	    else if(fMapADC[indSig][2]==25 || fMapADC[indSig][2]==49) fSector[1]=5;
	  }
	  else if((fMapADC[indSig][2]>=17 && fMapADC[indSig][2]<=21) || (fMapADC[indSig][2]>=41 && fMapADC[indSig][2]<=45)){
	    fSector[0] = 2; 
	    //
	    if(fMapADC[indSig][2]==17 || fMapADC[indSig][2]==41) fSector[1]=0;
	    else if(fMapADC[indSig][2]==18 || fMapADC[indSig][2]==42) fSector[1]=1;
	    else if(fMapADC[indSig][2]==19 || fMapADC[indSig][2]==43) fSector[1]=2;
	    else if(fMapADC[indSig][2]==20 || fMapADC[indSig][2]==44) fSector[1]=3;
	    else if(fMapADC[indSig][2]==21 || fMapADC[indSig][2]==45) fSector[1]=4;
	  }
	  else if(fMapADC[indSig][2]==22 || fMapADC[indSig][2]==23 || fMapADC[indSig][2]==46 || fMapADC[indSig][2]==47){
	    fSector[0] = 3; 
	    //
	    if(fMapADC[indSig][2]==22 || fMapADC[indSig][2]==46) fSector[1]=1;
	    else if(fMapADC[indSig][2]==23 || fMapADC[indSig][2]==49) fSector[1]=2;
	  }
	  //
	  printf("\t Signal %d -> fSector[0] %d, fSector[1] %d\n", fMapADC[indSig][2], fSector[0], fSector[1]);
	  
	  if(fADCModule<0 || fADCModule>3){
            AliWarning(Form("	 AliZDCRawStream -> No valid ADC module: %d\n",fADCModule));
            fRawReader->AddMajorErrorLog(kInvalidADCModule);
          }

        }//No underflow nor overflow	
      }//ADC data word
      // *** ADC EOB
      else if((fBuffer & 0x07000000) == 0x04000000){
        fIsADCEOB = kTRUE;
    	//printf("  AliZDCRawStream -> EOB --------------------------\n");
      }
     }//ADC module
        
    
  }
  fPosition++;

  return kTRUE;
}
