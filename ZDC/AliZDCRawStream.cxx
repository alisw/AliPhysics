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

#include <TSystem.h>
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
  fSODReading(kFALSE),
  fIsMapRead(kFALSE),
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
  fIsUnderflow(kFALSE),
  fIsOverflow(kFALSE),
  fScGeo(0),	  
  fScNWords(0),	  
  fScTriggerSource(0),	  
  fScTriggerNumber(0),
  fIsScEventGood(kTRUE),
  fIsScHeaderRead(kFALSE),
  fScStartCounter(0),
  fScEvCounter(0),
  fNChannelsOn(0),
  fNConnCh(-1),
  fCabledSignal(-1)
{
  // Create an object to read ZDC raw digits
  fRawReader->Reset();
  fRawReader->Select("ZDC");
  //
  const int kNch = 48;
  for(Int_t i=0; i<kNch; i++){
    for(Int_t j=0; j<5; j++) fMapADC[i][j]=-1;
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
  fSODReading(stream.fSODReading),
  fIsMapRead(stream.fIsMapRead),
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
  fIsUnderflow(stream.fIsUnderflow),
  fIsOverflow(stream.fIsOverflow),
  fScGeo(stream.GetScGeo()),	  
  fScNWords(stream.GetScNWords()),	  
  fScTriggerSource(stream.GetScTriggerSource()),	  
  fScTriggerNumber(stream.fScTriggerNumber),
  fIsScEventGood(stream.fIsScEventGood),
  fIsScHeaderRead(stream.fIsScHeaderRead),
  fScStartCounter(stream.fScStartCounter),
  fScEvCounter(stream.fScEvCounter),
  fNChannelsOn(stream.fNChannelsOn),
  fNConnCh(stream.fNConnCh),
  fCabledSignal(stream.GetCabledSignal())
{
  // Copy constructor
  const int kNch = 48;
  for(Int_t j=0; j<2; j++) fSector[j] = stream.GetSector(j);	 
  for(Int_t i=0; i<kNch; i++){
    for(Int_t j=0; j<5; j++) fMapADC[i][j] = stream.fMapADC[i][j];
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
void AliZDCRawStream::ReadChMap()
{
  // Reading channel map
  const int kNch = 48;
  //printf("\t Reading ZDC ADC mapping from OCDB\n");
  AliZDCChMap * chMap = GetChMap();
  //chMap->Print("");
  for(Int_t i=0; i<kNch; i++){
    fMapADC[i][0] = chMap->GetADCModule(i);
    fMapADC[i][1] = chMap->GetADCChannel(i);
    fMapADC[i][2] = -1;
    fMapADC[i][3] = chMap->GetDetector(i);
    fMapADC[i][4] = chMap->GetSector(i);
  }
  fIsMapRead = kTRUE;
}

//_____________________________________________________________________________
void AliZDCRawStream::ReadCDHHeader()
{
  // Reading CDH 
  const AliRawDataHeader* header = fRawReader->GetDataHeader();
  if(!header) {
      AliError(" No CDH in raw data streaming");
      fRawReader->AddMajorErrorLog(kCDHError);
      return;
  }
  else{
    //printf("\t AliZDCRawStream::ReadCDHHeader -> Data Size = %d\n",fRawReader->GetDataSize());

    fDARCEvBlockLenght = header->fSize;
    //printf("\t AliZDCRawStream::ReadCDHHeader -> fDARCEvBlockLenght = %d\n",fDARCEvBlockLenght);
    
    UChar_t message = header->GetAttributes();
    //printf("\t AliZDCRawStream::ReadCDHHeader -> Attributes %x\n",message);
    
    if(message & 0x0){ // PHYSICS RUN
       //printf("\t PHYSICS RUN raw data found\n");
    }
    else if(message & 0x10){ // COSMIC RUN
       //printf("\t STANDALONE_COSMIC RUN raw data found\n");
    }
    else if(message & 0x20){ // PEDESTAL RUN
       //printf("\t STANDALONE_PEDESTAL RUN raw data found\n");
    }
    else if(message & 0x30){ // LASER RUN
       //printf("\t STANDALONE_LASER RUN raw data found\n");
    }
    else if(message & 0x40){ // CALIBRATION_CENTRAL RUN
       //printf("\t CALIBRATION_CENTRAL RUN raw data found\n");
    }
    else if(message & 0x50){ // CALIBRATION_SEMICENTRAL
       //printf("\t CALIBRATION_SEMICENTRAL RUN raw data found\n");
    }
    else if(message & 0x60){ // CALIBRATION_MB
       //printf("\t CALIBRATION_MB RUN raw data found\n");
    }
    else if(message & 0x70){ // CALIBRATION_EMD
       //printf("\t CALIBRATION_EMD RUN raw data found\n");
    }
    
    if(header->GetL1TriggerMessage() & 0x1){ // Calibration bit set in CDH
      fIsCalib = kTRUE;
    }
    //printf("\t AliZDCRawStream::ReadCDHHeader -> L1TriggerMessage %x\n",header->GetL1TriggerMessage());
    //printf("\t AliZDCRawStream::ReadCDHHeader -> Calibration bit = %d\n",fIsCalib);    
    
    UInt_t status = header->GetStatus();
    //printf("\t AliZDCRawStream::ReadCDHHeader -> status = %d\n",status);
    if((status & 0x000f) == 0x0001){
      AliDebug(2,"CDH -> DARC trg0 overlap error");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if((status & 0x000f) == 0x0002){
      AliDebug(2,"CDH -> DARC trg0 missing error");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if((status & 0x000f) == 0x0004){
      AliDebug(2,"CDH -> DARC data parity error");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if((status & 0x000f) == 0x0008){
      AliDebug(2,"CDH -> DARC ctrl parity error");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    //
    if((status & 0x00f0) == 0x0010){
      AliDebug(2,"CDH -> DARC trg unavailable");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if((status & 0x00f0) == 0x0020){
      AliWarning("CDH -> DARC FEE error");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    //
    if((status & 0x0f00) == 0x0200){
      AliDebug(2,"CDH -> DARC L1 time violation");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if((status & 0x0f00) == 0x0400){
      AliWarning("CDH -> DARC L2 time-out");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    if((status & 0x0f00) == 0x0800){
      AliWarning("CDH -> DARC prepulse time violation");
      fRawReader->AddMajorErrorLog(kDARCError);
    }
    //
    if((status & 0xf000) == 0x1000){
      AliDebug(2,"CDH -> DARC other error");
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
  const int kNch = 48;
  Bool_t readScaler = kFALSE;
  fIsChMapping = kFALSE; fIsADCHeader = kFALSE; 
  fIsADCDataWord = kFALSE; fIsADCEOB = kFALSE;
  fIsUnderflow = kFALSE;
  fIsOverflow = kFALSE; 
  fSector[0] = fSector[1] = -1;
  
  fEvType = fRawReader->GetType();
  // CH. debug
  //printf("\n\t AliZDCRawStream::Next() -> ev. type %d\n",fEvType);
  //printf("\n  AliZDCRawStream::Next() - fBuffer[%d] = %x\n",fPosition, fBuffer);
  
  if(fPosition==0){
    //if(fEvType==7 || fEvType ==8){ //Physics or calibration event
      //ReadEventHeader();
      ReadCDHHeader();
    //}
    fNConnCh=0;
  }
  
  // *** End of ZDC event
  if(fBuffer == 0xcafefade){
    //printf("\n-> AliZDCRawStream::Next() ***** End of ZDC event *****\n\n");
    return kFALSE;
  }
  
  // -------------------------------------------
  // --- DARC header
  // -------------------------------------------
  // If the CDH has been read then 
  // the DARC header must follow
  if(fIsDARCHeader){
    //printf("\t ---- DARC header ----\n");
    if(fIsCalib){
      fDeadfaceOffset = 9;
      fDeadbeefOffset = 25;
      readScaler = kTRUE;
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
  // In the SOD event ADC ch. mapping is written
  if(fEvType==10 && fSODReading){
    //printf("\n-> AliZDCRawStream::Next() - fBuffer[%d] = %x\n",fPosition, fBuffer);
    
    if(fPosition>fDataOffset){
      if((fBuffer&0xff000000) == 0xff000000){
        if(fPosition==(fDataOffset+1)){ 
	   printf("\n\n\t AliZDCRawStream -> Reading mapping from StartOfData event\n");
	   fNConnCh=0;	
        }
	else{
	  //printf("\n\t AliZDCRawStream -> End of ZDC StartOfData event\n\n");
          //printf("AliZDCRawStream: fSODReading after SOD reading set to %d\n", fSODReading);
	  return kFALSE;
	}
      }
      else if((fBuffer&0x80000000)>>31 == 1){
        // Mapping identification
	fADCModule = ((fBuffer & 0x7f000000)>>24);
	fModType = ((fBuffer & 0x7ff00)>>8);
	fADCNChannels = (fBuffer & 0xff);
	//
	printf("  ******** GEO %d, mod. type %d, #ch. %d\n",fADCModule,fModType,fADCNChannels);
      }
      else if(fModType==1 && (fBuffer&0x80000000)>>31 == 0){
        // Channel signal
	if((fBuffer&0x40000000)>>30==0){ // high range chain ADC
	  fIsChMapping = kTRUE;
	  fADCChannel = ((fBuffer & 0x3fff0000)>>16);
	  fCabledSignal = (fBuffer&0xffff);
      	  //printf("\tADC ch. %d, signal %d\n",fADCChannel,fCabledSignal);
	  //
	  fMapADC[fNConnCh][0] = fADCModule;
	  fMapADC[fNConnCh][1] = fADCChannel;
	  fMapADC[fNConnCh][2] = fCabledSignal;
	  //
	  // Determining detector and sector
	  // -----------------------------------------
	  //  For the decoding of the following lines 
	  //  look the enum in AliZDCRawStream.h file
	  // -----------------------------------------
	  if((fCabledSignal>=2 && fCabledSignal<=6) || (fCabledSignal>=26 && fCabledSignal<=30)
	     || fCabledSignal==24 || fCabledSignal==48){
	    fMapADC[fNConnCh][3] = 4; //ZNA
	    //
	    if(fCabledSignal==2 || fCabledSignal==26)       fMapADC[fNConnCh][4]=0;
	    else if(fCabledSignal==3 || fCabledSignal==27)  fMapADC[fNConnCh][4]=1;
	    else if(fCabledSignal==4 || fCabledSignal==28)  fMapADC[fNConnCh][4]=2;
	    else if(fCabledSignal==5 || fCabledSignal==29)  fMapADC[fNConnCh][4]=3;
	    else if(fCabledSignal==6 || fCabledSignal==30)  fMapADC[fNConnCh][4]=4;
	    else if(fCabledSignal==24 || fCabledSignal==48) fMapADC[fNConnCh][4]=5;
	  }
	  else if((fCabledSignal>=7 && fCabledSignal<=11) || (fCabledSignal>=31 && fCabledSignal<=35)){
	    fMapADC[fNConnCh][3] = 5; //ZPA
	    //
	    if(fCabledSignal==7 || fCabledSignal==31) 	    fMapADC[fNConnCh][4]=0;
	    else if(fCabledSignal==8 || fCabledSignal==32)  fMapADC[fNConnCh][4]=1;
	    else if(fCabledSignal==9 || fCabledSignal==33)  fMapADC[fNConnCh][4]=2;
	    else if(fCabledSignal==10 || fCabledSignal==34) fMapADC[fNConnCh][4]=3;
	    else if(fCabledSignal==11 || fCabledSignal==35) fMapADC[fNConnCh][4]=4;
	  }
	  else if((fCabledSignal>=12 && fCabledSignal<=16) || (fCabledSignal>=36 && fCabledSignal<=40)
	     || fCabledSignal==25 || fCabledSignal==49){
	    fMapADC[fNConnCh][3] = 1; //ZNC
	    //
	    if(fCabledSignal==12 || fCabledSignal==36)      fMapADC[fNConnCh][4]=0;
	    else if(fCabledSignal==13 || fCabledSignal==37) fMapADC[fNConnCh][4]=1;
	    else if(fCabledSignal==14 || fCabledSignal==38) fMapADC[fNConnCh][4]=2;
	    else if(fCabledSignal==15 || fCabledSignal==39) fMapADC[fNConnCh][4]=3;
	    else if(fCabledSignal==16 || fCabledSignal==40) fMapADC[fNConnCh][4]=4;
	    else if(fCabledSignal==25 || fCabledSignal==49) fMapADC[fNConnCh][4]=5;
	  }
	  else if((fCabledSignal>=17 && fCabledSignal<=21) || (fCabledSignal>=41 && fCabledSignal<=45)){
	    fMapADC[fNConnCh][3] = 2; //ZPC
	    //
	    if(fCabledSignal==17 || fCabledSignal==41) 	    fMapADC[fNConnCh][4]=0;
	    else if(fCabledSignal==18 || fCabledSignal==42) fMapADC[fNConnCh][4]=1;
	    else if(fCabledSignal==19 || fCabledSignal==43) fMapADC[fNConnCh][4]=2;
	    else if(fCabledSignal==20 || fCabledSignal==44) fMapADC[fNConnCh][4]=3;
	    else if(fCabledSignal==21 || fCabledSignal==45) fMapADC[fNConnCh][4]=4;
	  }
	  else if(fCabledSignal==22 || fCabledSignal==23 || fCabledSignal==46 || fCabledSignal==47){
	    fMapADC[fNConnCh][3] = 3;
	    //
	    if(fCabledSignal==22 || fCabledSignal==46)      fMapADC[fNConnCh][4]=1;
	    else if(fCabledSignal==23 || fCabledSignal==47) fMapADC[fNConnCh][4]=2;
	  }
	  //
	  fNConnCh++;
	  if(fNConnCh>=kNch){
	    // Protection manually set since it returned:
	    // RawData48 mod. 3 ch. 2048 signal 515
	    // to be checked with Pietro!!!!!!!!!!!!!!!!!!!!!!!
	    AliDebug(2," No. of cabled channels > kNch !!!");
            fPosition++;
	    return kTRUE;
	  }
	}
      }// ModType=1 (ADC mapping)
      //else if(fModType!=1){
        
      //}
    }
    fPosition++;
    return kTRUE;
  } // SOD event
  
  // -------------------------------------------
  // --- DARC data
  // -------------------------------------------
  if(fPosition<fDeadfaceOffset){
    fPosition++;
    return kTRUE;
  }
  else if(fPosition==fDeadfaceOffset){
    if(fBuffer != 0xdeadface){
      AliWarning(" NO deadface after DARC data");
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
      AliWarning(" NO deadbeef after DARC global data");
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
    
    //printf("AliZDCRawStream: fSODReading = %d\n", fSODReading);
    if(!fSODReading && !fIsMapRead) ReadChMap();
    
    // Not valid datum before the event 
    // there MUST be a NOT valid datum before the event!!!
    if(fPosition==fDataOffset){
      //printf("\t **** ZDC data begin ****\n");
      if((fBuffer & 0x07000000) != 0x06000000){
        fRawReader->AddMajorErrorLog(kZDCDataError);
      }
      /*else{
        //printf("    AliZDCRawStream -> Not valid datum in ADC %d,"
        //       "position %d in word data buffer\n",fADCModule,fPosition);
      }*/
    }
    
    // If the not valid datum isn't followed by the 1st ADC header
    // the event is corrupted (i.e., 2 gates arrived before trigger)
    else if(fPosition==fDataOffset+1){
      if((fBuffer & 0x07000000) != 0x02000000){
        AliWarning("ZDC ADC -> The not valid datum is NOT followed by an ADC header!");
        fRawReader->AddMajorErrorLog(kZDCDataError);
      }
    }
     
    // Get geo address of current word to determine if:
    // - it is an ADC word (geo address <= 3)
    // - it is a scaler word (geo address == kScalerAddress)
    fADCModule = (Int_t) ((fBuffer & 0xf8000000)>>27);
    //printf("  AliZDCRawStream -> fADCModule %d\n",fADCModule);
    
    // ****** ADC MODULES ******
    if(fADCModule>=0 && fADCModule<=3 && fIsScHeaderRead==kFALSE){//ADC modules (0,1,2,3)
      // *** ADC header
      if((fBuffer & 0x07000000) == 0x02000000){
        fIsADCHeader = kTRUE;
    	fADCNChannels = ((fBuffer & 0x00003f00)>>8);
	if(fADCModule==0) fNChannelsOn = fADCNChannels;
	else fNChannelsOn += fADCNChannels;
    	//printf("  AliZDCRawStream -> ADC HEADER: mod.%d has %d ch. \n",fADCModule,fADCNChannels);
      }
      // *** ADC data word
      else if((fBuffer & 0x07000000) == 0x00000000){
        fIsADCDataWord = kTRUE;
        fADCChannel = ((fBuffer & 0x1e0000) >> 17);
        fADCGain = ((fBuffer & 0x10000) >> 16);       
        fADCValue = (fBuffer & 0xfff);  
    	
	//printf("  AliZDCRawStream -> DATA: ADC mod. %d ch. %d gain %d value %d\n",
	//  fADCModule,fADCChannel,fADCGain,fADCValue);
	
	// Checking if the channel map for the ADCs has been provided/read
	if(fMapADC[0][0]==-1){
	  printf("\t ATTENTION!!! No ADC mapping has been found/provided!!!\n");
	  return kFALSE;
	}
	//
	/*for(Int_t ci=0; ci<kNch; ci++){
	  printf("  %d mod %d ch %d det %d sec %d\n",ci,fMapADC[ci][0],
          fMapADC[ci][1], fMapADC[ci][3], fMapADC[ci][4]);
	}*/
		
	// Scan of the map to assign the correct volumes
	Int_t foundMapEntry = kFALSE;
	for(Int_t k=0; k<kNch; k++){
	   if(fADCModule==fMapADC[k][0] && fADCChannel==fMapADC[k][1]){
	     fSector[0] = fMapADC[k][3];
	     fSector[1] = fMapADC[k][4];
	     foundMapEntry = kTRUE;
	     break;
	   } 
	}
	if(foundMapEntry==kFALSE){
	  AliWarning(Form(" No valid entry found in ADC mapping for raw data %d ADCmod. %d ch. %d gain %d\n",
	      fPosition,fADCModule,fADCChannel,fADCGain));
	}
	//
	//printf("AliZDCRawStream -> ADCmod. %d ch. %d gain %d -> det %d sec %d\n",
	//  fADCModule,fADCChannel,fADCGain,fSector[0],fSector[1]);
	
	// Final checks
	if(foundMapEntry==kTRUE){
	  if(fSector[0]<1 || fSector[0]>5){
            AliWarning(Form(" No valid detector assignment: %d",fSector[0]));
            fRawReader->AddMajorErrorLog(kInvalidSector);
	  }
	  //
	  if(fSector[1]<0 || fSector[1]>5){
            AliWarning(Form(" No valid sector assignment: %d",fSector[1]));
            fRawReader->AddMajorErrorLog(kInvalidSector);
	  }
	  //
	  if(fADCModule<0 || fADCModule>3){
            AliError(Form(" No valid ADC module: %d",fADCModule));
            fRawReader->AddMajorErrorLog(kInvalidADCModule);
          }
	}

	// Checking the underflow and overflow bits
        if(fBuffer & 0x1000)       fIsUnderflow = kTRUE;
	else if (fBuffer & 0x2000) fIsOverflow = kTRUE; 
        	
      }//ADC data word
      // *** ADC EOB
      else if((fBuffer & 0x07000000) == 0x04000000){
        fIsADCEOB = kTRUE;
    	//printf("  AliZDCRawStream -> EOB --------------------------\n");
      }
    }//ADC module
    // *** DECODING SCALER
    else if(fADCModule == 8){
      if(fBuffer & 0x04000000 && fIsScHeaderRead==kFALSE){ // *** Scaler header
        fScGeo = (fBuffer & 0xf8000000)>>27;	   
        fScNWords = (fBuffer & 0x00fc0000)>>18;	   
        fScTriggerSource = (fBuffer & 0x00030000)>>16;	   
        fScTriggerNumber = (fBuffer & 0x0000ffff);
        fIsScHeaderRead = kTRUE; 
	fScStartCounter = (Int_t) (fPosition);
        //Ch. debug
        //printf("  AliZDCRawStream -> SCALER HEADER: geo %d Nwords %d TS %d TN %d\n",
        //   fScGeo,fScNWords,fScTriggerSource,fScTriggerNumber);
      } 
      else if(!(fBuffer & 0x04000000)){
        fIsScEventGood = kFALSE;
      }
    }    
    if(fIsScHeaderRead && fPosition>=fScStartCounter+1){ // *** Scaler word
      fScEvCounter = fBuffer;
      Int_t nWords = (Int_t) (fScNWords);
      if(fPosition == fScStartCounter+nWords) fIsScHeaderRead = kFALSE;
      //Ch. debug
      //printf("  AliZDCRawStream -> scaler datum %d", fScEvCounter);
    }
    
  }
  fPosition++;

  return kTRUE;
}

//_____________________________________________________________________________
AliCDBStorage* AliZDCRawStream::SetStorage(const char *uri) 
{
  // Setting the storage
  
  AliCDBStorage *storage = AliCDBManager::Instance()->GetStorage(uri); 

  return storage; 
}


//_____________________________________________________________________________
AliZDCChMap* AliZDCRawStream::GetChMap() const
{

  // Getting calibration object for ZDC

  AliCDBEntry  *entry = AliCDBManager::Instance()->Get("ZDC/Calib/ChMap");
  if(!entry) AliFatal("No calibration data loaded!");  

  AliZDCChMap *calibdata = dynamic_cast<AliZDCChMap*> (entry->GetObject());
  if(!calibdata) AliFatal("Wrong calibration object in calibration  file!");

  return calibdata;
}
