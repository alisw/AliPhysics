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

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to CTP DDL raw data.
///
/// The raw data format is taken form the trigger TDR.
/// The meaning of the trigger class and cluster masks
/// are given in the trigger description file (in /data)
/// and in the AliCentralTrigger class.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliCTPRawStream.h"
#include "AliRawReader.h"
#include "AliLog.h"
#include "AliDAQ.h"
#include "AliTriggerIR.h"

ClassImp(AliCTPRawStream)

//_____________________________________________________________________________
AliCTPRawStream::AliCTPRawStream(AliRawReader* rawReader) :
  fIRArray("AliTriggerIR",3),
  fOrbit(0),
  fBC(0),
  fL0TriggerInputs(0),
  fL1TriggerInputs(0),
  fL2TriggerInputs(0),
  fClassMask(0),
  fClassMaskNext50(0),
  fClusterMask(0),
  fRawReader(rawReader)
{
  // create an object to read CTP raw data
  //
  // select the raw data corresponding to
  // the CTP detector id
  fRawReader->Reset();
  AliDebug(1,Form("Selecting raw data for detector %d",AliDAQ::DetectorID("TRG")));
  fRawReader->Select("TRG");
}

//_____________________________________________________________________________
AliCTPRawStream::AliCTPRawStream(const AliCTPRawStream& stream) :
  TObject(stream),
  fIRArray("AliTriggerIR",3),
  fOrbit(0),
  fBC(0),
  fL0TriggerInputs(0),
  fL1TriggerInputs(0),
  fL2TriggerInputs(0),
  fClassMask(0),
  fClassMaskNext50(0),
  fClusterMask(0),
  fRawReader(NULL)
{
  // Copy constructor
  AliFatal("Copy constructor not implemented");
}

//_____________________________________________________________________________
AliCTPRawStream& AliCTPRawStream::operator = (const AliCTPRawStream& 
					      /* stream */)
{
  AliFatal("Assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliCTPRawStream::~AliCTPRawStream()
{
  // destructor
  fIRArray.Delete();
}

//_____________________________________________________________________________
void AliCTPRawStream::Reset()
{
  // reset raw stream params
  fIRArray.Clear();

  fClassMask = fClusterMask = 0;
  fClassMaskNext50 = 0;

  if (fRawReader) fRawReader->Reset();
}

//_____________________________________________________________________________
Bool_t AliCTPRawStream::Next()
{
  // read the whole CTP raw data stream
  // return kFALSE in case of error

  UChar_t *data = NULL;

  // CTP raw data does not contain CDH
  fRawReader->RequireHeader(kFALSE);

  if (!fRawReader->ReadNextData(data)) {
    fRawReader->RequireHeader(kTRUE);
    return kFALSE;
  }

  if ((fRawReader->GetDataSize()) < 32 ||
      ((fRawReader->GetDataSize() % 4) != 0)) {
    AliError(Form("Wrong CTP raw data size: %d",fRawReader->GetDataSize()));
    fRawReader->RequireHeader(kTRUE);
    return kFALSE;
  }

  fBC = data[0];
  fBC |= (data[1] & 0xF) << 8;

  fOrbit = data[4] << 12;
  fOrbit |= (data[5] & 0xF) << 20;
  fOrbit |= data[8];
  fOrbit |= (data[9] & 0xF) << 8;

  // 4th word bit 12= version, check here.
  UInt_t version=data[13] & 0x10;

  if(version == 0){
   return GetPayloadRun1(data);
  } 
  else{
   return GetPayloadRun2(data);
  } 
} 
Bool_t AliCTPRawStream::GetPayloadRun2(UChar_t *data)
{
  fClusterMask = data[12]&0xef;

  fClassMaskNext50 =  ((ULong64_t)data[16] & 0xf) << 46;   // 100..97

  fClassMaskNext50 |= ((ULong64_t)data[21] & 0xF) << 42;   // 96..93
  fClassMaskNext50 |= (ULong64_t)data[20] << 34;           // 92..85

  fClassMaskNext50 |= ((ULong64_t)data[25] & 0xF) << 30;   //84..81
  fClassMaskNext50 |= (ULong64_t)data[24] << 22;           //80..73

  fClassMaskNext50 |= ((ULong64_t)data[29] & 0xF) << 18;   //72..69
  fClassMaskNext50 |= (ULong64_t)data[28] << 10;           //68..61   

  fClassMaskNext50 |= ((ULong64_t)data[33] & 0xf ) << 6;   //60..57
  fClassMaskNext50 |= ((ULong64_t)data[32] & 0xfc )>> 2;   //56..51
  fClassMask        = ((ULong64_t)data[32] & 0x3 ) << 48;  //50..49

  fClassMask |= ((ULong64_t)data[37] & 0xf) << 44;         //48..45
  fClassMask |= (ULong64_t)data[36]  << 36;                //44..37

  fClassMask |= ((ULong64_t)data[41] & 0xf) << 32;         //36..33
  fClassMask |= (ULong64_t)data[40]  << 24;                //32..25

  fClassMask |= ((ULong64_t)data[45] & 0xf) << 20;         //24..21
  fClassMask |= (ULong64_t)data[44]  << 12;                //20..13

  fClassMask |= ((ULong64_t)data[49] & 0xf) << 8;         //12..9
  fClassMask |= (ULong64_t)data[48]  << 0;                //8..1

  if (fRawReader->GetDataSize() == 52) {
    AliDebug(1,"No trigger input and interaction records found");
    fRawReader->RequireHeader(kTRUE);
    return kTRUE;
  }

  // Read detector trigger inputs
  if (fRawReader->GetDataSize() < 72) {
    AliError(Form("Wrong CTP raw data size: %d",fRawReader->GetDataSize()));
    fRawReader->RequireHeader(kTRUE);
    return kFALSE;
  }

  fL0TriggerInputs = data[52] << 12;
  fL0TriggerInputs |= (data[53] & 0xF) << 20;
  fL0TriggerInputs |= data[56];
  fL0TriggerInputs |= (data[57] & 0xF) << 8;

  fL1TriggerInputs = data[60] << 12;
  fL1TriggerInputs |= (data[61] & 0xF) << 20;
  fL1TriggerInputs |= data[64];
  fL1TriggerInputs |= (data[65] & 0xF) << 8;

  fL2TriggerInputs = data[68] << 12;
  fL2TriggerInputs |= (data[69] & 0xF) << 20;

  if (fRawReader->GetDataSize() == 72) {
    AliDebug(1,"No interaction records found");
    fRawReader->RequireHeader(kTRUE);
    return kTRUE;
  }

  // Read IRs
  Int_t iword = 72;
  //Int_t iii=72;
  //while(iii<fRawReader->GetDataSize()){
  // printf("%i %02x %02x %02x %02x \n",iii,data[iii+3],data[iii+2],data[iii+1],data[iii]);
  // iii+=4;
  //}
  UInt_t ddl2=data[iword+3];
  UInt_t* atr=(UInt_t*)fRawReader->GetAttributes();
  if(atr){
    AliDebug(1,Form("USER ATTRIBUTES 0x%x 0x%x 0x%x\n",atr[0],atr[1],atr[2]));
  }else AliDebug(1,"NO USER ATTRIBUTE found.");
  //printf("ddl2= 0x%x \n",ddl2);
  Int_t ret=0;
  if(ddl2==0){ // DDL1 first
   ret = GetDDL1Data(data,iword);
   if(ret){ // Also DDL2 present
    Int_t ret2=GetDDL2Data(data,iword);
    if(ret2){
     AliWarning(Form("Invalid interaction record: DDL1-DDL2-DDL1 %i",iword));
    }
   }
  }else{
   // DDL2 first
   ret=GetDDL2Data(data,iword);
   if(ret){
    AliWarning(Form("Invalid interaction record: DDL2-DDL1 %i",iword));
   }
  }
  // Restore the raw-reader state!!
  fRawReader->RequireHeader(kTRUE);
 return kTRUE;
}
Int_t AliCTPRawStream::GetDDL1Data(UChar_t* data,Int_t& iword)
{
  AliDebug(1,"Doing DDL1");
  UChar_t level = 0;
  UInt_t *irdata = NULL;
  UInt_t irsize = 0;
  UInt_t orbit = 0;
  UInt_t norbits = 0;
  Bool_t incomplete = kFALSE, transerr = kFALSE;
  while ((iword < fRawReader->GetDataSize()) && (data[iword+3]==0)) {
    if (data[iword+1] & 0x80) {
      UChar_t flag = ((data[iword+1] >> 4) & 0x3);
      if (flag == 0) {
        norbits++;
	if (irdata) {
	  new (fIRArray[fIRArray.GetEntriesFast()])
	    AliTriggerIR(orbit,irsize,irdata,incomplete,transerr);
	  irdata = NULL; irsize = 0;
	}
	level = 1;
	orbit = data[iword] << 12;
	orbit |= (data[iword+1] & 0xF) << 20;
	transerr = ((data[iword+1] >> 6) & 0x1);
	iword += 4;
	continue;
      }
      else if (flag == 3) {
	if (level == 1) {
	  level = 2;
	  orbit |= data[iword];
	  orbit |= ((data[iword+1] & 0xF) << 8);
	  iword += 4;
	  AliDebug(1,Form("Orbit=0x%x\n",orbit));
	  continue;
	}
      }
      UShort_t bc = data[iword];
      bc |= ((data[iword+1] & 0xF) << 8);
      AliDebug(1,Form("BC=0x%x\n",bc));
      if (bc == 0xFFF) {
	incomplete = kTRUE;
      }
      else {
	if (level == 2) {
	  level = 3;
	  irdata = (UInt_t *)&data[iword];
	  irsize = 0;
	}
	if (level == 3) {
	  irsize++;
	}
      }
    }
    else
      AliWarning(Form("Invalid interaction record (%d %d)",iword,fRawReader->GetDataSize()));

    iword += 4;
  }
  if (irdata) {
    new (fIRArray[fIRArray.GetEntriesFast()])
      AliTriggerIR(orbit,irsize,irdata,incomplete,transerr);
    irdata = NULL; irsize = 0;
  }
 if(norbits != 5)AliWarning(Form("Number of orbits != %i",norbits));
 if((iword>=fRawReader->GetDataSize())) return 0;
 return 1;
}
Int_t AliCTPRawStream::GetDDL2Data(UChar_t* data,Int_t& iword)
{
   AliDebug(1,"Doing DDL2");
   ULong64_t irdata[3564];
   UInt_t irsize = 0;
   UInt_t orbit = 0x0;
   UInt_t norbits=0;
   Bool_t incomplete = kFALSE, transerr = kFALSE;
   while ((iword < fRawReader->GetDataSize()) && ((data[iword+3]+data[iword+2]) !=0)) {
      if((data[iword+3]&0xc0) == 0xc0){   // DDL2
        norbits++;
        if(irsize != 0){
	 new (fIRArray[fIRArray.GetEntriesFast()]) AliTriggerIR(orbit,irsize,irdata,incomplete,transerr);
         irsize=0;
	}
        orbit = data[iword]+(data[iword+1]<<8)+(data[iword+2]<<16);
      }  
      else if((data[iword+3]&0x80)==0x80){
        irdata[irsize]=data[iword]+(data[iword+1]<<8)+(data[iword+2]<<16)+((data[iword+3]&0x3f)<<24);
      }
      else if((data[iword+3]&0x40)==0x40){
        ULong64_t word2=data[iword]+(data[iword+1]<<8)+(data[iword+2]<<16)+((data[iword+3]&0x3f)<<24);
	irdata[irsize] += word2<<30;
	irsize++;
      }
      else{
        ULong64_t word=data[iword]+(data[iword+1]<<8)+(data[iword+2]<<16)+((data[iword+3]&0x3f)<<24);
        if(word==0x1fffffff)incomplete=kTRUE;
	else if(word==0x3fffffff)transerr=kTRUE;
	else{
          AliWarning(Form("Invalid interaction record (%d %d)",iword,fRawReader->GetDataSize()));
	}
      }
      iword += 4;
    }
    if (irsize != 0) {
       new (fIRArray[fIRArray.GetEntriesFast()]) AliTriggerIR(orbit,irsize,irdata,incomplete,transerr);
    }
 if(norbits != 5)AliWarning(Form("Number of orbits != %i",norbits));
 if((iword>=fRawReader->GetDataSize())) return 0;
 return 1;
}
Bool_t AliCTPRawStream::GetPayloadRun1(UChar_t *data)
{
  fClusterMask = data[12] >> 2;

  fClassMask =  ((ULong64_t)data[12] & 0x3) << 48;

  fClassMask |= (ULong64_t)data[16] << 36;
  fClassMask |= ((ULong64_t)data[17] & 0xF) << 44;

  fClassMask |= (ULong64_t)data[20] << 24;
  fClassMask |= ((ULong64_t)data[21] & 0xF) << 32;

  fClassMask |= (ULong64_t)data[24] << 12;
  fClassMask |= ((ULong64_t)data[25] & 0xF) << 20;

  fClassMask |= (ULong64_t)data[28];
  fClassMask |= ((ULong64_t)data[29] & 0xF) << 8;

  if (fRawReader->GetDataSize() == 32) {
    AliDebug(1,"No trigger input and interaction records found");
    fRawReader->RequireHeader(kTRUE);
    return kTRUE;
  }

  // Read detector trigger inputs
  if (fRawReader->GetDataSize() < 52) {
    AliError(Form("Wrong CTP raw data size: %d",fRawReader->GetDataSize()));
    fRawReader->RequireHeader(kTRUE);
    return kFALSE;
  }

  fL0TriggerInputs = data[32] << 12;
  fL0TriggerInputs |= (data[33] & 0xF) << 20;
  fL0TriggerInputs |= data[36];
  fL0TriggerInputs |= (data[37] & 0xF) << 8;

  fL1TriggerInputs = data[40] << 12;
  fL1TriggerInputs |= (data[41] & 0xF) << 20;
  fL1TriggerInputs |= data[44];
  fL1TriggerInputs |= (data[45] & 0xF) << 8;

  fL2TriggerInputs = data[48] << 12;
  fL2TriggerInputs |= (data[49] & 0xF) << 20;

  if (fRawReader->GetDataSize() == 52) {
    AliDebug(1,"No interaction records found");
    fRawReader->RequireHeader(kTRUE);
    return kTRUE;
  }

  // Read IRs
  Int_t iword = 52;
  UChar_t level = 0;
  UInt_t *irdata = NULL;
  UInt_t irsize = 0;
  UInt_t orbit = 0;
  Bool_t incomplete = kFALSE, transerr = kFALSE;
  while (iword < fRawReader->GetDataSize()) {
    if (data[iword+1] & 0x80) {
      UChar_t flag = ((data[iword+1] >> 4) & 0x3);
      if (flag == 0) {
	if (irdata) {
	  new (fIRArray[fIRArray.GetEntriesFast()])
	    AliTriggerIR(orbit,irsize,irdata,incomplete,transerr);
	  irdata = NULL; irsize = 0;
	}
	level = 1;
	orbit = data[iword] << 12;
	orbit |= (data[iword+1] & 0xF) << 20;
	transerr = ((data[iword+1] >> 6) & 0x1);
	iword += 4;
	continue;
      }
      else if (flag == 3) {
	if (level == 1) {
	  level = 2;
	  orbit |= data[iword];
	  orbit |= ((data[iword+1] & 0xF) << 8);
	  iword += 4;
	  AliDebug(1,Form("Orbit=0x%x\n",orbit));
	  continue;
	}
      }
      UShort_t bc = data[iword];
      bc |= ((data[iword+1] & 0xF) << 8);
      AliDebug(1,Form("BC=0x%x\n",bc));
      if (bc == 0xFFF) {
	incomplete = kTRUE;
      }
      else {
	if (level == 2) {
	  level = 3;
	  irdata = (UInt_t *)&data[iword];
	  irsize = 0;
	}
	if (level == 3) {
	  irsize++;
	}
      }
    }
    else
      AliWarning(Form("Invalid interaction record (%d %d)",iword,fRawReader->GetDataSize()));

    iword += 4;
  }

  if (irdata) {
    new (fIRArray[fIRArray.GetEntriesFast()])
      AliTriggerIR(orbit,irsize,irdata,incomplete,transerr);
    irdata = NULL; irsize = 0;
  }

  // Restore the raw-reader state!!
  fRawReader->RequireHeader(kTRUE);

  return kTRUE;
}

