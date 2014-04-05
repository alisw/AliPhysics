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

  // 4th word bit 11= version, check here.
  UInt_t version=data[13] & 0x4;

  if(version == 0){
   return GetPayloadRun1(data);
  } 
  else{
   return GetPayloadRun2(data);
  } 
} 
Bool_t AliCTPRawStream::GetPayloadRun2(UChar_t *data)
{
  fClusterMask = data[12];

  fClassMaskNext50 =  ((ULong64_t)data[16] & 0xf) << 46;   // 100..97

  fClassMaskNext50 |= (ULong64_t)data[20] << 38;           // 96..89
  fClassMaskNext50 |= ((ULong64_t)data[21] & 0xF) << 34;   // 88..85

  fClassMaskNext50 |= (ULong64_t)data[24] << 26;           //84..77
  fClassMaskNext50 |= ((ULong64_t)data[25] & 0xF) << 22;   //76..73

  fClassMaskNext50 |= (ULong64_t)data[28] << 14;           //72..65   
  fClassMaskNext50 |= ((ULong64_t)data[29] & 0xF) << 10;   //64..61

  fClassMaskNext50 |= (ULong64_t)data[32] << 2;            //60..53
  fClassMaskNext50 |= ((ULong64_t)data[33] & 0xc) >> 2;    //52..51
  fClassMask       = ((ULong64_t)data[33] & 0x3) << 48;    //50..49

  fClassMask |= (ULong64_t)data[36]  << 40;                //48..41
  fClassMask |= ((ULong64_t)data[37] & 0xf) << 36;         //40..37

  fClassMask |= (ULong64_t)data[40]  << 28;                //36..29
  fClassMask |= ((ULong64_t)data[41] & 0xf) << 24;         //28..25

  fClassMask |= (ULong64_t)data[44]  << 16;                //24..17
  fClassMask |= ((ULong64_t)data[45] & 0xf) << 12;         //16..13

  fClassMask |= (ULong64_t)data[48]  << 4;                //12..5
  fClassMask |= ((ULong64_t)data[49] & 0xf) << 0;         //4..1

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

