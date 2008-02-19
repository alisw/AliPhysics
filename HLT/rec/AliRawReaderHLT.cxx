// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliRawReaderHLT.cxx
    @author Matthias Richter
    @date   
    @brief  AliRawReader implementation which replaces original input of
            detectors with the appropriate HLT output.                    */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliRawReaderHLT.h"
#include "AliHLTOUTRawReader.h"
#include "AliHLTModuleAgent.h"
#include "AliHLTOUTHandler.h"
#include "AliHLTOUTHandlerEquId.h"
#include "AliLog.h"
#include "AliDAQ.h"            // RAW, for detector names and equipment ids
#include "TObjString.h"
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliRawReaderHLT)

AliRawReaderHLT::AliRawReaderHLT(AliRawReader* pRawreader, const char* options)
  :
  AliRawReader(),
  fpParentReader(pRawreader),
  fOptions(),
  fpData(NULL),
  fDataSize(0),
  fOffset(0),
  fEquipmentId(-1),
  fbHaveHLTData(false),
  fDetectors(),
  fpHLTOUT(NULL),
  fpDataHandler(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  fOptions=options;
  ScanOptions(options);
}

AliRawReaderHLT::~AliRawReaderHLT()
{
  // see header file for class documentation
  if (fpHLTOUT) {
    if (fpDataHandler) fpDataHandler->ReleaseProcessedData(fpData, fDataSize);
    else fpHLTOUT->ReleaseDataBuffer(fpData);
    fpDataHandler=NULL;
    delete fpHLTOUT;
    fpHLTOUT=NULL;
  }
}

UInt_t AliRawReaderHLT::GetType() const
{
  // see header file for class documentation
  return fpParentReader->GetType();
}

UInt_t AliRawReaderHLT::GetRunNumber() const
{
  // see header file for class documentation
  return fpParentReader->GetRunNumber();
}

const UInt_t* AliRawReaderHLT::GetEventId() const
{
  // see header file for class documentation
  return fpParentReader->GetEventId();
}

const UInt_t* AliRawReaderHLT::GetTriggerPattern() const
{
  // see header file for class documentation
  return fpParentReader->GetTriggerPattern();
}

const UInt_t* AliRawReaderHLT::GetDetectorPattern() const
{
  // see header file for class documentation
  return fpParentReader->GetDetectorPattern();
}

const UInt_t* AliRawReaderHLT::GetAttributes() const
{
  // see header file for class documentation
  return fpParentReader->GetAttributes();
}

const UInt_t* AliRawReaderHLT::GetSubEventAttributes() const
{
  // see header file for class documentation
  return fpParentReader->GetSubEventAttributes();
}

UInt_t AliRawReaderHLT::GetLDCId() const
{
  // see header file for class documentation
  return fpParentReader->GetLDCId();
}

UInt_t AliRawReaderHLT::GetGDCId() const
{
  // see header file for class documentation
  return fpParentReader->GetGDCId();
}

UInt_t AliRawReaderHLT::GetTimestamp() const
{
  // see header file for class documentation
  return fpParentReader->GetTimestamp();
}

const UInt_t* AliRawReaderHLT::GetEquipmentAttributes() const
{
  // see header file for class documentation
  return fpParentReader->GetEquipmentAttributes();
}

Int_t    AliRawReaderHLT::GetEquipmentElementSize() const
{
  // see header file for class documentation
  // don't know what it really means, bu the AliRawReaderFile
  // just sets it to 0
  // do the same if we have a valid equipment data set from
  // the HLT stream
  if (fEquipmentId>=0) return 0;
  return fpParentReader->GetEquipmentElementSize();
}

Int_t    AliRawReaderHLT::GetEquipmentHeaderSize() const
{
  // see header file for class documentation

  // equipment header means the additional data header?
  // if we have a valid equipment data set from the HLT stream
  // there is no additional header
  if (fEquipmentId>=0) return 0;
  return fpParentReader->GetEquipmentHeaderSize();
}

Int_t    AliRawReaderHLT::GetEquipmentSize() const
{
  // see header file for class documentation
  if (fEquipmentId>=0) return fDataSize+sizeof(AliRawDataHeader);
  return fpParentReader->GetEquipmentSize();
}

Int_t    AliRawReaderHLT::GetEquipmentType() const
{
  // see header file for class documentation
  return fpParentReader->GetEquipmentType();
}

Int_t    AliRawReaderHLT::GetEquipmentId() const
{
  // see header file for class documentation
  Int_t id=-1;
  if (fEquipmentId>=0) id=fEquipmentId;
  else id=fpParentReader->GetEquipmentId();
  return id;
}

Bool_t   AliRawReaderHLT::ReadHeader()
{
  // see header file for class documentation
  Bool_t result=fpParentReader->ReadHeader();
  fHeader=const_cast<AliRawDataHeader*>(fpParentReader->GetDataHeader());
  return result;
}

Bool_t   AliRawReaderHLT::ReadNextData(UChar_t*& data)
{
  // see header file for class documentation

  // this function is the backbone of the ReadNext functions, it gets the
  // whole data block either from the HLT stream or the parent raw reader.
  // Each call of ReadNextData directly jumps to the next data set.
  Bool_t result=kFALSE;
  if (fbHaveHLTData&=ReadNextHLTData()) {
    // all internal data variables set
    assert(fpData!=NULL);
    data=const_cast<AliHLTUInt8_t*>(fpData);
    result=kTRUE;
  }
  if (!result) {
    // no data in the HLT stream, read real data
    //AliInfo(Form("read from parent reader: min=%d max=%d", fSelectMinEquipmentId, fSelectMaxEquipmentId));

    // first set the selection back to the original one
    fpParentReader->SelectEquipment(fSelectEquipmentType, fSelectMinEquipmentId, fSelectMaxEquipmentId);

    // read data
    while (result=fpParentReader->ReadNextData(data)) {
      // continue if the Equipment Id is supposed to be replaced by the HLT stream
      // in that case we do not want to read it from the parent raw reader
      if (!IsHLTInput(fpParentReader->GetEquipmentId())) break;
    }

    // set the header of this reader from the parent reader.
    // This is necessary because of a few base class methods working directly
    // on the header
    fHeader=const_cast<AliRawDataHeader*>(fpParentReader->GetDataHeader());
    if (result) {
      fpData=data;
      fDataSize=fpParentReader->GetDataSize();
    } else {
      fpData=NULL;
      fDataSize=0;
    }
    fOffset=0;
    fEquipmentId=-1;
  }
  return result;
}

Bool_t   AliRawReaderHLT::ReadNextInt(UInt_t& data)
{
  // see header file for class documentation
  int iCopy=sizeof(UInt_t);
  UChar_t* dummy=NULL;
  do {
    if (fpData && (fDataSize-fOffset)>=iCopy) {
      data=*reinterpret_cast<const UInt_t*>(fpData+fOffset);
      fOffset+=iCopy;
      return kTRUE;
    }
  } while (ReadNextData(dummy));
  return kFALSE;
}

Bool_t   AliRawReaderHLT::ReadNextShort(UShort_t& data)
{
  // see header file for class documentation
  int iCopy=sizeof(UShort_t);
  UChar_t* dummy=NULL;
  do {
    if (fpData && (fDataSize-fOffset)>=iCopy) {
      data=*reinterpret_cast<const UShort_t*>(fpData+fOffset);
      fOffset+=iCopy;
      return kTRUE;
    }
  } while (ReadNextData(dummy));
  return kFALSE;
}

Bool_t   AliRawReaderHLT::ReadNextChar(UChar_t& data)
{
  // see header file for class documentation
  int iCopy=sizeof(UChar_t);
  UChar_t* dummy=NULL;
  do {
    if (fpData && (fDataSize-fOffset)>=iCopy) {
      data=*reinterpret_cast<const UChar_t*>(fpData+fOffset);
      fOffset+=iCopy;
      return kTRUE;
    }
  } while (ReadNextData(dummy));
  return kFALSE;
}

Bool_t   AliRawReaderHLT::ReadNext(UChar_t* data, Int_t size)
{
  // see header file for class documentation
  UChar_t* dummy=NULL;
  do {
    if (fpData && (fDataSize-fOffset)>=size) {
      // copy remaining data
      int iCopy=fDataSize-fOffset;
      if (iCopy>size) iCopy=size;
      memcpy(data, fpData+fOffset, iCopy);
      fOffset+=iCopy;
      return kTRUE;
    }
  } while (ReadNextData(dummy));
  return kFALSE;
}

Bool_t   AliRawReaderHLT::Reset()
{
  // see header file for class documentation
  Bool_t result=fpParentReader->Reset();
  fpData=NULL;
  fDataSize=0;
  fOffset=0;
  fEquipmentId=-1;
  if (fbHaveHLTData=(fDetectors.size()>0)) {
    vector<int>::iterator detector=fDetectors.begin();
    for (; detector!=fDetectors.end(); detector++) {
      int ddlOffset=AliDAQ::DdlIDOffset(*detector);
      int nofDDLs=AliDAQ::NumberOfDdls(*detector);
      if ((fSelectMinEquipmentId>=0 && fSelectMinEquipmentId>ddlOffset+nofDDLs) ||
	  (fSelectMinEquipmentId>=0 && fSelectMaxEquipmentId<ddlOffset))
	continue;
      break;
    }
    fbHaveHLTData=detector!=fDetectors.end();
  }

  if (fpHLTOUT) {
    if (fpDataHandler) fpDataHandler->ReleaseProcessedData(fpData, fDataSize);
    else fpHLTOUT->ReleaseDataBuffer(fpData);
    fpDataHandler=NULL;
    delete fpHLTOUT;
    fpHLTOUT=NULL;
  }

  return result;
}

Bool_t   AliRawReaderHLT::NextEvent()
{
  // see header file for class documentation
  Bool_t result=fpParentReader->NextEvent();
  if (result) {
    fEventNumber++;
    Reset();
  }
  return result;
}

Bool_t   AliRawReaderHLT::RewindEvents()
{
  // see header file for class documentation
  fEventNumber=-1;
  Reset();
  return fpParentReader->RewindEvents();
}

void AliRawReaderHLT::Select(Int_t detectorID, Int_t minDDLID, Int_t maxDDLID)
{
  // see header file for class documentation
  AliRawReader::Select(detectorID, minDDLID, maxDDLID);
  fpParentReader->Select(detectorID, minDDLID, maxDDLID);
}

// most likely we do not need this method since the base class directly forwards
// to this method
// void AliRawReaderHLT::Select(const char *detectorName, Int_t minDDLID, Int_t maxDDLID)
// {
//   AliInfo(Form("detectorName=%s, minDDLID=%d, maxDDLID=%d", detectorName, minDDLID, maxDDLID));
//   AliRawReader::Select(detectorName, minDDLID, maxDDLID);
//   fpParentReader->Select(detectorName, minDDLID, maxDDLID);
// }

void AliRawReaderHLT::SelectEquipment(Int_t equipmentType, Int_t minEquipmentId, Int_t maxEquipmentId)
{
  // see header file for class documentation

  //AliInfo(Form("equipmentType=%d, minEquipmentId=%d, maxEquipmentId=%d", equipmentType, minEquipmentId, maxEquipmentId));
  AliRawReader::Select(equipmentType, minEquipmentId, maxEquipmentId);
  fpParentReader->Select(equipmentType, minEquipmentId, maxEquipmentId);
}

void AliRawReaderHLT::SkipInvalid(Bool_t skip)
{
  // see header file for class documentation

  AliRawReader::SkipInvalid(skip);
  fpParentReader->SkipInvalid(skip);
}

void AliRawReaderHLT::SelectEvents(Int_t type)
{
  // see header file for class documentation

  //AliInfo(Form("type=%d", type));
  AliRawReader::SelectEvents(type);
  fpParentReader->SelectEvents(type);
}

int AliRawReaderHLT::ScanOptions(const char* options)
{
  // see header file for class documentation
  int iResult=0;
  TString optString(options);
  TString argument;
  TString parameter;
  TObjArray* pTokens=optString.Tokenize(" ");
  if (pTokens) {
    int iEntries=pTokens->GetEntries();
    for (int i =0; i<iEntries; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      // first scan all the other options
      // no other options for the moment

      // it must be a detector name
      int detId=AliDAQ::DetectorID(argument.Data());
      if (detId>=0) {
	fDetectors.push_back(detId);
      }
    }
    delete pTokens;
  }

  return iResult;
}

Bool_t   AliRawReaderHLT::ReadNextHLTData()
{
  // see header file for class documentation
  bool result=kTRUE;
  if (!fpHLTOUT) {
    fpHLTOUT=new AliHLTOUTRawReader(fpParentReader);
    if (result=(fpHLTOUT!=NULL)) {
      if (result=(fpHLTOUT->Init()>=0)) {
	result=fpHLTOUT->SelectFirstDataBlock(kAliHLTAnyDataType, kAliHLTVoidDataSpec,
					      AliHLTModuleAgent::kRawReader)>=0;
      }
    }
  } else {
    // first release the data buffer
    if (fpDataHandler) fpDataHandler->ReleaseProcessedData(fpData, fDataSize);
    else fpHLTOUT->ReleaseDataBuffer(fpData);
    fpDataHandler=NULL;
    if (!(result=fpHLTOUT->SelectNextDataBlock()>=0)) {
      delete fpHLTOUT;
      fpHLTOUT=NULL;
    }
  }
  if (result) {
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    fpHLTOUT->GetDataBlockDescription(dt, spec);
    AliHLTUInt32_t size=0;
    AliHLTOUTHandler* pHandler=fpHLTOUT->GetHandler();
    if (pHandler) {
      if (dynamic_cast<AliHLTOUTHandlerEquId*>(pHandler)!=NULL) {
	AliHLTOUT::AliHLTOUTLockGuard g(fpHLTOUT);
	fEquipmentId=pHandler->ProcessData(fpHLTOUT);
	fpData=NULL;
	fDataSize=pHandler->GetProcessedData(fpData);
	if (!fpData) {
	  result=fpHLTOUT->GetDataBuffer(fpData, size)>=0;
	  AliDebug(AliLog::kDebug, Form("forward data block from HLTOUT stream to equipment %d", fEquipmentId));
	  fDataSize=(int)size;
	} else {
	  // remember the current handler in order to properly release the data buffer
	  fpDataHandler=pHandler;
	  AliDebug(AliLog::kDebug, Form("forward decoded data block provided by handler to equipment %d", fEquipmentId));
	}
      } else {
	AliError(Form("handler is not of type AliHLTOUTHandlerEquId for block %x data type %s spec %#x; data block skipped",
		      fpHLTOUT->GetDataBlockIndex(), AliHLTComponent::DataType2Text(dt).c_str(), spec));
      }
    } else {
      AliWarning(Form("no data handler found for block %x data type %s spec %#x; data block skipped",
		      fpHLTOUT->GetDataBlockIndex(), AliHLTComponent::DataType2Text(dt).c_str(), spec));
    }
  } else {
    fpData=NULL;
    fDataSize=0;
    fOffset=0;
    fEquipmentId=-1;
  }
  return false;
}

Bool_t AliRawReaderHLT::IsHLTInput(int ddlid)
{
  // see header file for class documentation
  vector<int>::iterator detector=fDetectors.begin();
  for (; detector!=fDetectors.end(); detector++) {
    int ddlOffset=AliDAQ::DdlIDOffset(*detector);
    int nofDDLs=AliDAQ::NumberOfDdls(*detector);
    if (ddlid>=ddlOffset && ddlid<ddlOffset+nofDDLs)
      return kTRUE;
  }
  return kFALSE;
}

AliRawReader* AliRawReaderHLTCreateInstance(AliRawReader* pParentReader, const char* options)
{
  // see header file for class documentation
  if (!pParentReader) return NULL;
  return new AliRawReaderHLT(pParentReader, options);
}
