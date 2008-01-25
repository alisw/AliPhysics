// $Id: AliRawReaderHLT.cxx,v 1.3 2007/11/15 18:12:44 szostak Exp $

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
#include "AliLog.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliRawReaderHLT)

AliRawReaderHLT::AliRawReaderHLT(AliRawReader* pRawreader, const char* options)
  :
  AliRawReader(),
  fpParentReader(pRawreader),
  fOptions()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  fOptions=options;
}

AliRawReaderHLT::~AliRawReaderHLT()
{
  // see header file for class documentation
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
  return fpParentReader->GetEquipmentElementSize();
}

Int_t    AliRawReaderHLT::GetEquipmentHeaderSize() const
{
  // see header file for class documentation
  return fpParentReader->GetEquipmentHeaderSize();
}

Int_t    AliRawReaderHLT::GetEquipmentSize() const
{
  // see header file for class documentation
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
  Int_t id=fpParentReader->GetEquipmentId();
  //AliInfo(Form("id=%d",id));
  return id;
}

Bool_t   AliRawReaderHLT::ReadHeader()
{
  // see header file for class documentation
  Bool_t result=fpParentReader->ReadHeader();
  fHeader=const_cast<AliRawDataHeader*>(GetDataHeader());
  return result;
}

Bool_t   AliRawReaderHLT::ReadNextData(UChar_t*& data)
{
  // see header file for class documentation
  Bool_t result=fpParentReader->ReadNextData(data);
  fHeader=const_cast<AliRawDataHeader*>(GetDataHeader());
  return result;
}

Bool_t   AliRawReaderHLT::ReadNextInt(UInt_t& data)
{
  // see header file for class documentation
  Bool_t result=fpParentReader->ReadNextInt(data);
  fHeader=const_cast<AliRawDataHeader*>(GetDataHeader());
  return result;
}

Bool_t   AliRawReaderHLT::ReadNextShort(UShort_t& data)
{
  // see header file for class documentation
  Bool_t result=fpParentReader->ReadNextShort(data);
  fHeader=const_cast<AliRawDataHeader*>(GetDataHeader());
  return result;
}

Bool_t   AliRawReaderHLT::ReadNextChar(UChar_t& data)
{
  // see header file for class documentation
  Bool_t result=fpParentReader->ReadNextChar(data);
  fHeader=const_cast<AliRawDataHeader*>(GetDataHeader());
  return result;
}

Bool_t   AliRawReaderHLT::ReadNext(UChar_t* data, Int_t size)
{
  // see header file for class documentation
  Bool_t result=fpParentReader->ReadNext(data, size);
  fHeader=const_cast<AliRawDataHeader*>(GetDataHeader());
  return result;
}

Bool_t   AliRawReaderHLT::Reset()
{
  // see header file for class documentation
  return fpParentReader->Reset();
}

Bool_t   AliRawReaderHLT::NextEvent()
{
  // see header file for class documentation
  //AliInfo(Form("SelectEquipment: type=%d min=%d max=%d", fSelectEquipmentType, fSelectMinEquipmentId, fSelectMaxEquipmentId));
  //fpParentReader->SelectEquipment(fSelectEquipmentType, fSelectMinEquipmentId, fSelectMaxEquipmentId);
  Bool_t result=fpParentReader->NextEvent();
  if (result) fEventNumber++;
  AliInfo(Form("event %d", fEventNumber));
  return result;
}

Bool_t   AliRawReaderHLT::RewindEvents()
{
  // see header file for class documentation
  fEventNumber=-1;
  return fpParentReader->RewindEvents();
}

void AliRawReaderHLT::Select(Int_t detectorID, Int_t minDDLID, Int_t maxDDLID)
{
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
  AliInfo(Form("equipmentType=%d, minEquipmentId=%d, maxEquipmentId=%d", equipmentType, minEquipmentId, maxEquipmentId));
  AliRawReader::Select(equipmentType, minEquipmentId, maxEquipmentId);
  fpParentReader->Select(equipmentType, minEquipmentId, maxEquipmentId);
}

void AliRawReaderHLT::SkipInvalid(Bool_t skip)
{
  AliRawReader::SkipInvalid(skip);
  fpParentReader->SkipInvalid(skip);
}

void AliRawReaderHLT::SelectEvents(Int_t type)
{
  AliInfo(Form("type=%d", type));
  AliRawReader::SelectEvents(type);
  fpParentReader->SelectEvents(type);
}

AliRawReader* AliRawReaderHLTCreateInstance(AliRawReader* pParentReader, const char* options)
{
  // see header file for class documentation
  if (!pParentReader) return NULL;
  return new AliRawReaderHLT(pParentReader, options);
}
