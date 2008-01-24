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
  return fpParentReader->GetEquipmentId();
}

Bool_t   AliRawReaderHLT::ReadHeader()
{
  // see header file for class documentation
  return fpParentReader->ReadHeader();
}

Bool_t   AliRawReaderHLT::ReadNextData(UChar_t*& data)
{
  // see header file for class documentation
  return fpParentReader->ReadNextData(data);
}

Bool_t   AliRawReaderHLT::ReadNextInt(UInt_t& data)
{
  // see header file for class documentation
  return fpParentReader->ReadNextInt(data);
}

Bool_t   AliRawReaderHLT::ReadNextShort(UShort_t& data)
{
  // see header file for class documentation
  return fpParentReader->ReadNextShort(data);

}

Bool_t   AliRawReaderHLT::ReadNextChar(UChar_t& data)
{
  // see header file for class documentation
  return fpParentReader->ReadNextChar(data);
}

Bool_t   AliRawReaderHLT::ReadNext(UChar_t* data, Int_t size)
{
  // see header file for class documentation
  return fpParentReader->ReadNext(data, size);
}

Bool_t   AliRawReaderHLT::Reset()
{
  // see header file for class documentation
  return fpParentReader->Reset();
}

Bool_t   AliRawReaderHLT::NextEvent()
{
  // see header file for class documentation
  fpParentReader-NextEvent();
}

Bool_t   AliRawReaderHLT::RewindEvents()
{
  // see header file for class documentation
  return fpParentReader->RewindEvents();
}

AliRawReader* AliRawReaderHLTCreateInstance(AliRawReader* pParentReader, const char* options)
{
  // see header file for class documentation
  return new AliRawReaderHLT(pParentReader, options);
}
