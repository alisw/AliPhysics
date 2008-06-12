// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTOUTRawReader.cxx
    @author Matthias Richter
    @date   
    @brief  HLTOUT data wrapper for AliRawReader.                         
*/

#include "AliHLTOUTRawReader.h"
#include "AliRawReader.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTRawReader)

AliHLTOUTRawReader::AliHLTOUTRawReader()
  :
  AliHLTOUTHomerCollection(),
  fpRawreader(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUTRawReader::AliHLTOUTRawReader(AliRawReader* pRawreader, int event, AliHLTEsdManager* pEsdManager)
  :
  AliHLTOUTHomerCollection(event, pEsdManager),
  fpRawreader(pRawreader)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUTRawReader::~AliHLTOUTRawReader()
{
  // see header file for class documentation
}

Bool_t AliHLTOUTRawReader::ReadNextData(UChar_t*& data)
{
  // see header file for class documentation
  if (!fpRawreader) return kFALSE;
  return fpRawreader->ReadNextData(data);
}

int AliHLTOUTRawReader::Reset()
{
  // see header file for class documentation
  if (fpRawreader) return fpRawreader->Reset();
  return 0;
}

int AliHLTOUTRawReader::GetDataSize()
{
  // see header file for class documentation
  if (fpRawreader) return fpRawreader->GetDataSize();
  return 0;
}

const AliRawDataHeader* AliHLTOUTRawReader::GetDataHeader()
{
  // see header file for class documentation
  if (fpRawreader) return fpRawreader->GetDataHeader();
  return NULL;
}

void AliHLTOUTRawReader::SelectEquipment(int equipmentType, int minEquipmentId, int maxEquipmentId)
{
  // see header file for class documentation
  if (fpRawreader) fpRawreader->SelectEquipment(equipmentType, minEquipmentId, maxEquipmentId);
}

int AliHLTOUTRawReader::GetEquipmentId()
{
  // see header file for class documentation
  if (fpRawreader) return fpRawreader->GetEquipmentId();
  return -1;
}
