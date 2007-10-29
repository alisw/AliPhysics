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

/** @file   AliHLTOUTRawReader.cxx
    @author Matthias Richter
    @date   
    @brief  HLTOUT data wrapper for AliRawReader.                         */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTOUTRawReader.h"
#include "AliHLTHOMERLibManager.h"
#include "AliRawReader.h"
#include "AliHLTHOMERReader.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTRawReader)

AliHLTOUTRawReader::AliHLTOUTRawReader(AliRawReader* pRawreader)
  :
  AliHLTOUTHomerBuffer(NULL, 0),
  fpRawreader(pRawreader),
  fpCurrent(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const int AliHLTOUTRawReader::fgkIdShift=16;

AliHLTOUTRawReader::~AliHLTOUTRawReader()
{
  // see header file for class documentation
  if (fpManager) {
    if (fpCurrent) fpManager->DeleteReader(fpCurrent);
    fpCurrent=NULL;
  }
}

int AliHLTOUTRawReader::GenerateIndex()
{
  // see header file for class documentation
  // step through all HLT ddls, create HOMER reader and
  // scan data block
  int iResult=0;
  if (fpRawreader && fpManager) {
    fpRawreader->Reset();
    fpRawreader->Select("HLT");
    UChar_t* pSrc=NULL;
    while (fpRawreader->ReadNextData(pSrc) && pSrc!=NULL && iResult>=0) {
      AliHLTUInt32_t id=(fpRawreader->GetEquipmentId())<<fgkIdShift;
      int size=fpRawreader->GetDataSize();
      AliHLTHOMERReader* pReader=fpManager->OpenReader(pSrc, size);
      if (pReader) {
	iResult=ScanReader(pReader, id);
	fpManager->DeleteReader(pReader);
      }
    }
  } else {
    iResult=-ENODEV;
  }
  return iResult;
}

int AliHLTOUTRawReader::GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
					AliHLTUInt32_t& size)
{
  // see header file for class documentation
  int iResult=0;
  if (fpManager) {
    AliHLTUInt32_t id=index>>fgkIdShift;
    AliHLTUInt32_t blockNo=index&((0x1<<fgkIdShift)-1);

    // block from the same ddl requested?
    if (fpCurrent && fpRawreader->GetEquipmentId()!=id) {
      fpManager->DeleteReader(fpCurrent);
      fpCurrent=NULL;
    }

    // open ddl for equipment id and create HOMER reader
    if (!fpCurrent) {
      fpRawreader->Reset();
      fpRawreader->SelectEquipment(-1, id, id);
      UChar_t* pSrc=NULL;
      if (fpRawreader->ReadNextData(pSrc) && pSrc!=NULL) {
	fpCurrent=fpManager->OpenReader(pSrc, size);
      } else {
	iResult=-ENOSYS;
      }
    }

    // get data
    if (fpCurrent) {
      if ((pBuffer=static_cast<const AliHLTUInt8_t*>(fpCurrent->GetBlockData(blockNo)))!=NULL) {
	size=fpCurrent->GetBlockDataLength(blockNo);
      } else {
	iResult=-ENOENT;
      }
    }
  } else {
    iResult=-ENODEV;
  }
  return iResult;
}
