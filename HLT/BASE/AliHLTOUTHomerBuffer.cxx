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

/** @file   AliHLTOUTHomerBuffer.cxx
    @author Matthias Richter
    @date   
    @brief  HLTOUT data wrapper for buffer in HOMER format.               */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cerrno>
#include <cassert>
#include "AliHLTOUTHomerBuffer.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTHOMERLibManager.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTHomerBuffer)

AliHLTOUTHomerBuffer::AliHLTOUTHomerBuffer(const AliHLTUInt8_t* pBuffer, int size)
  :
  AliHLTOUT(),
  fpBuffer(pBuffer),
  fSize(size),
  fpReader(NULL),
  fpManager(new AliHLTHOMERLibManager)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  assert(sizeof(homer_uint64)==kAliHLTComponentDataTypefIDsize);
  assert(sizeof(homer_uint32)==kAliHLTComponentDataTypefOriginSize);
  assert(fpManager);
}

AliHLTOUTHomerBuffer::~AliHLTOUTHomerBuffer()
{
  // see header file for class documentation
  if (fpManager) {
    if (fpReader) fpManager->DeleteReader(fpReader);
    delete fpManager;
    fpManager=NULL;
    fpReader=NULL;
  }
}

int AliHLTOUTHomerBuffer::GenerateIndex()
{
  // see header file for class documentation
  int iResult=0;
  if (!fpReader) {
    if (fpManager) {
      fpReader=fpManager->OpenReader(fpBuffer, fSize);
    }
  }
  if (fpReader) {
    iResult=ScanReader(fpReader);
  } else {
    iResult=-ENODEV;
  }
  return iResult;
}

int AliHLTOUTHomerBuffer::GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
					AliHLTUInt32_t& size)
{
  // see header file for class documentation
  int iResult=0;
  if (fpReader) {
    if ((pBuffer=static_cast<const AliHLTUInt8_t*>(fpReader->GetBlockData(index)))!=NULL) {
      size=fpReader->GetBlockDataLength(index);
    } else {
      iResult=-ENOENT;
    }
  } else {
    iResult=-ENODEV;
  }
  return iResult;
}

AliHLTOUT::AliHLTOUTByteOrder_t AliHLTOUTHomerBuffer::CheckBlockByteOrder(AliHLTUInt32_t index)
{
  if (fpReader) {
    return static_cast<AliHLTOUTByteOrder_t>(fpReader->GetBlockByteOrder(index));
  }
  return kInvalidByteOrder;
}

int AliHLTOUTHomerBuffer::CheckBlockAlignment(AliHLTUInt32_t index, AliHLTOUT::AliHLTOUTDataType_t type)
{
  if (fpReader) {
    return fpReader->GetBlockTypeAlignment(index, static_cast<homer_uint8>(type));
  }
  return -ENODATA;
}

int AliHLTOUTHomerBuffer::ScanReader(AliHLTHOMERReader* pReader, AliHLTUInt32_t offset)
{
  // see header file for class documentation
  int iResult=0;
  if (pReader) {
    AliHLTUInt32_t nofBlocks=pReader->GetBlockCnt();
    AliHLTUInt32_t tmp1=0x1;
    AliHLTUInt32_t tmp2=offset;

    // first check if the offset allows to add all data blocks without exceeding the
    // range
    while (nofBlocks<tmp1 && tmp2>0) {
      if (tmp2&0x1) {
	HLTError("index range %#x exceeded for %d data blocks", nofBlocks, offset);
	iResult=-ERANGE;
      }
      tmp2>>1;
      tmp1<<1;
    }

    // loop over data blocks
    HLTDebug("generating index for %d data blocks of reader with offset %#x", nofBlocks, offset);
    for (AliHLTUInt32_t i=0; i<nofBlocks && iResult>=0; i++) {
      homer_uint64 id=pReader->GetBlockDataType( i );
      homer_uint32 origin=pReader->GetBlockDataOrigin( i );
      homer_uint32 spec=pReader->GetBlockDataSpec( i );
      AliHLTComponentDataType dt;
      memcpy(&dt.fID, &id, kAliHLTComponentDataTypefIDsize);
      memcpy(&dt.fOrigin, &origin, kAliHLTComponentDataTypefOriginSize);
      AliHLTOUTBlockDescriptor desc(dt, spec, offset|i);
      iResult=AddBlockDescriptor(desc);
    }
  } else {
    iResult=-ENODEV;
  }
  return iResult;
}

