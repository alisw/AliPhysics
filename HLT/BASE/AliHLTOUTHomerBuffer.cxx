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
    @brief  HLTOUT data wrapper for AliRawReader.                         */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cerrno>
#include "AliHLTOUTHomerBuffer.h"
#include "AliHLTHOMERReader.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTHomerBuffer)

AliHLTOUTHomerBuffer::AliHLTOUTHomerBuffer(const AliHLTUInt8_t* pBuffer)
  :
  AliHLTOUT(),
  fpBuffer(pBuffer),
  fpReader(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUTHomerBuffer::~AliHLTOUTHomerBuffer()
{
  // see header file for class documentation
}

int AliHLTOUTHomerBuffer::GenerateIndex()
{
  // see header file for class documentation
  int iResult=0;
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
