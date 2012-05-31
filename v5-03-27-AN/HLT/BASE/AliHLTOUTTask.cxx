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

/** @file   AliHLTOUTTask.cxx
    @author Matthias Richter
    @date   
    @brief  A special HLTOUT sibling working as a data sink in chains
*/

#include "AliHLTOUTTask.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTTask)

AliHLTOUTTask::AliHLTOUTTask(const char* chains)
  :
  AliHLTOUT(),
  AliHLTDumpTask(chains)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUTTask::~AliHLTOUTTask()
{
  // see header file for class documentation
}

int AliHLTOUTTask::GenerateIndex()
{
  // see header file for class documentation
  int iResult=0;
  const AliHLTComponentBlockDataList& list=GetDataBlocks();
  for (unsigned int i=0; i<list.size(); i++) {
    AliHLTOUTBlockDescriptor desc(list[i].fDataType, list[i].fSpecification, i, this);
    //HLTDebug("adding block %d: %s %#x", i, AliHLTComponent::DataType2Text(list[i].fDataType).c_str(), list[i].fSpecification);
    iResult=AddBlockDescriptor(desc);
  }
  return iResult;
}

int AliHLTOUTTask::GetDataBuffer(AliHLTUInt32_t index, const AliHLTUInt8_t* &pBuffer, 
				 AliHLTUInt32_t& size)
{
  // see header file for class documentation
  int iResult=0;
  const AliHLTComponentBlockDataList& list=GetDataBlocks();
  if (index>=list.size()) return -ENOENT;
  pBuffer=reinterpret_cast<AliHLTUInt8_t*>(list[index].fPtr);
  size=list[index].fSize;
  return iResult;
}

AliHLTOUT::AliHLTOUTByteOrder AliHLTOUTTask::CheckBlockByteOrder(AliHLTUInt32_t /*index*/)
{
  // see header file for class documentation
  return kInvalidByteOrder;
}

int AliHLTOUTTask::CheckBlockAlignment(AliHLTUInt32_t /*index*/, AliHLTOUT::AliHLTOUTDataType /*type*/)
{
  // see header file for class documentation
  int iResult=0;
  return iResult;
}

int AliHLTOUTTask::ResetInput()
{
  // see header file for class documentation
  return ReleaseDataBlocks();
}
