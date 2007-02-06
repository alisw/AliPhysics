// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTDataSink.cxx
    @author Matthias Richter
    @date   
    @brief  Base class implementation for HLT data source components. */

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTDataSink.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataSink)

AliHLTDataSink::AliHLTDataSink()
{ 
}

AliHLTDataSink::~AliHLTDataSink()
{ 
}

AliHLTComponentDataType AliHLTDataSink::GetOutputDataType()
{
  AliHLTComponentDataType dt =
    {sizeof(AliHLTComponentDataType),
     kAliHLTVoidDataTypeID,
     kAliHLTVoidDataOrigin};
  return dt;
}

void AliHLTDataSink::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  constBase=0;
  inputMultiplier=0;
}

int AliHLTDataSink::DoProcessing( const AliHLTComponentEventData& evtData,
				    const AliHLTComponentBlockData* blocks, 
				    AliHLTComponentTriggerData& trigData,
				    AliHLTUInt8_t* outputPtr, 
				    AliHLTUInt32_t& size,
				    AliHLTUInt32_t& outputBlockCnt, 
				    AliHLTComponentBlockData*& outputBlocks,
				    AliHLTComponentEventDoneData*& edd )
{
  int iResult=0;
  if (outputPtr==NULL
      && size==0 
      && outputBlockCnt==0 
      && outputBlocks==NULL
      && edd==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  vector<AliHLTComponentBlockData> blockData;
  iResult=DumpEvent(evtData, blocks, trigData);
  return iResult;
}
