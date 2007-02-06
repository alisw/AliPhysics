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

/** @file   AliHLTDataSource.cxx
    @author Matthias Richter
    @date   
    @brief  Base class implementation for HLT data source components. */

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTDataSource.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataSource)

AliHLTDataSource::AliHLTDataSource()
{ 
}

AliHLTDataSource::~AliHLTDataSource()
{ 
}

void AliHLTDataSource::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  list.clear(); // there are no input data types
}


int AliHLTDataSource::DoProcessing( const AliHLTComponentEventData& evtData,
				    const AliHLTComponentBlockData* blocks, 
				    AliHLTComponentTriggerData& trigData,
				    AliHLTUInt8_t* outputPtr, 
				    AliHLTUInt32_t& size,
				    AliHLTUInt32_t& outputBlockCnt, 
				    AliHLTComponentBlockData*& outputBlocks,
				    AliHLTComponentEventDoneData*& edd )
{
  int iResult=0;
  if (blocks) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  vector<AliHLTComponentBlockData> blockData;
  if (evtData.fBlockCnt > 0) {
    HLTWarning("Data source component skips input data blocks");
  }
  iResult=GetEvent(evtData, trigData, outputPtr, size, blockData);
  HLTDebug("component %s (%p) GetEvent finished (%d)", GetComponentID(), this, iResult);
  if (iResult>=0) {
    iResult=MakeOutputDataBlockList(blockData, &outputBlockCnt, &outputBlocks);
    if (iResult<0) {
      HLTFatal("component %s (%p): can not convert output block descriptor list", GetComponentID(), this);
    }
  }
  edd = NULL;
  return iResult;
}
