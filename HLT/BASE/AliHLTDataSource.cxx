// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
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

/// @file   AliHLTDataSource.cxx
/// @author Matthias Richter
/// @date   
/// @brief  Base class implementation for HLT data source components.
///

#include "AliHLTDataSource.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataSource)

AliHLTDataSource::AliHLTDataSource()
{
  // Base class of HLT data source components.
  // The class provides a common interface for the implementation of HLT data
  // source components.
  // Source components do not consume any input consequently the processing
  // function is called 'GetEvent'.
}

AliHLTDataSource::~AliHLTDataSource()
{ 
  // destructor
}

void AliHLTDataSource::GetInputDataTypes( AliHLTComponentDataTypeList& list)
{
  // default method as source components do not consume input
  list.clear(); // there are no input data types
}


int AliHLTDataSource::DoProcessing( const AliHLTComponentEventData& evtData,
				    const AliHLTComponentBlockData* blocks, 
				    AliHLTComponentTriggerData& trigData,
				    AliHLTUInt8_t* outputPtr, 
				    AliHLTUInt32_t& size,
				    AliHLTComponentBlockDataList& outputBlocks,
				    AliHLTComponentEventDoneData*& edd )
{
  // Processing method, calls child's GetEvent
  int iResult=0;
  if (evtData.fBlockCnt > 0) {
    int unknown=-1;
    for (unsigned int block=0; block<evtData.fBlockCnt; block++) {
      if (blocks[block].fDataType==kAliHLTDataTypeSOR ||
	  blocks[block].fDataType==kAliHLTDataTypeEOR ||
	  blocks[block].fDataType==kAliHLTDataTypeEvent ||
	  blocks[block].fDataType==kAliHLTDataTypeRunType ||
	  blocks[block].fDataType==kAliHLTDataTypeComponentStatistics ||
	  blocks[block].fDataType==kAliHLTDataTypeComponentTable ||
	  blocks[block].fDataType==kAliHLTDataTypeECSParam) {
	continue;
      }
      unknown=block;
      break;
    }
    static int warningCount=0;
    if (unknown>=0 && warningCount++<5) {
      HLTWarning("Data source component skips input data blocks: first unknown block %s",
		 DataType2Text(blocks[unknown].fDataType).c_str());
    }
  }
  iResult=GetEvent(evtData, trigData, outputPtr, size, outputBlocks);
  HLTDebug("component %s (%p) GetEvent finished (%d)", GetComponentID(), this, iResult);
  edd = NULL;
  return iResult;
}

int AliHLTDataSource::GetEvent( const AliHLTComponentEventData& evtData,
				AliHLTComponentTriggerData& trigData,
				AliHLTUInt8_t* /*outputPtr*/, 
				AliHLTUInt32_t& size,
				AliHLTComponentBlockDataList& /*outputBlocks*/ )
{
  // we just forward to the high level method, all other parameters already
  // have been stored internally
  size=0;
  return GetEvent(evtData, trigData);
}

int AliHLTDataSource::GetEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // default method: one of GetEvent methods must be implemented
  HLTFatal("no processing method implemented");
  return -ENOSYS;
}
