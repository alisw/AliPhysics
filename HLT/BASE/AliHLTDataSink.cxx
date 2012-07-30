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

/// @file   AliHLTDataSink.cxx
/// @author Matthias Richter
/// @date   
/// @brief  Base class implementation for HLT data source components.
///

#include "AliHLTDataSink.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDataSink)

AliHLTDataSink::AliHLTDataSink()
{ 
  // Base class of HLT data sink components.
  // The class provides a common interface for the implementation of HLT data
  // sink components.
  // Sink components do not produce any output consequently the processing
  // function is called 'DumpEvent'.
}

AliHLTDataSink::~AliHLTDataSink()
{ 
  // destructor
}

AliHLTComponentDataType AliHLTDataSink::GetOutputDataType()
{
  // default method as sink components do not produce output
  AliHLTComponentDataType dt =
    {sizeof(AliHLTComponentDataType),
     kAliHLTVoidDataTypeID,
     kAliHLTVoidDataOrigin};
  return dt;
}

void AliHLTDataSink::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // default method as sink components do not produce output
  constBase=0;
  inputMultiplier=0;
}

int AliHLTDataSink::DoProcessing( const AliHLTComponentEventData& evtData,
				  const AliHLTComponentBlockData* blocks, 
				  AliHLTComponentTriggerData& trigData,
				  AliHLTUInt8_t* outputPtr, 
				  AliHLTUInt32_t& size,
				  AliHLTComponentBlockDataList& outputBlocks,
				  AliHLTComponentEventDoneData*& edd )
{
  // Processing method, calls child's DumpEvent
  int iResult=0;
  if (outputPtr==NULL
      && size==0 
      && edd==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  outputBlocks.clear();
  iResult=DumpEvent(evtData, blocks, trigData);
  return iResult;
}

int AliHLTDataSink::DumpEvent( const AliHLTComponentEventData& evtData,
			       const AliHLTComponentBlockData* /*blocks*/, 
			       AliHLTComponentTriggerData& trigData )
{
  // we just forward to the high level method, all other parameters already
  // have been stored internally
  return DumpEvent(evtData, trigData);
}

int AliHLTDataSink::DumpEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // default method: one of DumpEvent methods must be implemented
  HLTFatal("no processing method implemented");
  return -ENOSYS;
}
