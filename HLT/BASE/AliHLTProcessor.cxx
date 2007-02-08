// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
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

/** @file   AliHLTProcessor.cxx
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Base class implementation for HLT analysis components. */

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTProcessor.h"
#include <string.h>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTProcessor)

AliHLTProcessor::AliHLTProcessor()
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTProcessor::~AliHLTProcessor()
{ 
  // see header file for class documentation
}

int AliHLTProcessor::Init( AliHLTComponentEnvironment* environ, void* environ_param, int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  iResult=AliHLTComponent::Init(environ, environ_param, argc, argv);
  return iResult;
}

int AliHLTProcessor::Deinit()
{
  // see header file for class documentation
  int iResult=0;
  iResult=AliHLTComponent::Deinit();
  return iResult;
}

int AliHLTProcessor::DoProcessing( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
			    AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
			    AliHLTUInt32_t& size, AliHLTUInt32_t& outputBlockCnt, 
			    AliHLTComponentBlockData*& outputBlocks,
			    AliHLTComponentEventDoneData*& edd )
{
  // see header file for class documentation
  int iResult=0;
  vector<AliHLTComponentBlockData> blockData;
  iResult=DoEvent(evtData, blocks, trigData, outputPtr, size, blockData);
  if (iResult>=0) {
    iResult=MakeOutputDataBlockList(blockData, &outputBlockCnt, &outputBlocks);
  }
  edd = NULL;
  return iResult;
}
