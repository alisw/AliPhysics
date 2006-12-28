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

/** @file   AliHLTSampleComponent1.cxx
    @author Matthias Richter, Timm M. Steinbeck
    @date   
    @brief  A sample processing component for the HLT. */

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTSampleComponent1.h"

/**
 * The global object used for automatic component registration, 
 * @note DO NOT use this component for calculation!
 */
AliHLTSampleComponent1 gAliHLTSampleComponent1;

ClassImp(AliHLTSampleComponent1)

const AliHLTComponentDataType AliHLTSampleComponent1::inputDataTypes[]={{0,0}, {0,0}}; //'zero' terminated array
const AliHLTComponentDataType AliHLTSampleComponent1::outputDataType={0,0};

AliHLTSampleComponent1::AliHLTSampleComponent1()
{
}

AliHLTSampleComponent1::~AliHLTSampleComponent1()
{
}

int AliHLTSampleComponent1::DoInit( int argc, const char** argv ){
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component1, DoInit");
  return 0;
}

int AliHLTSampleComponent1::DoDeinit(){
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component1, DoDeinit");
  return 0;
}

int AliHLTSampleComponent1::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
				      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ) {
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component1, DoEvent");
  return 0;
}
