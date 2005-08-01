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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// a sample processing component for the HLT                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__== 3
using namespace std;
#endif

#include "AliL3StandardIncludes.h"
#include "AliHLTSampleComponent2.h"

// this is a global object used for automatic component registration, do not use this
AliHLTSampleComponent2 gAliHLTSampleComponent2;

ClassImp(AliHLTSampleComponent2)

const AliHLTComponent_DataType AliHLTSampleComponent2::inputDataTypes[]={{0,0}, {0,0}}; //'zero' terminated array
const AliHLTComponent_DataType AliHLTSampleComponent2::outputDataType={0,0};

AliHLTSampleComponent2::AliHLTSampleComponent2()
{
}

AliHLTSampleComponent2::~AliHLTSampleComponent2()
{
}

int AliHLTSampleComponent2::DoInit( int argc, const char** argv ){
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component2, DoInit");
  return 0;
}

int AliHLTSampleComponent2::DoDeinit(){
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component2, DoDeinit");
  return 0;
}

int AliHLTSampleComponent2::DoEvent( AliHLTComponent_EventData evtData, AliHLTComponent_BlockData* blocks, 
				      AliHLTComponent_TriggerData trigData, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t* size, vector<AliHLTComponent_BlockData>& outputBlocks ) {
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component2, DoEvent");
  return 0;
}
