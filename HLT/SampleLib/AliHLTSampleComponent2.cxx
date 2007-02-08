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

/** @file   AliHLTSampleComponent2.cxx
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Another sample processing component for the HLT. */

#if __GNUC__== 3
using namespace std;
#endif

#include "AliHLTSampleComponent2.h"

/**
 * The global object used for automatic component registration, 
 * @note DO NOT use this component for calculation!
 */
AliHLTSampleComponent2 gAliHLTSampleComponent2;

ClassImp(AliHLTSampleComponent2)

const AliHLTComponentDataType AliHLTSampleComponent2::fgInputDataTypes[]={kAliHLTVoidDataType,
									{0,"",""}}; //'zero' terminated array

AliHLTSampleComponent2::AliHLTSampleComponent2()
{
  // see header file for class documentation
}

AliHLTSampleComponent2::~AliHLTSampleComponent2()
{
  // see header file for class documentation
}

int AliHLTSampleComponent2::DoInit( int argc, const char** argv ){
  // see header file for class documentation
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component2, DoInit");
  if (argc==0 && argv==NULL) {
    // this is just to get rid of the warning "unused parameter"
  }
  return 0;
}

int AliHLTSampleComponent2::DoDeinit(){
  // see header file for class documentation
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component2, DoDeinit");
  return 0;
}

int AliHLTSampleComponent2::DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
				      AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
				      AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks ) {
  // see header file for class documentation
  Logging(kHLTLogInfo, "HLT", "Sample", "Sample component2, DoEvent");
  if (evtData.fStructSize==0 && blocks==NULL && trigData.fStructSize==0 &&
      outputPtr==0 && size==0)
  {
    outputBlocks.clear();
    // this is just to get rid of the warning "unused parameter"
  }
  return 0;
}
