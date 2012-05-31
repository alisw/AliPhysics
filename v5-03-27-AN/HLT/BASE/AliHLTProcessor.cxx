// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

/** @file   AliHLTProcessor.cxx
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Base class implementation for HLT analysis components. */

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTProcessor.h"
#include <string.h>
#include "TDatime.h"
#include "TString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTProcessor)

AliHLTProcessor::AliHLTProcessor()
  : AliHLTComponent()
  , fpDebugCounters(NULL)
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
  if (fpDebugCounters) delete fpDebugCounters;
  fpDebugCounters=NULL;
}

int AliHLTProcessor::DoProcessing( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
			    AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
			    AliHLTUInt32_t& size,
			    vector<AliHLTComponentBlockData>& outputBlocks,
			    AliHLTComponentEventDoneData*& edd )
{
  // see header file for class documentation
  int iResult=0;
  ReleaseEventDoneData();

  iResult=DoEvent(evtData, blocks, trigData, outputPtr, size, outputBlocks);

  edd = NULL;
  AliHLTComponentEventDoneData* eddTmp = GetCurrentEventDoneData();
  if (eddTmp) {
    int ret = GetEventDoneData(eddTmp->fDataSize, &edd);
    if (ret) {
      HLTError( "Cannot get event done data of %u bytes for event %lu: error (%d)", 
		eddTmp->fDataSize, evtData.fEventID, ret );
      return -ENOMEM;
    }
    edd->fStructSize = sizeof(AliHLTComponentEventDoneData);
    edd->fDataSize = eddTmp->fDataSize;
    edd->fData = reinterpret_cast<AliHLTUInt8_t*>(edd)+edd->fStructSize;
    memcpy( edd->fData, eddTmp->fData, eddTmp->fDataSize );

    // 2009-12-07 want to make this switchable, but this first needs some
    // extension in the online framework to change the log level settings
    // in the component while running
    if (false/*CheckFilter(kHLTLogDebug)*/) {
      if (!fpDebugCounters) {
	fpDebugCounters=new AliHLTProcessorCounters;
      }
      if (fpDebugCounters) {
	int wordCnt=edd->fDataSize/4;
	AliHLTUInt32_t* buffer=reinterpret_cast<AliHLTUInt32_t*>(edd->fData);
	int word=0;
	while (word<wordCnt) {
	  switch (buffer[word]) {
	  case 3: 
	    fpDebugCounters->fReadoutFilter++; 
	    word+=1+buffer[word+1]*4;
	    break;
	  case 4:
	    fpDebugCounters->fMonitoringFilter++; 
	    word+=1+buffer[word+1]*4;
	    break;
	  case 5:
	    fpDebugCounters->fMonitoringEvent++; 
	    break;
	  default:
	    fpDebugCounters->fMismatch++;
	    break;
	  }
	  word++;
	}

	static UInt_t lastTime=0;
	TDatime time;
	if (time.Get()-lastTime>1) {
	  lastTime=time.Get();
	  HLTImportant("EventDoneData size %d: readout %d, monitoring filter %d, monitoring event %d, format error %d", 
		       edd->fDataSize, fpDebugCounters->fReadoutFilter, fpDebugCounters->fMonitoringFilter, fpDebugCounters->fMonitoringEvent, fpDebugCounters->fMismatch);
	  for (int i=0; i< wordCnt; ) {
	    TString message;
	    for (int j=0; j<4 && i<wordCnt; j++) {
	      TString number; number.Form("0x%08x ", buffer[i++]);
	      message+=number;
	    }
	    HLTImportant("   %s", message.Data());
	  }
	}
      }
    }
  }
  return iResult;
}

int AliHLTProcessor::DoEvent( const AliHLTComponentEventData& evtData,
			      const AliHLTComponentBlockData* /*blocks*/, 
			      AliHLTComponentTriggerData& trigData,
			      AliHLTUInt8_t* /*outputPtr*/, 
			      AliHLTUInt32_t& size,
			      vector<AliHLTComponentBlockData>& /*outputBlocks*/ )
{
  // we just forward to the high level method, all other parameters already
  // have been stored internally
  size=0;
  return DoEvent(evtData, trigData);
}

int AliHLTProcessor::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  HLTFatal("no processing method implemented");
  return -ENOSYS;
}
