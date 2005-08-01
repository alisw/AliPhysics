// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          Artur Szostak <artursz@iafrica.com>                           *
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
// global HLT module management                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__== 3
using namespace std;
#endif

#include <errno.h>
#include <string.h>
#include "AliL3StandardIncludes.h"
#include "AliHLTSystem.h"
#include "AliHLTComponentHandler.h"
#include "AliHLTComponent.h"

ClassImp(AliHLTSystem)

char AliHLTSystem::fLogBuffer[LOG_BUFFER_SIZE]="";

AliHLTSystem::AliHLTSystem()
{
  fpComponentHandler=new AliHLTComponentHandler();
  if (fpComponentHandler) {
    AliHLTComponentEnvironment env;
    memset(&env, 0, sizeof(AliHLTComponentEnvironment));
    env.fLoggingFunc=AliHLTSystem::Logging;
    fpComponentHandler->SetEnvironment(&env);
  }
}


AliHLTSystem::~AliHLTSystem()
{
}

int AliHLTSystem::Logging(void *param, AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message) {
  int iResult=0;
  const char* strSeverity="";
  switch (severity) {
  case kHLTLogBenchmark: 
    strSeverity="benchmark";
    break;
  case kHLTLogDebug:
    strSeverity="debug";
    break;
  case kHLTLogInfo:
    strSeverity="info";
    break;
  case kHLTLogWarning:
    strSeverity="warning";
    break;
  case kHLTLogError:
    strSeverity="error";
    break;
  case kHLTLogFatal:
    strSeverity="fatal";
    break;
  default:
    break;
  }
  cout << "HLT Log " << strSeverity << ": " << origin << " (" << keyword << ") " << message << endl;
  return iResult;
}

const char* AliHLTSystem::BuildLogString(const char *format, va_list ap) {
  int tgtLen=0;
  int iBufferSize=LOG_BUFFER_SIZE;
  char* tgtBuffer=fLogBuffer;
  tgtBuffer[tgtLen]=0;

  tgtLen = snprintf(tgtBuffer, iBufferSize, LOG_PREFIX); // add logging prefix
  if (tgtLen>=0) {
    tgtBuffer+=tgtLen; iBufferSize-=tgtLen;
    tgtLen = vsnprintf(tgtBuffer, iBufferSize, format, ap);
    if (tgtLen>0) {
      tgtBuffer+=tgtLen;
//       if (tgtLen<LOG_BUFFER_SIZE-1) {
// 	*tgtBuffer++='\n'; // add newline if space in buffer
//      }
      *tgtBuffer=0; // terminate the buffer
    }
  }
  return fLogBuffer;
}


