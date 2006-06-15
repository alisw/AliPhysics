// $Id: 

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
// HLT logging tools                                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if __GNUC__>= 3
using namespace std;
#endif

#include <cerrno>
#include "AliL3StandardIncludes.h"
#include "AliHLTLogging.h"
#include <cstdarg>
#include <string>

ClassImp(AliHLTLogging)

char AliHLTLogging::fLogBuffer[LOG_BUFFER_SIZE]="";
char AliHLTLogging::fOriginBuffer[LOG_BUFFER_SIZE]="";

AliHLTComponent_LogSeverity AliHLTLogging::fGlobalLogFilter=kHLTLogAll;
AliHLTfctLogging AliHLTLogging::fLoggingFunc=NULL;

AliHLTLogging::AliHLTLogging()
{
  fpDefaultKeyword=NULL;
  fpCurrentKeyword=NULL;
  //fLocalLogFilter=kHLTLogDefault;
  fLocalLogFilter=kHLTLogAll;
}


AliHLTLogging::~AliHLTLogging()
{
}

int AliHLTLogging::Init(AliHLTfctLogging pFun) 
{ 
  if (fLoggingFunc!=NULL && fLoggingFunc!=pFun) {
    (*fLoggingFunc)(NULL/*fParam*/, kHLTLogWarning, "AliHLTLogging::Init", "no key", "overriding previously initialized logging function");    
  }
  fLoggingFunc=pFun; 
  return 0;
}

int AliHLTLogging::Message(void *param, AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message) {
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
  cout << "HLT Log " << strSeverity << ": " << origin << " " << message;
  if (strcmp(keyword, HLT_DEFAULT_LOG_KEYWORD)!=0)
    cout << " (" << keyword << ")";
  cout << endl;
  return iResult;
}

const char* AliHLTLogging::BuildLogString(const char *format, va_list ap) {
  int tgtLen=0;
  int iBufferSize=LOG_BUFFER_SIZE;
  char* tgtBuffer=fLogBuffer;
  tgtBuffer[tgtLen]=0;

#if (defined LOG_PREFIX)
  tgtLen = snprintf(tgtBuffer, iBufferSize, LOG_PREFIX); // add logging prefix
#endif
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

int AliHLTLogging::Logging(AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* format, ... ) {
  int iResult=CheckFilter(severity);
  if (iResult>0) {
    va_list args;
    va_start(args, format);
    if (fLoggingFunc) {
      iResult = (*fLoggingFunc)(NULL/*fParam*/, severity, origin, keyword, AliHLTLogging::BuildLogString(format, args ));
    } else {
      iResult = Message(NULL/*fParam*/, severity, origin, keyword, AliHLTLogging::BuildLogString(format, args ));
    }
  }
  return iResult;
}

int AliHLTLogging::LoggingVarargs( AliHLTComponent_LogSeverity severity, const char* origin_class, const char* origin_func,  ... )
{
  int iResult=CheckFilter(severity);
  if (iResult>0) {
    int iMaxSize=LOG_BUFFER_SIZE-1;
    int iPos=0;
    const char* separator="";
    fOriginBuffer[iPos]=0;
    if (origin_class) {
      if ((int)strlen(origin_class)<iMaxSize-iPos) {
	strcpy(&fOriginBuffer[iPos], origin_class);
	iPos+=strlen(origin_class);
	separator="::";
      }
    }
    if (origin_func) {
      if ((int)strlen(origin_func)+(int)strlen(separator)<iMaxSize-iPos) {
	strcpy(&fOriginBuffer[iPos], separator);
	iPos+=strlen(separator);
	strcpy(&fOriginBuffer[iPos], origin_func);
	iPos+=strlen(origin_func);
      }
    }
    va_list args;
    va_start(args, origin_func);
    const char* format = va_arg(args, const char*);

    const char* message=format;
    char* qualifier=NULL;
    if ((qualifier=strchr(format, '%'))!=NULL) {
      message=AliHLTLogging::BuildLogString(format, args);
    }
    if (fLoggingFunc) {
      iResult=(*fLoggingFunc)(NULL/*fParam*/, severity, fOriginBuffer, GetKeyword(), message);
    } else {
      iResult=Message(NULL/*fParam*/, severity, fOriginBuffer, GetKeyword(), message);
    }
    va_end(args);
  }
  return iResult;
}

int AliHLTLogging::CheckFilter(AliHLTComponent_LogSeverity severity)
{
  int iResult=severity==kHLTLogNone || (severity&fGlobalLogFilter)>0 && (severity&fLocalLogFilter)>0;
  return iResult;
}
