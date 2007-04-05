// @(#) $Id$

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

/** @file   AliHLTLogging.cxx
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Implementation of HLT logging primitives.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#ifndef NOALIROOT_LOGGING
#include "AliLog.h"
#endif
#include "AliHLTStdIncludes.h"
#include "AliHLTLogging.h"
#include "TString.h"
#include "TArrayC.h"
#include "Varargs.h"
#include <string>
#include <sstream>
#include <iostream>

/** target stream for AliRoot logging methods */
ostringstream gLogstr;

#ifndef NOALIROOT_LOGGING
/**
 * Notification callback for AliRoot logging methods
 */
void LogNotification(AliLog::EType_t level, const char* message)
{
  AliHLTLogging hltlog;
  hltlog.SwitchAliLog(0);
  hltlog.Logging(kHLTLogInfo, "NotificationHandler", "AliLog", gLogstr.str().c_str());
  gLogstr.clear();
  string empty("");
  gLogstr.str(empty);
}
#endif

/**
 * The global logging buffer.
 * The buffer is created with an initial size and grown dynamically on
 * demand.
 */
TArrayC gAliHLTLoggingTarget(200);

/** the maximum size of the buffer */
const int gALIHLTLOGGING_MAXBUFFERSIZE=10000;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTLogging)

AliHLTLogging::AliHLTLogging()
  :
  //fLocalLogFilter(kHLTLogDefault),
  fLocalLogFilter(kHLTLogAll),
  fpDefaultKeyword(NULL),
  fpCurrentKeyword(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTLogging::AliHLTLogging(const AliHLTLogging&)
  :
  fLocalLogFilter(kHLTLogAll),
  fpDefaultKeyword(NULL),
  fpCurrentKeyword(NULL)
{
  // see header file for class documentation
  HLTFatal("copy constructor untested");
}

AliHLTLogging& AliHLTLogging::operator=(const AliHLTLogging&)
{ 
  // see header file for class documentation
  HLTFatal("assignment operator untested");
  return *this;
}

AliHLTComponentLogSeverity AliHLTLogging::fGlobalLogFilter=kHLTLogAll;
AliHLTfctLogging AliHLTLogging::fLoggingFunc=NULL;
int AliHLTLogging::fgUseAliLog=1;

AliHLTLogging::~AliHLTLogging()
{
  // see header file for class documentation
}

int AliHLTLogging::Init(AliHLTfctLogging pFun) 
{
  // see header file for class documentation
  if (fLoggingFunc!=NULL && fLoggingFunc!=pFun) {
    (*fLoggingFunc)(NULL/*fParam*/, kHLTLogWarning, "AliHLTLogging::Init", "no key", "overriding previously initialized logging function");    
  }
  fLoggingFunc=pFun;
  // older versions of AliLog does not support the notification callback and
  // stringstreams, but they support the logging macros in general
#ifndef NOALIROOT_LOGGING
#ifndef NO_ALILOG_NOTIFICATION
  AliLog* log=new AliLog;
  log->SetLogNotification(LogNotification);
  log->SetStreamOutput(&gLogstr);
#endif // NO_ALILOG_NOTIFICATION
#endif // NOALIROOT_LOGGING
  
  return 0;
}

int AliHLTLogging::Message(void *param, AliHLTComponentLogSeverity severity,
			   const char* origin, const char* keyword,
			   const char* message) 
{
  // see header file for class documentation
  int iResult=0;
  if (param==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }

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
  TString out="HLT Log ";
  out+=strSeverity;
  if (origin && origin[0]!=0) {out+=": <"; out+=origin; out+="> ";}
  out+=" "; out+=message;
  if (keyword!=NULL && strcmp(keyword, HLT_DEFAULT_LOG_KEYWORD)!=0) {
    out+=" ("; out+=keyword; out +=")";
  }
  cout << out.Data() << endl;
  return iResult;
}

#ifndef NOALIROOT_LOGGING
int AliHLTLogging::AliMessage(AliHLTComponentLogSeverity severity, 
			      const char* origin_class, const char* origin_func,
			      const char* file, int line, const char* message) 
{
  // see header file for class documentation

  switch (severity) {
  case kHLTLogBenchmark: 
    AliLog::Message(AliLog::kInfo, message, "HLT", origin_class, origin_func, file, line);
    break;
  case kHLTLogDebug:
    AliLog::Message(AliLog::kDebug, message, "HLT", origin_class, origin_func, file, line);
    break;
  case kHLTLogInfo:
    AliLog::Message(AliLog::kInfo, message, "HLT", origin_class, origin_func, file, line);
    break;
  case kHLTLogWarning:
    AliLog::Message(AliLog::kWarning, message, "HLT", origin_class, origin_func, file, line);
    break;
  case kHLTLogError:
    AliLog::Message(AliLog::kError, message, "HLT", origin_class, origin_func, file, line);
    break;
  case kHLTLogFatal:
    AliLog::Message(AliLog::kWarning, message, "HLT", origin_class, origin_func, file, line);
    break;
  default:
    break;
  }
  return 0;
}
#endif

const char* AliHLTLogging::BuildLogString(const char *format, va_list ap) 
{
  // see header file for class documentation

  int iResult=0;
#ifdef R__VA_COPY
  va_list bap;
  R__VA_COPY(bap, ap);
#endif //R__VA_COPY

  // take the first argument from the list as format string if no
  // format was given
  const char* fmt = format;
  if (fmt==NULL) fmt=va_arg(ap, const char*);

  gAliHLTLoggingTarget[0]=0;
  while (fmt!=NULL) {
    iResult=vsnprintf(gAliHLTLoggingTarget.GetArray(), gAliHLTLoggingTarget.GetSize(), fmt, ap);
    if (iResult==-1)
      // for compatibility with older version of vsnprintf
      iResult=gAliHLTLoggingTarget.GetSize()*2;
    else if (iResult<gAliHLTLoggingTarget.GetSize())
      break;

    // terminate if buffer is already at the limit
    if (gAliHLTLoggingTarget.GetSize()>=gALIHLTLOGGING_MAXBUFFERSIZE) 
    {
      gAliHLTLoggingTarget[gAliHLTLoggingTarget.GetSize()-1]=0;
      break;
    }

    // check limitation and grow the buffer
    if (iResult>gALIHLTLOGGING_MAXBUFFERSIZE) iResult=gALIHLTLOGGING_MAXBUFFERSIZE;
    gAliHLTLoggingTarget.Set(iResult+1);

    // copy the original list and skip the first argument if this was the format string
#ifdef R__VA_COPY
    va_end(ap);
    R__VA_COPY(ap, bap);
#else
    gAliHLTLoggingTarget[gAliHLTLoggingTarget.GetSize()-1]=0;
    break;
#endif //R__VA_COPY
    if (format==NULL) va_arg(ap, const char*);
  }     
#ifdef R__VA_COPY
  va_end(bap);
#endif //R__VA_COPY

  return gAliHLTLoggingTarget.GetArray();
}

int AliHLTLogging::Logging(AliHLTComponentLogSeverity severity,
			   const char* origin, const char* keyword,
			   const char* format, ... ) 
{
  // see header file for class documentation
  int iResult=CheckFilter(severity);
  if (iResult>0) {
    va_list args;
    va_start(args, format);
    if (fLoggingFunc) {
      iResult = (*fLoggingFunc)(NULL/*fParam*/, severity, origin, keyword, AliHLTLogging::BuildLogString(format, args ));
    } else {
#ifndef NOALIROOT_LOGGING
      if (fgUseAliLog)
	iResult=AliMessage(severity, NULL, origin, NULL, 0, AliHLTLogging::BuildLogString(format, args ));
      else
#endif
        iResult=Message(NULL/*fParam*/, severity, origin, keyword, AliHLTLogging::BuildLogString(format, args ));
    }
    va_end(args);
  }
  return iResult;
}

int AliHLTLogging::LoggingVarargs(AliHLTComponentLogSeverity severity, 
				  const char* origin_class, const char* origin_func,
				  const char* file, int line,  ... ) const
{
  // see header file for class documentation

  if (file==NULL && line==0) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  int iResult=CheckFilter(severity);
  if (iResult>0) {
    const char* separator="";
    TString origin;
    if (origin_class) {
	origin+=origin_class;
	separator="::";
    }
    if (origin_func) {
	origin+=separator;
	origin+=origin_func;
    }
    va_list args;
    va_start(args, line);

    if (fLoggingFunc) {
      iResult=(*fLoggingFunc)(NULL/*fParam*/, severity, origin.Data(), GetKeyword(), AliHLTLogging::BuildLogString(NULL, args ));
    } else {
#ifndef NOALIROOT_LOGGING
      if (fgUseAliLog)
	iResult=AliMessage(severity, origin_class, origin_func, file, line, AliHLTLogging::BuildLogString(NULL, args ));
      else
#endif
	iResult=Message(NULL/*fParam*/, severity, origin.Data(), GetKeyword(), AliHLTLogging::BuildLogString(NULL, args ));
    }
    va_end(args);
  }
  return iResult;
}

int AliHLTLogging::CheckFilter(AliHLTComponentLogSeverity severity) const
{
  // see header file for class documentation

  int iResult=severity==kHLTLogNone || (severity&fGlobalLogFilter)>0 && (severity&fLocalLogFilter)>0;
  return iResult;
}

void AliHLTLogging::SetGlobalLoggingLevel(AliHLTComponentLogSeverity level)
{
  // see header file for class documentation

  fGlobalLogFilter=level;
}

void AliHLTLogging::SetLocalLoggingLevel(AliHLTComponentLogSeverity level)
{
  // see header file for class documentation

  fLocalLogFilter=level;
}
