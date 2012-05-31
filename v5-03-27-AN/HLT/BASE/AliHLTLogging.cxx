// @(#) $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTLogging.cxx
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Implementation of HLT logging primitives.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTStdIncludes.h"
#include "AliHLTLogging.h"
#include "AliHLTComponentHandler.h"
#include "TString.h"
#include "Varargs.h"
#include <string>
#include <sstream>
#include <iostream>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTLogging);

AliHLTLogging::AliHLTLogging()
  :
  fLocalLogFilter(fgLocalLogDefault),
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

ostringstream AliHLTLogging::fgLogstr;
AliHLTComponentLogSeverity AliHLTLogging::fgGlobalLogFilter=kHLTLogAll;
AliHLTComponentLogSeverity AliHLTLogging::fgLocalLogDefault=kHLTLogAll;
AliHLTfctLogging AliHLTLogging::fgLoggingFunc=NULL;
AliHLTLogging::AliHLTDynamicMessage AliHLTLogging::fgAliLoggingFunc=NULL;
int AliHLTLogging::fgUseAliLog=1;

TString AliHLTLogging::fgBlackList="";
TString AliHLTLogging::fgWhiteList="";

AliHLTLogging::~AliHLTLogging()
{
  // see header file for class documentation
}

// the array will be grown dynamically, this is just an initial size 
TArrayC AliHLTLogging::fgAliHLTLoggingTarget(200);
// the maximum size of the array
const int AliHLTLogging::fgkALIHLTLOGGINGMAXBUFFERSIZE=10000;

int AliHLTLogging::Init(AliHLTfctLogging pFun) 
{
  // see header file for class documentation
  if (fgLoggingFunc!=NULL && fgLoggingFunc!=pFun) {
    (*fgLoggingFunc)(NULL/*fParam*/, kHLTLogWarning, "AliHLTLogging::Init", "no key", "overriding previously initialized logging function");    
  }
  fgLoggingFunc=pFun;
  
  return 0;
}

int AliHLTLogging::InitAliLogTrap(AliHLTComponentHandler* pHandler)
{
  // see header file for class documentation
  // init the AliRoot logging trap
  // AliLog messages are redirected to PubSub,
  int iResult=0;
  if (pHandler) {
    // set temporary loglevel of component handler
    AliHLTComponentLogSeverity loglevel=pHandler->GetLocalLoggingLevel();
    pHandler->SetLocalLoggingLevel(kHLTLogError);

    // load library containing AliRoot dependencies and initialization handler
    pHandler->LoadLibrary(ALILOG_WRAPPER_LIBRARY, 0/* do not activate agents */);

    // restore loglevel
    pHandler->SetLocalLoggingLevel(loglevel);

    // find the symbol
    InitAliDynamicMessageCallback pFunc=(InitAliDynamicMessageCallback)pHandler->FindSymbol(ALILOG_WRAPPER_LIBRARY, "InitAliDynamicMessageCallback"); 
    if (pFunc) {
      iResult=(*pFunc)();
    } else {
      Message(NULL, kHLTLogError, "AliHLTLogging::InitAliLogTrap", "init logging",
	      "can not initialize AliLog callback");
      iResult=-ENOSYS;
    }
  } else {
    iResult=-EINVAL;
  }
  
  return iResult;
}

int AliHLTLogging::InitAliLogFunc(AliHLTComponentHandler* pHandler)
{
  // see header file for class documentation
  int iResult=0;
  if (pHandler) {
    // set temporary loglevel of component handler
    AliHLTComponentLogSeverity loglevel=pHandler->GetLocalLoggingLevel();
    pHandler->SetLocalLoggingLevel(kHLTLogError);

    // load library containing AliRoot dependencies and initialization handler
    pHandler->LoadLibrary(ALILOG_WRAPPER_LIBRARY, 0/* do not activate agents */);

    // restore loglevel
    pHandler->SetLocalLoggingLevel(loglevel);

    // find the symbol
    fgAliLoggingFunc=(AliHLTLogging::AliHLTDynamicMessage)pHandler->FindSymbol(ALILOG_WRAPPER_LIBRARY, "AliDynamicMessage");
    if (fgAliLoggingFunc==NULL) {
      Message(NULL, kHLTLogError, "AliHLTLogging::InitAliLogFunc", "init logging",
	      "symbol lookup failure: can not find AliDynamicMessage, switching to HLT logging system");
      iResult=-ENOSYS;
    }
  } else {
    iResult=-EINVAL;
  }
  
  return iResult;
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
  case kHLTLogImportant:
    strSeverity="notify";
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

#if 0
int AliHLTLogging::AliMessage(AliHLTComponentLogSeverity severity, 
			      const char* originClass, const char* originFunc,
			      const char* file, int line, const char* message) 
{
  // see header file for class documentation

  switch (severity) {
  case kHLTLogBenchmark: 
    AliLog::Message(AliLog::kInfo, message, "HLT", originClass, originFunc, file, line);
    break;
  case kHLTLogDebug:
    AliLog::Message(AliLog::kDebug, message, "HLT", originClass, originFunc, file, line);
    break;
  case kHLTLogInfo:
    AliLog::Message(AliLog::kInfo, message, "HLT", originClass, originFunc, file, line);
    break;
  case kHLTLogWarning:
    AliLog::Message(AliLog::kWarning, message, "HLT", originClass, originFunc, file, line);
    break;
  case kHLTLogError:
    AliLog::Message(AliLog::kError, message, "HLT", originClass, originFunc, file, line);
    break;
  case kHLTLogFatal:
    AliLog::Message(AliLog::kWarning, message, "HLT", originClass, originFunc, file, line);
    break;
  default:
    break;
  }
  return 0;
}
#endif

const char* AliHLTLogging::BuildLogString(const char *format, va_list &ap, bool bAppend) 
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

  unsigned int iOffset=0;
  if (bAppend) {
    iOffset=strlen(fgAliHLTLoggingTarget.GetArray());
  } else {
    fgAliHLTLoggingTarget[0]=0;
  }
  while (fmt!=NULL) {
    iResult=vsnprintf(fgAliHLTLoggingTarget.GetArray()+iOffset, fgAliHLTLoggingTarget.GetSize()-iOffset, fmt, ap);
    if (iResult==-1)
      // for compatibility with older version of vsnprintf
      iResult=fgAliHLTLoggingTarget.GetSize()*2;
    else
      iResult+=iOffset;

    if (iResult<fgAliHLTLoggingTarget.GetSize())
      // everything in the limit
      break;

    // terminate if buffer is already at the limit
    if (fgAliHLTLoggingTarget.GetSize()>=fgkALIHLTLOGGINGMAXBUFFERSIZE) 
    {
      fgAliHLTLoggingTarget[fgAliHLTLoggingTarget.GetSize()-1]=0;
      break;
    }

    // check limitation and grow the buffer
    if (iResult>fgkALIHLTLOGGINGMAXBUFFERSIZE) iResult=fgkALIHLTLOGGINGMAXBUFFERSIZE;
    fgAliHLTLoggingTarget.Set(iResult+1);

    // copy the original list and skip the first argument if this was the format string
#ifdef R__VA_COPY
    va_end(ap);
    R__VA_COPY(ap, bap);
#else
    fgAliHLTLoggingTarget[fgAliHLTLoggingTarget.GetSize()-1]=0;
    break;
#endif //R__VA_COPY
    if (format==NULL) va_arg(ap, const char*);
  }     
#ifdef R__VA_COPY
  va_end(bap);
#endif //R__VA_COPY

  return fgAliHLTLoggingTarget.GetArray();
}

const char* AliHLTLogging::SetLogString(const void* p, const char* pfmt, const char *format, ...)
{
  // see header file for class documentation
  if (!p || !pfmt) return NULL;
  TString formatstr=format;
  TString pstr;
#ifdef __DEBUG
  pstr.Form(pfmt, p);
#endif
  formatstr.ReplaceAll("_pfmt_", pstr);
  va_list args;
  va_start(args, format);

  const char* message=BuildLogString(formatstr.Data(), args);
  va_end(args);

  return message;
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
    if (fgLoggingFunc) {
      iResult = (*fgLoggingFunc)(NULL/*fParam*/, severity, origin, keyword, AliHLTLogging::BuildLogString(format, args ));
    } else {
      if (fgUseAliLog!=0 && fgAliLoggingFunc!=NULL)
	iResult=(*fgAliLoggingFunc)(severity, NULL, origin, NULL, 0, AliHLTLogging::BuildLogString(format, args ));
      else
        iResult=Message(NULL/*fParam*/, severity, origin, keyword, AliHLTLogging::BuildLogString(format, args ));
    }
    va_end(args);
  }
  return iResult;
}

int AliHLTLogging::LoggingVarargs(AliHLTComponentLogSeverity severity, 
				  const char* originClass, const char* originFunc,
				  const char* file, int line,  ... ) const
{
  // see header file for class documentation
  int iResult=0;

  va_list args;
  va_start(args, line);

  iResult=SendMessage(severity, originClass, originFunc, file, line, AliHLTLogging::BuildLogString(NULL, args ));
  va_end(args);

  return iResult;
}

int AliHLTLogging::SendMessage(AliHLTComponentLogSeverity severity, 
			       const char* originClass, const char* originFunc,
			       const char* file, int line,
			       const char* message) const
{
  // see header file for class documentation
  int iResult=0;
  const char* separator="";
  TString origin;
  if (originClass) {
    origin+=originClass;
    separator="::";
  }
  if (originFunc) {
    origin+=separator;
    origin+=originFunc;
  }

  if (fgLoggingFunc) {
    iResult=(*fgLoggingFunc)(NULL/*fParam*/, severity, origin.Data(), GetKeyword(), message);
  } else {
    if (fgUseAliLog!=0 && fgAliLoggingFunc!=NULL)
      iResult=(*fgAliLoggingFunc)(severity, originClass, originFunc, file, line, message);
    else
      iResult=Message(NULL/*fParam*/, severity, origin.Data(), GetKeyword(), message);
  }
  return iResult;
}

int AliHLTLogging::CheckFilter(AliHLTComponentLogSeverity severity) const
{
  // see header file for class documentation

  int iResult=severity==kHLTLogNone || ((severity&fgGlobalLogFilter)>0 && (severity&fLocalLogFilter)>0);
  return iResult;
}

void AliHLTLogging::SetGlobalLoggingLevel(AliHLTComponentLogSeverity level)
{
  // see header file for class documentation

  fgGlobalLogFilter=level;
}

AliHLTComponentLogSeverity AliHLTLogging::GetGlobalLoggingLevel()
{
  // see header file for class documentation

  return fgGlobalLogFilter;
}

void AliHLTLogging::SetLocalLoggingLevel(AliHLTComponentLogSeverity level)
{
  // see header file for class documentation

  fLocalLogFilter=level;
}


AliHLTComponentLogSeverity AliHLTLogging::GetLocalLoggingLevel()
{
  // see header file for class documentation

  return fLocalLogFilter;
}

void AliHLTLogging::SetLocalLoggingDefault(AliHLTComponentLogSeverity level)
{
  // see header file for class documentation
  fgLocalLogDefault=level;
}

int AliHLTLogging::CheckGroup(const char* /*originClass*/) const
{
  // see header file for class documentation

  return 1;
}

int AliHLTLogging::SetBlackList(const char* classnames)
{
  // see header file for class documentation

  if (classnames)
    fgBlackList=classnames;
  return 0;
}

int AliHLTLogging::SetWhiteList(const char* classnames)
{
  // see header file for class documentation

  if (classnames)
    fgWhiteList=classnames;
  return 0;
}
