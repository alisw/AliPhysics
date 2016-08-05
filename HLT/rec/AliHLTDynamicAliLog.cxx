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

/** @file   AliHLTDynamicAliLog.cxx
    @author Matthias Richter
    @date   
    @brief  Implementation of dynamically loaded AliLog functionality
*/

#include <sstream>
#include <iostream>
#include "TString.h"
#include "AliLog.h"
#include "AliHLTLogging.h"
#include "AliHLTDataTypes.h"

/**
 * Notification callback for AliRoot logging methods
 */
void LogNotification(AliLog::EType_t /*level*/, const char* /*message*/)
{
  // Notification callback for AliRoot logging methods

  AliHLTLogging hltlog;
  // in case of the initialized callback we never want to redirect
  // HLT logging messages to AliLog (that would be a circular function call)
  hltlog.SwitchAliLog(0);
  AliHLTComponentLogSeverity level=kHLTLogNone;
  int offset=2;
  TString logstring(AliHLTLogging::fgLogstr.str().c_str());
  if (logstring.Length()<2) return;
  switch (logstring[0]) {
  case 'D':
    level=kHLTLogDebug;
    break;
  case 'I':
    level=kHLTLogInfo;
    break;
  case 'W':
    level=kHLTLogWarning;
    break;
  case 'E':
    level=kHLTLogError;
    break;
  case 'F':
    level=kHLTLogFatal;
    break;
  default:
    level=kHLTLogInfo;
    offset=0;
  }
  
  TString origin=&logstring[offset];
  TString message=origin;
  int blank=origin.First(' ');
  if (blank>0 && origin[blank-1]==':') {
    origin.Remove(blank-1, origin.Length());
    message.Remove(0, blank+1);
  } else {
    origin="";
  }
  message=message.Strip(TString::kTrailing, '\n');

  hltlog.Logging(level, origin.Data(), "AliLog", message.Data());
  AliHLTLogging::fgLogstr.clear();
  string empty("");
  AliHLTLogging::fgLogstr.str(empty);
}

/**
 * This is the entry point for AliLog messages.
 * The function pointer is fetched by the AliLogging class after libAliHLTUtil
 * was loaded dynamically. By that we can keep libHLTbase free of AliRoot
 * libraries.
 */
extern "C" int AliDynamicMessage(AliHLTComponentLogSeverity severity, 
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
    AliLog::Message(AliLog::kFatal, message, "HLT", originClass, originFunc, file, line);
    break;
  case kHLTLogImportant:
    AliLog::Message(AliLog::kInfo, message, "HLT", originClass, originFunc, file, line);
    break;
  default:
    AliLog::Message(AliLog::kInfo, message, "HLT", originClass?originClass:"AliHLT", originFunc, file, line);
  }
  return 0;
}

/**
 * Init the AliLog callback.
 * If libHLTbase is used within AliRoot, no message callback is initialized since
 * all logging happens through AliRoot. If externally used by other frameworks (e.g.
 * PubSub), all messages from components to AliLog must be trapped and redirected
 * to the external callback.
 */
extern "C" int InitAliDynamicMessageCallback() 
{
  // older versions of AliLog does not support the notification callback and
  // stringstreams, but they support the logging macros in general
#ifndef NO_ALILOG_NOTIFICATION
#ifndef NO_ALILOG_GETROOTLOGGER
  AliLog* log = AliLog::GetRootLogger();
#else
  AliLog* log = new AliLog;
#endif //NO_ALILOG_GETROOTLOGGER
  log->SetLogNotification(LogNotification);
  log->SetStreamOutput(&AliHLTLogging::fgLogstr);
  log->SetPrintScope(true);
  return 0;
#endif // NO_ALILOG_NOTIFICATION
  return -ENOSYS;
}
