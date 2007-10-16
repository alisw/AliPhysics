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
  hltlog.Logging(kHLTLogInfo, "NotificationHandler", "AliLog", AliHLTLogging::fgLogstr.str().c_str());
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
    AliLog::Message(AliLog::kWarning, message, "HLT", originClass, originFunc, file, line);
    break;
  default:
    break;
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
  AliLog* log=new AliLog;
  log->SetLogNotification(LogNotification);
  log->SetStreamOutput(&AliHLTLogging::fgLogstr);
  return 0;
#endif // NO_ALILOG_NOTIFICATION
  return -ENOSYS;
}
