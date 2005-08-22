// @(#) $Id: 

#ifndef ALIHLTLOGGING_H
#define ALIHLTLOGGING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTLogging
   the HLT logging methods
 */


#include "AliHLTDataTypes.h"
#include <TObject.h>
#include <stdio.h>

#define LOG_BUFFER_SIZE 100 // global logging buffer
#define LOG_PREFIX ""       // logging prefix, for later extensions

#define DebugMsg( ... ) LoggingVarargs(kHLTLogDebug, this->Class_Name() , __func__ ,  __VA_ARGS__ )

class AliHLTLogging {
public:
  AliHLTLogging();
  virtual ~AliHLTLogging();

  /* logging filter for all objects
   */
  static AliHLTComponent_LogSeverity SetGlobalLogLevel(AliHLTComponent_LogSeverity iLogFilter) {fGlobalLogFilter=iLogFilter; return fGlobalLogFilter;}

  /* logging filter for individual object
   */
  AliHLTComponent_LogSeverity SetLocalLogLevel(AliHLTComponent_LogSeverity iLogFilter) {fLocalLogFilter=iLogFilter; return fLocalLogFilter;}
  
  static int Init(AliHLTfctLogging pFun) { fLoggingFunc=pFun; return 0;}

  // genaral logging function
  //
  static int Logging( AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message, ... );

  // logging function with two origin parameters, used by the log macros
  //
  int LoggingVarargs( AliHLTComponent_LogSeverity severity, const char* origin_class, const char* origin_func,  ... );

  // apply filter, return 1 if message should pass
  //
  int CheckFilter(AliHLTComponent_LogSeverity severity);

  static int Message(void * param, AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message);

  static const char* BuildLogString(const char *format, va_list ap);

  virtual void* GetParameter() {return NULL;}
protected:

private:
  static  AliHLTComponent_LogSeverity fGlobalLogFilter;
  AliHLTComponent_LogSeverity fLocalLogFilter;
  static char fLogBuffer[LOG_BUFFER_SIZE];
  static char fOriginBuffer[LOG_BUFFER_SIZE];
  static AliHLTfctLogging fLoggingFunc;

  ClassDef(AliHLTLogging, 0)
};
#endif

