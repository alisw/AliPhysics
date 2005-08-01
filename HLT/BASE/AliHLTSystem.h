// @(#) $Id$

#ifndef ALIHLTSYSTEM_H
#define ALIHLTSYSTEM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTComponentHandler
   global HLT module management 
 */


#include "AliL3RootTypes.h"
#include "AliHLTDataTypes.h"
#include <stdarg.h>
class AliHLTComponentHandler;

#define LOG_BUFFER_SIZE 100 // global logging buffer
#define LOG_PREFIX " "       // logging prefix, for later extensions

class AliHLTSystem {
public:
  AliHLTSystem();
  virtual ~AliHLTSystem();

  Int_t SetLogLevel(Int_t iLogLevel) {return fLogLevel;}
  static int Logging(void * param, AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message);
  static const char* BuildLogString(const char *format, va_list ap);
  AliHLTComponentHandler* fpComponentHandler;
protected:

private:
  Int_t fLogLevel;
  static char fLogBuffer[LOG_BUFFER_SIZE];

  ClassDef(AliHLTSystem, 0)
};
#endif

