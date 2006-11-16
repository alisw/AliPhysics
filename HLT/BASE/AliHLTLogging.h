// @(#) $Id: 

#ifndef ALIHLTLOGGING_H
#define ALIHLTLOGGING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTLogging.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  HLT module logging primitives.
*/

#include "AliHLTDataTypes.h"
#include <TObject.h>
#include "AliHLTStdIncludes.h"

//#define LOG_PREFIX ""       // logging prefix, for later extensions


/* the logging macros can be used inside methods of classes which inherit from 
 * AliHLTLogging
 */
// HLTMessage is not filtered
#define HLTMessage( ... )   LoggingVarargs(kHLTLogNone,      NULL , NULL ,  __VA_ARGS__ )

// the following macros are filtered by the Global and Local Log Filter
#define HLTBenchmark( ... ) LoggingVarargs(kHLTLogBenchmark, this->Class_Name() , __func__ ,  __VA_ARGS__ )
#define HLTDebug( ... )     LoggingVarargs(kHLTLogDebug,     this->Class_Name() , __func__ ,  __VA_ARGS__ )
#define HLTInfo( ... )      LoggingVarargs(kHLTLogInfo,      this->Class_Name() , __func__ ,  __VA_ARGS__ )
#define HLTWarning( ... )   LoggingVarargs(kHLTLogWarning,   this->Class_Name() , __func__ ,  __VA_ARGS__ )
#define HLTError( ... )     LoggingVarargs(kHLTLogError,     this->Class_Name() , __func__ ,  __VA_ARGS__ )
#define HLTFatal( ... )     LoggingVarargs(kHLTLogFatal,     this->Class_Name() , __func__ ,  __VA_ARGS__ )

// helper macro to set the keyword
#define HLTLogKeyword(a)    AliHLTKeyword __hltlog_tmpkey__LINE__(this, a)

#define HLT_DEFAULT_LOG_KEYWORD "no key"

class AliHLTLogging {
public:
  AliHLTLogging();
  AliHLTLogging(const AliHLTLogging&);
  AliHLTLogging& operator=(const AliHLTLogging&);
  virtual ~AliHLTLogging();

  // logging filter for all objects
  //
  static AliHLTComponent_LogSeverity SetGlobalLogLevel(AliHLTComponent_LogSeverity iLogFilter) {fGlobalLogFilter=iLogFilter; return fGlobalLogFilter;}

  // logging filter for individual object
  //
  AliHLTComponent_LogSeverity SetLocalLogLevel(AliHLTComponent_LogSeverity iLogFilter) {fLocalLogFilter=iLogFilter; return fLocalLogFilter;}

  // set the default key word
  // the keyword is intended to simplify the use of logging macros
  // 
  void SetDefaultKeyword(const char* keyword) { fpDefaultKeyword=keyword; }

  // set a temporary keyword
  // returns the old key value
  const char* SetKeyword(const char* keyword) 
    { 
      const char* currentKeyword=fpCurrentKeyword;
      fpCurrentKeyword=keyword;
      return currentKeyword; 
    }

  // get the current keyword
  //
  const char* GetKeyword() const
    {
      if (fpCurrentKeyword) return fpCurrentKeyword;
      else if (fpDefaultKeyword) return fpDefaultKeyword;
      return HLT_DEFAULT_LOG_KEYWORD;
    }
  
  static int Init(AliHLTfctLogging pFun);

  // genaral logging function
  //
  int Logging( AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message, ... );

  // logging function with two origin parameters, used by the log macros
  //
  int LoggingVarargs( AliHLTComponent_LogSeverity severity, const char* origin_class, const char* origin_func,  ... ) const;

  // apply filter, return 1 if message should pass
  //
  int CheckFilter(AliHLTComponent_LogSeverity severity) const;

  static int Message(void * param, AliHLTComponent_LogSeverity severity, const char* origin, const char* keyword, const char* message);

  static const char* BuildLogString(const char *format, va_list ap);

  virtual void* GetParameter() {return NULL;}
protected:

private:
  static  AliHLTComponent_LogSeverity fGlobalLogFilter;
  AliHLTComponent_LogSeverity fLocalLogFilter;
  static AliHLTfctLogging fLoggingFunc;
  const char* fpDefaultKeyword;
  const char* fpCurrentKeyword;

  ClassDef(AliHLTLogging, 0)
};

/* the class AliHLTKeyword is a simple helper class used by the HLTLogKeyword macro
 * HLTLogKeyword("a keyword") creates an object of AliHLTKeyword which sets the keyword for the logging class
 * the object is destroyed automatically when the current scope is left and so the keyword is set
 * to the original value
 */
class AliHLTKeyword {
 public:
  AliHLTKeyword()
    :
    fpParent(NULL),
    fpOriginal(NULL)
    {
    }

  AliHLTKeyword(AliHLTLogging* parent, const char* keyword)
    :
    fpParent(parent),
    fpOriginal(NULL)
    {
      if (parent) {
	fpOriginal=fpParent->SetKeyword(keyword);
      }
    }

  AliHLTKeyword(const AliHLTKeyword& kw)
    :
    fpParent(kw.fpParent),
    fpOriginal(kw.fpOriginal)
    {
    }

  AliHLTKeyword& operator=(const AliHLTKeyword& kw)
    { 
      fpParent=kw.fpParent;
      fpOriginal=kw.fpOriginal;
      return *this;
    }

  ~AliHLTKeyword()
    {
      if (fpParent) {
	fpParent->SetKeyword(fpOriginal);
      }
    }

 private:
  AliHLTLogging* fpParent;
  const char* fpOriginal;
};
#endif

