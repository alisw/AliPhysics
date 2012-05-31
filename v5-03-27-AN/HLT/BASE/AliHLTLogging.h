//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTLOGGING_H
#define ALIHLTLOGGING_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTLogging.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  HLT module logging primitives.
*/

#include "AliHLTDataTypes.h"
#include "AliHLTStdIncludes.h"
#include "TString.h"
#include "TObject.h"
#include "TArrayC.h"

class AliHLTComponentHandler;
//#define LOG_PREFIX ""       // logging prefix, for later extensions

#define ALILOG_WRAPPER_LIBRARY "libHLTrec.so"

/* the logging macros can be used inside methods of classes which inherit from 
 * AliHLTLogging
 */
// function name
#if defined(__GNUC__) || defined(__ICC) || defined(__ECC) || defined(__APPLE__)
#define FUNCTIONNAME() __FUNCTION__
#else
#define FUNCTIONNAME() "???"
#endif

#ifndef ALIHLTLOGGINGVARIADICFREE_H
// HLTMessage is not filtered
#define HLTMessage( ... )   LoggingVarargs(kHLTLogNone,      NULL , NULL , __FILE__ , __LINE__ , __VA_ARGS__ )

// the following macros are filtered by the Global and Local Log Filter
#define HLTLog( level, ... ) if (CheckFilter(level))         LoggingVarargs(level, Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , __VA_ARGS__ )
#define HLTBenchmark( ... ) LoggingVarargs(kHLTLogBenchmark, Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , __VA_ARGS__ )
#ifdef __DEBUG
#define HLTDebug( ... )     if (CheckFilter(kHLTLogDebug) && CheckGroup(Class_Name())) LoggingVarargs(kHLTLogDebug,     Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , __VA_ARGS__ )
#else
#define HLTDebug( ... )
#endif
#define HLTInfo( ... )      if (CheckFilter(kHLTLogInfo))    LoggingVarargs(kHLTLogInfo,      Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , __VA_ARGS__ )
#define HLTWarning( ... )   if (CheckFilter(kHLTLogWarning)) LoggingVarargs(kHLTLogWarning,   Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , __VA_ARGS__ )
#define HLTError( ... )     if (CheckFilter(kHLTLogError))   LoggingVarargs(kHLTLogError,     Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , __VA_ARGS__ )
#define HLTFatal( ... )     if (CheckFilter(kHLTLogFatal))   LoggingVarargs(kHLTLogFatal,     Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , __VA_ARGS__ )
#define HLTImportant( ... ) if (CheckFilter(kHLTLogImportant))LoggingVarargs(kHLTLogImportant,Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , __VA_ARGS__ )

// the same macros are defined variadic free, in that case the message must be complete
// include AliHLTLoggingVariadicFree.h
#else //ALIHLTLOGGINGVARIADICFREE_H
#define HLTMessage( message )   LoggingVarargs(kHLTLogNone,      NULL , NULL , __FILE__ , __LINE__ , message )
#define HLTLog( level, message) if (CheckFilter(level))          LoggingVarargs(level, Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , message )
#define HLTBenchmark( message ) LoggingVarargs(kHLTLogBenchmark, Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , message )
#ifdef __DEBUG
#define HLTDebug( message )     if (CheckFilter(kHLTLogDebug) && CheckGroup(Class_Name())) LoggingVarargs(kHLTLogDebug,     Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , message )
#else
#define HLTDebug( message )
#endif
#define HLTInfo( message )      if (CheckFilter(kHLTLogInfo))    LoggingVarargs(kHLTLogInfo,      Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , message )
#define HLTWarning( message )   if (CheckFilter(kHLTLogWarning)) LoggingVarargs(kHLTLogWarning,   Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , message )
#define HLTError( message )     if (CheckFilter(kHLTLogError))   LoggingVarargs(kHLTLogError,     Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , message )
#define HLTFatal( message )     if (CheckFilter(kHLTLogFatal))   LoggingVarargs(kHLTLogFatal,     Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , message )
#define HLTImportant( message ) if (CheckFilter(kHLTLogImportant))LoggingVarargs(kHLTLogImportant,Class_Name() , FUNCTIONNAME() , __FILE__ , __LINE__ , message )
#endif //ALIHLTLOGGINGVARIADICFREE_H

// helper macro to set the keyword
#define HLTLogKeyword(a)    AliHLTKeyword hltlogTmpkey(this, a)

#define HLT_DEFAULT_LOG_KEYWORD "no key"

/**
 * @class AliHLTLogging
 * Basic logging class. All classes inherit the besic HLT logging functionality.
 * Logging levels are controlled by a global logging filter and a local logging
 * filter.
 * 
 * @section alihlt_logging_levels Logging Levels
 * Logging levels are switched by a bit pattern,  AliHLTComponentLogSeverity {
 * - ::kHLTLogNone no logging (0)
 * - ::kHLTLogBenchmark benchmark messages (0x1)
 * - ::kHLTLogDebug debug messages (0x2)
 * - ::kHLTLogInfo info messages (0x4)
 * - ::kHLTLogWarning warning messages (0x8)
 * - ::kHLTLogError error messages (0x10)
 * - ::kHLTLogFatal fatal error messages (0x20)
 * - ::kHLTLogImportant few important messages not to be filtered out (0x40)
 * - ::kHLTLogAll special value to enable all messages (0x7f)
 * - ::kHLTLogDefault the default logging level: Warning, Error, Fatal, Important (0x79)
 *
 * @section alihlt_logging_filter Logging Filters
 * The class provides a global and a local logging filter, the AND beween both
 * defines whether a message is printed or not.
 *
 * The global filter is by default set to ::kHLTLogAll. Please note that AliHLTSystem
 * changes the global logging level to ::kHLTLogDefault. The global filter can be
 * adjusted by means of SetGlobalLoggingLevel().
 *
 * The local filter is set to ::kHLTLogAll and can be adjusted by
 * SetLocalLoggingLevel(). The default can be changed for all objects by
 * SetLocalLoggingDefault(). Please note that a change of the default level only
 * applies to objects generated after the change of the default.
 *
 * @section alihlt_logging_external Redirection
 * - external logging function
 * - keyword
 *
 * @section alihlt_logging_aliroot AliRoot Redirection
 * - switching of redirection
 * - logging options in AliSimulation/AliReconstruction
 *
 * @ingroup alihlt_component
 */
class AliHLTLogging {
public:
  AliHLTLogging();
  AliHLTLogging(const AliHLTLogging&);
  AliHLTLogging& operator=(const AliHLTLogging&);
  virtual ~AliHLTLogging();

  /** set the default key word
   * the keyword is intended to simplify the use of logging macros
   */ 
  void SetDefaultKeyword(const char* keyword) { fpDefaultKeyword=keyword; }

  /**
   * Set a temporary keyword
   * Keywords need to be static const strings, the class handles only
   * pointers and assumes the strings to be persistent.
   * returns the old key value
   */
  const char* SetKeyword(const char* keyword) 
    { 
      const char* currentKeyword=fpCurrentKeyword;
      fpCurrentKeyword=keyword;
      return currentKeyword; 
    }

  /**
   * Get the current keyword
   */
  const char* GetKeyword() const
    {
      if (fpCurrentKeyword) return fpCurrentKeyword;
      else if (fpDefaultKeyword) return fpDefaultKeyword;
      return HLT_DEFAULT_LOG_KEYWORD;
    }
  
  /**
   * Init the AliLogging class for use from external package.
   * This initializes the logging callback. <br>
   * Only deployed by external users of the C wrapper interface, not used
   * when running in AliRoot
   */
  static int Init(AliHLTfctLogging pFun);

  /**
   * Init the message trap in AliLog.
   * This initializes the AliLog trap, the AliLog class is the logging
   * mechanism of AliRoot. The trap can fetch log messages written to
   * AliLog, components and detector algorithms can use the AliLog
   * mechanism to be as close as possible to Offline habits. <br>
   * Only used with external users of the C wrapper interface, not used
   * when running in AliRoot
   */
  static int InitAliLogTrap(AliHLTComponentHandler* pHandler);

  /**
   * Init the AliRoot logging function.
   * All log messages are redirected to AliLog when running in AliRoot.
   * Note: when running in PubSub, AliLog messages are redirected to PubSub,
   * see AliHLTLogging::InitAliLogTrap
   */
  static int InitAliLogFunc(AliHLTComponentHandler* pHandler);

  /**
   * Genaral logging function
   */
  int Logging( AliHLTComponentLogSeverity severity, const char* origin, const char* keyword, const char* message, ... );

  /**
   * Logging function with two origin parameters, used by the log macros
   */
  virtual int LoggingVarargs(AliHLTComponentLogSeverity severity, 
			     const char* originClass, const char* originFunc,
			     const char* file, int line, ... ) const;

  /**
   * Send formatted string to the different channels
   */
  int SendMessage(AliHLTComponentLogSeverity severity, 
		  const char* originClass, const char* originFunc,
		  const char* file, int line,
		  const char* message) const;

  /**
   * Evaluate the group of the debug message from the class name.
   * @return 1 if message should be printed
   */
  int CheckGroup(const char* originClass) const;

  /**
   * Set the black list of classes.
   * If the white list is set, debug messages are skipped for
   * all classes matching one of the regular expressions in the string.
   */
  static int SetBlackList(const char* classnames);

  /**
   * Set the white list of classes.
   * If the white list is set, debug messages are only printed for
   * classes matching one of the regular expressions in the string.
   */
  static int SetWhiteList(const char* classnames);

  /**
   * Apply filter
   * @return 1 if message should pass
   */
  int CheckFilter(AliHLTComponentLogSeverity severity) const;

  /**
   * Set global logging level
   * logging filter for all objects
   */
  static void SetGlobalLoggingLevel(AliHLTComponentLogSeverity level);

  /**
   * Get global logging level
   * logging filter for all objects
   */
  static AliHLTComponentLogSeverity GetGlobalLoggingLevel();

  /**
   * Set local logging level
   * logging filter for individual object
   */
  virtual void SetLocalLoggingLevel(AliHLTComponentLogSeverity level);

  /**
   * Set local logging default
   * Default logging filter for individual objects.
   */
  static void SetLocalLoggingDefault(AliHLTComponentLogSeverity level);

  /**
   * Get default setting for local logging filter for individual objects.
   */
  static AliHLTComponentLogSeverity GetLocalLoggingDefault() { return fgLocalLogDefault; }

  /**
   * Get local logging level
   * logging filter for individual object
   */
  AliHLTComponentLogSeverity GetLocalLoggingLevel();

  /**
   * Print message to stdout
   */
  static int Message(void * param, AliHLTComponentLogSeverity severity, const char* origin, const char* keyword, const char* message);

  /**
   * Build the log string from format specifier and variadac arguments
   * @param format     format string of printf style
   * @param ap         opened and initialized argument list
   * @param bAppend    append to current content
   * @return const char string with the formatted message 
   */
  static const char* BuildLogString(const char *format, va_list &ap, bool bAppend=false);

  /**
   * Set the log string from format specifier and from variable arguments.
   * @param format     format string of printf style
   * @return const char string with the formatted message 
   */
  static const char* SetLogString(const void* p, const char* pfmt, const char *format, ... );

  /**
   * Get parameter given by the external caller.
   * This functionality is not yet implemented. It is intended
   * to pass the parameter pointer given to the component at
   * initialization back to the caller.
   */
  virtual void* GetParameter() const {return NULL;}

  /**
   * Switch logging through AliLog on or off
   * @param sw          1 = logging through AliLog
   */
  void SwitchAliLog(int sw) {fgUseAliLog=(sw!=0);}

  /** target stream for AliRoot logging methods */
  static ostringstream fgLogstr;                                   //! transient

  /** 
   * The message function for dynamic use.
   * In order to avoid dependencies on AliRoot libraries, libHLTbase loads
   * the library dynamically and looks for the symbol.
   */
  typedef int (*AliHLTDynamicMessage)(AliHLTComponentLogSeverity severity, 
				      const char* originClass, 
				      const char* originFunc,
				      const char* file, int line, 
				      const char* message); 

  /**
   * The init function of the message callback for dynamic use.
   * In order to avoid dependencies on AliRoot libraries, libHLTbase loads
   * the library dynamically and looks for the symbol.
   */
  typedef int (*InitAliDynamicMessageCallback)();
  
protected:
  /** the AliRoot logging function */
  static AliHLTDynamicMessage fgAliLoggingFunc;                    //! transient

private:
  /** the global logging filter */
  static  AliHLTComponentLogSeverity fgGlobalLogFilter;            // see above
  /** the local logging filter for one class */
  AliHLTComponentLogSeverity fLocalLogFilter;                      // see above
  /** the global logging filter */
  static  AliHLTComponentLogSeverity fgLocalLogDefault;            // see above
  /** logging callback from the framework */
  static AliHLTfctLogging fgLoggingFunc;                           // see above
  /** default keyword */
  const char* fpDefaultKeyword;                                    //! transient
  /** current keyword */
  const char* fpCurrentKeyword;                                    //! transient
  /** switch for logging through AliLog, default on */
  static int fgUseAliLog;                                          // see above
  /**
   * The global logging buffer.
   * The buffer is created with an initial size and grown dynamically on
   * demand.
   */
  static TArrayC fgAliHLTLoggingTarget;                            //! transient
  
  /** the maximum size of the buffer */
  static const int fgkALIHLTLOGGINGMAXBUFFERSIZE;                  //! transient

  /** groups of classes not to print debug messages */
  static TString fgBlackList;                                      //! transient
  
  /** groups of classes not to print debug messages */
  static TString fgWhiteList;                                      //! transient
  
  ClassDef(AliHLTLogging, 0)
};

/* the class AliHLTKeyword is a simple helper class used by the HLTLogKeyword macro
 * HLTLogKeyword("a keyword") creates an object of AliHLTKeyword which sets the keyword for the logging class
 * the object is destroyed automatically when the current scope is left and so the keyword is set
 * to the original value. Please note that all keywords need to be static strings, only pointyers
 * are handled and the strings required to ber persistent.
 */
class AliHLTKeyword {
 public:
  AliHLTKeyword()
    :
    fpParent(NULL),
    fpOriginal(NULL)
    {
    }

  AliHLTKeyword(const AliHLTLogging* parent, const char* keyword)
    :
    fpParent(const_cast<AliHLTLogging*>(parent)),
    fpOriginal(NULL)
    {
      // the const cast is on purpose in order to be allowed to use
      // HLTLogKeyword from const member functions
      if (fpParent) {
	fpOriginal=fpParent->SetKeyword(keyword);
      }
    }

  ~AliHLTKeyword()
    {
      if (fpParent) {
	fpParent->SetKeyword(fpOriginal);
      }
    }

 private:
  /// copy constructor prohibited
  AliHLTKeyword(const AliHLTKeyword& kw);
  /// assignment operator prohibited
  AliHLTKeyword& operator=(const AliHLTKeyword& kw);

  AliHLTLogging* fpParent;                                         //! transient
  const char* fpOriginal;                                          //! transient
};
#endif

