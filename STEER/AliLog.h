#ifndef ALILOG_H
#define ALILOG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///
/// class for logging debug, info and error messages
///

#include <TObject.h>
#include <TObjArray.h>
#include <TString.h>


class AliLog: public TObject {
 public:
  AliLog();
  virtual ~AliLog();
  static AliLog* Instance() {return fgInstance;}

  enum EType {kFatal = 0, kError, kWarning, kInfo, kDebug, kMaxType};

  static void  EnableDebug(Bool_t enabled);
  static void  SetGlobalLogLevel(EType type);
  static Int_t GetGlobalLogLevel();
  static void  SetGlobalDebugLevel(Int_t level);
  static Int_t GetGlobalDebugLevel();
  static void  SetModuleDebugLevel(const char* module, Int_t level);
  static void  ClearModuleDebugLevel(const char* module);
  static void  SetClassDebugLevel(const char* className, Int_t level);
  static void  ClearClassDebugLevel(const char* className);

  static void  SetStandardOutput();
  static void  SetStandardOutput(EType type);
  static void  SetErrorOutput();
  static void  SetErrorOutput(EType type);
  static void  SetFileOutput(const char* fileName);
  static void  SetFileOutput(EType type, const char* fileName);
  static void  Flush();

  static void  SetHandleRootMessages(Bool_t on);

  static void  SetPrintType(Bool_t on);
  static void  SetPrintType(EType type, Bool_t on);
  static void  SetPrintModule(Bool_t on);
  static void  SetPrintModule(EType type, Bool_t on);
  static void  SetPrintScope(Bool_t on);
  static void  SetPrintScope(EType type, Bool_t on);
  static void  SetPrintLocation(Bool_t on);
  static void  SetPrintLocation(EType type, Bool_t on);

  static void  SetPrintRepetitions(Bool_t on);

  static void  Write(const char* name, Int_t option = 0);

  // the following public methods are used by the preprocessor macros 
  // and should not be called directly
  static Bool_t IsDebugEnabled() {return fgDebugEnabled;}
  static Int_t GetDebugLevel(const char* module, const char* className);
  static void  Message(UInt_t level, const char* message, 
                       const char* module, const char* className,
                       const char* function, const char* file, Int_t line);
  static void  Debug(UInt_t level, const char* message, 
                     const char* module, const char* className,
                     const char* function, const char* file, Int_t line);

  static Int_t RedirectStdoutTo(EType type, UInt_t level, const char* module, 
                                const char* className, const char* function,
                                const char* file, Int_t line, Bool_t print);
  static Int_t RedirectStderrTo(EType type, UInt_t level, const char* module, 
                                const char* className, const char* function,
                                const char* file, Int_t line, Bool_t print);
  static void  RestoreStdout(Int_t original);
  static void  RestoreStderr(Int_t original);

  static ostream& Stream(EType type, UInt_t level,
                         const char* module, const char* className,
                         const char* function, const char* file, Int_t line);

 private:
  AliLog(const AliLog& log);
  AliLog& operator = (const AliLog& log);

  void           ReadEnvSettings();

  static void    RootErrorHandler(Int_t level, Bool_t abort, 
				  const char* location, const char* message);

  void           CloseFile(Int_t type);
  FILE*          GetOutputStream(Int_t type);

  UInt_t         GetLogLevel(const char* module, const char* className) const;
  void           PrintMessage(UInt_t type, const char* message, 
                              const char* module, const char* className,
                              const char* function, 
                              const char* file, Int_t line);
  void           PrintRepetitions();

  Int_t          RedirectTo(FILE* stream, EType type, UInt_t level,
                            const char* module, const char* className,
                            const char* function,
                            const char* file, Int_t line, Bool_t print);

  ostream&       GetStream(EType type, UInt_t level,
                           const char* module, const char* className,
                           const char* function, const char* file, Int_t line);

  enum {kDebugOffset = kDebug-1};

  static AliLog* fgInstance;                 //! pointer to current instance

  static Bool_t  fgDebugEnabled;             // flag for debug en-/disabling

  UInt_t         fGlobalLogLevel;            // global logging level
  TObjArray      fModuleDebugLevels;         // debug levels for modules
  TObjArray      fClassDebugLevels;          // debug levels for classes

  Int_t          fOutputTypes[kMaxType];     // types of output streams
  TString        fFileNames[kMaxType];       // file names
  FILE*          fOutputFiles[kMaxType];     //! log output files
  ofstream*      fOutputStreams[kMaxType];   //! log output streams

  Bool_t         fPrintType[kMaxType];       // print type on/off
  Bool_t         fPrintModule[kMaxType];     // print module on/off
  Bool_t         fPrintScope[kMaxType];      // print scope/class name on/off
  Bool_t         fPrintLocation[kMaxType];   // print file and line on/off

  Bool_t         fPrintRepetitions;          // print number of repetitions instead of repeated message on/off

  Int_t          fRepetitions;               //! counter of repetitions
  UInt_t         fLastType;                  //! type of last message
  TString        fLastMessage;               //! last message
  TString        fLastModule;                //! module name of last message
  TString        fLastClassName;             //! class name of last message
  TString        fLastFunction;              //! function name of last message
  TString        fLastFile;                  //! file name of last message
  Int_t          fLastLine;                  //! line number of last message

  ClassDef(AliLog, 1)   // class for logging debug, info and error messages
};


// module name
#ifdef __MODULE__
#define MODULENAME() __MODULE__
#else
#define MODULENAME() "NoModule"
#endif

// function name
#if defined(__GNUC__) || defined(__ICC) || defined(__ECC) || defined(__APPLE__)
#define FUNCTIONNAME() __FUNCTION__
// #elif defined(__HP_aCC) || defined(__alpha) || defined(__DECCXX)
// #define FUNCTIONNAME() __FUNC__
#else
#define FUNCTIONNAME() "???"
#endif

// redirection
#define REDIRECTSTDOUT(type, level, scope, whatever) {Int_t originalStdout = AliLog::RedirectStdoutTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; AliLog::RestoreStdout(originalStdout);}
#define REDIRECTSTDERR(type, level, scope, whatever) {Int_t originalStderr = AliLog::RedirectStderrTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; AliLog::RestoreStderr(originalStderr);}
#define REDIRECTSTDOUTANDSTDERR(type, level, scope, whatever) {Int_t originalStdout = AliLog::RedirectStdoutTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); Int_t originalStderr = AliLog::RedirectStderrTo(type, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; AliLog::RestoreStderr(originalStderr); AliLog::RestoreStdout(originalStdout);}


// debug level
#ifdef LOG_NO_DEBUG
#define AliDebugLevel() -1
#define AliDebugLevelClass() -1
#define AliDebugLevelGeneral(scope) -1
#else
#define AliDebugLevel() ((AliLog::IsDebugEnabled()) ? AliLog::GetDebugLevel(MODULENAME(), ClassName()) : -1)
#define AliDebugLevelClass() ((AliLog::IsDebugEnabled()) ? AliLog::GetDebugLevel(MODULENAME(), Class()->GetName()) : -1)
#define AliDebugLevelGeneral(scope) ((AliLog::IsDebugEnabled()) ? AliLog::GetDebugLevel(MODULENAME(), scope) : -1)
#endif

// debug messages
#ifdef LOG_NO_DEBUG
#define AliDebug(level, message)
#define AliDebugClass(level, message)
#define AliDebugGeneral(scope, level, message)
#else
#define AliDebug(level, message) {if (AliLog::IsDebugEnabled()) AliLog::Debug(level, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliDebugClass(level, message) {if (AliLog::IsDebugEnabled()) AliLog::Debug(level, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliDebugGeneral(scope, level, message) {if (AliLog::IsDebugEnabled()) AliLog::Debug(level, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}
#endif

// redirection to debug
#define StdoutToAliDebug(level, whatever) REDIRECTSTDOUT(AliLog::kDebug, level, ClassName(), whatever)
#define StderrToAliDebug(level, whatever) REDIRECTSTDERR(AliLog::kDebug, level, ClassName(), whatever)
#define ToAliDebug(level, whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kDebug, level, ClassName(), whatever)
#define StdoutToAliDebugClass(level, whatever) REDIRECTSTDOUT(AliLog::kDebug, level, Class()->GetName(), whatever)
#define StderrToAliDebugClass(level, whatever) REDIRECTSTDERR(AliLog::kDebug, level, Class()->GetName(), whatever)
#define ToAliDebugClass(level, whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kDebug, level, Class()->GetName(), whatever)
#define StdoutToAliDebugGeneral(scope, level, whatever) REDIRECTSTDOUT(AliLog::kDebug, level, scope, whatever)
#define StderrToAliDebugGeneral(scope, level, whatever) REDIRECTSTDERR(AliLog::kDebug, level, scope, whatever)
#define ToAliDebugGeneral(scope, level, whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kDebug, level, scope, whatever)

// debug stream objects
#define AliDebugStream(level) AliLog::Stream(AliLog::kDebug, level, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define AliDebugClassStream(level) AliLog::Stream(AliLog::kDebug, level, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define AliDebugGeneralStream(scope, level) AliLog::Stream(AliLog::kDebug, level, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// info messages
#ifdef LOG_NO_INFO
#define AliInfo(message)
#define AliInfoClass(message)
#define AliInfoGeneral(scope, message)
#else
#define AliInfo(message) {AliLog::Message(AliLog::kInfo, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliInfoClass(message) {AliLog::Message(AliLog::kInfo, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliInfoGeneral(scope, message) {AliLog::Message(AliLog::kInfo, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}
#endif

// redirection to info
#define StdoutToAliInfo(whatever) REDIRECTSTDOUT(AliLog::kInfo, 0, ClassName(), whatever)
#define StderrToAliInfo(whatever) REDIRECTSTDERR(AliLog::kInfo, 0, ClassName(), whatever)
#define ToAliInfo(whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kInfo, 0, ClassName(), whatever)
#define StdoutToAliInfoClass(whatever) REDIRECTSTDOUT(AliLog::kInfo, 0, Class()->GetName(), whatever)
#define StderrToAliInfoClass(whatever) REDIRECTSTDERR(AliLog::kInfo, 0, Class()->GetName(), whatever)
#define ToAliInfoClass(whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kInfo, 0, Class()->GetName(), whatever)
#define StdoutToAliInfoGeneral(scope, whatever) REDIRECTSTDOUT(AliLog::kInfo, 0, scope, whatever)
#define StderrToAliInfoGeneral(scope, whatever) REDIRECTSTDERR(AliLog::kInfo, 0, scope, whatever)
#define ToAliInfoGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kInfo, 0, scope, whatever)

// info stream objects
#define AliInfoStream() AliLog::Stream(AliLog::kInfo, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define AliInfoClassStream() AliLog::Stream(AliLog::kInfo, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define AliInfoGeneralStream(scope) AliLog::Stream(AliLog::kInfo, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// warning messages
#ifdef LOG_NO_WARNING
#define AliWarning(message)
#define AliWarningClass(message)
#define AliWarningGeneral(scope, message)
#else
#define AliWarning(message) {AliLog::Message(AliLog::kWarning, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliWarningClass(message) {AliLog::Message(AliLog::kWarning, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliWarningGeneral(scope, message) {AliLog::Message(AliLog::kWarning, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}
#endif

// redirection to warning
#define StdoutToAliWarning(whatever) REDIRECTSTDOUT(AliLog::kWarning, 0, ClassName(), whatever)
#define StderrToAliWarning(whatever) REDIRECTSTDERR(AliLog::kWarning, 0, ClassName(), whatever)
#define ToAliWarning(whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kWarning, 0, ClassName(), whatever)
#define StdoutToAliWarningClass(whatever) REDIRECTSTDOUT(AliLog::kWarning, 0, Class()->GetName(), whatever)
#define StderrToAliWarningClass(whatever) REDIRECTSTDERR(AliLog::kWarning, 0, Class()->GetName(), whatever)
#define ToAliWarningClass(whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kWarning, 0, Class()->GetName(), whatever)
#define StdoutToAliWarningGeneral(scope, whatever) REDIRECTSTDOUT(AliLog::kWarning, 0, scope, whatever)
#define StderrToAliWarningGeneral(scope, whatever) REDIRECTSTDERR(AliLog::kWarning, 0, scope, whatever)
#define ToAliWarningGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kWarning, 0, scope, whatever)

// warning stream objects
#define AliWarningStream() AliLog::Stream(AliLog::kWarning, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define AliWarningClassStream() AliLog::Stream(AliLog::kWarning, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define AliWarningGeneralStream(scope) AliLog::Stream(AliLog::kWarning, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// error messages
#define AliError(message) {AliLog::Message(AliLog::kError, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliErrorClass(message) {AliLog::Message(AliLog::kError, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliErrorGeneral(scope, message) {AliLog::Message(AliLog::kError, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}

// redirection to error
#define StdoutToAliError(whatever) REDIRECTSTDOUT(AliLog::kError, 0, ClassName(), whatever)
#define StderrToAliError(whatever) REDIRECTSTDERR(AliLog::kError, 0, ClassName(), whatever)
#define ToAliError(whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kError, 0, ClassName(), whatever)
#define StdoutToAliErrorClass(whatever) REDIRECTSTDOUT(AliLog::kError, 0, Class()->GetName(), whatever)
#define StderrToAliErrorClass(whatever) REDIRECTSTDERR(AliLog::kError, 0, Class()->GetName(), whatever)
#define ToAliErrorClass(whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kError, 0, Class()->GetName(), whatever)
#define StdoutToAliErrorGeneral(scope, whatever) REDIRECTSTDOUT(AliLog::kError, 0, scope, whatever)
#define StderrToAliErrorGeneral(scope, whatever) REDIRECTSTDERR(AliLog::kError, 0, scope, whatever)
#define ToAliErrorGeneral(scope, whatever) REDIRECTSTDOUTANDSTDERR(AliLog::kError, 0, scope, whatever)

// error stream objects
#define AliErrorStream() AliLog::Stream(AliLog::kError, 0, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define AliErrorClassStream() AliLog::Stream(AliLog::kError, 0, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__)
#define AliErrorGeneralStream(scope) AliLog::Stream(AliLog::kError, 0, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__)


// fatal messages
#define AliFatal(message) {AliLog::Message(AliLog::kFatal, message, MODULENAME(), ClassName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliFatalClass(message) {AliLog::Message(AliLog::kFatal, message, MODULENAME(), Class()->GetName(), FUNCTIONNAME(), __FILE__, __LINE__);}
#define AliFatalGeneral(scope, message) {AliLog::Message(AliLog::kFatal, message, MODULENAME(), scope, FUNCTIONNAME(), __FILE__, __LINE__);}

#endif
