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

 private:
  AliLog(const AliLog& log);
  AliLog& operator = (const AliLog& log);

  static void    RootErrorHandler(Int_t level, Bool_t abort, 
				  const char* location, const char* message);

  void           CloseFile(Int_t type);
  FILE*          GetOutputStream(Int_t type);

  UInt_t         GetLogLevel(const char* module, const char* className) const;

  enum {kDebugOffset = kDebug-1};

  static AliLog* fgInstance;                 //! pointer to current instance

  static Bool_t  fgDebugEnabled;             // flag for debug en-/disabling

  UInt_t         fGlobalLogLevel;            // global logging level
  TObjArray      fModuleDebugLevels;         // debug levels for modules
  TObjArray      fClassDebugLevels;          // debug levels for classes

  Int_t          fOutputTypes[kMaxType];     // types of output streams
  TString        fFileNames[kMaxType];       // file names
  FILE*          fOutputFiles[kMaxType];     //! log output files

  Bool_t         fPrintType[kMaxType];       // print type on/off
  Bool_t         fPrintModule[kMaxType];     // print module on/off
  Bool_t         fPrintScope[kMaxType];      // print scope/class name on/off
  Bool_t         fPrintLocation[kMaxType];   // print file and line on/off

  ClassDef(AliLog, 1)   // class for logging debug, info and error messages
};

#ifndef __GNUC__
#ifndef __APPLE__
#define __FUNCTION__ "???"
#endif
#endif

#ifdef LOG_NO_DEBUG
#define AliDebugLevel() -1
#define AliDebugLevelClass() -1
#define AliDebugLevelGeneral(scope) -1
#else
#define AliDebugLevel() ((AliLog::IsDebugEnabled()) ? AliLog::GetDebugLevel(__MODULE__, ClassName()) : -1)
#define AliDebugLevelClass() ((AliLog::IsDebugEnabled()) ? AliLog::GetDebugLevel(__MODULE__, Class()->GetName()) : -1)
#define AliDebugLevelGeneral(scope) ((AliLog::IsDebugEnabled()) ? AliLog::GetDebugLevel(__MODULE__, scope) : -1)
#endif

#ifdef LOG_NO_DEBUG
#define AliDebug(level, message)
#define AliDebugClass(level, message)
#define AliDebugGeneral(scope, level, message)
#else
#define AliDebug(level, message) {if (AliLog::IsDebugEnabled()) AliLog::Debug(level, message, __MODULE__, ClassName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliDebugClass(level, message) {if (AliLog::IsDebugEnabled()) AliLog::Debug(level, message, __MODULE__, Class()->GetName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliDebugGeneral(scope, level, message) {if (AliLog::IsDebugEnabled()) AliLog::Debug(level, message, __MODULE__, scope, __FUNCTION__, __FILE__, __LINE__);}
#endif

#ifdef LOG_NO_INFO
#define AliInfo(message)
#define AliInfoClass(message)
#define AliInfoGeneral(scope, message)
#else
#define AliInfo(message) {AliLog::Message(AliLog::kInfo, message, __MODULE__, ClassName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliInfoClass(message) {AliLog::Message(AliLog::kInfo, message, __MODULE__, Class()->GetName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliInfoGeneral(scope, message) {AliLog::Message(AliLog::kInfo, message, __MODULE__, scope, __FUNCTION__, __FILE__, __LINE__);}
#endif

#ifdef LOG_NO_WARNING
#define AliWarning(message)
#define AliWarningClass(message)
#define AliWarningGeneral(scope, message)
#else
#define AliWarning(message) {AliLog::Message(AliLog::kWarning, message, __MODULE__, ClassName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliWarningClass(message) {AliLog::Message(AliLog::kWarning, message, __MODULE__, Class()->GetName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliWarningGeneral(scope, message) {AliLog::Message(AliLog::kWarning, message, __MODULE__, scope, __FUNCTION__, __FILE__, __LINE__);}
#endif

#define AliError(message) {AliLog::Message(AliLog::kError, message, __MODULE__, ClassName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliErrorClass(message) {AliLog::Message(AliLog::kError, message, __MODULE__, Class()->GetName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliErrorGeneral(scope, message) {AliLog::Message(AliLog::kError, message, __MODULE__, scope, __FUNCTION__, __FILE__, __LINE__);}

#define AliFatal(message) {AliLog::Message(AliLog::kFatal, message, __MODULE__, ClassName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliFatalClass(message) {AliLog::Message(AliLog::kFatal, message, __MODULE__, Class()->GetName(), __FUNCTION__, __FILE__, __LINE__);}
#define AliFatalGeneral(scope, message) {AliLog::Message(AliLog::kFatal, message, __MODULE__, scope, __FUNCTION__, __FILE__, __LINE__);}

#endif
