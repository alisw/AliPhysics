/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for logging debug, info and error messages                          //
//                                                                           //
// The AliLog class is a singleton class. It allows to steer the output      //
// level and output streams for different types of messages via static       //
// methods.                                                                  //
//                                                                           //
// It also handles the messages produces by the preprocessor macros defined  //
// in the header file: AliDebug, AliInfo, AliWarning, AliError, AliFatal.    //
//                                                                           //
// More details about the message logging can be found on the ALICE Offline  //
// web page.                                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TError.h>
#include <TNamed.h>
#include <TSystem.h>

#include "AliLog.h"

ClassImp(AliLog)


AliLog* AliLog::fgInstance = NULL;

Bool_t AliLog::fgDebugEnabled = kTRUE;


//_____________________________________________________________________________
AliLog::AliLog() :
  TObject(),
  fGlobalLogLevel(kInfo),
  fModuleDebugLevels(),
  fClassDebugLevels()
{
// default constructor: set default values

  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    fOutputTypes[iType] = 0;
    fFileNames[iType] = "";
    fOutputFiles[iType] = NULL;

    fPrintType[iType] = kTRUE;
    fPrintModule[iType] = kFALSE;
    fPrintScope[iType] = kTRUE;
    fPrintLocation[iType] = (iType == kDebug);  
  }

  SetHandleRootMessages(kTRUE);

  // replace the previous instance by this one
  if (fgInstance) delete fgInstance;
  fgInstance = this;
}

//_____________________________________________________________________________
AliLog::~AliLog()
{
// destructor: clean up and reset instance pointer

  for (Int_t i = 0; i < fModuleDebugLevels.GetEntriesFast(); i++) {
    if (fModuleDebugLevels[i]) fModuleDebugLevels[i]->Delete();
  }
  fClassDebugLevels.Delete();
  for (Int_t i = 0; i < fClassDebugLevels.GetEntriesFast(); i++) {
    if (fClassDebugLevels[i]) fClassDebugLevels[i]->Delete();
  }
  fClassDebugLevels.Delete();

  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    CloseFile(iType);
  }
  fflush(stderr);
  fflush(stdout);

  fgInstance = NULL;
}

//_____________________________________________________________________________
AliLog::AliLog(const AliLog& log) :
  TObject(log)
{
// copy constructor

  Fatal("AliLog", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliLog& AliLog::operator = (const AliLog& /*log*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}


//_____________________________________________________________________________
void AliLog::RootErrorHandler(Int_t level, Bool_t abort, 
			      const char* location, const char* message)
{
// new error handler for messages from root

  switch (level) {
  case ::kFatal    : level = kFatal; break;
  case ::kSysError :
    DefaultErrorHandler(level, abort, location, message);
    return;
  case ::kBreak    :
    DefaultErrorHandler(level, abort, location, message);
    return;
  case ::kError    : level = kError; break;
  case ::kWarning  : level = kWarning; break;
  case ::kInfo     : level = kInfo; break;
  default          : level = kDebug; break;
  }
  AliLog::Message(level, message, "ROOT", NULL, location, NULL, 0);
}


//_____________________________________________________________________________
void AliLog::EnableDebug(Bool_t enabled)
{
// enable or disable debug output

  fgDebugEnabled = enabled;
}

//_____________________________________________________________________________
void AliLog::SetGlobalLogLevel(EType type)
{
// set the global debug level

  if (!fgInstance) new AliLog; 
  fgInstance->fGlobalLogLevel = type;
}

//_____________________________________________________________________________
Int_t AliLog::GetGlobalLogLevel()
{
// get the global debug level

  if (!fgInstance) new AliLog;
  return fgInstance->fGlobalLogLevel;
}

//_____________________________________________________________________________
void AliLog::SetGlobalDebugLevel(Int_t level)
{
// set the global debug level

  if (!fgInstance) new AliLog;
  if (level < -kDebugOffset) level = -kDebugOffset;
  fgInstance->fGlobalLogLevel = kDebugOffset + level;
}

//_____________________________________________________________________________
Int_t AliLog::GetGlobalDebugLevel()
{
// get the global debug level

  if (!fgInstance) new AliLog;
  return fgInstance->fGlobalLogLevel - kDebugOffset;
}

//_____________________________________________________________________________
void AliLog::SetModuleDebugLevel(const char* module, Int_t level)
{
// set the debug level for the given module

  if (!module) return;
  if (!fgInstance) new AliLog;
  TObject* obj = fgInstance->fModuleDebugLevels.FindObject(module);
  if (!obj) {
    obj = new TNamed(module, module);
    fgInstance->fModuleDebugLevels.Add(obj);
  }
  level += kDebugOffset;
  if (level < kFatal) level = kFatal;
  obj->SetUniqueID(level);
}

//_____________________________________________________________________________
void AliLog::ClearModuleDebugLevel(const char* module)
{
// remove the setting of the debug level for the given module

  if (!module) return;
  if (!fgInstance) new AliLog;
  TObject* obj = fgInstance->fModuleDebugLevels.FindObject(module);
  if (obj) delete fgInstance->fModuleDebugLevels.Remove(obj);
}

//_____________________________________________________________________________
void AliLog::SetClassDebugLevel(const char* className, Int_t level)
{
// set the debug level for the given class

  if (!className) return;
  if (!fgInstance) new AliLog;
  TObject* obj = fgInstance->fClassDebugLevels.FindObject(className);
  if (!obj) {
    obj = new TNamed(className, className);
    fgInstance->fClassDebugLevels.Add(obj);
  }
  level += kDebugOffset;
  if (level < kFatal) level = kFatal;
  obj->SetUniqueID(level);
}

//_____________________________________________________________________________
void AliLog::ClearClassDebugLevel(const char* className)
{
// remove the setting of the debug level for the given class

  if (!className) return;
  if (!fgInstance) new AliLog;
  TObject* obj = fgInstance->fClassDebugLevels.FindObject(className);
  if (obj) delete fgInstance->fClassDebugLevels.Remove(obj);
}


//_____________________________________________________________________________
void AliLog::SetStandardOutput()
{
// write all log messages to the standard output (stdout)

  if (!fgInstance) new AliLog;
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    fgInstance->CloseFile(iType);
    fgInstance->fOutputTypes[iType] = 0;
  }
}

//_____________________________________________________________________________
void AliLog::SetStandardOutput(EType type)
{
// write log messages of the given type to the standard output (stdout)

  if ((type < kFatal) || (type >= kMaxType)) return;
  if (!fgInstance) new AliLog;
  fgInstance->CloseFile(type);
  fgInstance->fOutputTypes[type] = 0;
}

//_____________________________________________________________________________
void AliLog::SetErrorOutput()
{
// write all log messages to the error output (stderr)

  if (!fgInstance) new AliLog;
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    fgInstance->CloseFile(iType);
    fgInstance->fOutputTypes[iType] = 1;
  }
}

//_____________________________________________________________________________
void AliLog::SetErrorOutput(EType type)
{
// write log messages of the given type to the error output (stderr)

  if ((type < kFatal) || (type >= kMaxType)) return;
  if (!fgInstance) new AliLog;
  fgInstance->CloseFile(type);
  fgInstance->fOutputTypes[type] = 1;
}

//_____________________________________________________________________________
void AliLog::SetFileOutput(const char* fileName)
{
// write all log messages to the given file

  if (!fgInstance) new AliLog;
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    if ((fgInstance->fOutputTypes[iType] = 2) && 
	(fgInstance->fFileNames[iType].CompareTo(fileName) != 0)) {
      fgInstance->CloseFile(iType);
    }
    fgInstance->fOutputTypes[iType] = 2;
    fgInstance->fFileNames[iType] = fileName;
    fgInstance->fOutputFiles[iType] = NULL;
  }
}

//_____________________________________________________________________________
void AliLog::SetFileOutput(EType type, const char* fileName)
{
// write log messages of the given type to the given file

  if ((type < kFatal) || (type >= kMaxType)) return;
  if (!fgInstance) new AliLog;
  if ((fgInstance->fOutputTypes[type] = 2) && 
      (fgInstance->fFileNames[type].CompareTo(fileName) != 0)) {
    fgInstance->CloseFile(type);
  }
  fgInstance->fOutputTypes[type] = 2;
  fgInstance->fFileNames[type] = fileName;
  fgInstance->fOutputFiles[type] = NULL;
}

//_____________________________________________________________________________
void AliLog::CloseFile(Int_t type)
{
// close the file for the given type if needed

  if ((fOutputTypes[type] == 2) && fOutputFiles[type]) {
    Bool_t closeFile = kTRUE;
    for (Int_t iType = kFatal; iType < kMaxType; iType++) {
      if ((iType != type) && (fOutputFiles[iType] == fOutputFiles[type])) {
	closeFile = kFALSE;
      }
    }
    if (closeFile) fclose(fOutputFiles[type]);
  }
  fOutputFiles[type] = NULL;
  fFileNames[type] = "";
  fOutputTypes[type] = 0;
}

//_____________________________________________________________________________
FILE* AliLog::GetOutputStream(Int_t type)
{
// get the output stream for the given type of messages

  if (fOutputTypes[type] == 0) return stdout;
  else if (fOutputTypes[type] == 1) return stderr;
  else if (fOutputTypes[type] == 2) {
    if (!fOutputFiles[type]) {
      FILE* file = NULL;
      if (!fFileNames[type].IsNull()) {
	for (Int_t iType = kFatal; iType < kMaxType; iType++) {
	  if ((iType != type) && 
	      (fFileNames[iType].CompareTo(fFileNames[type]) == 0) &&
	      fOutputFiles[iType]) {
	    file = fOutputFiles[iType];
	    break;
	  }
	}
	if (!file) file = fopen(fFileNames[type], "a");
      }
      fOutputFiles[type] = file;
      if (!file) CloseFile(type);
    }
    if (fOutputFiles[type]) return fOutputFiles[type];
  }

  return stdout;
}

//_____________________________________________________________________________
void AliLog::Flush()
{
// flush the output streams

  if (!fgInstance) new AliLog;
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    if (fgInstance->fOutputFiles[iType]) {
      fflush(fgInstance->fOutputFiles[iType]);
    }
  }
  fflush(stderr);
  fflush(stdout);
}


//_____________________________________________________________________________
void AliLog::SetHandleRootMessages(Bool_t on)
{
// enable or disable the handling of messages form root

  if (on) {
    SetErrorHandler(RootErrorHandler);
  } else {
    SetErrorHandler(DefaultErrorHandler);
  }
}


//_____________________________________________________________________________
void AliLog::SetPrintType(Bool_t on)
{
// switch on or off the printing of the message type for all message types

  if (!fgInstance) new AliLog;
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    fgInstance->fPrintType[iType] = on;
  }
}

//_____________________________________________________________________________
void AliLog::SetPrintType(EType type, Bool_t on)
{
// switch on or off the printing of the message type for the given message type

  if ((type < kFatal) || (type >= kMaxType)) return;
  if (!fgInstance) new AliLog;
  fgInstance->fPrintType[type] = on;
}

//_____________________________________________________________________________
void AliLog::SetPrintModule(Bool_t on)
{
// switch on or off the printing of the module for all message types

  if (!fgInstance) new AliLog;
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    fgInstance->fPrintModule[iType] = on;
  }
}

//_____________________________________________________________________________
void AliLog::SetPrintModule(EType type, Bool_t on)
{
// switch on or off the printing of the module for the given message type

  if ((type < kFatal) || (type >= kMaxType)) return;
  if (!fgInstance) new AliLog;
  fgInstance->fPrintModule[type] = on;
}

//_____________________________________________________________________________
void AliLog::SetPrintScope(Bool_t on)
{
// switch on or off the printing of the scope/class name for all message types

  if (!fgInstance) new AliLog;
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    fgInstance->fPrintScope[iType] = on;
  }
}

//_____________________________________________________________________________
void AliLog::SetPrintScope(EType type, Bool_t on)
{
// switch on or off the printing of the scope/class name
// for the given message type

  if ((type < kFatal) || (type >= kMaxType)) return;
  if (!fgInstance) new AliLog;
  fgInstance->fPrintScope[type] = on;
}

//_____________________________________________________________________________
void AliLog::SetPrintLocation(Bool_t on)
{
// switch on or off the printing of the file name and line number
// for all message types

  if (!fgInstance) new AliLog;
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    fgInstance->fPrintLocation[iType] = on;
  }
}

//_____________________________________________________________________________
void AliLog::SetPrintLocation(EType type, Bool_t on)
{
// switch on or off the printing of the file name and line number 
// for the given message type

  if ((type < kFatal) || (type >= kMaxType)) return;
  if (!fgInstance) new AliLog;
  fgInstance->fPrintLocation[type] = on;
}


//_____________________________________________________________________________
void AliLog::Write(const char* name, Int_t option)
{
// write the log object with the given name and option to the current file

  if (!fgInstance) new AliLog;
  fgInstance->TObject::Write(name, option);
}


//_____________________________________________________________________________
UInt_t AliLog::GetLogLevel(const char* module, const char* className) const
{
// get the logging level for the given module and class

  if (!fgInstance) new AliLog;
  TObject* obj = fgInstance->fClassDebugLevels.FindObject(className);
  if (obj) return obj->GetUniqueID();
  obj = fgInstance->fModuleDebugLevels.FindObject(module);
  if (obj) return obj->GetUniqueID();
  return fgInstance->fGlobalLogLevel;
}

//_____________________________________________________________________________
Int_t AliLog::GetDebugLevel(const char* module, const char* className)
{
// get the debug level for the given module and class

  if (!fgInstance) new AliLog;
  return fgInstance->GetLogLevel(module, className) - kDebugOffset;
}

//_____________________________________________________________________________
void AliLog::Message(UInt_t level, const char* message, 
		     const char* module, const char* className,
		     const char* function, const char* file, Int_t line)
{
// print a log message

  if (!fgInstance) new AliLog;

  // get the message type
  static const char* typeNames[kMaxType] = 
    {"Fatal", "Error", "Warning", "Info", "Debug"};
  UInt_t type = level;
  if (type >= kMaxType) type = kMaxType - 1;

  // print the message if the debug level allows
  if (level <= fgInstance->GetLogLevel(module, className)) {
    if (fgInstance->fPrintType[type]) {
      fprintf(fgInstance->GetOutputStream(type), "%s in ", typeNames[type]);
    }
    fprintf(fgInstance->GetOutputStream(type), "<");
    if (fgInstance->fPrintModule[type] && module) {
      fprintf(fgInstance->GetOutputStream(type), "%s/", module);
    }
    if (fgInstance->fPrintScope[type] && className) {
      fprintf(fgInstance->GetOutputStream(type), "%s::", className);
    }
    fprintf(fgInstance->GetOutputStream(type), "%s>: %s", function, message);
    if (fgInstance->fPrintLocation[type] && file) {
      fprintf(fgInstance->GetOutputStream(type), " (%s:%.0d)", file, line);
    }
    fprintf(fgInstance->GetOutputStream(type), "\n");
  }

  // abort in case of a fatal message
  if (type == kFatal) {
    delete fgInstance;
    if (gSystem) {
      gSystem->StackTrace();
      gSystem->Abort();
    } else {
      ::abort();
    }
  }
}

//_____________________________________________________________________________
void AliLog::Debug(UInt_t level, const char* message, 
		   const char* module, const char* className,
		   const char* function, const char* file, Int_t line)
{
// print a debug message

  if (level == 0) level = 1;
  level += kDebugOffset;
  Message(level, message, module, className, function, file, line);
}
