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

#include <strings.h>
#include <Riostream.h>
#include <TError.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TEnv.h>

#include "AliLog.h"

ClassImp(AliLog)


AliLog* AliLog::fgInstance = NULL;

Bool_t AliLog::fgDebugEnabled = kTRUE;


//_____________________________________________________________________________
AliLog::AliLog() :
  TObject(),
  fGlobalLogLevel(kInfo),
  fModuleDebugLevels(),
  fClassDebugLevels(),
  fPrintRepetitions(kTRUE),
  fRepetitions(0),
  fLastType(0),
  fLastMessage(),
  fLastModule(),
  fLastClassName(),
  fLastFunction(),
  fLastFile(),
  fLastLine(0)
{
// default constructor: set default values

  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    fOutputTypes[iType] = 0;
    fFileNames[iType] = "";
    fOutputFiles[iType] = NULL;
    fOutputStreams[iType] = NULL;

    fPrintType[iType] = kTRUE;
    fPrintModule[iType] = kFALSE;
    fPrintScope[iType] = kTRUE;
    fPrintLocation[iType] = (iType == kDebug);  
  }

  // replace the previous instance by this one
  if (fgInstance) delete fgInstance;
  fgInstance = this;

  SetHandleRootMessages(kTRUE);

  // read the .rootrc settings
  ReadEnvSettings();
}

//_____________________________________________________________________________
AliLog::~AliLog()
{
// destructor: clean up and reset instance pointer

  if (fRepetitions > 0) PrintRepetitions();

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
void AliLog::ReadEnvSettings()
{
// load settings from the root configuration file (.rootrc)
// and from environment variables

  static const char* typeNames[kMaxType] = 
    {"kFatal", "kError", "kWarning", "kInfo", "kDebug"};

  // debug en- or disabling
  if (gSystem->Getenv("LOG_NO_DEBUG")) {
    fgDebugEnabled = kFALSE;
  } else if (gEnv->Defined("AliRoot.AliLog.EnableDebug")) {
    fgDebugEnabled = gEnv->GetValue("AliRoot.AliLog.EnableDebug", 
                                    fgDebugEnabled);
    AliInfo(Form("debug %sabled", ((fgDebugEnabled) ? "en" : "dis")));
  }

  // global log level
  if (gEnv->Defined("AliRoot.AliLog.GlobalLogLevel")) {
    const char* type = gEnv->GetValue("AliRoot.AliLog.GlobalLogLevel", "");
    for (Int_t iType = kFatal; iType < kMaxType; iType++) {
      if (strcmp(type, typeNames[iType]) == 0) fGlobalLogLevel = iType;
    }
    AliDebug(3, Form("global log level set to %d", fGlobalLogLevel));
  }

  // global debug level
  if (gEnv->Defined("AliRoot.AliLog.GlobalDebugLevel")) {
    Int_t level = gEnv->GetValue("AliRoot.AliLog.GlobalDebugLevel",
                                 Int_t(fGlobalLogLevel - kDebugOffset));
    if (level < -kDebugOffset) level = kDebugOffset;
    fGlobalLogLevel = kDebugOffset + level;
    AliDebug(3, Form("global debug level set to %d",
        	     fGlobalLogLevel - kDebugOffset));
  }

  // module debug level
  if (gEnv->Defined("AliRoot.AliLog.ModuleDebugLevel")) {
    TString levels = gEnv->GetValue("AliRoot.AliLog.ModuleDebugLevel", "");
    char* p = const_cast<char*>(levels.Data());
    while (const char* module = strtok(p, " ")) {
      p = NULL;
      char* pos = index(module, ':');
      if (!pos) continue;
      *(pos++) = '\0';
      Int_t level = atoi(pos);
      SetModuleDebugLevel(module, level);
      AliDebug(3, Form("debug level for module %s set to %d", module, level));
    }
  }

  // class debug level
  if (gEnv->Defined("AliRoot.AliLog.ClassDebugLevel")) {
    TString levels = gEnv->GetValue("AliRoot.AliLog.ClassDebugLevel", "");
    char* p = const_cast<char*>(levels.Data());
    while (const char* className = strtok(p, " ")) {
      p = NULL;
      char* pos = index(className, ':');
      if (!pos) continue;
      *(pos++) = '\0';
      Int_t level = atoi(pos);
      SetClassDebugLevel(className, level);
      AliDebug(3, Form("debug level for class %s set to %d", 
                       className, level));
    }
  }

  // general output stream
  if (gEnv->Defined("AliRoot.AliLog.Output")) {
    TString stream = gEnv->GetValue("AliRoot.AliLog.Output", "Standard");
    if (stream.CompareTo("standard", TString::kIgnoreCase) == 0) {
      SetStandardOutput();
      AliDebug(3, "output stream set to standard output for all types");
    } else if (stream.CompareTo("error", TString::kIgnoreCase) == 0) {
      SetErrorOutput();
      AliDebug(3, "output stream set to error output for all types");
    } else if (!stream.IsNull()) {
      SetFileOutput(stream);
      AliDebug(3, Form("output stream set to file %s for all types", 
                       stream.Data()));
    }
  }

  // individual output streams
  for (Int_t iType = kFatal; iType < kMaxType; iType++) {
    TString name("AliRoot.AliLog.Output.");
    name += &typeNames[iType][1];
    if (gEnv->Defined(name)) {
      TString stream = gEnv->GetValue(name, "Standard");
      if (stream.CompareTo("standard", TString::kIgnoreCase) == 0) {
        SetStandardOutput(EType(iType));
        AliDebug(3, Form("output stream set to standard output for type %s",
                         typeNames[iType]));
      } else if (stream.CompareTo("error", TString::kIgnoreCase) == 0) {
        SetErrorOutput(EType(iType));
        AliDebug(3, Form("output stream set to error output for type %s",
                         typeNames[iType]));
      } else if (!stream.IsNull()) {
        SetFileOutput(EType(iType), stream);
        AliDebug(3, Form("output stream set to file %s for type %s", 
                         stream.Data(), typeNames[iType]));
      }
    }
  }

  // handling of root error messages
  if (gEnv->Defined("AliRoot.AliLog.HandleRootMessages")) {
    Bool_t on = gEnv->GetValue("AliRoot.AliLog.HandleRootMessages", kTRUE);
    SetHandleRootMessages(on);
    AliDebug(3, Form("handling of root messages %sabled",
                     ((on) ? "en" : "dis")));
  }

  // printout settings
  static const char* settingNames[4] = 
    {"Type", "Module", "Scope", "Location"};
  Bool_t* settings[] = 
    {fPrintType, fPrintModule, fPrintScope, fPrintLocation};
  for (Int_t iSetting = 0; iSetting < 4; iSetting++) {
    TString name("AliRoot.AliLog.Print");
    name += settingNames[iSetting];
    if (gEnv->Defined(name)) {
      Bool_t on = gEnv->GetValue(name, settings[iSetting][0]);
      for (Int_t iType = kFatal; iType < kMaxType; iType++) {
        settings[iSetting][iType] = on;
      }
      AliDebug(3, Form("printing of %s %sabled for all types",
                       settingNames[iSetting], ((on) ? "en" : "dis")));
    }

    for (Int_t iType = kFatal; iType < kMaxType; iType++) {
      TString nameType = name + "." + &typeNames[iType][1];
      if (gEnv->Defined(nameType)) {
        Bool_t on = gEnv->GetValue(nameType, settings[iSetting][iType]);
        settings[iSetting][iType] = on;
        AliDebug(3, Form("printing of %s %sabled for type %s",
                         settingNames[iSetting], ((on) ? "en" : "dis"),
                         typeNames[iType]));
      }
    }
  }

  // repetition of messages
  if (gEnv->Defined("AliRoot.AliLog.PrintRepetitions")) {
    Bool_t on = gEnv->GetValue("AliRoot.AliLog.PrintRepetitions", kTRUE);
    fPrintRepetitions = on;
    AliDebug(3, Form("printing of message repetitions %sabled",
                     ((on) ? "en" : "dis")));
  }
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
    if ((fgInstance->fOutputTypes[iType] == 2) && 
	(fgInstance->fFileNames[iType].CompareTo(fileName) != 0)) {
      fgInstance->CloseFile(iType);
    }
    fgInstance->fOutputTypes[iType] = 2;
    fgInstance->fFileNames[iType] = fileName;
    fgInstance->fOutputFiles[iType] = NULL;
    fgInstance->fOutputStreams[iType] = NULL;
  }
}

//_____________________________________________________________________________
void AliLog::SetFileOutput(EType type, const char* fileName)
{
// write log messages of the given type to the given file

  if ((type < kFatal) || (type >= kMaxType)) return;
  if (!fgInstance) new AliLog;
  if ((fgInstance->fOutputTypes[type] == 2) && 
      (fgInstance->fFileNames[type].CompareTo(fileName) != 0)) {
    fgInstance->CloseFile(type);
  }
  fgInstance->fOutputTypes[type] = 2;
  fgInstance->fFileNames[type] = fileName;
  fgInstance->fOutputFiles[type] = NULL;
  fgInstance->fOutputStreams[type] = NULL;
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
    if (closeFile) {
      fclose(fOutputFiles[type]);
      fOutputStreams[type]->close();
      delete fOutputStreams[type];
    }
  }
  fOutputFiles[type] = NULL;
  fOutputStreams[type] = NULL;
  fFileNames[type] = "";
  fOutputTypes[type] = 0;
}

//_____________________________________________________________________________
FILE* AliLog::GetOutputStream(Int_t type)
{
// get the output stream for the given type of messages

  if (type > kDebug) type = kDebug;
  if (fOutputTypes[type] == 0) return stdout;
  else if (fOutputTypes[type] == 1) return stderr;
  else if (fOutputTypes[type] == 2) {
    if (!fOutputFiles[type]) {
      FILE* file = NULL;
      ofstream* stream = NULL;
      if (!fFileNames[type].IsNull()) {
	for (Int_t iType = kFatal; iType < kMaxType; iType++) {
	  if ((iType != type) && 
	      (fFileNames[iType].CompareTo(fFileNames[type]) == 0) &&
	      fOutputFiles[iType]) {
	    file = fOutputFiles[iType];
	    stream = fOutputStreams[iType];
	    break;
	  }
	}
	if (!file) {
	  file = fopen(fFileNames[type], "a");
	  stream = new ofstream(fFileNames[type], ios::app);
	}
      }
      fOutputFiles[type] = file;
      fOutputStreams[type] = stream;
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
      fgInstance->fOutputStreams[iType]->flush();
    }
  }
  fflush(stderr);
  fflush(stdout);
}


//_____________________________________________________________________________
void AliLog::SetHandleRootMessages(Bool_t on)
{
// enable or disable the handling of messages form root

  if (!fgInstance) new AliLog;
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
void AliLog::SetPrintRepetitions(Bool_t on)
{
// switch on or off the printing of the number of repetitions of a message
// instead of repeating the same message

  if (!fgInstance) new AliLog;
  if (!on && (fgInstance->fRepetitions > 0)) fgInstance->PrintRepetitions();
  fgInstance->fPrintRepetitions = on;
}


//_____________________________________________________________________________
void AliLog::WriteToFile(const char* name, Int_t option)
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
  if (className) {
    TObject* obj = fgInstance->fClassDebugLevels.FindObject(className);
    if (obj) return obj->GetUniqueID();
  }
  if (module) {
    TObject* obj = fgInstance->fModuleDebugLevels.FindObject(module);
    if (obj) return obj->GetUniqueID();
  }
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
void AliLog::PrintMessage(UInt_t type, const char* message, 
                          const char* module, const char* className,
                          const char* function, const char* file, Int_t line)
{
// print the given message

  // don't print the message if it is repeated
  if (fPrintRepetitions &&
      (fLastType == type) && 
      (message && (fLastMessage.CompareTo(message) == 0)) &&
      ((module && (fLastModule.CompareTo(module) == 0)) ||
       (!module && fLastModule.IsNull())) &&
      ((className && (fLastClassName.CompareTo(className) == 0)) ||
       (!className && fLastClassName.IsNull())) &&
      ((function && (fLastFunction.CompareTo(function) == 0)) ||
       (!function && fLastFunction.IsNull()))&&
      ((file && (fLastFile.CompareTo(file) == 0)) ||
       (!file && fLastFile.IsNull())) &&
      (fLastLine == line)) {
    fRepetitions++;
    return;
  }

  // print number of repetitions
  if (fRepetitions > 0) PrintRepetitions();

  // remember this message
  fRepetitions = 0;
  fLastType = type;
  fLastMessage = message;
  fLastModule = module;
  fLastClassName = className;
  fLastFunction = function;
  fLastFile = file;
  fLastLine = line;

  // print the message
  FILE* stream = GetOutputStream(type);
  static const char* typeNames[kMaxType] = 
    {"Fatal", "Error", "Warning", "Info", "Debug"};

  if (fPrintType[type]) {
    fprintf(stream, "%c-", typeNames[type][0]);
  }
  if (fPrintModule[type] && module) {
    fprintf(stream, "%s/", module);
  }
  if (fPrintScope[type] && className) {
    fprintf(stream, "%s::", className);
  }
  if (message) {
    fprintf(stream, "%s: %s", function, message);
  } else {
    fprintf(stream, "%s", function);
  }
  if (fPrintLocation[type] && file) {
    fprintf(stream, " (%s:%.0d)", file, line);
  }
  if (message) {
    fprintf(stream, "\n");
  } else {
    fprintf(stream, ": ");
  }
}

//_____________________________________________________________________________
void AliLog::PrintRepetitions()
{
// print number of repetitions

  fprintf(GetOutputStream(fLastType), " <message repeated %d time%s>\n", 
          fRepetitions, (fRepetitions > 1) ? "s" : "");
}

//_____________________________________________________________________________
void AliLog::Message(UInt_t level, const char* message, 
		     const char* module, const char* className,
		     const char* function, const char* file, Int_t line)
{
// print a log message

  if (!fgInstance) new AliLog;

  // get the message type
  UInt_t type = level;
  if (type >= kMaxType) type = kMaxType - 1;

  // print the message if the debug level allows
  if (level <= fgInstance->GetLogLevel(module, className)) {
    fgInstance->PrintMessage(type, message, 
                             module, className, function, file, line);
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


//_____________________________________________________________________________
Int_t AliLog::RedirectStdoutTo(EType type, UInt_t level, const char* module, 
                               const char* className, const char* function,
                               const char* file, Int_t line, Bool_t print)
{
// redirect the standard output to the stream of the given type

  if (!fgInstance) new AliLog;
  return fgInstance->RedirectTo(stdout, type, level, module, className, 
                                function, file, line, print);
}

//_____________________________________________________________________________
Int_t AliLog::RedirectStderrTo(EType type, UInt_t level, const char* module, 
                               const char* className, const char* function,
                               const char* file, Int_t line, Bool_t print)
{
// redirect the standard error output to the stream of the given type

  if (!fgInstance) new AliLog;
  return fgInstance->RedirectTo(stderr, type, level, module, className, 
                                function, file, line, print);
}

//_____________________________________________________________________________
Int_t AliLog::RedirectTo(FILE* stream, EType type, UInt_t level, 
                         const char* module, const char* className,
                         const char* function, const char* file, Int_t line,
			 Bool_t print)
{
// redirect the standard (error) output stream to the stream of the given type

  // get the original file descriptor to be able to restore it later
  Int_t original = dup(fileno(stream));
  fflush(stream);

  // flush the stream of the selected type
  FILE* newStream = GetOutputStream(type);
  fflush(newStream);

  // redirect stream
  if ((type == kDebug) && (level > 0)) level--;
  if (type + level > GetLogLevel(module, className)) { // /dev/null
    freopen("/dev/null", "a", stream);
  } else if (fOutputTypes[type] == 0) {         // stdout
    if (stream != stdout) dup2(fileno(stdout), fileno(stream));
  } else if (fOutputTypes[type] == 1) {         // stderr
    if (stream != stderr) dup2(fileno(stderr), fileno(stream));
  } else if (fOutputTypes[type] == 2) {         // file
    freopen(fFileNames[type], "a", stream);
  }

  // print information
  if (print) {
    PrintMessage(type, NULL, module, className, function, file, line);
    fflush(newStream);
  }

  return original;
}

//_____________________________________________________________________________
void AliLog::RestoreStdout(Int_t original)
{
// restore the standard output

  fflush(stdout);
  dup2(original, fileno(stdout));  
  close(original);
}

//_____________________________________________________________________________
void AliLog::RestoreStderr(Int_t original)
{
// restore the standard error output

  fflush(stderr);
  dup2(original, fileno(stderr));  
  close(original);
}


//_____________________________________________________________________________
ostream& AliLog::Stream(EType type, UInt_t level,
                        const char* module, const char* className,
                        const char* function, const char* file, Int_t line)
{
// get the stream object for the given output type

  if (!fgInstance) new AliLog;
  return fgInstance->GetStream(type, level, module, className, 
                               function, file, line);
}

//_____________________________________________________________________________
ostream& AliLog::GetStream(EType type, UInt_t level,
                           const char* module, const char* className,
                           const char* function, const char* file, Int_t line)
{
// get the stream object for the given output type

  if ((type == kDebug) && (level > 0)) level--;
  Bool_t noOutput = (type + level > GetLogLevel(module, className));

  if (!noOutput) {
    PrintMessage(type, NULL, module, className, function, file, line);
  }
  fflush(GetOutputStream(type));

  static ofstream nullStream("/dev/null");
  if (noOutput) {
    return nullStream;
  } else if (fOutputTypes[type] == 0) {
    return cout;
  } else if (fOutputTypes[type] == 1) {
    return cerr;
  } else if (fOutputTypes[type] == 2) {
    return *fOutputStreams[type];
  }

  return nullStream;
}

