//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTERRORGUARD_H
#define ALIHLTERRORGUARD_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTErrorGuard.h
/// @author Matthias Richter
/// @date   01.07.2010
/// @brief  Helper class for suppression of error floods.

#include "AliHLTLogging.h"
#include "TString.h"
#include "Varargs.h"

/**
 * @class AliHLTErrorGuard
 * Helper class for suppression of error message floods caused by error messages
 * occurring rather frequently, e.g. for every event. The class suppresses the
 * printout after an adjustable number of occurences and prints an error summary
 * when the instance gets out of scope.
 *
 * Examples:
 * <pre>
 * if (nastyerror) {
 *   ALIHLTERRORGUARD(5, "nasty error, first occurence in event %d", event);
 * }
 * </pre>
 * <pre>
 * if (nastyerror) {
 *   static AliHLTErrorGuard g("classname", "functionname", "message");
 *   g.Throw(5);
 * }
 * </pre>
 * Both examples will throw the error for the first 5 occurrences. The macro
 * ALIHLTERRORGUARD handles also class and function name, source file and line
 * number, and supports variable messages through variadic macros.
 *
 * The second example illustrates usage of the class directly. The 'static'
 * attribute causes the object not to be destroyed at run time, only when the
 * program is terminated the object is deleted. This will print the error summary
 * at the very end of the program execution.
 *
 * @ingroup alihlt_base
 */
class AliHLTErrorGuard : public AliHLTLogging {
 public:
  /// constructor
 AliHLTErrorGuard(const char* classname, const char* functionname, const char* message, const char* file=NULL, int line=0)
   : fClass(classname), fFunction(functionname), fFile(file?file:""), fMessage(message), fLine(line), fOccurrence(0) {}

  /// set variable message
  void SetMessage( int dummy, ... )
  {
    va_list args;
    va_start(args, dummy);
    fMessage=AliHLTLogging::BuildLogString(NULL, args );
    va_end(args);
  }

  /// destructor
  virtual ~AliHLTErrorGuard() {
    Throw(-1, "Postponed message: %s - %d time(s)");
  }

  /// prefix increment operator
  AliHLTErrorGuard& operator++() {fOccurrence++; return *this;}

  int GetOccurrence() const {return fOccurrence;}

  void Throw(int maxoccurrence=1, const char* format="%s") {
    if (fOccurrence<=maxoccurrence || maxoccurrence<0) 
      LoggingVarargs(kHLTLogError, fClass.Data(), fFunction.Data(), fFile.Data(), fLine, format, fMessage.Data(), GetOccurrence());
  }

 protected:

 private:
  /** standard constructor prohibited */
  AliHLTErrorGuard();
  /** copy constructor prohibited */
  AliHLTErrorGuard(const AliHLTErrorGuard&);
  /** assignment operator prohibited */
  AliHLTErrorGuard& operator=(const AliHLTErrorGuard&);

  TString fClass; //! transient
  TString fFunction; //! transient
  TString fFile; //! transient
  TString fMessage; //! transient
  int fLine; //!transient
  int fOccurrence; //!transient

  ClassDef(AliHLTErrorGuard, 0)
};

#define ALIHLTERRORGUARD( max, ... )  {                                                \
    static AliHLTErrorGuard g(Class_Name() , FUNCTIONNAME() , "", __FILE__, __LINE__); \
  if (g.GetOccurrence()==0) g.SetMessage( 0, __VA_ARGS__ );                            \
  (++g).Throw(max);				                                       \
}

#endif
