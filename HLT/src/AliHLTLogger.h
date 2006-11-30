// @(#) $Id$

#ifndef ALIL3LOGGER_H
#define ALIL3LOGGER_H

#include "AliHLTRootTypes.h"
#include "AliHLTStandardIncludes.h"

class MLUCLogServer;

class AliHLTLogger {

 public:
 
  AliHLTLogger();
  virtual ~AliHLTLogger();

  void Set(Int_t l);
  void UnSet(Int_t l);
  void UseDevNull();
  void UseStdout();
  void UseStderr();
  void UseStream(char *name="AliLevel3.log");
  void NotUseDevNull();
  void NotUseStdout();
  void NotUseStderr();
  void NotUseStream();

 protected:

  static Int_t fgAll;   //level all
  static Int_t fgDebug; //level debug
  static Int_t fgInformational; //level info
  static Int_t fgWarning; //level warning
  static Int_t fgError; //level error
  static Int_t fgFatal; //level fatal

 private:

  MLUCLogServer *fdn; //!
  MLUCLogServer *fso; //!
  MLUCLogServer *fse; //!
  MLUCLogServer *fsm; //!
#if __GNUC__ == 3
  std::ofstream *fof; //!
#else  
  ofstream *fof; //!
#endif

  ClassDef(AliHLTLogger,0)
};

typedef AliHLTLogger AliL3Logger; // for backward compatibility

#endif

