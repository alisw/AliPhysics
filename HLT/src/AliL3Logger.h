// @(#) $Id$

#ifndef ALIL3LOGGER_H
#define ALIL3LOGGER_H

#include "AliL3RootTypes.h"
#include "AliL3StandardIncludes.h"

class MLUCLogServer;

class AliL3Logger{
  public:
  static Int_t kAll;
  static Int_t kDebug;
  static Int_t kInformational;
  static Int_t kWarning;
  static Int_t kError;
  static Int_t kFatal;
  AliL3Logger();
  virtual ~AliL3Logger();
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
  private:
  MLUCLogServer *dn; //!
  MLUCLogServer *so; //!
  MLUCLogServer *se; //!
  MLUCLogServer *sm; //!
#if __GNUC__ == 3
  std::ofstream *of; //!
#else  
  ofstream *of; //!
#endif

  ClassDef(AliL3Logger,1)
};

#endif

