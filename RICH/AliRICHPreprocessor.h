#ifndef AliRICHPreprocessor_h
#define AliRICHPreprocessor_h

#include <AliPreprocessor.h> //base class

class TMap;

class AliRICHPreprocessor : public AliPreprocessor
{
public:
           AliRICHPreprocessor(AliShuttleInterface* pShuttle):AliPreprocessor("RICH",pShuttle)  {}
  virtual ~AliRICHPreprocessor(                             )                                   {}
  static  void      Test(); 
  static  TMap*     SimulateDcsMap(); 
protected:
  virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process   (TMap* pDcsMap                              );
  
  static const Int_t fgkNalias=2;
  static const char *fgkAliasName[fgkNalias];         //array of DCS alias names
  ClassDef(AliRICHPreprocessor, 0);
};

#endif
