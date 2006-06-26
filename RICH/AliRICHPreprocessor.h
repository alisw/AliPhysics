#ifndef AliRICHPreprocessor_h
#define AliRICHPreprocessor_h

#include <AliPreprocessor.h> //base class

// test preprocessor that writes data to AliTestDataDCS

class AliTestDataDCS;

class AliRICHPreprocessor : public AliPreprocessor
{
public:
           AliRICHPreprocessor(AliShuttleInterface* pShuttle):AliPreprocessor("RICH",pShuttle)  {}
  virtual ~AliRICHPreprocessor(                             )                                   {}

protected:
  virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process   (TMap* pDcsMap                              );
  ClassDef(AliRICHPreprocessor, 0);
};

#endif
