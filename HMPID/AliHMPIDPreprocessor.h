#ifndef AliHMPIDPreprocessor_h
#define AliHMPIDPreprocessor_h

#include <AliPreprocessor.h> //base class

class TMap;

class AliHMPIDPreprocessor : public AliPreprocessor
{
public:
           AliHMPIDPreprocessor(AliShuttleInterface* pShuttle):AliPreprocessor("HMPID",pShuttle)  {}
  virtual ~AliHMPIDPreprocessor(                             )                                   {}
protected:
  virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process   (TMap* pDcsMap                              );
  ClassDef(AliHMPIDPreprocessor, 0);
};

typedef AliHMPIDPreprocessor AliRICHPreprocessor; // for backward compatibility

#endif
