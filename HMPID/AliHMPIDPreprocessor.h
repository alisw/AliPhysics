#ifndef AliHMPIDPreprocessor_h
#define AliHMPIDPreprocessor_h

#include <AliPreprocessor.h> //base class

class TMap;

class AliHMPIDPreprocessor : public AliPreprocessor
{
public:
           AliHMPIDPreprocessor(AliShuttleInterface* pShuttle):AliPreprocessor("HMP",pShuttle) {}
  virtual ~AliHMPIDPreprocessor(                             )                                 {}
  static char    *fP;     // Name of the aliases provided by the DCS
  static char    *fHV;     // Name of the aliases provided by the DCS
  static char    *fT1; // Name of the aliases provided by the DCS
  static char    *fT2; // Name of the aliases provided by the DCS
protected:
  virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process   (TMap* pDcsMap                              );
  ClassDef(AliHMPIDPreprocessor, 0);
};

typedef AliHMPIDPreprocessor AliRICHPreprocessor; // for backward compatibility

#endif
