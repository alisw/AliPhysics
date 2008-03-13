#ifndef AliHMPIDPreprocessor_h
#define AliHMPIDPreprocessor_h

#include <AliPreprocessor.h> //base class

//.
//.
//HMPID Preprocessor base class 
//.
//.
class TMap;

class AliHMPIDPreprocessor : public AliPreprocessor
{
public:
           AliHMPIDPreprocessor(AliShuttleInterface* pShuttle):AliPreprocessor("HMP",pShuttle) {}
  virtual ~AliHMPIDPreprocessor(                             )                                 {}
protected:
  virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime); //
  virtual UInt_t Process   (TMap* pDcsMap                              ); //process everthing
          Bool_t ProcDcs   (TMap* pDcsMap                              ); //process DCS data points
          Bool_t ProcPed   (                                           ); //process pedestal files
  ClassDef(AliHMPIDPreprocessor, 0);
};

#endif
