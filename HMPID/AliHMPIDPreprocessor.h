#ifndef AliHMPIDPreprocessor_h
#define AliHMPIDPreprocessor_h

#include <AliPreprocessor.h> //base class

//.
//HMPID Preprocessor base class 
//.
class TMap;

class AliHMPIDPreprocessor : public AliPreprocessor
{
public:
           AliHMPIDPreprocessor(AliShuttleInterface* pShuttle):AliPreprocessor("HMP",pShuttle) {}
  virtual ~AliHMPIDPreprocessor(                             )                                 {}
  
protected:
  static char    *fgP;     // Name of the aliases provided by the DCS
  static char    *fgHV;    // Name of the aliases provided by the DCS
  static char    *fgT1;    // Name of the aliases provided by the DCS
  static char    *fgT2;    // Name of the aliases provided by the DCS
  virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
  virtual UInt_t Process   (TMap* pDcsMap                              );
  ClassDef(AliHMPIDPreprocessor, 0);
};

#endif
