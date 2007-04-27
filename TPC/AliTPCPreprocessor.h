#ifndef ALI_TPC_PREPROCESSOR_H
#define ALI_TPC_PREPROCESSOR_H

#include "AliPreprocessor.h"


// test preprocessor that writes data to AliTestDataDCS

class AliTestDataDCS;
class AliTPCSensorTempArray;

class AliTPCPreprocessor : public AliPreprocessor
{
  public:
    AliTPCPreprocessor(AliShuttleInterface* shuttle);
//    AliTPCPreprocessor(const AliTPCPreprocessor &org);
    virtual ~AliTPCPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);
    UInt_t  MapTemperature(TMap* dcsAliasMap);
    AliTPCPreprocessor& operator = (const AliTPCPreprocessor& rhs);

  private:
    AliTPCSensorTempArray  *fTemp;    // CDB class for temperature sensors

    ClassDef(AliTPCPreprocessor, 1)
};

#endif
