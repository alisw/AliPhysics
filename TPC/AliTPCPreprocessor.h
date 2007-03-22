#ifndef ALI_TPC_PREPROCESSOR_H
#define ALI_TPC_PREPROCESSOR_H

#include "AliTPCPreprocessor.h"

// test preprocessor that writes data to AliTestDataDCS

class AliTestDataDCS;

class AliTPCPreprocessor : public AliPreprocessor
{
  public:
    AliTPCPreprocessor(const char* detector, AliShuttleInterface* shuttle);
    virtual ~AliTPCPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);
    UInt_t AliTPCPreprocessor::MapTemperature(TMap* dcsAliasMap);


  private:
    AliTPCSensorTempArray  *fTemp;    // CDB class for temperature sensors

    ClassDef(AliTPCPreprocessor, 0);
};

#endif
