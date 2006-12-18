#ifndef ALI_TEST_PREPROCESSOR_H
#define ALI_TEST_PREPROCESSOR_H

#include "AliPreprocessor.h"

// test preprocessor that writes data to AliTestDataDCS

class AliTestDataDCS;

class AliTestPreprocessor : public AliPreprocessor
{
  public:
    AliTestPreprocessor(AliShuttleInterface* shuttle);
    virtual ~AliTestPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);

  private:
    AliTestDataDCS *fData;    // CDB class that stores the data

    ClassDef(AliTestPreprocessor, 0);
};

#endif
