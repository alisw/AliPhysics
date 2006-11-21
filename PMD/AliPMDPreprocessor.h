#ifndef ALI_PMD_PREPROCESSOR_H
#define ALI_PMD_PREPROCESSOR_H

#include "AliPreprocessor.h"

// test preprocessor that writes data to AliPMDDataDAQ  


class AliPMDPreprocessor : public AliPreprocessor
{
  public:
    AliPMDPreprocessor(const char* detector, AliShuttleInterface* shuttle);
    virtual ~AliPMDPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* pdaqAliasMap);

  private:
//    AliPMDDataDAQ *fData;    // CDB class that stores the data

    ClassDef(AliPMDPreprocessor, 0);
};

#endif
