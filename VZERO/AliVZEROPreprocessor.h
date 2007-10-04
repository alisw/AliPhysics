#ifndef ALI_VZERO_PREPROCESSOR_H
#define ALI_VZERO_PREPROCESSOR_H

#include "AliPreprocessor.h"
#include "AliVZERODataDCS.h"
#include "AliShuttleInterface.h"

// VZERO PreProcessor  header 


class AliVZEROPreprocessor : public AliPreprocessor
{
  public:
    AliVZEROPreprocessor(AliShuttleInterface* shuttle);
    virtual ~AliVZEROPreprocessor();
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

  protected:
    virtual UInt_t Process(TMap* dcsAliasMap);

    AliVZERODataDCS *fData;    // CDB class that stores the data
    
    ClassDef(AliVZEROPreprocessor, 0);
};

#endif
