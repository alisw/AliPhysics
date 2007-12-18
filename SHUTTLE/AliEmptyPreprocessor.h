#ifndef ALIEMPTYPREPROCESSOR_H
#define ALIEMPTYPREPROCESSOR_H

#include "AliPreprocessor.h"

// This preprocessor is used as a placeholder for non-existing preprocessors
// during the FDR. Its task is just to fail, so that the run does not stay
// in processing state forever.

class AliEmptyPreprocessor : public AliPreprocessor
{
  public:
    AliEmptyPreprocessor(AliShuttleInterface* shuttle, const char* detector);
    virtual ~AliEmptyPreprocessor();

  protected:
    virtual UInt_t Process(TMap*) { Printf("Dummy preprocessor. FAILING..."); return 1; }
    virtual Bool_t ProcessDCS() { return kFALSE; }

  private:
    ClassDef(AliEmptyPreprocessor, 0);
};

#endif
