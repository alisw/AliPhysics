#ifndef ALI_ITS_PREPROCESSOR_SPD_H
#define ALI_ITS_PREPROCESSOR_SPD_H

///////////////////////////////////////////////
//  Author: Henrik Tydesjo                   //
//  Preprocessor Class for the SPD           //
//                                           //
///////////////////////////////////////////////

#include "AliPreprocessor.h"

class AliITSPreprocessorSPD : public AliPreprocessor
{
  public:
    AliITSPreprocessorSPD(AliShuttleInterface* shuttle);
    virtual ~AliITSPreprocessorSPD();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);

  private:

    ClassDef(AliITSPreprocessorSPD, 0);
};

#endif
