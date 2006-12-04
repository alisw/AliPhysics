#ifndef ALI_TRD_PREPROCESSOR_H
#define ALI_TRD_PREPROCESSOR_H

#include "AliPreprocessor.h"

/////////////////////////////////////////////////////
//
// TRD preprocessor
//
//////////////////////////////////////////////////

class AliTRDPreprocessor : public AliPreprocessor
{

  public:

    AliTRDPreprocessor(const Char_t *detector, AliShuttleInterface *shuttle);
    virtual ~AliTRDPreprocessor();

  protected:

    virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* /*dcsAliasMap*/);

  private:
    
    ClassDef(AliTRDPreprocessor,0);

};

#endif
