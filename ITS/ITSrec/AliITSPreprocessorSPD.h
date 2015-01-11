#ifndef ALI_ITS_PREPROCESSOR_SPD_H
#define ALI_ITS_PREPROCESSOR_SPD_H

///////////////////////////////////////////////
//  Author: Henrik Tydesjo                   //
//  Preprocessor Class for the SPD           //
//                                           //
///////////////////////////////////////////////

/* $Id$ */

#include "AliPreprocessor.h"
#include <TList.h>

class AliITSPreprocessorSPD : public AliPreprocessor
{
  public:
    AliITSPreprocessorSPD(AliShuttleInterface* shuttle);
    virtual ~AliITSPreprocessorSPD();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);

  private:
    TList  fIdList; // list of ids for files that should be retrieved from FXS
    Bool_t RemoveIdFromList(const Char_t *id);
    Bool_t StoreRefForIdStartingWith(const Char_t *idStart);
    Bool_t StoreRefFromTarForId(const Char_t *id);

    ClassDef(AliITSPreprocessorSPD, 0);
};

#endif
