#ifndef ALI_ITS_PREPROCESSOR_SSD_H
#define ALI_ITS_PREPROCESSOR_SSD_H

#include "AliPreprocessor.h"

//
// Author: Enrico Fragiacomo
// Date: 13/10/2006
// 
// SHUTTLE preprocessing class for SSD calibration files

/* $Id$ */

class AliITSPreprocessorSSD : public AliPreprocessor
{
  public:
  AliITSPreprocessorSSD(AliShuttleInterface* shuttle):
    AliPreprocessor("SSD",shuttle) {}
  virtual ~AliITSPreprocessorSSD() {;}
    enum {kDDLperLDC = 4};      // number of DDLs in LDC

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);

  private:
    ClassDef(AliITSPreprocessorSSD, 0);
};

#endif
