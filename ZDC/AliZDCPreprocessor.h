#ifndef ALI_ZDC_PREPROCESSOR_H
#define ALI_ZDC_PREPROCESSOR_H

#include "AliPreprocessor.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  			Zero Degree Calorimeter			             //
// ZDC Preprocessor -> writes data to AliZDCDataDCS
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class AliZDCDataDCS;

class AliZDCPreprocessor : public AliPreprocessor
{
  public:
    AliZDCPreprocessor(const char* detector, AliShuttleInterface* shuttle);
    virtual ~AliZDCPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);
    AliZDCPreprocessor(const AliZDCPreprocessor&);
    AliZDCPreprocessor& operator=(const AliZDCPreprocessor&);

  private:
    AliZDCDataDCS *fData;    // CDB class that stores the data

    ClassDef(AliZDCPreprocessor, 0);
};

    

#endif
