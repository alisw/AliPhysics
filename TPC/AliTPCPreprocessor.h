#ifndef ALI_TPC_PREPROCESSOR_H
#define ALI_TPC_PREPROCESSOR_H

#include "AliPreprocessor.h"


// test preprocessor that writes data to AliTestDataDCS

class AliTestDataDCS;
class AliTPCSensorTempArray;
class AliDCSSensorArray;
class AliTPCROC;
class TEnv;

class AliTPCPreprocessor : public AliPreprocessor
{
  public:
    AliTPCPreprocessor(AliShuttleInterface* shuttle);
//    AliTPCPreprocessor(const AliTPCPreprocessor &org);
    virtual ~AliTPCPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);
    UInt_t  MapTemperature(TMap* dcsAliasMap);
    UInt_t  MapPressure(TMap* dcsAliasMap);
    UInt_t  ExtractPedestals(Int_t sourceFXS);
    UInt_t  ExtractPulser(Int_t sourceFXS);
    UInt_t  ExtractCE(Int_t sourceFXS);
    AliTPCPreprocessor& operator = (const AliTPCPreprocessor& rhs);

  private:
    TEnv                   *fConfEnv;  // Preprocessor configuration map
    AliTPCSensorTempArray  *fTemp;     // CDB class for temperature sensors
    Bool_t                 fConfigOK;  // Identify succesful reading of OCDB Config
    AliTPCROC              *fROC;      // TPC Read-Out configuration

    ClassDef(AliTPCPreprocessor, 2)
};

#endif
