#ifndef ALI_AD_PREPROCESSOR_H
#define ALI_AD_PREPROCESSOR_H

#include "AliPreprocessor.h"

class  AliShuttleInterface;
class  AliADDataDCS;

//  AD Preprocessor  header 

//  Calibration object from DCS and DAQ is written into  OCDB/AD/Calib/Data
//  Calibration object from DAQ is written into OCDB/AD/Calib/TimeSlewing

class AliADPreprocessor : public AliPreprocessor
{
  public:
    AliADPreprocessor(AliShuttleInterface* shuttle);
    virtual ~AliADPreprocessor();
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

  protected:
    virtual UInt_t Process(TMap* dcsAliasMap);
    UInt_t ProcessTimeSlewing();
    UInt_t ProcessTrendings();

    AliADDataDCS *fDCSData;    // CDB class that stores the calibration data

 private:
    AliADPreprocessor(const AliADPreprocessor&); // Not implemented
    AliADPreprocessor& operator=(const AliADPreprocessor&); // Not implemented
    
    ClassDef(AliADPreprocessor, 1);
};

#endif
