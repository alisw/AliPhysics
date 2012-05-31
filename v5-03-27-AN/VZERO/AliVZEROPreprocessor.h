#ifndef ALI_VZERO_PREPROCESSOR_H
#define ALI_VZERO_PREPROCESSOR_H

#include "AliPreprocessor.h"

class  AliShuttleInterface;
class  AliVZERODataFEE;
class  AliVZERODataDCS;

//  VZERO Preprocessor  header 

//  1 calibration object from DCS and DAQ is written into  OCDB/VZERO/Calib/Data
//  1 trigger-simulation object from DCS is written into OCDB/VZERO/Trigger/Data

class AliVZEROPreprocessor : public AliPreprocessor
{
  public:
    AliVZEROPreprocessor(AliShuttleInterface* shuttle);
    virtual ~AliVZEROPreprocessor();
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);

  protected:
    virtual UInt_t Process(TMap* dcsAliasMap);

    AliVZERODataDCS *fData;    // CDB class that stores the data for HV
    AliVZERODataFEE *fFEEData; // CDB class that stores the data for FEE


 private:
    AliVZEROPreprocessor(const AliVZEROPreprocessor&); // Not implemented
    AliVZEROPreprocessor& operator=(const AliVZEROPreprocessor&); // Not implemented
    
    ClassDef(AliVZEROPreprocessor, 1);
};

#endif
