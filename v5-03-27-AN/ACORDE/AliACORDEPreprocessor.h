//#ifndef ALI_ACORDE_PREPROCESSOR_H
//#define ALI_ACORDE_PREPROCESSOR_H
#ifndef ALIACORDEPREPROCESSOR_H
#define AliACORDEPREPROCESSOR_H

#include "AliPreprocessor.h"

// test preprocessor that writes data to AliACORDECalibModule

class AliACORDECalibData;
class AliACORDEDataDCS;

class AliACORDEPreprocessor : public AliPreprocessor
{
  public:
    enum{kNModules=60};
    AliACORDEPreprocessor(AliShuttleInterface* shuttle);
    virtual ~AliACORDEPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);
    void CreateTableofReference();
    //virtual Bool_t ProcessDCS();

  private:

    AliACORDEPreprocessor(const AliACORDEPreprocessor &proc); //copy constructor
    AliACORDEPreprocessor& operator = (const AliACORDEPreprocessor & proc);
    AliACORDECalibData *fCalData;    // CDB class that stores the data
    AliACORDEDataDCS   *fDataDCS;   // ACORDE Data DCS  

    ClassDef(AliACORDEPreprocessor, 0);
};

#endif

