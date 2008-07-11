#ifndef ALI_ZDC_PREPROCESSOR_H
#define ALI_ZDC_PREPROCESSOR_H

#include "AliPreprocessor.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  			Zero Degree Calorimeter			             //
// ZDC Preprocessor -> DCS data are passed to AliZDCDataDCS class to be      //
//      processed, DAQ output files are processed according to Run Type      //
// 	1 alignment object with DCS data is written to OCDB                  //
// 	1 calibration object with DAQ data is written to OCDB                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class AliZDCDataDCS;

class AliZDCPreprocessor : public AliPreprocessor
{
  public:
    AliZDCPreprocessor(AliShuttleInterface* shuttle);
    virtual ~AliZDCPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);
    virtual UInt_t ProcessChMap(TString runType);
    AliZDCPreprocessor(const AliZDCPreprocessor&);
    AliZDCPreprocessor& operator=(const AliZDCPreprocessor&);

  private:
    AliZDCDataDCS *fData;    // CDB class that stores the data

    ClassDef(AliZDCPreprocessor, 0);
};

    

#endif
