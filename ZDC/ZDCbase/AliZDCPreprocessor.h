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
    virtual Bool_t ProcessDCS();
    AliZDCPreprocessor(const AliZDCPreprocessor&);
    AliZDCPreprocessor& operator=(const AliZDCPreprocessor&);
    void SetPedSubMethFlag(Bool_t ifv) {fPedSubMethFlag = ifv;}

  private:
    UInt_t ProcessDCSData(TMap* dcsAliasMap);
    UInt_t ProcessChMap();
    UInt_t ProcessppData();
    UInt_t ProcessCalibData(Float_t beamEnergy);
    UInt_t ProcessPedestalData();
    UInt_t ProcessLaserData();
    UInt_t ProcessMBCalibData();

    AliZDCDataDCS *fData;    // OCDB class that stores DCS data
    Bool_t         fPedSubMethFlag; //flag for pedestal subtraction mode (from RUN2)

    ClassDef(AliZDCPreprocessor, 1);
};

    

#endif
