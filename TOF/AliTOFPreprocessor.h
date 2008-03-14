#ifndef ALI_TOF_PREPROCESSOR_H
#define ALI_TOF_PREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliPreprocessor.h"

// TOF preprocessor. It takes care of both  
// DCS Data Points
// and DAQ histograms to compute online calibration constants

class AliTOFDataDCS;
class TObjArray;
class TH2S;

class AliTOFPreprocessor : public AliPreprocessor
{
  public:
    AliTOFPreprocessor(AliShuttleInterface* shuttle);
    virtual ~AliTOFPreprocessor();
    void   SetStoreRefData(Bool_t in){fStoreRefData=in;};
    Bool_t GetStoreRefData() const {return fStoreRefData;}
    void SetFDRFlag(Bool_t flag) {fFDRFlag = flag;}
    Bool_t GetFDRFlag() const {return fFDRFlag;}

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);
    virtual Bool_t ProcessDCS();

  private:
    AliTOFPreprocessor(const AliTOFPreprocessor & proc); // copy constructor
    AliTOFPreprocessor& operator=(const AliTOFPreprocessor & proc);
    UInt_t ProcessDCSDataPoints(TMap* dcsAliasMap);
    UInt_t ProcessOnlineDelays();
    UInt_t ProcessPulserData();
    UInt_t ProcessNoiseData();
    UInt_t ProcessFEEData(); // dummy, for the time being

    static const Int_t fgkBinRangeAve;    // number of bins where to 
                                          // calculate the mean
    static const Double_t fgkIntegralThr; // min number of entries per channel 
                                          // to perform calculation of delay
    static const Double_t fgkThrPar;      // parameter used to trigger the 
                                          // calculation of the delay
    AliTOFDataDCS *fData;                 // CDB class that stores the data
    TH2S *fh2;                            // TH2S from DAQ for histograms 
                                          // for delays  
    TObjArray *fCal;                      // TOF Calibration object
    TObjArray *fCalStatus;                // TOF Calibration object from pulser/noise
    TObjArray *fFEEStatus;                // TOF Calibration object from FEE
    Int_t fNChannels;                     // number of TOF channels
    Bool_t fStoreRefData;                 // Flag to decide storage of Ref Data
    Bool_t fFDRFlag;                      // Flag for FDR runs 
    ClassDef(AliTOFPreprocessor, 4);
};
#endif
