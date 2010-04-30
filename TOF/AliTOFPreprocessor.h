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
class AliTOFLvHvDataPoints;
class AliTOFChannelOnlineStatusArray;
class AliTOFChannelOnlineArray;
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
    virtual UInt_t Process(TMap *dcsAliasMap);
    virtual Bool_t ProcessDCS();

  private:
    AliTOFPreprocessor(const AliTOFPreprocessor & proc); // copy constructor
    AliTOFPreprocessor& operator=(const AliTOFPreprocessor & proc);
    UInt_t ProcessDCSDataPoints(TMap *dcsAliasMap);
    UInt_t ProcessHVandLVdps(TMap *dcsAliasMap);
    UInt_t ProcessOnlineDelays();
    UInt_t ProcessPulserData();
    UInt_t ProcessNoiseData();
    UInt_t ProcessFEEData(); // dummy, for the time being
    UInt_t ProcessT0Fill();
    UInt_t ProcessNoiseCalibTrg();

    void FillWithCosmicCalibration(AliTOFChannelOnlineArray *cal); // fill with cosmic calibration
    void FillWithCableLengthMap(AliTOFChannelOnlineArray *cal); // fill with cable-lenght map

    static const Int_t fgkBinRangeAve;       // number of bins where to 
                                             // calculate the mean
    static const Double_t fgkIntegralThr;    // min number of entries per channel 
                                             // to perform calculation of delay
    static const Double_t fgkThrPar;         // parameter used to trigger the 
                                             // calculation of the delay
    AliTOFDataDCS *fData;                    // CDB class that stores the data
    AliTOFLvHvDataPoints *fHVLVmaps;         // HV and LV status maps
    AliTOFChannelOnlineArray *fCal;          // TOF Calibration object
    Int_t fNChannels;                        // number of TOF channels
    Bool_t fStoreRefData;                    // Flag to decide storage of Ref Data
    Bool_t fFDRFlag;                         // Flag for FDR runs 
    AliTOFChannelOnlineStatusArray *fStatus; // Array with TOF channels' status
    Int_t *fMatchingWindow;                  //[fNChannels]
                                             // Array of matching windows (one per channel) - to be used in noise runs
    Int_t *fLatencyWindow;                   //[fNChannels]
                                             // Array of latency windows (one per channel)
    Bool_t fIsStatusMapChanged;              // flag to check is the status map OCDB has to be updated

    ClassDef(AliTOFPreprocessor, 1);

};
#endif
