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
    Bool_t GetStoreRefData() const {return fStoreRefData;};

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);

  private:
    AliTOFPreprocessor(const AliTOFPreprocessor & proc); // copy constructor
    AliTOFPreprocessor& operator=(const AliTOFPreprocessor & proc);

    static const Int_t fgkBinRangeAve;   // number of bins where to 
                                         // calculate the mean
    static const Double_t fgkThrPar;     // parameter used to trigger the 
                                         // calculation of the delay
    AliTOFDataDCS *fData;                // CDB class that stores the data
    TH2S *fh2;                           // TH2S from DAQ for histograms 
                                         // for delays  
    TObjArray *fCal;                     // TOF Calibration object
    Int_t fNChannels;                    // number of TOF channels
    Bool_t fStoreRefData;                // Flag to decide storage of Ref Data
    ClassDef(AliTOFPreprocessor, 2);
};
#endif
