#ifndef ALITOFPREPROCESSORFDR_H
#define ALITOFPREPROCESSORFDR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* 
$Id$
*/

#include "AliPreprocessor.h"

// TOF preprocessor. It takes care of both  
// DCS Data Points
// and DAQ histograms to compute online calibration constants

class TObjArray;
class TH2S;

class AliTOFPreprocessorFDR : public AliPreprocessor
{
  public:
    AliTOFPreprocessorFDR(AliShuttleInterface* shuttle);
    virtual ~AliTOFPreprocessorFDR();
    void   SetStoreRefData(Bool_t in){fStoreRefData=in;};
    Bool_t GetStoreRefData() const {return fStoreRefData;};

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);

  private:
    AliTOFPreprocessorFDR(const AliTOFPreprocessorFDR & proc); // copy constructor
    AliTOFPreprocessorFDR& operator=(const AliTOFPreprocessorFDR & proc);
    UInt_t ProcessDCSDataPoints(TMap* dcsAliasMap);

    Bool_t fStoreRefData;                // Flag to decide storage of Ref Data
    ClassDef(AliTOFPreprocessorFDR, 2);
};
#endif
