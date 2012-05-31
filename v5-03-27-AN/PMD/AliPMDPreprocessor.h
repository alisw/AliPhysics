#ifndef ALI_PMD_PREPROCESSOR_H
#define ALI_PMD_PREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/**********************************
 *
 * Pre-Processor CODE for PMD 
 *
 **********************************/ 


#include "AliPreprocessor.h"

class TTimeStamp;
class TSystem;

class AliPMDPreprocessor : public AliPreprocessor
{
  public:
    AliPMDPreprocessor(AliShuttleInterface* shuttle);
    virtual ~AliPMDPreprocessor();

  protected:
    virtual void Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* pdaqAliasMap);
    virtual Bool_t ProcessDAQ();
  
 private:

    Bool_t StorePmdPED();             // PMD PEDESTAL data
    Bool_t StorePmdGAIN();            // PMD GAIN Data
    Bool_t StorePmdHOT();             // PMD HOT data 
    Bool_t StorePmdMEAN();            // PMD SM MEAN data not used now 
    Bool_t StorePmdDCS(TMap *sdaqAM); // PMD DCS data points
    
    ClassDef(AliPMDPreprocessor, 3);
};

#endif
