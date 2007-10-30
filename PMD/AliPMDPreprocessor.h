#ifndef ALI_PMD_PREPROCESSOR_H
#define ALI_PMD_PREPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Header File : AliPMDPreprocessor.h                 //
//                                                     //
// test preprocessor that writes data to AliPMDDataDAQ //
//-----------------------------------------------------//
// Author - A. Ahmed

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

  private:
//    AliPMDDataDAQ *fData;    // CDB class that stores the data

    ClassDef(AliPMDPreprocessor, 1);
};

#endif
