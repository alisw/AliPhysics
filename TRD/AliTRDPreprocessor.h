#ifndef ALI_TRD_PREPROCESSOR_H
#define ALI_TRD_PREPROCESSOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TRD preprocessor for the database SHUTTLE                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliPreprocessor.h"

class TMap;

class AliTRDPreprocessor : public AliPreprocessor
{

  public:

    AliTRDPreprocessor(AliShuttleInterface *shuttle);
    virtual ~AliTRDPreprocessor();

    enum {
      kEExtractDCS    =   1   // error in case of failure by extracting DCS variables
      ,kEStore        =   2   // Store or StoreReferenceData
      ,kEGetFileHLT   =   4   // GetFileSources and GetFile HLT
      ,kEEmptyListHLT =   8   // Empty list HLT
      ,kEGetFileDAQ   =  16   // GetFileSources and GetFile DAQ
      ,kEEmptyListDAQ =  32   // Empty list DAQ
    };

  protected:

    virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);


    void    ExtractPedestals();
    void    ExtractDriftVelocityDAQ();
    void    ExtractHLT();
    void    ProcessDCS(TMap*dcsAliasMap);

  private:
    

    UInt_t  fResult;                // result preprocessor
    Bool_t  fVdriftHLT;             // HLT Vdrift
    ClassDef(AliTRDPreprocessor,1) // The SHUTTLE preprocessor for TRD

};

#endif
