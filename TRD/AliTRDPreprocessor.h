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

class AliTRDPreprocessor : public AliPreprocessor
{

  public:

    AliTRDPreprocessor(AliShuttleInterface *shuttle);
    virtual ~AliTRDPreprocessor();

    enum {
       EExtractDCS    =   1  // error in case of failure by extracting DCS variables
      ,EStoreRefDCS   =   2  // error in case of failure by storing DCS variables references
      ,EFitDCS        =   4  // error in case of failure by fitting DCS variables
      ,EStoreDCS      =   8  // error in case of failure by storing DCS variables fit results
      ,EListFileHLT   =  16  // error in case of failure by taking the listof HLT files
      ,EOpenFileHLT   =  32  // error in case of failure by opening the HLTfile
      ,ETakeHistoHLT  =  64  // error in case of failure by taking the histos HLT
      ,EStoreHistoHLT = 128  // error in case of failure by storing the reference data HLT
      ,EFitHistoHLT   = 256  // error in case of failure by fitting the histos HLT
      ,EStoreCalHLT   = 512  // error in case of failure by storing the HLTcal objects
    };

  protected:

    virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);

  private:
    
    ClassDef(AliTRDPreprocessor,1) // The SHUTTLE preprocessor for TRD

};

#endif
