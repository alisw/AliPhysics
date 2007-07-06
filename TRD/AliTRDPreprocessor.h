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
       kEExtractDCS    =   1   // error in case of failure by extracting DCS variables
      ,kEStoreRefDCS   =   2   // error in case of failure by storing DCS variables references
      ,kEFitDCS        =   4   // error in case of failure by fitting DCS variables
      ,kEStoreDCS      =   8   // error in case of failure by storing DCS variables fit results
      ,kEListFileHLT   =  16   // error in case of failure by taking the listof HLT files
      ,kEOpenFileHLT   =  32   // error in case of failure by opening the HLTfile
      ,kETakeHistoHLT  =  64   // error in case of failure by taking the histos HLT
      ,kEStoreHistoHLT = 128   // error in case of failure by storing the reference data HLT
      ,kEFitHistoHLT   = 256   // error in case of failure by fitting the histos HLT
      ,kEStoreCalHLT   = 512   // error in case of failure by storing the HLTcal objects
      ,kEListFileDAQ   = 1024  // error in case of failure by taking the list of DAQ files 
      ,kEOpenFileDAQ   = 2048  // error in case of failure by opening the DAQfile
      ,kETakeObjectDAQ = 4096  // error in case of failure by taking the objects DAQ
      ,kEStoreRefDAQ   = 8192  // error in case of failure by storing the reference data DAQ
      ,kEFitObjectDAQ  = 16384 // error in case of failure by fitting the DAQ objects
      ,kEStoreCalDAQ   = 32768 // error in case of failure by storing the DAQcal objects   

    };

  protected:

    virtual void   Initialize(Int_t run, UInt_t startTime, UInt_t endTime);
    virtual UInt_t Process(TMap* dcsAliasMap);

  private:
    
    ClassDef(AliTRDPreprocessor,1) // The SHUTTLE preprocessor for TRD

};

#endif
