#ifndef ALIMUONRAWSTREAMTRACKER_H
#define ALIMUONRAWSTREAMTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup base
/// \class AliMUONRawStreamTracker
/// \brief Class for reading MUON raw digits
///
/// \author Christian Finck
///
///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to MUON digits in raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>
#include "AliMUONPayloadTracker.h"

class AliRawReader;
class AliMUONDDLTracker;


class AliMUONRawStreamTracker: public TObject {
  public :
    AliMUONRawStreamTracker();
    AliMUONRawStreamTracker(AliRawReader* rawReader);
    virtual ~AliMUONRawStreamTracker();

    virtual Bool_t   Next();
    virtual Bool_t   NextDDL();

    Int_t GetMaxDDL()   const {return fMaxDDL;}
    Int_t GetMaxBlock() const {return  fPayload->GetMaxBlock();}
    Int_t GetMaxDsp()   const {return  fPayload->GetMaxDsp();}
    Int_t GetMaxBus()   const {return  fPayload->GetMaxBus();}

    // check input before assigment
    void SetMaxDDL(Int_t ddl);
    void SetMaxBlock(Int_t blk);

    // does not check, done via BusPatchManager
    void SetMaxDsp(Int_t dsp) {fPayload->SetMaxDsp(dsp);}
    void SetMaxBus(Int_t bus) {fPayload->SetMaxBus(bus);}


    void SetReader(AliRawReader* rawReader) {fRawReader = rawReader;}

    AliMUONDDLTracker*      GetDDLTracker()   const {return fPayload->GetDDLTracker();}
    Int_t                   GetDDL()          const {return fDDL - 1;}

  private :

    AliRawReader*    fRawReader;    ///< object for reading the raw data
 
    Int_t  fDDL;          ///< number of DDL
    Int_t  fBusPatchId;   ///< entry of buspatch structure
    Int_t  fDspId;        ///< entry of Dsp header
    Int_t  fBlkId;        ///< entry of Block header

    Bool_t fNextDDL;      ///< flag for next DDL to be read

    Int_t  fMaxDDL;       ///< maximum number of DDL in DATE file

    AliMUONPayloadTracker* fPayload; ///< pointer to payload decoder

    AliMUONRawStreamTracker(const AliMUONRawStreamTracker& stream);
    AliMUONRawStreamTracker& operator = (const AliMUONRawStreamTracker& stream);

    ClassDef(AliMUONRawStreamTracker, 2)    // base class for reading MUON raw digits
};

#endif
