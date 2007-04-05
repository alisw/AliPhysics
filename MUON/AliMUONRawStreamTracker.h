#ifndef ALIMUONRAWSTREAMTRACKER_H
#define ALIMUONRAWSTREAMTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup base
/// \class AliMUONRawStreamTracker
/// \brief Class for reading MUON raw digits
///
//  Author: Christian Finck

#include <TObject.h>
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

    /// Return maximum number of DDL in DATE file
    Int_t GetMaxDDL()   const {return fMaxDDL;}
    /// Return maximum number of block per DDL in DATE file
    Int_t GetMaxBlock() const {return  fPayload->GetMaxBlock();}
    /// Return maximum number of Dsp per block in DATE file
    Int_t GetMaxDsp()   const {return  fPayload->GetMaxDsp();}
    /// Return maximum number of Buspatch per Dsp in DATE file
    Int_t GetMaxBus()   const {return  fPayload->GetMaxBus();}

    // check input before assigment
    void SetMaxDDL(Int_t ddl);
    void SetMaxBlock(Int_t blk);

    /// Set maximum number of Dsp per block in DATE file
    /// does not check, done via BusPatchManager
    void SetMaxDsp(Int_t dsp) {fPayload->SetMaxDsp(dsp);}
    /// Set maximum number of Buspatch per Dsp in DATE file
    /// does not check, done via BusPatchManager
    void SetMaxBus(Int_t bus) {fPayload->SetMaxBus(bus);}

    /// Set object for reading the raw data
    void SetReader(AliRawReader* rawReader) {fRawReader = rawReader;}

    /// Return pointer for DDL
    AliMUONDDLTracker*      GetDDLTracker() const {return fPayload->GetDDLTracker();}

    /// Return pointer for payload
    AliMUONPayloadTracker*  GetPayLoad()    const {return fPayload;}

    /// Return number of DDL
    Int_t                   GetDDL()        const {return fDDL - 1;}

  private :
    /// Not implemented
    AliMUONRawStreamTracker(const AliMUONRawStreamTracker& stream);
    /// Not implemented
    AliMUONRawStreamTracker& operator = (const AliMUONRawStreamTracker& stream);

    AliRawReader*    fRawReader;    ///< object for reading the raw data
 
    Int_t  fDDL;          ///< number of DDL
    Int_t  fBusPatchId;   ///< entry of buspatch structure
    Int_t  fDspId;        ///< entry of Dsp header
    Int_t  fBlkId;        ///< entry of Block header

    Bool_t fNextDDL;      ///< flag for next DDL to be read

    Int_t  fMaxDDL;       ///< maximum number of DDL in DATE file

    AliMUONPayloadTracker* fPayload; ///< pointer to payload decoder

    ClassDef(AliMUONRawStreamTracker, 2)    // base class for reading MUON raw digits
};

#endif
