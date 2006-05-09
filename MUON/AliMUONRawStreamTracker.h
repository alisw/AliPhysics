#ifndef ALIMUONRAWSTREAMTRACKER_H
#define ALIMUONRAWSTREAMTRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
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
#include "AliMpBusPatch.h"

class AliRawReader;
class AliMUONDDLTracker;
class AliMUONBusStruct;
class AliMUONDspHeader;
class AliMUONBlockHeader;

class AliMUONRawStreamTracker: public TObject {
  public :
    AliMUONRawStreamTracker();
    AliMUONRawStreamTracker(AliRawReader* rawReader);
    AliMUONRawStreamTracker(const AliMUONRawStreamTracker& stream);
    AliMUONRawStreamTracker& operator = (const AliMUONRawStreamTracker& stream);
    virtual ~AliMUONRawStreamTracker();

    virtual Bool_t   Next();
    virtual Bool_t   NextDDL();
    virtual void     ResetDDL();

    Int_t GetMaxDDL()   const {return fMaxDDL;}
    Int_t GetMaxBlock() const {return fMaxBlock;}
    Int_t GetMaxDsp()   const {return fMaxDsp;}
    Int_t GetMaxBus()   const {return fMaxBus;}

    // check input before assigment
    void SetMaxDDL(Int_t ddl);
    void SetMaxBlock(Int_t blk);

    // does not check, done via BusPatchManager
    void SetMaxDsp(Int_t dsp) {fMaxDsp = dsp;}
    void SetMaxBus(Int_t bus) {fMaxBus = bus;}


    void SetReader(AliRawReader* rawReader) {fRawReader = rawReader;}

    AliMUONBusStruct*       GetBusPatchInfo() const {return fBusStructPtr;}
    AliMUONDDLTracker*      GetDDLTracker()   const {return fDDLTracker;}
    Int_t                   GetDDL()          const {return fDDL - 1;}

  protected :

    AliRawReader*    fRawReader;    // object for reading the raw data
 
    Int_t  fDDL;          // number of DDL
    Int_t  fBusPatchId;   // entry of buspatch structure
    Int_t  fDspId;        // entry of Dsp header
    Int_t  fBlkId;        // entry of Block header

    Bool_t fNextDDL;      // flag for next DDL to be read

    Int_t fMaxDDL;        // maximum number of DDL in DATE file
    Int_t fMaxBlock;      // maximum number of block per DDL in DATE file
    Int_t fMaxDsp;        // maximum number of Dsp per block in DATE file
    Int_t fMaxBus;        // maximum number of Buspatch per Dsp in DATE file


    AliMpBusPatch* fBusPatchManager; //! buspatch versus DE's & DDL

    AliMUONDDLTracker*      fDDLTracker;      //! pointer for buspatch structure
    AliMUONBusStruct*       fBusStruct;       //! pointer for local structure
    AliMUONBlockHeader*     fBlockHeader;     //! pointer for block structure 
    AliMUONDspHeader*       fDspHeader;       //! pointer for dsp structure 

    AliMUONBusStruct*       fBusStructPtr;       //! pointer for local structure

    ClassDef(AliMUONRawStreamTracker, 1)    // base class for reading MUON raw digits
};

#endif
