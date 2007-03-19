#ifndef ALIMUONPAYLOADTRACKER_H
#define ALIMUONPAYLOADTRACKER_H 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONPayloadTracker
/// \brief Class for decoding the payload for tracker raw data 
///
//  Author Christian Finck

#include <TObject.h>

class AliMUONDDLTracker;
class AliMUONBusStruct;
class AliMUONDspHeader;
class AliMUONBlockHeader;

class AliMUONPayloadTracker: public TObject {
  public :
    AliMUONPayloadTracker();
    virtual ~AliMUONPayloadTracker();

    /// Return maximum number of block per DDL in DATE file
    Int_t GetMaxBlock() const {return fMaxBlock;}
    /// Return maximum number of Dsp per block in DATE file
    Int_t GetMaxDsp()   const {return fMaxDsp;}
    /// Return maximum number of Buspatch per Dsp in DATE file
    Int_t GetMaxBus()   const {return fMaxBus;}

    // check input before assigment
    void SetMaxBlock(Int_t blk);

    /// \brief Set maximum number of Dsp per block in DATE file
    /// does not check, done via BusPatchManager
    void SetMaxDsp(Int_t dsp) {fMaxDsp = dsp;}
    /// \brief Set maximum number of Buspatch per Dsp in DATE file
    /// does not check, done via BusPatchManager
    void SetMaxBus(Int_t bus) {fMaxBus = bus;}

    void ResetDDL();

    Bool_t Decode(UInt_t* buffer, Int_t datasize);

    /// Return pointer for local structure
    AliMUONBusStruct*       GetBusPatchInfo() const {return fBusStruct;}
    /// Return pointer for buspatch structure
    AliMUONDDLTracker*      GetDDLTracker()   const {return fDDLTracker;}

  private :
    /// Not implemented
    AliMUONPayloadTracker(const AliMUONPayloadTracker& stream);
    /// Not implemented
    AliMUONPayloadTracker& operator = (const AliMUONPayloadTracker& stream);

    Bool_t CheckDataParity();

    Int_t  fBusPatchId;   ///< entry of buspatch structure
    Int_t  fDspId;        ///< entry of Dsp header
    Int_t  fBlkId;        ///< entry of Block header

    Int_t fMaxDDL;        ///< maximum number of DDL in DATE file
    Int_t fMaxBlock;      ///< maximum number of block per DDL in DATE file
    Int_t fMaxDsp;        ///< maximum number of Dsp per block in DATE file
    Int_t fMaxBus;        ///< maximum number of Buspatch per Dsp in DATE file

    AliMUONDDLTracker*      fDDLTracker;      //!< pointer for buspatch structure
    AliMUONBusStruct*       fBusStruct;       //!< pointer for local structure
    AliMUONBlockHeader*     fBlockHeader;     //!< pointer for block structure 
    AliMUONDspHeader*       fDspHeader;       //!< pointer for dsp structure 

    ClassDef(AliMUONPayloadTracker, 1)    // base class for reading MUON raw digits
};

#endif
