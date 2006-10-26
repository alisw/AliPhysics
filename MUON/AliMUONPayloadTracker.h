#ifndef ALIMUONPAYLOADTRACKER_H
#define ALIMUONPAYLOADTRACKER_H 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/

/// \ingroup raw
/// \class AliMUONPayloadTracker
/// \brief Class for decoding trackerrawdata 
///
/// \author Christian Finck
///
///////////////////////////////////////////////////////////////////////////////
///
/// This class decode the payload for tracker raw data.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliMUONDDLTracker;
class AliMUONBusStruct;
class AliMUONDspHeader;
class AliMUONBlockHeader;

class AliMUONPayloadTracker: public TObject {
  public :
    AliMUONPayloadTracker();
    virtual ~AliMUONPayloadTracker();

    Int_t GetMaxBlock() const {return fMaxBlock;}
    Int_t GetMaxDsp()   const {return fMaxDsp;}
    Int_t GetMaxBus()   const {return fMaxBus;}

    // check input before assigment
    void SetMaxBlock(Int_t blk);

    // does not check, done via BusPatchManager
    void SetMaxDsp(Int_t dsp) {fMaxDsp = dsp;}
    void SetMaxBus(Int_t bus) {fMaxBus = bus;}

    void ResetDDL();

    Bool_t Decode(UInt_t* buffer, Int_t datasize);

    AliMUONBusStruct*       GetBusPatchInfo() const {return fBusStruct;}
    AliMUONDDLTracker*      GetDDLTracker()   const {return fDDLTracker;}

  private :

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

    AliMUONPayloadTracker(const AliMUONPayloadTracker& stream);
    AliMUONPayloadTracker& operator = (const AliMUONPayloadTracker& stream);

    Bool_t CheckDataParity();

    ClassDef(AliMUONPayloadTracker, 1)    // base class for reading MUON raw digits
};

#endif
