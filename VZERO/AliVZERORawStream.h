#ifndef ALIVZERORAWSTREAM_H
#define ALIVZERORAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading the VZERO DDL raw data
/// The format of the raw data corresponds to the one
/// implemented in AliVZEROBuffer class.
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliRawReader;

class AliVZERORawStream: public TObject {
  public :
    AliVZERORawStream(AliRawReader* rawReader);
    virtual ~AliVZERORawStream();

    virtual void             Reset();
    virtual Bool_t           Next();

    UShort_t                 GetADC(Int_t channel) const
      { return fADC[channel][kNEvOfInt/2]; }
    UShort_t                 GetPedestal(Int_t channel, Int_t event) const
      { return fADC[channel][event]; }
    UInt_t                   GetTime(Int_t channel) const
      { return fTime[channel]; }

    enum EVZERORawDataParams {
      kNChannels = 64,
      kNEvOfInt  = 21, // Number of events of interest
      kNScalers  = 16,
      kNBunches  = 10
    };

    enum EVZERORawStreamError {
      kRawDataSizeErr = 1
    };

  private:

    AliVZERORawStream(const AliVZERORawStream& stream);
    AliVZERORawStream& operator = (const AliVZERORawStream& stream);

    UInt_t GetNextWord();
    UShort_t GetNextShort();

    ULong64_t       fBBScalers[kNChannels];       // 'Beam-Beam' scalers for all channels
    ULong64_t       fBGScalers[kNChannels];       // 'Beam-Gas' scalers for all channels
    UInt_t          fScalers[kNScalers];          // Trigger scalers
    UInt_t          fBunchNumbers[kNBunches];     // Bunch numbers for the previous 10 MB events
    UShort_t        fChargeMB[kNChannels][kNBunches]; // ADC counts for all channels for the previous 10 MB events
    Bool_t          fIsIntMB[kNChannels][kNBunches];  // 'Integrator' flag for all channels for the previous 10 MB events
    Bool_t          fIsBBMB[kNChannels][kNBunches];   // 'Beam-Beam' flag for all channels for the previous 10 MB events
    Bool_t          fIsBGMB[kNChannels][kNBunches];   // 'Beam-Gas' for all channels for the previous 10 MB events

    UShort_t        fADC[kNChannels][kNEvOfInt];   // ADC counts for all channels and all events of interest
    Bool_t          fIsInt[kNChannels][kNEvOfInt]; // 'Integrator' flag for all channels
    Bool_t          fIsBB[kNChannels][kNEvOfInt];  // 'Beam-Beam' flag for all channels
    Bool_t          fIsBG[kNChannels][kNEvOfInt];  // 'Beam-Gas' flag for all channels
    Int_t           fTime[kNChannels];             // leading time for all channels
    Int_t           fWidth[kNChannels];            // signal width for all channels

    UShort_t        fTrigger;      // VZERO trigger inputs
    UShort_t        fTriggerMask;  // VZERO trigger inputs mask

    Int_t           fPosition;     // current position in the raw-data payload

    AliRawReader*   fRawReader;    // object for reading the raw data

    UChar_t*        fData;         // pointer to raw data payload

    ClassDef(AliVZERORawStream, 0) // class for reading VZERO DDL raw data
};

#endif
