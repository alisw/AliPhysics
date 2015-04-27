#ifndef ALIADRAWSTREAM_H
#define ALIADRAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class for reading the AD DDL raw data
/// The format of the raw data corresponds to the one
/// implemented in AliADBuffer class.
///
/// PLEASE NOTE that Int_t channel is here the FEE channel from 0 to 16 in the 
/// not yet defined order which may be not the same order as the order 
/// defined in aliroot by the naming and numbering conventions.
/// Therefore one has to go from FEE_Channel to AliRoot_Channel using 
/// GetOfflineChannel(Int_t channel)  when going from Online to Offline !!!
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TMath.h>

class AliRawReader;

class AliADRawStream: public TObject {
  public :
    AliADRawStream(AliRawReader* rawReader);
    virtual ~AliADRawStream();

    virtual void      Reset();
    virtual Bool_t    Next();

// Getters of various scalers and Minimum Bias flags :

   ULong64_t          GetBBScalers(Int_t channel) const  
      { return        fBBScalers[channel]; }
   ULong64_t          GetBGScalers(Int_t channel) const  
      { return        fBGScalers[channel]; }
   UInt_t             GetTriggerScalers(Int_t num_scaler) const 
      { return        fScalers[num_scaler]; }
   UInt_t             GetBunchNumbersMB(Int_t num_bunch) const 
      { return        fBunchNumbers[num_bunch]; }
   UShort_t           GetChargeMB(Int_t channel, Int_t num_bunch) const  
      { return        fChargeMB[channel][num_bunch]; } 
   Bool_t             GetIntMBFlag(Int_t channel, Int_t num_bunch) const   
      { return        fIsIntMB[channel][num_bunch]; } 
   Bool_t             GetBBMBFlag(Int_t channel, Int_t num_bunch) const   
      { return        fIsBBMB[channel][num_bunch]; }  
   Bool_t             GetBGMBFlag(Int_t channel, Int_t num_bunch) const   
      { return        fIsBGMB[channel][num_bunch]; }      
       
// Getters of ADC signals, ADC pedestals, time information and corresponding flags :

    Short_t           GetADC(Int_t channel) const
      { return TMath::MaxElement(kNEvOfInt, fADC[channel]); }    // maximum value instead of central clock
//    { return fADC[channel][kNEvOfInt/2]; }
            
    Short_t           GetPedestal(Int_t channel, Int_t event) const
      { return fADC[channel][20-event]; }
    Bool_t            GetIntegratorFlag(Int_t channel, Int_t event) const
      { return fIsInt[channel][20-event]; }
    Bool_t            GetBBFlag(Int_t channel, Int_t event) const
      { return fIsBB[channel][20-event]; } 
    Bool_t            GetBGFlag(Int_t channel, Int_t event) const
      { return fIsBG[channel][20-event]; }   
    Short_t           GetTime(Int_t channel) const
      { return fTime[channel]; }
    Short_t           GetWidth(Int_t channel) const
      { return fWidth[channel]; }

    UShort_t          GetTriggerInputs() const
      { return fTrigger; }
    UShort_t          GetTriggerInputsMask() const
      { return fTriggerMask; }

    enum EADRawDataParams {
      kNChannels = 16, // number of electronic channels in AD (FEE numbering)
      kNEvOfInt  = 21, // number of events of interest
      kNScalers  = 16, // number of scalers
      kNBunches  = 10  // number of bunches used in Minimum Bias information 
    };

    enum EADRawStreamError {
      kRawDataSizeErr = 1
    };

  private:

    AliADRawStream(const AliADRawStream& stream);
    AliADRawStream& operator = (const AliADRawStream& stream);

    UInt_t GetNextWord();
    UShort_t GetNextShort();

    ULong64_t     fBBScalers[kNChannels];        // 'Beam-Beam' scalers for all channels
    ULong64_t     fBGScalers[kNChannels];        // 'Beam-Gas' scalers for all channels
    UInt_t        fScalers[kNScalers];           // Trigger scalers
    UInt_t        fBunchNumbers[kNBunches];      // Bunch numbers for the previous 10 MB events
    UShort_t      fChargeMB[kNChannels][kNBunches]; // ADC counts for all channels for the previous 10 MB events
    Bool_t        fIsIntMB[kNChannels][kNBunches];  // 'Integrator' flag for all channels for the previous 10 MB events
    Bool_t        fIsBBMB[kNChannels][kNBunches];   // 'Beam-Beam' flag for all channels for the previous 10 MB events
    Bool_t        fIsBGMB[kNChannels][kNBunches];   // 'Beam-Gas' for all channels for the previous 10 MB events

    Short_t       fADC[kNChannels][kNEvOfInt];   // ADC counts for all channels and all events of interest
    Bool_t        fIsInt[kNChannels][kNEvOfInt]; // 'Integrator' flag for all channels 
    Bool_t        fIsBB[kNChannels][kNEvOfInt];  // 'Beam-Beam' flag for all channels
    Bool_t        fIsBG[kNChannels][kNEvOfInt];  // 'Beam-Gas' flag for all channels
    Short_t       fTime[kNChannels];             // leading time for all channels - from HPTDC - in HPTDC units
    Short_t       fWidth[kNChannels];            // pulse width for all channels - from HPTDC - in HPTDC units

    UShort_t      fTrigger;        // AD trigger inputs
    UShort_t      fTriggerMask;    // AD trigger inputs mask

    Int_t         fPosition;       // current position in the raw-data payload
    

    AliRawReader* fRawReader;      // object for reading the raw data

    UChar_t*      fData;           // pointer to raw data payload

    ClassDef(AliADRawStream, 0) // class for reading AD DDL raw data
};

#endif
