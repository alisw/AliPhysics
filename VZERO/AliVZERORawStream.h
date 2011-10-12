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
/// PLEASE NOTE that Int_t channel is here the FEE channel from 0 to 63 in the 
/// order defined by Yannick Zoccarato which is not the same order as the order 
/// defined in aliroot by the naming and numbering conventions.
/// Therefore one has to go from FEE_Channel to AliRoot_Channel using 
/// GetOfflineChannel(Int_t channel)  when going from Online to Offline !!!
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TMath.h>

class AliRawReader;
class AliVZEROTriggerData;

class AliVZERORawStream: public TObject {
  public :
    AliVZERORawStream(AliRawReader* rawReader);
    virtual ~AliVZERORawStream();

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
      { return fADC[channel][event]; }
    Bool_t            GetIntegratorFlag(Int_t channel, Int_t event) const
      { return fIsInt[channel][event]; }
    Bool_t            GetBBFlag(Int_t channel, Int_t event) const
      { return fIsBB[channel][event]; } 
    Bool_t            GetBGFlag(Int_t channel, Int_t event) const
      { return fIsBG[channel][event]; }   
    Short_t           GetTime(Int_t channel) const
      { return fTime[channel]; }
    Short_t           GetWidth(Int_t channel) const
      { return fWidth[channel]; }

    UShort_t          GetTriggerInputs() const
      { return fTrigger; }
    UShort_t          GetTriggerInputsMask() const
      { return fTriggerMask; }

    void              CalculateChargeForCentrTriggers(AliVZEROTriggerData *triggerData,
						      UShort_t &chargeA, UShort_t &chargeC) const;

// Getter of Offline Channel number as used in aliroot (defined by aliroot 
// numbering convention) from FEE channel (electronic channel number given 
// by the V0 electronics readout) - See comment above - 

    Int_t              GetOfflineChannel(Int_t channel)  const
      { Int_t  fOfflineChannel[64] = {39, 38, 37, 36, 35, 34, 33, 32, 
                                      47, 46, 45, 44, 43, 42, 41, 40, 
			              55, 54, 53, 52, 51, 50, 49, 48, 
			              63, 62, 61, 60, 59, 58, 57, 56,
			               7,  6,  5,  4,  3,  2,  1,  0, 
			              15, 14, 13, 12, 11, 10,  9,  8,
			              23, 22, 21, 20, 19, 18, 17, 16, 
			              31, 30, 29, 28, 27, 26, 25, 24};
               return fOfflineChannel[channel]; }	

    enum EVZERORawDataParams {
      kNChannels = 64, // number of electronic channels in V0 (FEE numbering)
      kNEvOfInt  = 21, // number of events of interest
      kNScalers  = 16, // number of scalers
      kNBunches  = 10  // number of bunches used in Minimum Bias information 
    };

    enum EVZERORawStreamError {
      kRawDataSizeErr = 1
    };

  private:

    AliVZERORawStream(const AliVZERORawStream& stream);
    AliVZERORawStream& operator = (const AliVZERORawStream& stream);

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

    UShort_t      fTrigger;        // VZERO trigger inputs
    UShort_t      fTriggerMask;    // VZERO trigger inputs mask

    Int_t         fPosition;       // current position in the raw-data payload
    

    AliRawReader* fRawReader;      // object for reading the raw data

    UChar_t*      fData;           // pointer to raw data payload

    ClassDef(AliVZERORawStream, 0) // class for reading VZERO DDL raw data
};

#endif
