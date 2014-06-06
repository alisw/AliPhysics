#ifndef ALIVVVZEROFRIEND_H
#define ALIVVVZEROFRIEND_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
///
/// This is a class with access methods for VZERO DDL raw data
/// It is written to the ESD-friend file
///
///////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliVVVZEROfriend {
  public :
    AliVVVZEROfriend();
    virtual ~AliVVVZEROfriend();

    AliVVVZEROfriend(const AliVVVZEROfriend& vzerofriend);
    AliVVVZEROfriend& operator = (const AliVVVZEROfriend& vzerofriend);

    virtual void Reset();

// Getters of various scalers and Minimum Bias flags :

   virtual ULong64_t          GetBBScalers(Int_t channel) const  
      { return        fBBScalers[channel]; }
   virtual ULong64_t          GetBGScalers(Int_t channel) const  
      { return        fBGScalers[channel]; }
   virtual UInt_t             GetTriggerScalers(Int_t num_scaler) const 
      { return        fScalers[num_scaler]; }
   virtual UInt_t             GetBunchNumbersMB(Int_t num_bunch) const 
      { return        fBunchNumbers[num_bunch]; }
   virtual UShort_t           GetChargeMB(Int_t channel, Int_t num_bunch) const  
      { return        fChargeMB[channel][num_bunch]; } 
   virtual Bool_t             GetIntMBFlag(Int_t channel, Int_t num_bunch) const   
      { return        fIsIntMB[channel][num_bunch]; } 
   virtual Bool_t             GetBBMBFlag(Int_t channel, Int_t num_bunch) const   
      { return        fIsBBMB[channel][num_bunch]; }  
   virtual Bool_t             GetBGMBFlag(Int_t channel, Int_t num_bunch) const   
      { return        fIsBGMB[channel][num_bunch]; }      
       
// Getters of ADC signals, ADC pedestals, time information and corresponding flags :

    virtual Float_t           GetADC(Int_t channel) const
      { return fADC[channel][kNEvOfInt/2]; }
    virtual Float_t           GetPedestal(Int_t channel, Int_t event) const
      { return fADC[channel][event]; }
    virtual Bool_t            GetIntegratorFlag(Int_t channel, Int_t event) const
      { return fIsInt[channel][event]; }
    virtual Bool_t            GetBBFlag(Int_t channel, Int_t event) const
      { return fIsBB[channel][event]; } 
    virtual Bool_t            GetBGFlag(Int_t channel, Int_t event) const
      { return fIsBG[channel][event]; }   
    virtual Float_t            GetTime(Int_t channel) const
      { return fTime[channel]; }
    virtual Float_t            GetWidth(Int_t channel) const
      { return fWidth[channel]; }

    // Setters
    virtual void              SetBBScalers(Int_t channel, ULong64_t scalers)
      { fBBScalers[channel] = scalers; }
    virtual void              SetBGScalers(Int_t channel, ULong64_t scalers)
      { fBGScalers[channel] = scalers; }
    virtual void              SetTriggerScalers(Int_t num_scaler, UInt_t scaler)
      { fScalers[num_scaler] = scaler; }
    virtual void              SetBunchNumbersMB(Int_t num_bunch, UInt_t bunch)
      { fBunchNumbers[num_bunch] = bunch; }
    virtual void              SetChargeMB(Int_t channel,Int_t num_bunch, UShort_t charge)
      { fChargeMB[channel][num_bunch] = charge; }
    virtual void              SetIntMBFlag(Int_t channel,Int_t num_bunch, Bool_t flag)
      { fIsIntMB[channel][num_bunch] = flag; }
    virtual void              SetBBMBFlag(Int_t channel,Int_t num_bunch, Bool_t flag)
      { fIsBBMB[channel][num_bunch] = flag; }
    virtual void              SetBGMBFlag(Int_t channel,Int_t num_bunch, Bool_t flag)
      { fIsBGMB[channel][num_bunch] = flag; }

    virtual void              SetPedestal(Int_t channel, Int_t event, Float_t adc)
      { fADC[channel][event] = adc; }
    virtual void              SetIntegratorFlag(Int_t channel, Int_t event, Bool_t flag)
      { fIsInt[channel][event] = flag; }
    virtual void              SetBBFlag(Int_t channel, Int_t event, Bool_t flag)
      { fIsBB[channel][event] = flag; }
    virtual void              SetBGFlag(Int_t channel, Int_t event, Bool_t flag)
      { fIsBG[channel][event] = flag; }
    virtual void              SetTime(Int_t channel, Float_t time)
      { fTime[channel] = time; }
    virtual void              SetWidth(Int_t channel, Float_t width)
      { fWidth[channel] = width; }

    virtual UShort_t          GetTriggerInputs() const
      { return fTrigger; }
    virtual UShort_t          GetTriggerInputsMask() const
      { return fTriggerMask; }
    virtual void              SetTriggerInputs(UShort_t inputs)
      { fTrigger = inputs; }
    virtual void              SetTriggerInputsMask(UShort_t mask)
      { fTriggerMask = mask; }

    enum EESDVZEROfriendParams {
      kNChannels = 64, // number of electronic channels in V0 (FEE numbering)
      kNEvOfInt  = 21, // number of events of interest
      kNScalers  = 16, // number of scalers
      kNBunches  = 10  // number of bunches used in Minimum Bias information 
    };
};

#endif
