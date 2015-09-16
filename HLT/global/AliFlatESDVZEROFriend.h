#ifndef ALIFLATESDVZEROFRIEND_H
#define ALIFLATESDVZEROFRIEND_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/**
 * Flat structure representing a ESD VZEROFriend friend 
 */

#include "Rtypes.h"
#include "AliVMisc.h"
#include "AliVVZEROfriend.h"
#include "AliESDVZEROfriend.h"

class AliFlatESDVZEROFriend: public AliVVZEROfriend 
{
 public:
  // -- Constructor / Destructors
  
  AliFlatESDVZEROFriend();
 ~AliFlatESDVZEROFriend() {}

  // constructor and method for reinitialisation of virtual table
  AliFlatESDVZEROFriend( AliVConstructorReinitialisationFlag );
  void Reinitialize() { new (this) AliFlatESDVZEROFriend ( AliVReinitialize ); } 

  void SetFromESDVZEROfriend(const AliESDVZEROfriend &v );
  void GetESDVZEROfriend( AliESDVZEROfriend &v ) const;

  static size_t GetSize() { return sizeof(AliFlatESDVZEROFriend); }

  void Reset();

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

    Float_t           GetADC(Int_t channel) const
      { return fADC[channel][kNEvOfInt/2]; }
    Float_t           GetPedestal(Int_t channel, Int_t event) const
      { return fADC[channel][event]; }
    Bool_t            GetIntegratorFlag(Int_t channel, Int_t event) const
      { return fIsInt[channel][event]; }
    Bool_t            GetBBFlag(Int_t channel, Int_t event) const
      { return fIsBB[channel][event]; } 
    Bool_t            GetBGFlag(Int_t channel, Int_t event) const
      { return fIsBG[channel][event]; }   
    Float_t            GetTime(Int_t channel) const
      { return fTime[channel]; }
    Float_t            GetWidth(Int_t channel) const
      { return fWidth[channel]; }

    // Setters

    void              SetBBScalers(Int_t channel, ULong64_t scalers)
      { fBBScalers[channel] = scalers; }
    void              SetBGScalers(Int_t channel, ULong64_t scalers)
      { fBGScalers[channel] = scalers; }
    void              SetTriggerScalers(Int_t num_scaler, UInt_t scaler)
      { fScalers[num_scaler] = scaler; }
    void              SetBunchNumbersMB(Int_t num_bunch, UInt_t bunch)
      { fBunchNumbers[num_bunch] = bunch; }
    void              SetChargeMB(Int_t channel,Int_t num_bunch, UShort_t charge)
      { fChargeMB[channel][num_bunch] = charge; }
    void              SetIntMBFlag(Int_t channel,Int_t num_bunch, Bool_t flag)
      { fIsIntMB[channel][num_bunch] = flag; }
    void              SetBBMBFlag(Int_t channel,Int_t num_bunch, Bool_t flag)
      { fIsBBMB[channel][num_bunch] = flag; }
    void              SetBGMBFlag(Int_t channel,Int_t num_bunch, Bool_t flag)
      { fIsBGMB[channel][num_bunch] = flag; }

    void              SetPedestal(Int_t channel, Int_t event, Float_t adc)
      { fADC[channel][event] = adc; }
    void              SetIntegratorFlag(Int_t channel, Int_t event, Bool_t flag)
      { fIsInt[channel][event] = flag; }
    void              SetBBFlag(Int_t channel, Int_t event, Bool_t flag)
      { fIsBB[channel][event] = flag; }
    void              SetBGFlag(Int_t channel, Int_t event, Bool_t flag)
      { fIsBG[channel][event] = flag; }
    void              SetTime(Int_t channel, Float_t time)
      { fTime[channel] = time; }
    void              SetWidth(Int_t channel, Float_t width)
      { fWidth[channel] = width; }

    UShort_t          GetTriggerInputs() const
      { return fTrigger; }
    UShort_t          GetTriggerInputsMask() const
      { return fTriggerMask; }
    void              SetTriggerInputs(UShort_t inputs)
      { fTrigger = inputs; }
    void              SetTriggerInputsMask(UShort_t mask)
      { fTriggerMask = mask; }

private: 

  ULong64_t     fBBScalers[kNChannels];        // 'Beam-Beam' scalers for all channels
  ULong64_t     fBGScalers[kNChannels];        // 'Beam-Gas' scalers for all channels
  UInt_t        fScalers[kNScalers];           // Trigger scalers
  UInt_t        fBunchNumbers[kNBunches];      // Bunch numbers for the previous 10 MB events
  UShort_t      fChargeMB[kNChannels][kNBunches]; // ADC counts for all channels for the previous 10 MB events
  Bool_t        fIsIntMB[kNChannels][kNBunches];  // 'Integrator' flag for all channels for the previous 10 MB events
  Bool_t        fIsBBMB[kNChannels][kNBunches];   // 'Beam-Beam' flag for all channels for the previous 10 MB events
  Bool_t        fIsBGMB[kNChannels][kNBunches];   // 'Beam-Gas' for all channels for the previous 10 MB events
  
  Float_t       fADC[kNChannels][kNEvOfInt];   // ADC counts for all channels and all events of interest
  Bool_t        fIsInt[kNChannels][kNEvOfInt]; // 'Integrator' flag for all channels 
  Bool_t        fIsBB[kNChannels][kNEvOfInt];  // 'Beam-Beam' flag for all channels
  Bool_t        fIsBG[kNChannels][kNEvOfInt];  // 'Beam-Gas' flag for all channels
  Float_t       fTime[kNChannels];             // leading time for all channels - from HPTDC - in nanoseconds
  Float_t       fWidth[kNChannels];            // pulse width for all channels - from HPTDC - in nanoseconds

  UShort_t      fTrigger;        // VZERO trigger inputs
  UShort_t      fTriggerMask;    // VZERO trigger inputs mask
};

inline AliFlatESDVZEROFriend::AliFlatESDVZEROFriend() :
  AliVVZEROfriend(),
  fTrigger(0),
  fTriggerMask(0)
{
  // default constructor
  Reset();
}

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatESDVZEROFriend::AliFlatESDVZEROFriend( AliVConstructorReinitialisationFlag ){}  // only sets virtual table pointer
#pragma GCC diagnostic warning "-Weffc++" 

inline void AliFlatESDVZEROFriend::Reset()
{
  // Reset the contents of the object
  fTrigger = 0;
  fTriggerMask = 0;

  for (Int_t iScaler = 0; iScaler < kNScalers; iScaler++)
    fScalers[iScaler] = 0;

  for (Int_t iBunch = 0; iBunch < kNBunches; iBunch++)
    fBunchNumbers[iBunch] = 0;

  for (Int_t iChannel = 0; iChannel < kNChannels; iChannel++) {
    fBBScalers[iChannel] = 0;
    fBGScalers[iChannel] = 0;
    for (Int_t iBunch = 0; iBunch < kNBunches; iBunch++) {
      fChargeMB[iChannel][iBunch] = 0;
      fIsIntMB[iChannel][iBunch]  = kFALSE;
      fIsBBMB[iChannel][iBunch]   = kFALSE;
      fIsBGMB[iChannel][iBunch]   = kFALSE;
    }
    for (Int_t iEv = 0; iEv < kNEvOfInt; iEv++) {
      fADC[iChannel][iEv]   = 0.0;
      fIsInt[iChannel][iEv] = kFALSE;
      fIsBB[iChannel][iEv]  = kFALSE;
      fIsBG[iChannel][iEv]  = kFALSE;
    }
    fTime[iChannel]  = 0.0;
    fWidth[iChannel] = 0.0;
  }  
}

inline void AliFlatESDVZEROFriend::SetFromESDVZEROfriend(const AliESDVZEROfriend &v )
{
  // set from ESD VZERO friend object

  fTrigger = v.GetTriggerInputs();
  fTriggerMask = v.GetTriggerInputsMask();

  for (Int_t iScaler = 0; iScaler < kNScalers; iScaler++)
    fScalers[iScaler] = v.GetTriggerScalers(iScaler);

  for (Int_t iBunch = 0; iBunch < kNBunches; iBunch++)
    fBunchNumbers[iBunch] = v.GetBunchNumbersMB(iBunch);

  for (Int_t iChannel = 0; iChannel < kNChannels; iChannel++) {
    fBBScalers[iChannel] = v.GetBBScalers(iChannel);
    fBGScalers[iChannel] = v.GetBGScalers(iChannel);
    for (Int_t iBunch = 0; iBunch < kNBunches; iBunch++) {
      fChargeMB[iChannel][iBunch] = v.GetChargeMB(iChannel, iBunch);
      fIsIntMB[iChannel][iBunch]  = v.GetIntMBFlag(iChannel, iBunch);
      fIsBBMB[iChannel][iBunch]   = v.GetBBMBFlag(iChannel, iBunch);
      fIsBGMB[iChannel][iBunch]   = v.GetBGMBFlag(iChannel, iBunch);
    }
    for (Int_t iEv = 0; iEv < kNEvOfInt; iEv++) {
      fADC[iChannel][iEv]   = v.GetPedestal(iChannel, iEv);
      fIsInt[iChannel][iEv] = v.GetIntegratorFlag(iChannel, iEv);
      fIsBB[iChannel][iEv]  = v.GetBBFlag(iChannel, iEv);
      fIsBG[iChannel][iEv]  = v.GetBBFlag(iChannel, iEv);
    }
    fTime[iChannel]  = v.GetTime(iChannel);
    fWidth[iChannel] = v.GetWidth(iChannel);
  }  
}

inline void AliFlatESDVZEROFriend::GetESDVZEROfriend( AliESDVZEROfriend &v ) const
{
  // copy content to ESD VZEROFriend object

  v.SetTriggerInputs( fTrigger );
  v.SetTriggerInputsMask( fTriggerMask );

  for (Int_t iScaler = 0; iScaler < kNScalers; iScaler++)
    v.SetTriggerScalers( iScaler, fScalers[iScaler] );

  for (Int_t iBunch = 0; iBunch < kNBunches; iBunch++)
    v.SetBunchNumbersMB( iBunch, fBunchNumbers[iBunch] );

  for (Int_t iChannel = 0; iChannel < kNChannels; iChannel++) {
    v.SetBBScalers(iChannel, fBBScalers[iChannel] );
    v.SetBGScalers(iChannel, fBGScalers[iChannel] );
    for (Int_t iBunch = 0; iBunch < kNBunches; iBunch++) {
      v.SetChargeMB( iChannel, iBunch, fChargeMB[iChannel][iBunch] );
      v.SetIntMBFlag( iChannel, iBunch, fIsIntMB[iChannel][iBunch] );
      v.SetBBMBFlag( iChannel, iBunch, fIsBBMB[iChannel][iBunch] );
      v.SetBGMBFlag( iChannel, iBunch, fIsBGMB[iChannel][iBunch] );
    }
    for (Int_t iEv = 0; iEv < kNEvOfInt; iEv++) {
      v.SetPedestal( iChannel, iEv, fADC[iChannel][iEv] );
      v.SetIntMBFlag( iChannel, iEv, fIsInt[iChannel][iEv] );
      v.SetBBMBFlag( iChannel, iEv, fIsBB[iChannel][iEv] );
      v.SetBGMBFlag( iChannel, iEv, fIsBG[iChannel][iEv] );
    }
    v.SetTime( iChannel, fTime[iChannel] );
    v.SetWidth( iChannel, fWidth[iChannel] );
  }
}

#endif
