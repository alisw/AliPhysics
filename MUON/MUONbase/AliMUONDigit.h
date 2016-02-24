#ifndef ALIMUONDIGIT_H
#define ALIMUONDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONDigit
/// \brief MUON digit

#ifndef ALIMUONVDIGIT_H
#  include "AliMUONVDigit.h"
#endif

class AliMUONDigit : public AliMUONVDigit 
{
 public:
  AliMUONDigit();
  AliMUONDigit(Int_t detElemId, Int_t manuId, Int_t manuChannel, Int_t cathode);
  
    AliMUONDigit(const AliMUONDigit& rhs);
    virtual ~AliMUONDigit();

    AliMUONDigit& operator=(const AliMUONDigit& rhs);
    
    /// Own clone methods (as the default TObject::Clone turned out to be pretty slow !)
    virtual TObject* Clone(const char* /*newname*/ = "") const { return new AliMUONDigit(*this); }
    
    virtual Bool_t HasMCInformation() const { return kTRUE; }
    
    virtual Int_t DetElemId()const     {return fDetElemId;}  ///< Return detection element ID  
    virtual Int_t PadX() const         {return fPadX;}       ///< Return pad number along x
    virtual Int_t PadY() const         {return fPadY;}       ///< Return pad number along y
    virtual Int_t Cathode() const      {return fCathode;}    ///< Return cathode number
    
    virtual Float_t Charge() const     {return fSignal;}     ///< Return signal amplitude
    
    virtual Int_t Hit() const          {return fHit;}        ///< Return MC hit number
    
    virtual Float_t Time() const       {return fTime;}       /// Return MC hit age

    virtual Int_t Ntracks() const { return fNtracks; }       ///< Return MC tracks making to this digit
    virtual void AddTrack(Int_t trackNumber, Float_t trackCharge);
    virtual Int_t Track(Int_t i) const;
    virtual Float_t TrackCharge(Int_t i) const;
    
    virtual Int_t ADC() const { return fADC; }                 ///< Return ADC value
    virtual Int_t ManuId() const { return fManuId; }           ///< Return Id of the MANU chip
    virtual Int_t ManuChannel() const { return fManuChannel; } ///< Return Channel within the MANU chip
    virtual Bool_t IsSaturated() const;
    virtual Bool_t IsNoiseOnly() const;
    virtual Bool_t IsEfficiencyApplied() const;
    virtual Bool_t IsConverted() const;
    virtual Bool_t IsChargeInFC() const;
    virtual UInt_t StatusMap() const { return fStatusMap; }    ///< Return Neighbouring pad status
    
    virtual void NoiseOnly(Bool_t value=kTRUE);
    virtual void Saturated(Bool_t saturated=kTRUE);
    virtual void EfficiencyApplied(Bool_t value=kTRUE);
    virtual void Converted(Bool_t value=kTRUE);
    virtual void ChargeInFC(Bool_t value=kTRUE);
  
    virtual void SetADC(Int_t adc)         {fADC=adc; }        ///< Set ADC value
    virtual void SetPadXY(Int_t padx, Int_t pady)        {fPadX = padx; fPadY=pady; }      ///< Set pad number along x
    virtual void SetCharge(Float_t q)        {fSignal = q;}    ///< Set charge
    virtual void SetHit(Int_t n)           {fHit = n;}         ///< Set MC hit number
    virtual void SetTime(Float_t t) {fTime = t;}               ///< Set MC hit age
    virtual void SetStatusMap(UInt_t statusMap) { fStatusMap = statusMap; } ///< Set status map
    
    virtual void Copy(TObject& digit) const;
    
    /** Delete the internal track arrays (which are dynamically allocated).
      * This is to insure we can put those digits in e.g. TClonesArray
      * w/o leaking memory.
      */
    virtual void Clear(Option_t*);
    
    // Add mask to the track numbers.
    virtual void PatchTracks(Int_t mask);
    
    virtual Bool_t MergeWith(const AliMUONVDigit& other);

    virtual Bool_t IsUsed() const;
    virtual void Used(Bool_t value);
    
    virtual Bool_t IsCalibrated() const;
    virtual void Calibrated(Bool_t value);
    
    virtual UInt_t GetUniqueID() const;
    
private:
    Int_t fDetElemId;     ///< Detection element ID
    Int_t fManuId;        ///< Id of the MANU chip.
    Int_t fManuChannel;   ///< Channel within the MANU chip.
    Float_t fSignal;        ///< Signal amplitude    
      
    Int_t fPadX;          ///< Pad number along x
    Int_t fPadY;          ///< Pad number along y
    Int_t fCathode;       ///< Cathode number
    Int_t fADC;           ///< ADC value
    UInt_t fFlags;        ///< Special flags (e.g. is the signal an overflow ?)
    
    Int_t fNtracks;       ///< MC tracks making to this digit.
    
    /// charges of MC track making this digit
    Float_t* fTcharges;     //[fNtracks]  charges of MC track making this digit

    /// primary MC tracks making this digit
    Int_t* fTracks;       //[fNtracks]  primary MC tracks making this digit

    Int_t fHit;           ///< MC hit number - temporary solution
    Float_t fTime;        ///< MC hit age
  
    UInt_t fStatusMap; ///< Neighbouring pad status (whether ped, lv, hv were ok or not)
    
    static const UInt_t fgkSaturatedMask = 0x1; ///< the mask (part of fFlags) to indicate this digit is saturated
    static const UInt_t fgkUsedMask = 0x10; ///< whether this digit is used by whatever other object (typically a cluster, though)
    static const UInt_t fgkCalibratedMask = 0x100; ///< whether this digits has been calibrated
    static const UInt_t fgkNoiseOnlyMask = 0x1000; ///< indicate a simulated digit due to noise only
    static const UInt_t fgkEfficiencyMask = 0x2000; ///< indicate chamber efficiency has been applied to a simulated digit
    static const UInt_t fgkConverted       = 0x4000; ///< has been converted from a real digit
    static const UInt_t fgkChargeInFC      = 0x8000; ///< charge unit are femto coulomb
  
    ClassDef(AliMUONDigit,13)  //Digits for MUON
};

#endif
