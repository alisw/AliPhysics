#ifndef ALIMUONDIGIT_H
#define ALIMUONDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONDigit
/// \brief MUON digit

#include <TObject.h>

class AliMUONDigit : public TObject 
{
 public:
    AliMUONDigit();
    AliMUONDigit(const AliMUONDigit& rhs);
    virtual ~AliMUONDigit();

    AliMUONDigit& operator=(const AliMUONDigit& rhs);
    
    virtual Bool_t IsSortable() const {return kTRUE;}        ///< Return true if sortable
    virtual Int_t Compare(const TObject *obj) const;

    virtual Int_t DetElemId()const     {return fDetElemId;}  ///< Return detection element ID  
    virtual Int_t PadX() const         {return fPadX;}       ///< Return pad number along x
    virtual Int_t PadY() const         {return fPadY;}       ///< Return pad number along y
    virtual Int_t Cathode() const      {return fCathode;}    ///< Return cathode number
    
    virtual Float_t Signal() const       {return fSignal;}     ///< Return signal amplitude
    
    virtual Float_t Physics() const      {return fPhysics;}    ///< Return MC physics contribution to signal
    
    virtual Int_t Hit() const          {return fHit;}        ///< Return MC hit number
    
    virtual Int_t Ntracks() const { return fNtracks; }       ///< Return MC tracks making to this digit
    virtual void AddTrack(Int_t trackNumber, Float_t trackCharge);
    virtual Int_t Track(Int_t i) const;
    virtual Float_t TrackCharge(Int_t i) const;
    
    virtual Int_t ADC() const { return fADC; }                 ///< Return ADC value
    virtual Int_t ManuId() const { return fManuId; }           ///< Return Id of the MANU chip
    virtual Int_t ManuChannel() const { return fManuChannel; } ///< Return Channel within the MANU chip
    virtual Bool_t IsSaturated() const;
    virtual Bool_t IsNoiseOnly() const;
    virtual UInt_t StatusMap() const { return fStatusMap; }
    
    virtual void NoiseOnly(Bool_t value=kTRUE);
    virtual void Saturated(Bool_t saturated=kTRUE);
    virtual void SetElectronics(Int_t manuId, Int_t manuChannel);
    virtual void SetADC(Int_t adc) { fADC=adc; }
    virtual void SetDetElemId(Int_t id)    {fDetElemId = id;}  ///< Set detection element ID
    virtual void SetPadX(Int_t pad)        {fPadX = pad;}      ///< Set pad number along x
    virtual void SetPadY(Int_t pad)        {fPadY = pad;}      ///< Set pad number along y
    virtual void SetSignal(Float_t q)        {fSignal = q;}      ///< Set signal amplitude
    virtual void AddSignal(Float_t q)        {fSignal += q;}     ///< Add signal amplitude
    virtual void AddPhysicsSignal(Float_t q) {fPhysics += q;}    ///< Add MC physics contribution to signal
    virtual void SetHit(Int_t n)           {fHit = n;}         ///< Set MC hit number
    virtual void SetCathode(Int_t c)       {fCathode = c;}     ///< Set cathode number
    virtual void SetPhysicsSignal(Float_t q) {fPhysics = q; }    ///< Set MC physics contribution to signal
    virtual void SetStatusMap(UInt_t statusMap) { fStatusMap = statusMap; }
    
    virtual void Print(Option_t* opt="") const;
    
    virtual void Copy(TObject& digit) const;
    
    /** Delete the internal track arrays (which are dynamically allocated).
      * This is to insure we can put those digits in e.g. TClonesArray
      * w/o leaking memory.
      */
    virtual void Clear(Option_t*);
    
    // Add mask to the track numbers.
    virtual void PatchTracks(Int_t mask);
    
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

    Float_t fPhysics;       ///< MC physics contribution to signal 
    Int_t fHit;           ///< MC hit number - temporary solution
  
    UInt_t fStatusMap; ///< Neighbouring pad status (whether ped, gains, hv were ok or not)
    
    static const UInt_t fgkSaturatedMask = 0x1; ///< the mask (part of fFlags) to indicate this digit is saturated
    static const UInt_t fgkNoiseOnlyMask = 0x1000; ///< indicate a simulated digit due to noise only
    
    ClassDef(AliMUONDigit,6)  //Digits for MUON
};
#endif
