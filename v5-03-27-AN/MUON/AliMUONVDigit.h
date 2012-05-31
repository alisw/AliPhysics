#ifndef ALIMUONVDIGIT_H
#define ALIMUONVDIGIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONVDigit
/// \brief ABC of a MUON digit
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVDigit : public TObject
{
public:
  AliMUONVDigit();
  AliMUONVDigit(Int_t detElemId, Int_t manuId, Int_t manuChannel, Int_t cathode);
  virtual ~AliMUONVDigit();
  
  virtual Bool_t IsEqual(const TObject* object) const;
  /// Advertise that we can be sorted in TCollections
  virtual Bool_t IsSortable() const { return kTRUE; }
  virtual Int_t Compare(const TObject* object) const;
  
  virtual const char* GetName() const;
  
  /// The detection element this digit belongs to
  virtual Int_t DetElemId() const=0;
  /// The x-index of this digit (>=0)
  virtual Int_t PadX() const=0;
  /// The y-index of this digit (>=0)
  virtual Int_t PadY() const=0;
  /// Cathode number this digit is on (0 or 1)
  virtual Int_t Cathode() const=0;
  
  /// The charge of this digit, calibrated or not depending on IsCalibrated()
  virtual Float_t Charge() const=0;
  
  /// Raw ADC value of this digit
  virtual Int_t ADC() const = 0;
  
  /// The electronic card id this digit belongs to (manuId for tracker, localboardId for trigger)
  virtual Int_t ManuId() const = 0;
  /// The channel within ManuId() this digit belongs to (manuChannel for tracker, localBoardChannel for trigger)
  virtual Int_t ManuChannel() const=0;
  
  /// Whether the ADC has saturated
  virtual Bool_t IsSaturated() const=0;
  /// Set the saturation status
  virtual void Saturated(Bool_t saturated=kTRUE)=0;
  
  /// Whether this (simulated) digit is purely noise
  virtual Bool_t IsNoiseOnly() const=0;
  /// Set the noiseOnly status
  virtual void NoiseOnly(Bool_t /*value*/=kTRUE) { }
  
  /// Whether this (simulated) digit got corrected by chamber efficiency
  virtual Bool_t IsEfficiencyApplied() const=0;
  /// Set the efficiencyApplied status
  virtual void EfficiencyApplied(Bool_t /*value*/=kTRUE) {}
  
  /// Whether this digit has been calibrated or not (see note 1 in AliMUONVDigit.cxx)
  virtual Bool_t IsCalibrated() const=0;
  /// Set the calibrated status (see note 1 in AliMUONVDigit.cxx)
  virtual void Calibrated(Bool_t value)=0;

  /// Whether this digit has charge in femto coulomb (see note 1 in AliMUONVDigit.cxx)
  virtual Bool_t IsChargeInFC() const { return kFALSE; }
  /// Set the unit value (see note 1 in AliMUONVDigit.cxx)
  virtual void ChargeInFC(Bool_t value=kTRUE)=0;

  /// Whether or not this digit was obtained from a conversion (e.g. real to simulated)
  virtual Bool_t IsConverted() const { return kFALSE; }

  /// Whether this digit is used somewhere (typically in a cluster)
  virtual Bool_t IsUsed() const = 0;
  /// Set the used status
  virtual void Used(Bool_t value) = 0;
  
  /// A word describing the status of the neighbours of this digit
  virtual UInt_t StatusMap() const=0;
  /// Set the statusMap
  virtual void SetStatusMap(UInt_t statusMap)=0;
  
  /// Set the ADC value
  virtual void SetADC(Int_t adc)=0;
  /// Set the ix and iy of this digit
  virtual void SetPadXY(Int_t padx, Int_t pady)=0;
  /// Set the charge of this digit
  virtual void SetCharge(Float_t q)=0;
  /// Add a charge
  virtual void AddCharge(Float_t q) { SetCharge(Charge()+q); }
  
  /// Merge this with other
  virtual Bool_t MergeWith(const AliMUONVDigit& other)=0;
  
  /// Whether this digit is a tracker digit (false if belongs to trigger)
  virtual Bool_t IsTracker() const { return !IsTrigger(); }
  /** FIXME: how to get this information w/o hard-coding, yet being efficient ?
    Use one fFlags that must be set when creating the digit for instance ?
    */
  virtual Bool_t IsTrigger() const { return DetElemId()>=1100; }
  
  virtual void Print(Option_t* opt="") const;
  
  /// Below are methods only relevant for MC digigts.
  
  /// Whether we implement MC methods.
  virtual Bool_t HasMCInformation() const = 0; 
  
  /// Hit number that contributed to this simulated digit
  virtual Int_t Hit() const { return 0; }
  /// Set the hit number
  virtual void SetHit(Int_t /*n*/) { }
  /// Hit age
    virtual Float_t Time() const       {return 0;}
  /// Set hit age
      virtual void SetTime(Float_t /*t*/) { } 
  /// Number of tracks contributing to this digit
  virtual Int_t Ntracks() const { return 0; }
  /// Add a track (and its charge) to the list of tracks we handle
  virtual void AddTrack(Int_t /*trackNumber*/, Float_t /*trackCharge*/) {}
  /// Return the i-th track number
  virtual Int_t Track(Int_t /*i*/) const { return 0; }
  /// Return the i-th track charge
  virtual Float_t TrackCharge(Int_t /*i*/) const { return 0; }
  /// Patch track with a mask
  virtual void PatchTracks(Int_t /*mask*/) {}
  
  static UInt_t BuildUniqueID(Int_t detElemId, Int_t manuId, 
                              Int_t manuChannel, Int_t cathode);
  
  static void DecodeUniqueID(UInt_t uniqueID,
                             Int_t& detElemId, Int_t& manuId, 
                             Int_t& manuChannel, Int_t& cathode);
  
  static Int_t DetElemId(UInt_t uniqueID);
  static Int_t ManuId(UInt_t uniqueID);
  static Int_t ManuChannel(UInt_t uniqueID);
  static Int_t Cathode(UInt_t uniqueID);

  /// Return the localBoardNumber from the uniqueID  
  static Int_t LocalBoardNumber(UInt_t uniqueID) { return ManuId(uniqueID); }
  /// Return the localBoardChannel from the uniqueID
  static Int_t LocalBoardChannel(UInt_t uniqueID) { return ManuChannel(uniqueID); }
  
  
  ClassDef(AliMUONVDigit,0) // ABC of a MUON Digit
};

#endif
