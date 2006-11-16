#ifndef ALIMUONCLUSTER_H
#define ALIMUONCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONCluster
/// \brief A group of adjacent pads
/// 
// Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif
#ifndef ALI_MP_AREA_H
#  include "AliMpArea.h"
#endif
#ifndef ALI_MP_DIRECTION_H
#  include "AliMpDirection.h"
#endif
#ifndef ALI_MP_INT_PAIR_H
#  include "AliMpIntPair.h"
#endif

class AliMUONPad;
class TObjArray;

class AliMUONCluster : public TObject
{
public:
  AliMUONCluster();
  AliMUONCluster(const AliMUONCluster& rhs);
  AliMUONCluster& operator=(const AliMUONCluster& rhs);
  
  virtual ~AliMUONCluster();

  void AddPad(const AliMUONPad& pad);

  /// Area that contains all the pads.
  AliMpArea Area() const;
    
  Float_t Charge() const;
  Float_t Charge(Int_t cathode) const;

  /// Return the cathode's charges asymmetry
  Float_t ChargeAsymmetry() const;

  Float_t Chi2() const { return fChi2; }

  virtual void Copy(TObject& obj) const;

  Bool_t HasPosition() const { return fHasPosition; }

  /// Whether we have at least one saturated pad in a given cathode 
  Bool_t IsSaturated(Int_t cathode) const { return fIsSaturated[cathode]; }
  
  /// Whether we have one saturated pad on *each* cathode
  Bool_t IsSaturated() const { return IsSaturated(0) && IsSaturated(1); }
  
  Int_t MaxChargeCathode() const { return Charge(0) > Charge(1) ? 0:1; }

  Int_t MaxRawChargeCathode() const { return RawCharge(0) > RawCharge(1) ? 0:1; }

  /// Return the smallest pad dimensions for a given cathode
  TVector2 MinPadDimensions(Int_t cathode, Int_t statusMask, Bool_t matchMask) const;
  
  /// Return the smallest pad dimensions
  TVector2 MinPadDimensions( Int_t statusMask, Bool_t matchMask) const;
  
  Int_t Multiplicity() const;
  Int_t Multiplicity(Int_t cathode) const;

  /// Compute number of pads in X and Y direction for a given cathode.  
  AliMpIntPair NofPads(Int_t cathode, Int_t statusMask, Bool_t matchMask) const;
  
  /// Number of pads in (X,Y) direction, whatever the cathode.
  AliMpIntPair NofPads(Int_t statusMask, Bool_t matchMask=kTRUE) const;
  
  AliMUONPad* Pad(Int_t index) const;
  
  virtual void Paint(Option_t* opt="");
  
  TVector2 Position() const { return fPosition; }
  TVector2 PositionError() const { return fPositionError; }

  virtual void Print(Option_t* opt="") const;
      
  /// By default, return the average of both cathode RawCharges.
  Float_t RawCharge() const;
  
  /// Returns the RawCharge on the given cathode.
  Float_t RawCharge(Int_t cathode) const;
    
  /// Return the cathode's raw charges asymmetry
  Float_t RawChargeAsymmetry() const;
  
  void RemovePad(AliMUONPad* pad);
  
  void SetCharge(Float_t chargeCath0, Float_t chargeCath1)
  { fHasCharge = kTRUE; fCharge[0]=chargeCath0; fCharge[1]=chargeCath1; }
  
  void SetChi2(Float_t chi2) { fChi2 = chi2; }

  void SetPosition(const TVector2& pos, const TVector2& errorOnPos) 
  { fHasPosition = kTRUE; fPosition = pos; fPositionError = errorOnPos; }
  
  void Sort();
  
private:
  TObjArray* fPads; ///< AliMUONPad(s) composing this cluster
  Bool_t fHasPosition; ///< false for pre-cluster (i.e. not yet computed)
  TVector2 fPosition; ///< (x,y) of that cluster (only valid if fHasPosition is kTRUE)
  TVector2 fPositionError; ///< errors on (x,y)
  Int_t fMultiplicity[2]; ///< number of pads in each cathode
  Float_t fRawCharge[2]; ///< cathode RawCharges
  Bool_t fHasCharge; ///< false if SetCharge has not been called
  Float_t fCharge[2]; ///< cathode (re)computed charges
  Float_t fChi2; ///< chi2 of the RawCharge fit (if any)
  Bool_t fIsSaturated[2]; ///< saturation status of cathodes
  Bool_t fIsSorted; ///< whether pads are sorted or not
  ClassDef(AliMUONCluster,1) // A cluster of AliMUONPad
};

#endif
