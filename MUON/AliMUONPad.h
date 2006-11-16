#ifndef ALIMUONPAD_H
#define ALIMUONPAD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONPad
/// \brief Combination of digit and mppad informations.
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif
#ifndef ALI_MP_AREA_H
#  include "AliMpArea.h"
#endif

class AliMUONPad : public TObject
{
public:
  AliMUONPad();
  AliMUONPad(Double_t x, Double_t y,
             Double_t dx, Double_t dy, Double_t charge);
  AliMUONPad(const TVector2& position, const TVector2& dimensions,
             Double_t charge);
  AliMUONPad(Int_t detElemId, Int_t cathode,
             Int_t ix, Int_t iy, Double_t x, Double_t y,
             Double_t dx, Double_t dy, Double_t charge);
  virtual ~AliMUONPad() {}

  void BackupCharge() { fChargeBackup = fCharge; }
  
  Int_t Cathode() const { return fCathode; }

  Double_t Coord(Int_t ixy) const;
  
  Double_t Charge() const { return fCharge; }
  Double_t ChargeBackup() const { return fChargeBackup; }

  virtual Int_t Compare(const TObject* obj) const;

  Int_t DetElemId() const { return fDetElemId; }
  
  Int_t DigitIndex() const { return fDigitIndex; }

  TVector2 Dimensions() const { return fDimensions; }
  
  Double_t DX() const { return fDimensions.X(); }
  Double_t DY() const { return fDimensions.Y(); }

  Bool_t IsReal() const { return fIsReal; }

  Bool_t IsSaturated() const { return fIsSaturated; }
  
  Bool_t IsSortable() const { return kTRUE; }
  
  Bool_t IsUsed() const { return fClusterId >= 0; }

  Int_t Ix() const { return fIx; }
  Int_t Iy() const { return fIy; }

  virtual void Paint(Option_t* opt="");

  TVector2 Position() const { return fPosition; }

  void Print(Option_t* opt = "") const;

  // Detach this pad from a cluster
  void Release() { fClusterId = -1; }
  
  void RevertCharge() { fCharge = fChargeBackup; }

  void SetCharge(Double_t charge) { fCharge = charge; }
  
  void SetChargeBackup(Double_t charge) { fChargeBackup = charge; }

  void SetCoord(Int_t ixy, Double_t Coord);

  void SetDigitIndex(Int_t digitIndex) { fDigitIndex = digitIndex; }

  void SetStatus(Int_t status) { fStatus = status; }
    
  void SetClusterId(Int_t id) { fClusterId = id; }  
  
  void SetSaturated(Bool_t val) { fIsSaturated = val; }
  
  void SetSize(Int_t ixy, Double_t Size);
  
  void SetReal(Bool_t val) { fIsReal = val; }

  void Shift(Int_t ixy, Double_t shift);
  
  Double_t Size(Int_t ixy) const;

  Int_t Status() const { return fStatus; }
  
  Double_t X() const { return fPosition.X(); }
  Double_t Y() const { return fPosition.Y(); }
  
  static AliMpArea Overlap(const AliMUONPad& d1, const AliMUONPad& d2);
  
  static Bool_t AreOverlapping(const AliMUONPad& d1, const AliMUONPad& d2,
                               const TVector2& precision,
                               AliMpArea& overlapArea);

  static Bool_t AreOverlapping(const AliMUONPad& d1, const AliMUONPad& d2,
                               const TVector2& precision);
  
  static Bool_t AreNeighbours(const AliMUONPad& d1, const AliMUONPad& d2);

  
private:
    void Init(Int_t detElemId, Int_t cathode,
              Int_t ix, Int_t iy,
              const TVector2& position,
              const TVector2& dimensions,
              Double_t charge);
  
private:
  Bool_t fIsSaturated; /// whether this pad is saturated or not
  Bool_t fIsReal; /// whether this is a real pad or a virtual one
  Int_t fClusterId; /// cluster id this pad belongs to (-1 if not attached to a cluster)
  Int_t fCathode; /// cathode number
  Int_t fDetElemId; /// detection element id
  Int_t fDigitIndex; /// corresponding digit index (to get back the digit from which this pad was constructed, if needed)
  Int_t fIx; /// x-index
  Int_t fIy; /// y-index
  Int_t fStatus; /// status word
  TVector2 fDimensions; /// half dimensions in x and y (cm)
  TVector2 fPosition; /// positions in x and y (cm)
  Double_t fCharge; /// pad charge
  Double_t fChargeBackup; /// backup charge (usefull if clustering somehow plays with the charge, this one is the "original" one)
  
  ClassDef(AliMUONPad,1) // A full grown pad 
};

#endif
