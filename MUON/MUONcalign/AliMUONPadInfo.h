#ifndef ALIMUONPADINFO_H
#define ALIMUONPADINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calign 
/// \class AliMUONPadInfo
/// \brief Class to summarize ESD data at pad
//  Author Philippe Pillot, Subatech


#include <TObject.h>

class AliMUONPadInfo : public TObject {
public:
  AliMUONPadInfo(); // Constructor
  virtual ~AliMUONPadInfo(); //< Destructor
  AliMUONPadInfo(const AliMUONPadInfo& pad);
  AliMUONPadInfo& operator=(const AliMUONPadInfo& pad);
  
  /// Clear method (used by TClonesArray)
  virtual void Clear(Option_t* = "") {}
  
  void     Print(Option_t * option = "") const;  
  
  
  // ------ pad info ------
  /// Set pad ID
  void     SetPadId(UInt_t padId) {fPadId = padId; SetUniqueID(fPadId);}
  /// Return pad ID
  UInt_t   GetPadId() const {return fPadId;}  
  
  /// Return detection element id, part of the uniqueID
  Int_t    GetDetElemId() const   {return (fPadId & 0x00000FFF);}
  /// Return electronic card id, part of the uniqueID
  Int_t    GetManuId() const      {return (fPadId & 0x00FFF000) >> 12;}
  /// Return the channel within ManuId(), part of the uniqueID
  Int_t    GetManuChannel() const {return (fPadId & 0x3F000000) >> 24;}
  /// Return the cathode number, part of the uniqueID
  Int_t    GetCathode() const     {return (fPadId & 0x40000000) >> 30;}
  /// Return the plane type 0=Bending 1=NonBending
  Int_t    GetPadPlaneType() const     {return fPadPlaneType;}
  /// Set the plane type 0=Bending 1=NonBending
  void     SetPadPlaneType(Int_t planeType) {fPadPlaneType = planeType;}
  
  /// Set pad coordinates (cm)
  void     SetPadXY(Double_t x, Double_t y) {fPadX = x; fPadY = y;}
  /// Return pad X-position (cm)
  Double_t GetPadX() const {return fPadX;}
  /// Return pad Y-position (cm)
  Double_t GetPadY() const {return fPadY;}
  
  /// Set pad dimension (cm)
  void     SetPadDimXY(Double_t dX, Double_t dY) {fPadDimX = dX; fPadDimY = dY;}
  /// Return pad X-dimension (cm)
  Double_t GetPadDimX() const {return fPadDimX;}
  /// Return pad Y-dimension (cm)
  Double_t GetPadDimY() const {return fPadDimY;}
    
  /// Set the calibrated charge
  void     SetPadCharge(Double_t charge) {fPadCharge = charge;}
  /// Return the calibrated charge
  Double_t GetPadCharge() const {return fPadCharge;}
  
  /// Set the raw charge
  void     SetPadADC(Int_t adc) {fPadADC = adc;}
  /// Return the raw charge
  Int_t    GetPadADC() const {return fPadADC;}
  
  /// Set the pad as being calibrated or not
  void     SetCalibrated(Bool_t calibrated = kTRUE) {fPadCalibrated = calibrated;}
  /// Return kTRUE if the pad is calibrated
  Bool_t   IsCalibrated() const {return fPadCalibrated;}
  /// Set the pad as being saturated or not
  void     SetSaturated(Bool_t saturated = kTRUE) {fPadSaturated = saturated;}
  /// Return kTRUE if the pad is saturated
  Bool_t   IsSaturated() const {return fPadSaturated;}
  
  
  // ------ calibration parameters ------
  /// Set pedestal parameters
  void     SetPedestal(Float_t mean, Float_t sigma) {fPedMean = mean; fPedSigma = sigma;}
  /// Return the mean value of the pedestal
  Float_t  GetPedMean() const {return fPedMean;}
  /// Return the sigma of the pedestal
  Float_t  GetPedSigma() const {return fPedSigma;}
  
protected:
    
  // pad info
  UInt_t     fPadId;         ///< pad ID
  Int_t      fPadPlaneType;   ///< pad plane tye (0=Bending; 1=NonBending)
  Double32_t fPadX;          ///< pad X position
  Double32_t fPadY;          ///< pad Y position
  Double32_t fPadDimX;       ///< pad X dimension
  Double32_t fPadDimY;       ///< pad Y dimension
  Double32_t fPadCharge;     ///< pad calibrated charge
  Int_t      fPadADC;        ///< pad raw charge
  Bool_t     fPadSaturated;  ///< pad saturation flag
  Bool_t     fPadCalibrated; ///< pad calibration flag
  
  // calibration parameters
  Float_t    fPedMean;       ///< mean value of pedestal
  Float_t    fPedSigma;      ///< sigma of pedestal
  
  ClassDef(AliMUONPadInfo, 3) //Class to summarize ESD data at pad
};

#endif
