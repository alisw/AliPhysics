#ifndef ALIMUONRAWCLUSTERV2_H
#define ALIMUONRAWCLUSTERV2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONRawClusterV2
/// \brief MUON raw cluster
///
//  Author Philippe Pillot, Subatech

#include "AliMUONVCluster.h"
#include <TMath.h>

class AliMUONRawClusterV2 : public AliMUONVCluster {

 public:
  AliMUONRawClusterV2();
  AliMUONRawClusterV2(Int_t chamberId, Int_t detElemId, Int_t clusterIndex);
  virtual ~AliMUONRawClusterV2();
  AliMUONRawClusterV2(const AliMUONRawClusterV2& cluster);
  AliMUONRawClusterV2 & operator=(const AliMUONRawClusterV2& cluster);
  
  virtual void Clear(Option_t* = "");
  
	   /// Create a copy of the current cluster
  virtual AliMUONRawClusterV2* Clone(const char* = "") const {return new AliMUONRawClusterV2(*this);}
  
           /// Set coordinates (cm)
  virtual void     SetXYZ(Double_t x, Double_t y, Double_t z) {fX = x; fY = y; fZ = z;}
	   /// Return coordinate X (cm)
  virtual Double_t GetX() const {return fX;}
	   /// Return coordinate Y (cm)
  virtual Double_t GetY() const {return fY;}
	   /// Return coordinate Z (cm)
  virtual Double_t GetZ() const {return fZ;}
  
	   /// Set resolution (cm) on coordinates (X,Y)
  virtual void     SetErrXY(Double_t errX, Double_t errY) {fErrX2 = errX * errX; fErrY2 = errY * errY;}
           /// Return resolution (cm) on coordinate X
  virtual Double_t GetErrX() const {return TMath::Sqrt(fErrX2);}
           /// Return resolution**2 (cm**2) on coordinate X
  virtual Double_t GetErrX2() const {return fErrX2;}
           /// Return resolution (cm) on coordinate Y
  virtual Double_t GetErrY() const {return TMath::Sqrt(fErrY2);}
           /// Return resolution**2 (cm**2) on coordinate Y
  virtual Double_t GetErrY2() const {return fErrY2;}
  
           /// Set the cluster charge
  virtual void     SetCharge(Double_t q) {fQ = q;}
           /// Set the cluster charge
  virtual Double_t GetCharge() const {return fQ;}
  
           /// Return chamber Id
  virtual Int_t    GetChamberId() const {return AliMUONVCluster::GetChamberId(GetUniqueID());}
           /// Return detection element id
  virtual Int_t    GetDetElemId() const {return AliMUONVCluster::GetDetElemId(GetUniqueID());}
  
  virtual void     SetDigitsId(Int_t nDigits, const UInt_t *digitsId);
           /// Add a digit Id to the array of associated digits
  virtual void     AddDigitId(UInt_t id);
           /// Return number of associated digits
  virtual Int_t    GetNDigits() const {return fNDigits;}
           /// Return Id of digits i
  virtual UInt_t   GetDigitId(Int_t i) const {return (i < fNDigits && fDigitsId) ? fDigitsId[i] : 0;}
  
           /// Set chi2 of cluster
  virtual void     SetChi2( Double_t chi2) {fChi2 = chi2;}
           /// Return chi2 of cluster
  virtual Double_t GetChi2() const {return fChi2;}
  
           /// Set the corresponding MC track number
  virtual void     SetMCLabel(Int_t label) {fMCLabel = label;}
           /// Return the corresponding MC track number
  virtual Int_t    GetMCLabel() const {return fMCLabel;}
  
           /// Return true as the function Compare() is implemented
  Bool_t       IsSortable() const {return kTRUE;}
  Int_t        Compare(const TObject *obj) const;
  
  
private:
  
  Double32_t fX;	///< X of cluster
  Double32_t fY;	///< Y of cluster
  Double32_t fZ;	///< Z of cluster
  
  Double32_t fErrX2;	///< X coordinate error square
  Double32_t fErrY2;	///< Y coordinate error square
  
  Double32_t fQ;	///< Q of cluster (in ADC counts)
  
  Double32_t fChi2;	///< Chi2 of cluster
  
  Int_t    fNDigits;	///< Number of digits attached to the cluster
  /// Indices of digits attached to the cluster
  UInt_t   *fDigitsId;	//[fNDigits] Indices of digits attached to the cluster
  
  Int_t fMCLabel;       ///< Point to the corresponding MC track
  
  
  ClassDef(AliMUONRawClusterV2,2)  //Cluster class for MUON
};

#endif
