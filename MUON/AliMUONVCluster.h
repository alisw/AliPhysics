#ifndef ALIMUONVCLUSTER_H
#define ALIMUONVCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONVCluster
/// \brief abstract base class for clusters
/// 
//  Author Philippe Pillot, Subatech


#include <TObject.h>

class AliMUONVCluster : public TObject {
 public:
  AliMUONVCluster(); // Constructor
  AliMUONVCluster(Int_t chamberId, Int_t detElemId, Int_t clusterIndex);
  virtual ~AliMUONVCluster(); // Destructor
  
           /// Clear method (used by TClonesArray)
  virtual void Clear(Option_t*) = 0;

           /// Set coordinates (cm)
  virtual void     SetXYZ(Double_t x, Double_t y, Double_t z) = 0;
           /// Return coordinate X (cm)
  virtual Double_t GetX() const = 0;
           /// Return coordinate Y (cm)
  virtual Double_t GetY() const = 0;
           /// Return coordinate Z (cm)
  virtual Double_t GetZ() const = 0;
  
           /// Set resolution (cm) on coordinates (X,Y)
  virtual void     SetErrXY(Double_t errX, Double_t errY) = 0;
           /// Return resolution (cm) on coordinate X
  virtual Double_t GetErrX() const = 0;
           /// Return resolution**2 (cm**2) on coordinate X
  virtual Double_t GetErrX2() const = 0;
           /// Return resolution (cm) on coordinate Y
  virtual Double_t GetErrY() const = 0;
           /// Return resolution**2 (cm**2) on coordinate Y
  virtual Double_t GetErrY2() const = 0;
  
           /// Set the cluster charge
  virtual void     SetCharge(Double_t charge) = 0;
           /// Set the cluster charge
  virtual Double_t GetCharge() const = 0;
  
           /// Build a single integer with id information
  static  UInt_t   BuildUniqueID(Int_t chamberId, Int_t detElemId, Int_t clusterIndex)
			{return (((chamberId & 0xF) << 28) | ((detElemId & 0x7FF) << 17) | (clusterIndex & 0x1FFFF));}
           /// Return chamber id (0..), part of the uniqueID
  static  Int_t    GetChamberId(UInt_t uniqueID)    {return (uniqueID & 0xF0000000) >> 28;}
           /// Return detection element id, part of the uniqueID
  static  Int_t    GetDetElemId(UInt_t uniqueID)    {return (uniqueID & 0x0FFE0000) >> 17;}
           /// The index of this cluster (0..), part of the uniqueID
  static  Int_t    GetClusterIndex(UInt_t uniqueID) {return (uniqueID & 0x0001FFFF);}
           /// Return chamber Id
  virtual Int_t    GetChamberId() const = 0;
           /// Return detection element Id
  virtual Int_t    GetDetElemId() const = 0;
  
           /// Set Id of associated digits
  virtual void     SetDigitsId(Int_t nDigits, const UInt_t *digitsId) = 0;
           /// Add a digit Id to the array of associated digits
  virtual void     AddDigitId(UInt_t id) = 0;
           /// Return number of associated digits
  virtual Int_t    GetNDigits() const = 0;
           /// Return Id of digits i
  virtual UInt_t   GetDigitId(Int_t i) const = 0;
  
           /// Set chi2 of cluster
  virtual void     SetChi2(Double_t chi2) = 0;
           /// Return chi2 of cluster
  virtual Double_t GetChi2() const = 0;
  
           /// Set the corresponding MC track number
  virtual void     SetMCLabel(Int_t label) = 0;
	   /// Return the corresponding MC track number
  virtual Int_t    GetMCLabel() const = 0;
  
  virtual void     Print(Option_t *option = "") const;
  
  
  ClassDef(AliMUONVCluster, 1) // abstract base class for cluster
};
	
#endif
