#ifndef ALIESDMUONCLUSTER_H
#define ALIESDMUONCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \class AliESDMuonCluster
/// \brief Class to describe the MUON clusters in the Event Summary Data
//  Author Philippe Pillot, Subatech


#include <TObject.h>

class AliESDMuonCluster : public TObject {
public:
  AliESDMuonCluster(); // Constructor
  virtual ~AliESDMuonCluster() {} ///< Destructor
  AliESDMuonCluster(const AliESDMuonCluster& cluster);
  AliESDMuonCluster& operator=(const AliESDMuonCluster& cluster);
  
  /// Clear method (used by TClonesArray)
  void Clear(Option_t*) {}
  
  /// Set coordinates (cm)
  void     SetXYZ(Double_t x, Double_t y, Double_t z) {fXYZ[0] = x; fXYZ[1] = y; fXYZ[2] = z;}
  /// Return X-position (cm)
  Double_t GetX() const {return fXYZ[0];}
  /// Return Y-position (cm)
  Double_t GetY() const {return fXYZ[1];}
  /// Return Z-position (cm)
  Double_t GetZ() const {return fXYZ[2];}
  
  /// Set (X,Y) resolution (cm)
  void     SetErrXY(Double_t errX, Double_t errY) {fErrXY[0] = errX; fErrXY[1] = errY;}
  /// Return X-resolution (cm)
  Double_t GetErrX() const  {return fErrXY[0];}
  /// Return X-resolution**2 (cm**2)
  Double_t GetErrX2() const {return fErrXY[0]*fErrXY[0];}
  /// Return Y-resolution (cm)
  Double_t GetErrY() const  {return fErrXY[1];}
  /// Return Y-resolution**2 (cm**2)
  Double_t GetErrY2() const {return fErrXY[1]*fErrXY[1];}
  
  /// Return chamber id (0..), part of the uniqueID
  Int_t    GetChamberId() const    {return (GetUniqueID() & 0xF0000000) >> 28;}
  /// Return detection element id, part of the uniqueID
  Int_t    GetDetElemId() const    {return (GetUniqueID() & 0x0FFE0000) >> 17;}
  /// Returnt the index of this cluster (0..), part of the uniqueID
  Int_t    GetClusterIndex() const {return (GetUniqueID() & 0x0001FFFF);}
  
  void     Print(Option_t */*option*/ = "") const;
  
  
protected:
  Double32_t fXYZ[3];   ///< cluster position
  Double32_t fErrXY[2]; ///< transverse position errors
  
  
  ClassDef(AliESDMuonCluster, 1) // MUON ESD cluster class
};

#endif
