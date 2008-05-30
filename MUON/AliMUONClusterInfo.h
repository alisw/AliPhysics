#ifndef ALIMUONCLUSTERINFO_H
#define ALIMUONCLUSTERINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup evaluation 
/// \class AliMUONClusterInfo
/// \brief Class to summarize ESD data at cluster
//  Author Philippe Pillot, Subatech


#include <TObject.h>
#include <TClonesArray.h>
#include <AliMUONPadInfo.h>

class AliMUONClusterInfo : public TObject {
public:
  AliMUONClusterInfo(); // Constructor
  virtual ~AliMUONClusterInfo(); //< Destructor
  AliMUONClusterInfo(const AliMUONClusterInfo& cluster);
  AliMUONClusterInfo& operator=(const AliMUONClusterInfo& cluster);
  
  virtual void Clear(Option_t* opt = "");
  
  void     Print(Option_t * option = "") const;  
  
  
  // ------ general info ------
  /// set event number
  void     SetEventId(Int_t eventId) {fEventId = eventId;}
  /// return event ID
  Int_t    GetEventId() const {return fEventId;}
  
  /// Set cluster/track Z-position (cm)
  void     SetZ(Double_t z) {fZ = z;}
  /// Return cluster/track Z-position (cm)
  Double_t GetZ() const {return fZ;}
  
  
  // ------ cluster info ------
  /// set cluster ID
  void     SetClusterId(UInt_t clusterId) {fClusterId = clusterId; SetUniqueID(clusterId);}
  /// return cluster ID
  UInt_t   GetClusterId() const    {return fClusterId;}
  /// Return chamber ID (0..), part of the cluster ID
  Int_t    GetChamberId() const    {return (fClusterId & 0xF0000000) >> 28;}
  /// Return detection element ID, part of the cluster ID
  Int_t    GetDetElemId() const    {return (fClusterId & 0x0FFE0000) >> 17;}
  /// Return the index of this cluster (0..), part of the cluster ID
  Int_t    GetClusterIndex() const {return (fClusterId & 0x0001FFFF);}
  
  /// Set cluster coordinates (cm)
  void     SetClusterXY(Double_t x, Double_t y) {fClusterX = x; fClusterY = y;}
  /// Return cluster X-position (cm)
  Double_t GetClusterX() const {return fClusterX;}
  /// Return cluster Y-position (cm)
  Double_t GetClusterY() const {return fClusterY;}
  
  /// Set cluster resolution (cm)
  void     SetClusterXYErr(Double_t xErr, Double_t yErr) {fClusterXErr = xErr; fClusterYErr = yErr;}
  /// Return cluster X-resolution (cm)
  Double_t GetClusterXErr() const {return fClusterXErr;}
  /// Return cluster Y-resolution (cm)
  Double_t GetClusterYErr() const {return fClusterYErr;}
  
  /// set cluster Chi2
  void     SetClusterChi2(Double_t clusterChi2) {fClusterChi2 = clusterChi2;}
  /// return cluster Chi2
  Double_t GetClusterChi2() const {return fClusterChi2;}

  /// Set the total cluster charge
  void     SetClusterCharge(Double_t charge) {fClusterCharge = charge;}
  /// Return the total cluster charge
  Double_t GetClusterCharge() const {return fClusterCharge;}
  
  
  // ------ track info ------
  /// set track ID
  void     SetTrackId(UInt_t trackId) {fTrackId = trackId;}
  /// return track ID
  UInt_t   GetTrackId() const {return fTrackId;}
  
  /// Set track coordinates (cm)
  void     SetTrackXY(Double_t x, Double_t y) {fTrackX = x; fTrackY = y;}
  /// Return track X-position (cm)
  Double_t GetTrackX() const {return fTrackX;}
  /// Return track Y-position (cm)
  Double_t GetTrackY() const {return fTrackY;}
  /// Set track angles (radian)
  void     SetTrackThetaXY(Double_t thetaX, Double_t thetaY) {fTrackThetaX = thetaX; fTrackThetaY = thetaY;}
  /// Return track ThetaX angle (radian)
  Double_t GetTrackThetaX() const {return fTrackThetaX;}
  /// Return track ThetaY angle (radian)
  Double_t GetTrackThetaY() const {return fTrackThetaY;}
  /// Set track momentum (MeV/c)
  void     SetTrackP(Double_t p) {fTrackP = p;}
  /// Return track momentum (MeV/c)
  Double_t GetTrackP() const {return fTrackP;}
  
  /// Set track resolution (cm)
  void     SetTrackXYErr(Double_t xErr, Double_t yErr) {fTrackXErr = xErr; fTrackYErr = yErr;}
  /// Return track X-resolution (cm)
  Double_t GetTrackXErr() const {return fTrackXErr;}
  /// Return track Y-resolution (cm)
  Double_t GetTrackYErr() const {return fTrackYErr;}

  /// set track Chi2
  void     SetTrackChi2(Double_t trackChi2) {fTrackChi2 = trackChi2;}
  /// return track Chi2
  Double_t GetTrackChi2() const {return fTrackChi2;}
  
  /// Set the muon charge
  void     SetTrackCharge(Short_t charge) {fTrackCharge = charge;}
  /// Return the muon charge
  Short_t  GetTrackCharge() const {return fTrackCharge;}
  
  
  // ------ pad info ------
  /// return the number of pads attached to the cluster
  Int_t    GetNPads() const {return fPads->GetEntriesFast();}
  /// return the array of pads attached to the cluster
  TClonesArray& GetPads() const {return *fPads;}
  /// attach a pad to the cluster
  void     AddPad(const AliMUONPadInfo &pad) {new ((*fPads)[fNPads++]) AliMUONPadInfo(pad);}

    
protected:
  
  // general info
  Int_t      fEventId;    ///< event number
  Double32_t fZ;          ///< track/cluster Z position
  
  // cluster info
  UInt_t     fClusterId;     ///< cluster ID
  Double32_t fClusterX;      ///< cluster X position
  Double32_t fClusterY;      ///< cluster Y position
  Double32_t fClusterXErr;   ///< cluster X resolution
  Double32_t fClusterYErr;   ///< cluster Y resolution
  Double32_t fClusterChi2;   ///< cluster chi2
  Double32_t fClusterCharge; ///< cluster charge
  
  // track info
  UInt_t     fTrackId;     ///< track ID
  Double32_t fTrackX;      ///< track X position
  Double32_t fTrackY;      ///< track Y position
  Double32_t fTrackThetaX; ///< track Theta_X angle
  Double32_t fTrackThetaY; ///< track Theta_Y angle
  Double32_t fTrackP;      ///< track momentum
  Double32_t fTrackXErr;   ///< track X resolution
  Double32_t fTrackYErr;   ///< track Y resolution
  Double32_t fTrackChi2;   ///< track normalized chi2
  Short_t    fTrackCharge; ///< track charge  
  
  Int_t         fNPads; ///< nPads  
  TClonesArray* fPads;  ///< Array of pads attached to the cluster
    
  ClassDef(AliMUONClusterInfo, 2)
};

#endif
