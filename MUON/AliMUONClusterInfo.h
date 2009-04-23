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
  /// set run number
  void     SetRunId(Int_t runId) {fRunId = runId;}
  /// return run ID
  Int_t    GetRunId() const {return fRunId;}

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
  /// Return the cluster charge for cathode iC
  Double_t GetClusterCharge(Int_t iC) const ;
  /// Return the bending cluster charge
  Double_t GetClusterChargeB() const {return GetClusterCharge(0);}
  /// Return the non bending cluster charge
  Double_t GetClusterChargeNB() const {return GetClusterCharge(1);}
  
  
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

  /// Get the total number of hits associated to the track leaving this cluster
  UChar_t  GetTrackNHits(void) const {return fTrackNHits;}
  /// Set the total number of hits associated to the track leaving this cluster
    void     SetTrackNHits(UInt_t NHits) {fTrackNHits = NHits;}

  /// Get the map of hit chambers
  UInt_t   GetTrackChamberHitMap() const {return fTrackChamberHitMap;}
  /// Set the map of hit chambers
  void     SetTrackChamberHitMap(UInt_t trackChamberHitMap) {fTrackChamberHitMap = trackChamberHitMap;}
  /// Is chamber hit by track
  Bool_t   IsChamberHit(Int_t chamber) const {return (Bool_t) ((fTrackChamberHitMap & BIT(chamber)) != 0);}
  
  // ------ pad info ------
  /// return the number of pads attached to the cluster
  Int_t    GetNPads() const {return fPads->GetEntriesFast();}
  /// return the number of pads attached to the cluster in cathode iC
  Int_t    GetNPads(Int_t iC) const ;
  /// return the number of bending pads attached to the cluster
  Int_t    GetNPadsB() const {return GetNPads(0);}
  /// return the number of non bending pads attached to the cluster
  Int_t    GetNPadsNB() const {return GetNPads(1);}
  /// return the number of pads attached to the cluster
  Int_t    GetNPadsX(Int_t iC) const ;
  /// return the number of pads attached to the cluster
  Int_t    GetNPadsXB() const {return GetNPadsX(0);}
  /// return the number of pads attached to the cluster
  Int_t    GetNPadsXNB() const {return GetNPadsX(1);}
  /// return the number of pads attached to the cluster
  Int_t    GetNPadsY(Int_t iC) const ;
  /// return the number of pads attached to the cluster
  Int_t    GetNPadsYB() const {return GetNPadsY(0);}
  /// return the number of pads attached to the cluster
  Int_t    GetNPadsYNB() const {return GetNPadsY(1);}
  /// return the array of pads attached to the cluster
  TClonesArray& GetPads() const {return *fPads;}
  /// attach a pad to the cluster
  void     AddPad(const AliMUONPadInfo &pad) {new ((*fPads)[fNPads++]) AliMUONPadInfo(pad);}

    
protected:
  
  // general info
  Int_t      fRunId;    ///< run number
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
  UChar_t    fTrackNHits;   ///< track number of hits
  UInt_t     fTrackChamberHitMap; ///< Map of clusters in tracking chambers

  Int_t         fNPads; ///< nPads  
  TClonesArray* fPads;  ///< Array of pads attached to the cluster
    
  ClassDef(AliMUONClusterInfo, 3)
};

#endif
