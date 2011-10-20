#ifndef AliMFTCluster_H
#define AliMFTCluster_H 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Class for the description of the clusters of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliMUONRawCluster.h"
#include "AliMUONVCluster.h"
#include "TObject.h"

//====================================================================================================================================================

class AliMFTCluster : public TObject {

public:

  AliMFTCluster();
//   AliMFTCluster(const AliMFTCluster& pt);
//   AliMFTCluster& operator=(const AliMFTCluster &source);

  virtual ~AliMFTCluster() {};  // destructor

  void SetXYZ(Double_t x, Double_t y, Double_t z) { fX=x; fY=y; fZ=z; }

  void SetX(Double_t x) { fX = x; }
  void SetY(Double_t y) { fY = y; }
  void SetZ(Double_t z) { fZ = z; }

  Double_t GetX() const { return fX; }
  Double_t GetY() const { return fY; }
  Double_t GetZ() const { return fZ; }
  
  void SetErrXYZ(Double_t errX, Double_t errY, Double_t errZ) { fErrX = errX; fErrY = errY; fErrZ = errZ; }
       
  void SetErrX(Double_t errX) { fErrX = errX; }
  void SetErrY(Double_t errY) { fErrY = errY; }
  void SetErrZ(Double_t errZ) { fErrZ = errZ; }

  Double_t GetErrX()  const { return fErrX; }
  Double_t GetErrY()  const { return fErrY; }
  Double_t GetErrZ()  const { return fErrZ; }
  Double_t GetErrX2() const { return fErrX*fErrX; }
  Double_t GetErrY2() const { return fErrY*fErrY; }
  Double_t GetErrZ2() const { return fErrZ*fErrZ; }
  
  void     SetNElectrons(Double_t nElectrons) { fNElectrons = nElectrons; }
  Double_t GetNElectrons() const { return fNElectrons; }
  
  void  AddMCLabel(Int_t label) { if (fNMCTracks==fNMaxMCTracks) return; else fMCLabel[fNMCTracks++]=label; }
  Int_t GetNMCTracks() const { return fNMCTracks; }
  Int_t GetMCLabel(Int_t track) const { if (track<fNMCTracks && track>=0) return fMCLabel[track]; else return -1; }

  void  SetPlane(Int_t plane) { fPlane = plane; }
  Int_t GetPlane() const { return fPlane; }

  void  SetSize(Int_t size) { fSize = size; }
  Int_t GetSize() const { return fSize; }

  void SetLocalChi2(Double_t chi2) { fLocalChi2 = chi2; }
  void SetTrackChi2(Double_t chi2) { fTrackChi2 = chi2; }

  Double_t GetLocalChi2() { return fLocalChi2; }
  Double_t GetTrackChi2() { return fTrackChi2; }

  AliMUONRawCluster* CreateMUONCluster();
  
private:

  static const Int_t fNMaxMCTracks = 30;

  Double_t fX, fY, fZ;   // cluster global coordinates
  Double_t fErrX, fErrY, fErrZ;

  Double_t fNElectrons;
  Int_t fNMCTracks;
  Int_t fPlane;
  Int_t fMCLabel[fNMaxMCTracks];

  Int_t fSize;   // the number of digits composing the cluster

  Double_t fTrackChi2; // Chi2 of the track when the associated cluster was attached
  Double_t fLocalChi2; // Local chi2 of the associated cluster with respect to the track
  
  ClassDef(AliMFTCluster, 1)

};

//====================================================================================================================================================
	
#endif

