#ifndef ALIESDPMDTRACK_H
#define ALIESDPMDTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Event Data Summary Class for pmd tracks
// This is part of the reconstructed ESD events
// for the PMD detector

#include "TObject.h"

class AliESDPmdTrack : public TObject {
 public:
  AliESDPmdTrack();
  virtual ~AliESDPmdTrack(){}
  AliESDPmdTrack (const AliESDPmdTrack &PMDTrack);  // copy constructor
  AliESDPmdTrack &operator=(const AliESDPmdTrack &PMDTrack); // assignment op
  
  void SetDetector(Int_t idet) {fDet = idet;}
  void SetTheta(Float_t theta) {fTheta = theta;}
  void SetPhi(Float_t phi) {fPhi = phi;}
  void SetClusterADC(Float_t cluadc) {fCluADC = cluadc;}
  void SetClusterPID(Float_t clupid) {fCluPID = clupid;}

  Int_t   GetDetector() const {return fDet;}
  Float_t GetTheta() const {return fTheta;}
  Float_t GetPhi() const {return fPhi;}
  Float_t GetClusterADC() const {return fCluADC;}
  Float_t GetClusterPID() const {return fCluPID;}
  
 protected:
  Int_t fDet;      // Detector, 0:PRE, 1:CPV
  Float_t fTheta;  // Theta of the Cluster in radian
  Float_t fPhi;    // Phi of the Cluster in radian
  Float_t fCluADC; // Cluster Energy in ADC
  Float_t fCluPID; // Cluster probability, 1: Photon, 0: Hadron

  ClassDef(AliESDPmdTrack,1)  //PMD ESD track class 
};

#endif 
