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

  void SetClusterX(Float_t xglobal) {fX = xglobal;}
  void SetClusterY(Float_t yglobal) {fY = yglobal;}
  void SetClusterZ(Float_t zglobal) {fZ = zglobal;}
  void SetClusterADC(Float_t cluadc) {fCluADC = cluadc;}
  void SetClusterCells(Float_t ncell) {fNcell = ncell;}
  void SetClusterPID(Float_t clupid) {fCluPID = clupid;}

  Int_t   GetDetector() const {return fDet;}
  Float_t GetClusterX() const {return fX;}
  Float_t GetClusterY() const {return fY;}
  Float_t GetClusterZ() const {return fZ;}
  Float_t GetClusterADC() const {return fCluADC;}
  Float_t GetClusterCells() const {return fNcell;}
  Float_t GetClusterPID() const {return fCluPID;}
  
 protected:
  Int_t fDet;      // Detector, 0:PRE, 1:CPV
  Float_t fX;      // Cluster X position
  Float_t fY;      // Cluster Y position
  Float_t fZ;      // Cluster Z position (vertex uncorrected)
  Float_t fCluADC; // Cluster Energy in ADC
  Float_t fNcell;  // Cluster cells
  Float_t fCluPID; // Cluster probability, 1: Photon, 0: Hadron

  ClassDef(AliESDPmdTrack,2)  //PMD ESD track class 
};

#endif 
