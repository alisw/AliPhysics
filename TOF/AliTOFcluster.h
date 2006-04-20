#ifndef ALITOFCLUSTER_H
#define ALITOFCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////
//                                          //
//     Class for TOF cluster definition     //
//                                          //
//////////////////////////////////////////////

#include "TMath.h"
#include "TObject.h"

class AliTOFcluster : public TObject {
 public:
  AliTOFcluster(); // default ctor
  AliTOFcluster(Double_t *h, Int_t *l, Int_t *ind, Int_t idx, Float_t ToT, Double_t TdcND); // ctor
  AliTOFcluster(Double_t *h, Int_t *ind); // new ctor
  AliTOFcluster(const AliTOFcluster & cluster); // copy ctor
  virtual ~AliTOFcluster(); // dtor

  Double_t GetR() const   {return fR;}   // Cluster Radius
  Double_t GetPhi() const {return fPhi;} // Cluster Phi
  Double_t GetZ()   const {return fZ;}   // Cluster Z
  Double_t GetTDC() const {return fTDC;} // Cluster ToF
  Double_t GetTDCND() const {return fTdcND;} // Cluster ToF
  Double_t GetADC() const {return TMath::Abs(fADC);}  // Cluster Charge
  Double_t GetToT() const {return fToT;}  // Cluster Charge
  Int_t    IsUsed() const {return (fADC<0) ? 1 : 0;}  // Flagging
  Int_t    GetLabel(Int_t n) const  {return fLab[n];} // Labels of tracks in Cluster
  Int_t    GetDetInd(Int_t n) const {return fdetIndex[n];} // Cluster Detector Indeces
  Int_t    GetIndex() const         {return fIdx;}         // Cluster Index
  void     Use() {fADC=-fADC;}
  Double_t GetQuality() const {return fQuality;}
  void     SetQuality(Double_t quality) {fQuality = quality;}
  void     SetToT(Float_t ToT) {fToT = ToT;}
  void     SetTDC(Float_t Tdc) {fTDC = Tdc;}
  void     SetTDCND(Float_t Tdc) {fTdcND = Tdc;}
 private:

  Int_t fLab[3];      // track labels
  Int_t fIdx;         // index of this cluster
  Int_t fdetIndex[5]; // Cluster detector Indeces (sector,plate,strip,padz,padx)
  Double_t fR;        // r-coordinate
  Double_t fPhi;      // phi-coordinate
  Double_t fZ;        // z-coordinate
  Double_t fTDC;      // TDC count
  Double_t fADC;      // ADC count
  Double_t fQuality;  // quality of the best track 
  Float_t  fToT;       // ToT
  Double_t fTdcND;      // TDC count

  ClassDef(AliTOFcluster, 2) // TOF cluster
};

#endif
