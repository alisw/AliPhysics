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
  AliTOFcluster(Double_t *h, Int_t *ind, Int_t *par, Bool_t status,Int_t *l, Int_t idx); // ctor
  AliTOFcluster(Double_t *h, Int_t *ind, Int_t *par); // ctor
  AliTOFcluster(const AliTOFcluster & cluster); // copy ctor
  virtual ~AliTOFcluster(); // dtor

  Double_t GetR() const   {return fR;}   // Cluster Radius
  Double_t GetPhi() const {return fPhi;} // Cluster Phi
  Double_t GetZ()   const {return fZ;}   // Cluster Z
  Int_t GetTDC() const {return fTDC;} // Cluster ToF
  Int_t GetTDCND() const {return fTdcND;} // Cluster ToF
  Int_t GetTDCRAW() const {return fTdcRAW;} // Cluster Raw time
  Int_t GetADC() const {return TMath::Abs(fADC);}  // Cluster Charge
  Int_t GetToT() const {return fToT;}  // Cluster Charge
  Int_t IsUsed() const {return (fADC<0) ? 1 : 0;}  // Flagging
  Int_t GetLabel(Int_t n) const  {return fLab[n];} // Labels of tracks in Cluster
  Int_t GetDetInd(Int_t n) const {return fdetIndex[n];} // Cluster Detector Indeces
  Int_t GetIndex() const         {return fIdx;}         // Cluster Index
  void     Use() {fADC=-fADC;}
  Double_t GetQuality() const {return fQuality;}
  void     SetQuality(Double_t quality) {fQuality = quality;}
  Bool_t   GetStatus() const {return fStatus;}
  void     SetStatus(Bool_t status) {fStatus = status;}
  void     SetToT(Int_t ToT) {fToT = ToT;}
  void     SetTDC(Int_t Tdc) {fTDC = Tdc;}
  void     SetTDCND(Int_t Tdc) {fTdcND = Tdc;}
  void     SetTDCRAW(Int_t Tdc) {fTdcRAW = Tdc;}

 private:

  Int_t fLab[3];      // track labels
  Int_t fIdx;         // index of this cluster
  Int_t fdetIndex[5]; // Cluster detector Indeces (sector,plate,strip,padz,padx)
  // Cluster Position
  Double_t fR;        // r-coordinate
  Double_t fPhi;      // phi-coordinate
  Double_t fZ;        // z-coordinate
  Double_t fQuality;  // quality of the best track 

  // TOF Signal parameters
  Int_t  fTDC;      // TDC count
  Int_t  fToT;       // ToT
  Int_t  fADC;      // ADC count
  Int_t  fTdcND;      // TDC count
  Int_t  fTdcRAW;      // RAW TDC count
  Bool_t fStatus;      // cluster online status 

  ClassDef(AliTOFcluster, 5) // TOF cluster
};

#endif
