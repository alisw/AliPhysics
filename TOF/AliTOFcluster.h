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
#include "AliCluster3D.h"

class AliTOFcluster : public AliCluster3D {
 public:
  AliTOFcluster(); // default ctor
  AliTOFcluster(UShort_t volId, 
     Float_t x,   Float_t y,   Float_t z,
     Float_t sx2, Float_t sxy, Float_t sxz,
                  Float_t sy2, Float_t syz,
                               Float_t sz2,
		Int_t *lab, Int_t * const ind, Int_t *par, 
		Bool_t status, Int_t idx); // ctor
  AliTOFcluster(const AliTOFcluster & cluster); // copy ctor
  AliTOFcluster& operator = (const AliTOFcluster& cls);
  virtual ~AliTOFcluster(); // dtor

  // Getters and Setters
  Double_t GetR() const   {return fR;}   // Cluster Radius
  Double_t GetPhi() const {return fPhi;} // Cluster Phi

  Double_t GetQuality() const {return fQuality;} // Cluster quality getter
  Bool_t   GetStatus()  const {return fStatus;}  // Cluster status getter
  Int_t GetToT()    const {return fToT;}    // Cluster Charge getter
  Int_t GetTDC()    const {return fTDC;}    // Cluster ToF getter
  Int_t GetTDCND()  const {return fTdcND;}  // Cluster ToF getter
  Int_t GetTDCRAW() const {return fTdcRAW;} // Cluster Raw time getter
  Int_t GetADC()    const {return TMath::Abs(fADC);} // Cluster Charge getter
  Int_t IsUsed()    const {return (fADC<0) ? 1 : 0;} // Flagging
  Int_t GetDetInd(Int_t n) const {return fdetIndex[n];} // Cluster Detector Indices getter
  Int_t GetIndex()  const {return fIdx;} // Digit Index getter
  Int_t GetDeltaBC() const {return fDeltaBC;}; // deltaBC
  Int_t GetL0L1Latency() const {return fL0L1Latency;}; // L0L1 latency

  void  Use(Int_t = 0) {fADC=-fADC;}        //  setter
  void  SetQuality(Double_t quality) {fQuality = quality;} // Cluster quality setter
  void  SetStatus(Bool_t status) {fStatus = status;}       // Cluster status setter
  void  SetToT(Int_t ToT) {fToT = ToT;}       // Cluster ToT setter
  void  SetTDC(Int_t Tdc) {fTDC = Tdc;}       // Cluster ToF setter
  void  SetTDCND(Int_t Tdc) {fTdcND = Tdc;}   // Cluster ToFnd setter
  void  SetTDCRAW(Int_t Tdc) {fTdcRAW = Tdc;} // Cluster ToF-raw setter
  void SetDeltaBC(Int_t value) {fDeltaBC = value;}; // deltaBC
  void SetL0L1Latency(Int_t value) {fL0L1Latency = value;}; // L0-L1 latency

 private:

  Int_t fIdx;         // index of the digit related to this cluster
  Int_t fdetIndex[5]; // Cluster detector indices
                      // (sector, plate, strip, padz, padx)
  // Cluster Quality
  Double_t fQuality;  // quality of the best track 

  // Cluster Global Position
  Double_t fR;        // r-coordinate
  Double_t fPhi;      // phi-coordinate

  // TOF Signal parameters
  Int_t  fTDC;        // TDC count
  Int_t  fToT;        // ToT
  Int_t  fADC;        // ADC count
  Int_t  fTdcND;      // TDC count
  Int_t  fTdcRAW;     // RAW TDC count
  Bool_t fStatus;     // cluster online status 
  Int_t fDeltaBC; // deltaBC
  Int_t fL0L1Latency; // L0L1 latency

  ClassDef(AliTOFcluster, 7) // TOF cluster
};

#endif
