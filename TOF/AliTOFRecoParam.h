#ifndef ALITOFRECOPARAM_H
#define ALITOFRECOPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with TOF reconstruction parameters                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliDetectorRecoParam.h"

#include "AliTOFGeometry.h"

class AliTOFRecoParam : public AliDetectorRecoParam
{
 public: 
  AliTOFRecoParam();
  AliTOFRecoParam(const AliTOFRecoParam &p); //Copy Ctor 
  AliTOFRecoParam& operator=(const AliTOFRecoParam &p); // ass. op.
  virtual ~AliTOFRecoParam();

  virtual void PrintParameters() const;

  //Getters
  Bool_t   UseTimeZero()           const  { return fTimeZero;}
  Bool_t   GetTimeZerofromT0()     const  { return fTimeZerofromT0;}
  Bool_t   GetTimeZerofromTOF()    const  { return fTimeZerofromTOF;}
  Bool_t   GetTimeWalkCorr()       const  { return fTimeWalkCorr;}
  Bool_t   GetApplyPbPbCuts()      const  { return fApplyPbPbCuts;}

  Float_t  GetWindowSizeMaxY()     const  { return fWindowSizeMaxY;}
  Float_t  GetWindowSizeMaxZ()     const  { return fWindowSizeMaxZ;}
  Float_t  GetWindowScaleFact()    const  { return fWindowScaleFact;}
  Float_t  GetDistanceCut()        const  { return fDistanceCut;}
  Float_t  GetSensRadius()         const  { return fSensRadius;}
  Float_t  GetStepSize()           const  { return fStepSize;}
  Double_t  GetMaxChi2()           const  { return fMaxChi2;}
  Double_t  GetMaxChi2TRD()        const  { return fMaxChi2TRD;}
  Double_t  GetTimeResolution()    const  { return fTimeResolution;}
  Double_t  GetTimeNSigma()        const  { return fTimeNSigma;}
  Int_t     GetMaxDeltaTime()      const  { return fMaxDeltaTime;}
  Int_t     GetStartTimeType()     const  { return fStartTimeType;}

  //Setters

  void   SetTimeZero( Bool_t flag)        {fTimeZero=flag;}
  void   SetTimeZerofromT0( Bool_t flag)  {fTimeZerofromT0=flag;}
  void   SetTimeZerofromTOF(Bool_t flag)  {fTimeZerofromTOF=flag;}
  void   SetTimeWalkCorr(Bool_t flag)     {fTimeWalkCorr=flag;}
  void   SetApplyPbPbCuts(Bool_t flag)    {fApplyPbPbCuts=flag;}

  void  SetWindowSizeMaxY(Float_t in)   {fWindowSizeMaxY=in;}
  void  SetWindowSizeMaxZ(Float_t in)   {fWindowSizeMaxZ=in;}
  void  SetWindowScaleFact(Float_t in) {fWindowScaleFact=in;}
  void  SetDistanceCut(Float_t in)  {fDistanceCut=in;}
  void  SetSensRadius(Float_t in)  {fSensRadius=in;}
  void  SetStepSize(Float_t in)  {fStepSize=in;}
  void  SetMaxChi2(Double_t in)  {fMaxChi2=in;}
  void  SetMaxChi2TRD(Double_t in)  {fMaxChi2TRD=in;}
  void  SetTimeResolution(Double_t in)  {fTimeResolution=in;}
  void  SetTimeNSigma(Double_t in)  {fTimeNSigma=in;}
  void  SetMaxDeltaTime(Int_t a)   {fMaxDeltaTime = a;}; // to set deltaTime [bin number]
  void  SetMaxDeltaTime(Float_t a) {fMaxDeltaTime = (Int_t)(a/AliTOFGeometry::TdcBinWidth());}; // to set deltaTime [ps]
  void  SetStartTimeType(Int_t a) {fStartTimeType = a;};

  static   AliTOFRecoParam *GetPbPbparam();       // reco param for PbPb.
  static   AliTOFRecoParam *GetPPparam();         // reco param for PP
 private:

  Bool_t fTimeZero; //use Time Zero info in Pid
  Bool_t fTimeZerofromT0; // Use Time Zero as determined by T0
  Bool_t fTimeZerofromTOF; //Use Time Zero as determined from TOF
  Bool_t fTimeWalkCorr; // Correct for signal time walk in z
  Bool_t fApplyPbPbCuts; //apply "high flux" cuts

  Float_t fWindowSizeMaxY;  // cluster search window size, Y (cm)
  Float_t fWindowSizeMaxZ;  // cluster search window size, Z (cm)
  Float_t fWindowScaleFact;  // cluster search window, scale factor
  Float_t fDistanceCut;  // cut on the closest approach distance
  Float_t fSensRadius;  // Average radius of sensitive volumes (cm)
  Float_t fStepSize;  // Propagation step size (cm)
  Double_t fMaxChi2;  // maximum X2 track-tof clusters
  Double_t fMaxChi2TRD;  // maximum X2 track-tof clusters (TRD)
  Double_t fTimeResolution;  // Time resolution for resp. function in PID (ps)
  Double_t fTimeNSigma;  // N-Sigma Range used for resp. function in PID 
  Int_t    fMaxDeltaTime; // max time difference in between two tof
                          // measurements for two neighbouring pads
                          // [number of bins]
  Int_t fStartTimeType; // flag to set start_time to be used for TOF PID

  ClassDef(AliTOFRecoParam, 5)
};

#endif
