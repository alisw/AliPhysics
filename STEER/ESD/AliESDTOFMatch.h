#ifndef ALIESDTOFMATCH_H
#define ALIESDTOFMATCH_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

//----------------------------------------------------------------------//
//                                                                      //
// AliESDTOFMatch Class                                                 //
//                                                                      //
//----------------------------------------------------------------------//

#include "AliVTOFMatch.h"

class AliESDTOFMatch : public AliVTOFMatch
{
 public:
  enum {
    kActive=0x1,
    kDeadZoneOK=0x2,
    kGainOK=0x4,
    kDriftOK=0x8,
    kHVROCOK=0x10,
    kHVDriftOK=0x20
  };
  AliESDTOFMatch();
  AliESDTOFMatch(Int_t i,Double_t inttimes[AliPID::kSPECIESC],Double_t dx,Double_t dy,Double_t dz,Double_t l);
  AliESDTOFMatch(AliESDTOFMatch &source);
  virtual ~AliESDTOFMatch() {}
  AliESDTOFMatch &operator=(const AliESDTOFMatch& source);
  virtual Float_t GetDx() const {return fDx;}
  virtual Float_t GetDy() const {return fDy;}
  virtual Float_t GetDz() const {return fDz;}
  virtual Float_t GetTrackLength() const {return fTrackLength;}
  virtual void SetDx(Double_t delta) {fDx = delta;}
  virtual void SetDy(Double_t delta) {fDy = delta;}
  virtual void SetDz(Double_t delta) {fDz = delta;}
  virtual void SetTrackLength(Double_t length) {fTrackLength = length;}
  //
  virtual Double_t GetIntegratedTimes(Int_t i) const {return fIntegratedTimes[i];}
  virtual void     GetIntegratedTimes(Double_t* t) const {for (int i=AliPID::kSPECIESC;i--;) t[i]=fIntegratedTimes[i];}
  virtual void     SetIntegratedTimes(Double_t* t)       {for (int i=AliPID::kSPECIESC;i--;) fIntegratedTimes[i]=t[i];}
  //
  virtual Int_t GetTrackIndex()           const {return GetUniqueID();}
  virtual void  SetTrackIndex(Int_t id)         {SetUniqueID(id);}
  void Print(const Option_t *opt=0) const;
  //
  Short_t GetX0Layer(Int_t layer){ return fX0Layer[layer];}
  Short_t GetRhoLayer(Int_t layer){ return fRhoLayer[layer];}
  Short_t GetTRDstatus(Int_t layer){ return fTRDStatus[layer];}
  Short_t GetTRDncls(Int_t layer){ return fTRDncls[layer];}
  Double_t GetChi2TPCTOFS(){return fChi2TPCTOFS;}
  Double_t GetChi2TPCTrDTOFS(){return fChi2TPCTRDTOFS;}
  void SetX0Layer(Int_t layer,Float_t value){ fX0Layer[layer]=value;}
  void SetRhoLayer(Int_t layer,Float_t value){ fRhoLayer[layer]=value;}
  void SetTRDstatus(Int_t layer,Short_t value){ fTRDStatus[layer]=value;}
  void SetTRDstatusBit(Int_t layer,Short_t bit){ fTRDStatus[layer]|=bit;}
  void SetTRDncls(Int_t layer,Short_t value){ fTRDncls[layer]=value;}
  void SetChi2TPCTOFS(Float_t value){fChi2TPCTOFS=value;}
  void SetChi2TPCTRDTOFS(Float_t value){fChi2TPCTRDTOFS=value;}

  //
 protected:
  Double32_t fDx;                // DeltaX residual
  Double32_t fDy;                //! DeltaY residual
  Double32_t fDz;                // DeltaZ residual
  Double32_t fTrackLength;       // track Length
  Double32_t fIntegratedTimes[AliPID::kSPECIESC]; // int timex
  Double32_t fX0Layer[6];        // [0,0,8]  X0 vector estimator
  Double32_t fRhoLayer[6];       // [0,0,8]  Rho vector estimator
  Short_t    fTRDStatus[6];      // TRD chamber status (active, dead zone, gain status, vdrift status, HV status)
  Short_t    fTRDncls[6];        // [0,60] - number of TRD clusters in road
  Double32_t fChi2TPCTOFS;       // [0,0,8] - TPC-TOF      sqrt(chi2)
  Double32_t fChi2TPCTRDTOFS;    // [0,0,8] - TPC-TRD-TOF  sqrt(chi2)
  //
  ClassDef(AliESDTOFMatch, 2) // TOF matchable hit
    //
};
#endif
