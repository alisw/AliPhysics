#ifndef ALITOFCALIBESD_H
#define ALITOFCALIBESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////
//  class for TOF calibration:: simulation of uncalibrated data //
//////////////////////////////////////////////////////////////////

#include "AliPID.h"
#include "AliESDtrack.h"
#include "TObject.h"

//class AliESDtrack;

class AliTOFcalibESD:public AliESDtrack{
public:
  AliTOFcalibESD();
  AliTOFcalibESD(const AliTOFcalibESD& UnCalib);
  ~AliTOFcalibESD();
  Float_t GetToT() const{return fToT;}         //Time Over Threshold
  Float_t GetTOFsignal() const{return fTOFtime;}
  Float_t GetTOFsignalND() const{return fTOFsignalND;}
  Float_t GetIntegratedLength() const{return fIntLen;}
  void GetExternalCovariance(Double_t cov[15]) const;
  void GetIntegratedTimes(Double_t exp[AliPID::kSPECIES]) const;
  Int_t GetCombID()const{return fCombID;}
  Float_t GetP()const{return fP;}
  Int_t GetSector() const {return fSector;}
  Int_t GetPlate() const {return fPlate;}
  Int_t GetStrip() const {return fStrip;}
  Int_t GetPadz() const {return fPadz;}
  Int_t GetPadx() const {return fPadz;}
  void SetToT(Float_t ToT) {fToT=ToT;}
  void SetTOFtime(Float_t TOFtime) {fTOFtime=TOFtime;}
  void SetTOFsignalND(Float_t TOFtimeND) {fTOFsignalND=TOFtimeND;}
  void SetP(Double_t p) {fP=p;}
  void SetIntegratedTime(const Double_t *tracktime);
  void SetCombID(Int_t ID){fCombID = ID;} // 0->pi, 1->K, 2->p
  void CopyFromAliESD(const AliESDtrack* track);
  Bool_t IsSortable() const {return kTRUE;}
  Int_t Compare(const TObject *uncobj) const;
  void SetSector(Int_t isector) {fSector=isector;}
  void SetPlate(Int_t iplate) {fPlate=iplate;}
  void SetStrip(Int_t istrip) {fStrip=istrip;}
  void SetPadz(Int_t ipadz) {fPadz=ipadz;}
  void SetPadx(Int_t ipadx) {fPadx=ipadx;}
private:
  Float_t  fToT;
  Float_t  fIntLen;
  Float_t  fTOFtime;
  Double_t fP;
  Double_t fTrackTime[AliPID::kSPECIES]; // TOFs estimated by the tracking
  Double_t fExtCov[15];
  Int_t    fCombID; 
  Int_t    fSector;
  Int_t    fPlate;
  Int_t    fStrip;
  Int_t    fPadz;
  Int_t    fPadx;
  Float_t  fTOFsignalND;

  ClassDef(AliTOFcalibESD,1);
};
#endif // AliTOFcalibESD_H
