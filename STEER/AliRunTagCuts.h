#ifndef ALIRUNTAGCUTS_H
#define ALIRUNTAGCUTS_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                       Class AliRunTagCuts
//              This is the class for the cuts in run tags
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>

class AliRunTag;

//___________________________________________________________________________
class AliRunTagCuts : public TObject {
 public:
  AliRunTagCuts();
  ~AliRunTagCuts();
  void Reset();
  
 //____________________________________________________//
  void SetRunId(Int_t Pid) {fAliceRunId = Pid; fAliceRunIdFlag = kTRUE;}
  void SetMagneticField(Float_t Pmag) {fAliceMagneticField = Pmag; fAliceMagneticFieldFlag = kTRUE;}
  void SetRunStartTimeRange(Int_t t0, Int_t t1) {fAliceRunStartTimeMin = t0; fAliceRunStartTimeMax = t1; fAliceRunStartTimeFlag = kTRUE;}
  void SetRunStopTimeRange(Int_t t0, Int_t t1) {fAliceRunStopTimeMin = t0; fAliceRunStopTimeMax = t1; fAliceRunStartTimeFlag = kTRUE;}
  void SetAlirootVersion(TString v) {fAlirootVersion = "VO_ALICE@AliRoot::"; fAlirootVersion += v; fAlirootVersionFlag = kTRUE;}
  void SetRootVersion(TString v) {fRootVersion = "VO_ALICE@ROOT::"; fRootVersion += v; fRootVersionFlag = kTRUE;}
  void SetGeant3Version(TString v) {fGeant3Version = "VO_ALICE@GEANT3::"; fGeant3Version += v; fGeant3VersionFlag = kTRUE;}
  void SetRunQuality(Int_t Pn) {fAliceRunQuality = Pn; fAliceRunQualityFlag = kTRUE;}
  void SetBeamEnergy(Float_t PE) {fAliceBeamEnergy = PE; fAliceBeamTypeFlag = kTRUE;}
  void SetBeamType(TString Ptype) {fAliceBeamType = Ptype; fAliceCalibrationVersionFlag = kTRUE;}
  void SetCalibVersion(Int_t Pn) {fAliceCalibrationVersion = Pn; fAliceCalibrationVersionFlag = kTRUE;}
  void SetDataType(Int_t i) {fAliceDataType = i; fAliceDataTypeFlag = kTRUE;}
 
  Bool_t IsAccepted(AliRunTag *RunTag) const;

  //____________________________________________________//
 private:
  Int_t   fAliceRunId;                  //the run id
  Bool_t  fAliceRunIdFlag;              //Shows whether this cut is used or not
  Float_t fAliceMagneticField;          //value of the magnetic field
  Bool_t  fAliceMagneticFieldFlag;      //Shows whether this cut is used or not
  Int_t   fAliceRunStartTimeMin;        //minimum run start date
  Int_t   fAliceRunStartTimeMax;        //maximum run start date
  Bool_t  fAliceRunStartTimeFlag;       //Shows whether this cut is used or not
  Int_t   fAliceRunStopTimeMin;         //minmum run stop date
  Int_t   fAliceRunStopTimeMax;         //maximum run stop date
  Bool_t  fAliceRunStopTimeFlag;        //Shows whether this cut is used or not
  TString fAlirootVersion;              //aliroot version
  Bool_t  fAlirootVersionFlag;          //Shows whether this cut is used or not
  TString fRootVersion;                 //root version
  Bool_t  fRootVersionFlag;             //Shows whether this cut is used or not
  TString fGeant3Version;               //geant3 version
  Bool_t  fGeant3VersionFlag;           //Shows whether this cut is used or not
  Bool_t  fAliceRunQuality;             //validation script
  Bool_t  fAliceRunQualityFlag;         //Shows whether this cut is used or not
  Float_t fAliceBeamEnergy;             //beam energy cm
  Bool_t  fAliceBeamEnergyFlag;         //Shows whether this cut is used or not
  TString fAliceBeamType;               //run type (pp, AA, pA)
  Bool_t  fAliceBeamTypeFlag;           //Shows whether this cut is used or not
  Int_t   fAliceCalibrationVersion;     //calibration version  
  Bool_t  fAliceCalibrationVersionFlag; //Shows whether this cut is used or not
  Int_t   fAliceDataType;               //0: simulation -- 1: data  
  Bool_t  fAliceDataTypeFlag;           //Shows whether this cut is used or not

  ClassDef(AliRunTagCuts, 1)
};

#endif
