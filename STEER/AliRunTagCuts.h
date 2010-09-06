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
  void SetDipoleField(Float_t Pmag) {fAliceDipoleField = Pmag; fAliceDipoleFieldFlag = kTRUE;}
  void SetRunStartTimeRange(Int_t t0, Int_t t1) {fAliceRunStartTimeMin = t0; fAliceRunStartTimeMax = t1; fAliceRunStartTimeFlag = kTRUE;}
  void SetRunStopTimeRange(Int_t t0, Int_t t1) {fAliceRunStopTimeMin = t0; fAliceRunStopTimeMax = t1; fAliceRunStartTimeFlag = kTRUE;}
  void SetAlirootVersion(TString v) {fAlirootVersion = "VO_ALICE@AliRoot::"; fAlirootVersion += v; fAlirootVersionFlag = kTRUE;}
  void SetRootVersion(TString v) {fRootVersion = "VO_ALICE@ROOT::"; fRootVersion += v; fRootVersionFlag = kTRUE;}
  void SetGeant3Version(TString v) {fGeant3Version = "VO_ALICE@GEANT3::"; fGeant3Version += v; fGeant3VersionFlag = kTRUE;}
  void SetLHCPeriod(TString v) {fLHCPeriod = v; fLHCPeriodFlag = kTRUE; }
  void SetReconstructionPass(TString v) {fRecPass = v; fRecPassFlag = kTRUE; }
  void SetProductionName(TString v) {fProdName = v; fProdNameFlag = kTRUE; }
  void SetRunValidation(Int_t Pn) {fAliceRunValidation = Pn; fAliceRunValidationFlag = kTRUE;}
  void AddRunQualityValue(Int_t qval);
  void SetBeamEnergy(Float_t PE) {fAliceBeamEnergy = PE; fAliceBeamTypeFlag = kTRUE;}
  void SetBeamType(TString Ptype) {fAliceBeamType = Ptype; fAliceCalibrationVersionFlag = kTRUE;}
  void SetCalibVersion(Int_t Pn) {fAliceCalibrationVersion = Pn; fAliceCalibrationVersionFlag = kTRUE;}
  void SetDataType(Int_t i) {fAliceDataType = i; fAliceDataTypeFlag = kTRUE;}
  void SetBeamTriggersRange(ULong_t tmin, ULong_t tmax) { fBeamTriggerRange[0] = tmin; fBeamTriggerRange[1] = tmax; fBeamTriggerFlag = kTRUE; }
  void SetCollisionTriggersRange(ULong_t tmin, ULong_t tmax) { fCollisionTriggerRange[0] = tmin; fCollisionTriggerRange[1] = tmax; fCollisionTriggerFlag = kTRUE; }
  void SetEmptyTriggersRange(ULong_t tmin, ULong_t tmax) { fEmptyTriggerRange[0] = tmin; fEmptyTriggerRange[1] = tmax; fEmptyTriggerFlag = kTRUE; }
  void SetASideTriggersRange(ULong_t tmin, ULong_t tmax) { fASideTriggerRange[0] = tmin; fASideTriggerRange[1] = tmax; fASideTriggerFlag = kTRUE; }
  void SetCSideTriggersRange(ULong_t tmin, ULong_t tmax) { fCSideTriggerRange[0] = tmin; fCSideTriggerRange[1] = tmax; fCSideTriggerFlag = kTRUE; }
  void SetHMTriggersRange(ULong_t tmin, ULong_t tmax) { fHMTriggerRange[0] = tmin; fHMTriggerRange[1] = tmax; fHMTriggerFlag = kTRUE; }
  void SetMuonTriggersRange(ULong_t tmin, ULong_t tmax) { fMuonTriggerRange[0] = tmin; fMuonTriggerRange[1] = tmax; fMuonTriggerFlag = kTRUE; }
  void SetCollisionRatesRange(ULong_t tmin, ULong_t tmax) { fCollisionRateRange[0] = tmin; fCollisionRateRange[1] = tmax; fCollisionRateFlag = kTRUE; }
  void SetMeanVertexsRange(ULong_t tmin, ULong_t tmax) { fMeanVertexRange[0] = tmin; fMeanVertexRange[1] = tmax; fMeanVertexFlag = kTRUE; }
  void SetVertexQualitysRange(ULong_t tmin, ULong_t tmax) { fVertexQualityRange[0] = tmin; fVertexQualityRange[1] = tmax; fVertexQualityFlag = kTRUE; }

  Bool_t IsAccepted(AliRunTag *RunTag) const;

  //____________________________________________________//
 private:
  Int_t   fAliceRunId;                  //the run id
  Bool_t  fAliceRunIdFlag;              //Shows whether this cut is used or not
  Float_t fAliceMagneticField;          //value of the magnetic field
  Bool_t  fAliceMagneticFieldFlag;      //Shows whether this cut is used or not
  Float_t fAliceDipoleField;            //value of the dipole field
  Bool_t  fAliceDipoleFieldFlag;        //Shows whether this cut is used or not
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
  TString fLHCPeriod;                   //LHC period version
  Bool_t  fLHCPeriodFlag;               //Shows whether this cut is used or not
  TString fRecPass;                     //Reconstruction pass
  Bool_t  fRecPassFlag;                 //Shows whether this cut is used or not
  TString fProdName;                    //Production Name
  Bool_t  fProdNameFlag;                //Shows whether this cut is used or not
  Bool_t  fAliceRunValidation;          //validation script
  Bool_t  fAliceRunValidationFlag;      //Shows whether this cut is used or not
  TString fAliceRunQualities;           //selected qualities
  Bool_t  fAliceRunQualitiesFlag;       //Shows whether this cut is used or not
  Float_t fAliceBeamEnergy;             //beam energy cm
  Bool_t  fAliceBeamEnergyFlag;         //Shows whether this cut is used or not
  TString fAliceBeamType;               //run type (pp, AA, pA)
  Bool_t  fAliceBeamTypeFlag;           //Shows whether this cut is used or not
  Int_t   fAliceCalibrationVersion;     //calibration version  
  Bool_t  fAliceCalibrationVersionFlag; //Shows whether this cut is used or not
  Int_t   fAliceDataType;               //0: simulation -- 1: data  
  Bool_t  fAliceDataTypeFlag;           //Shows whether this cut is used or not
  ULong_t fBeamTriggerRange[2];         //Beam trigger maximum and minimum
  Bool_t  fBeamTriggerFlag;             //Shows whether this cut is used or not
  ULong_t fCollisionTriggerRange[2];    //Collision trigger maximum and minimum
  Bool_t  fCollisionTriggerFlag;        //Shows whether this cut is used or not
  ULong_t fEmptyTriggerRange[2];        //Empty trigger maximum and minimum
  Bool_t  fEmptyTriggerFlag;            //Shows whether this cut is used or not
  ULong_t fASideTriggerRange[2];        //ASide trigger maximum and minimum
  Bool_t  fASideTriggerFlag;            //Shows whether this cut is used or not
  ULong_t fCSideTriggerRange[2];        //CSide trigger maximum and minimum
  Bool_t  fCSideTriggerFlag;            //Shows whether this cut is used or not
  ULong_t fHMTriggerRange[2];           //High Multiplicity trigger maximum and minimum
  Bool_t  fHMTriggerFlag;               //Shows whether this cut is used or not
  ULong_t fMuonTriggerRange[2];         //Muon trigger maximum and minimum
  Bool_t  fMuonTriggerFlag;             //Shows whether this cut is used or not
  Float_t fCollisionRateRange[2];       //Collision rate range
  Bool_t  fCollisionRateFlag;           //Shows whether this cut is used or not
  Float_t fMeanVertexRange[2];          //Mean Vertex Postion
  Bool_t  fMeanVertexFlag;              //Shows whether this cut is used or not
  Float_t fVertexQualityRange[2];       //Mean Vertex quality
  Bool_t  fVertexQualityFlag;           //Shows whether this cut is used or not

  ClassDef(AliRunTagCuts, 2)
};

#endif
