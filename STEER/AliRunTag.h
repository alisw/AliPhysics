#ifndef ALIRUNTAG_H
#define ALIRUNTAG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliRunTag
//   This is the class to deal with the tags for the run level
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TString.h>
#include <TClonesArray.h>
#include "AliLHCTag.h"
#include "AliDetectorTag.h"
#include "AliQA.h"

class AliEventTag;
//class AliDetectorTag;


//___________________________________________________________________________
class AliRunTag : public TObject {
 public:
  AliRunTag();
  virtual ~AliRunTag();
  AliRunTag(const AliRunTag& qa) ;   
  AliRunTag& operator = (const AliRunTag& tag) ;
  //____________________________________________________//
  void SetRunId(Int_t Pid) {fAliceRunId = Pid;}
  void SetMagneticField(Float_t Pmag) {fAliceMagneticField = Pmag;}
  void SetDipoleField(Float_t Pmag) {fAliceDipoleField = Pmag;}
  void SetRunStartTime(Int_t Pt0) {fAliceRunStartTime = Pt0;}
  void SetRunStopTime(Int_t Pt1) {fAliceRunStopTime = Pt1;}
  void SetAlirootVersion(TString v) {fAlirootVersion = v;}
  void SetRootVersion(TString v) {fRootVersion = v;}
  void SetGeant3Version(TString v) {fGeant3Version = v;}
  void SetLHCPeriod(TString v) {fLHCPeriod = v; }
  void SetReconstructionPass(TString v) {fRecPass = v; }
  void SetProductionName(TString v) {fProductionName = v; }
  void SetRunQuality(Int_t Pn) {fAliceRunQuality = Pn;}
  void SetBeamEnergy(Float_t PE) {fAliceBeamEnergy = PE;}
  void SetBeamType(TString Ptype) {fAliceBeamType = Ptype;}
  void SetCalibVersion(Int_t Pn) {fAliceCalibrationVersion = Pn;}
  void SetDataType(Int_t i) {fAliceDataType = i;}
  void SetNEvents(Int_t Pn) { fNumEvents = Pn; }
  void SetLHCTag(Float_t Plumin, TString type);
  void SetDetectorTag(UInt_t mask, UInt_t maskReco=0);
  void SetQA(const AliQA &qa) { fQA=qa; }  	
  void SetQAArray(ULong_t * qa, Int_t qalength) ; 
  void SetEventSpecies(Bool_t * es, Int_t eslength) ;
  void AddEventTag(const AliEventTag &t);
  void Clear(const char * opt = "");

  void CopyStandardContent(AliRunTag *oldtag);
  
  //____________________________________________________//
  Int_t       GetRunId() const {return fAliceRunId;}
  Float_t     GetMagneticField() const {return fAliceMagneticField;}
  Float_t     GetDipoleField() const {return fAliceDipoleField;}
  Int_t       GetRunStartTime() const {return fAliceRunStartTime;}
  Int_t       GetRunStopTime() const {return fAliceRunStopTime;}
  const char* GetAlirootVersion() const {return fAlirootVersion.Data();}
  const char* GetRootVersion() const {return fRootVersion.Data();}
  const char* GetGeant3Version() const {return fGeant3Version.Data();}
  const char* GetLHCPeriod() const {return fLHCPeriod.Data();}
  const char* GetReconstructionPass() const {return fRecPass.Data();}
  const char* GetProductionName() const {return fProductionName.Data();}
  Int_t       GetRunQuality() const {return fAliceRunQuality;}
  Float_t     GetBeamEnergy() const {return fAliceBeamEnergy;}
  const char *GetBeamType() const {return fAliceBeamType.Data();}
  Int_t       GetCalibVersion() const {return fAliceCalibrationVersion;}
  Int_t       GetDataType() const {return fAliceDataType;}
  Int_t       GetNEvents() const {return fNumEvents;}
  AliLHCTag  *GetLHCTag() {return &fLHCTag; } 
  AliDetectorTag *GetDetectorTags() {return &fDetectorTag;}
  const TClonesArray *GetEventTags() const {return &fEventTag;}
  const AliQA *GetQA() const {return &fQA;}
  ULong_t *  GetQAArray() const {return fQAArray;}	
  Int_t      GetQALength() const { return fQALength ; }
  Bool_t *   GetEventSpecies() const {return fEventSpecies;}	
  Int_t      GetESLength() const { return fESLength ; }
  
  //____________________________________________________//
 private:
  Int_t        fAliceRunId;              //the run id
  Float_t      fAliceMagneticField;      //value of the magnetic field
  Float_t      fAliceDipoleField;        //value of the magnetic field in dipole
  Int_t        fAliceRunStartTime;       //run start date
  Int_t        fAliceRunStopTime;        //run stop date
  TString      fAlirootVersion;          //aliroot version
  TString      fRootVersion;             //root version
  TString      fGeant3Version;           //geant3 version
  TString      fLHCPeriod;               //datataking period
  TString      fRecPass;                 //Reconstruction pass number
  TString      fProductionName;          //production name
  Bool_t       fAliceRunQuality;         //validation script
  Float_t      fAliceBeamEnergy;         //beam energy cm
  TString      fAliceBeamType;           //run type (pp, AA, pA)
  Int_t        fAliceCalibrationVersion; //calibration version  
  Int_t        fAliceDataType;           //0: simulation -- 1: data  
  Int_t        fNumEvents;               //number of events per file
  Int_t        fNumDetectors;            //number of detector configs per file
  TClonesArray fEventTag;                //array with all event tags
  AliDetectorTag fDetectorTag;           //array with all the detector tags
  AliLHCTag    fLHCTag;                  //LHC tag object
  AliQA        fQA;                      //needed for backward compaibility
  Int_t        fQALength;                // Length of the fQA array  
  ULong_t *    fQAArray ;                //[fQALength] QA objects's data	
  Int_t        fESLength;                // Length of the Event Specie Length
  Bool_t *     fEventSpecies;           //[fESLength] EventSpecies in this run	
  
  ClassDef(AliRunTag,7)  //(ClassName, ClassVersion)
};
//___________________________________________________________________________

#endif
