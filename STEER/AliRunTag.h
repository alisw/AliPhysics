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

class AliEventTag;
class AliDetectorTag;


//______________________________________________________________________________
class AliRunTag : public TObject
{
 public:
  AliRunTag();
  virtual ~AliRunTag();
  
  void          SetRunId(Int_t Pid) {fAliceRunId = Pid;}
  void          SetMagneticField(Float_t Pmag) {fAliceMagneticField = Pmag;}
  void          SetRunStartTime(Int_t Pt0) {fAliceRunStartTime = Pt0;}
  void          SetRunStopTime(Int_t Pt1) {fAliceRunStopTime = Pt1;}
  void          SetRecoVersion(Int_t Pn) {fAliceReconstructionVersion = Pn;}
  void          SetRunQuality(Int_t Pn) {fAliceRunQuality = Pn;}
  void          SetBeamEnergy(Float_t PE) {fAliceBeamEnergy = PE;}
  void          SetBeamType(const char *Ptype) {fAliceBeamType = Ptype;}
  void          SetCalibVersion(Int_t Pn) {fAliceCalibrationVersion = Pn;}
  
  void          SetDataType(Int_t i) {fAliceDataType = i;}
  
  void          SetNEvents(Int_t Pn) { fNumEvents = Pn; }
  
  void          SetLHCTag(Float_t Plumin, char *type);
  void          SetDetectorTag(const AliDetectorTag &t);
  void          AddEventTag(const AliEventTag &t);
  void          Clear(const char * opt = "");
  
  
  Int_t         GetRunId() const {return fAliceRunId;}
  Float_t       GetMagneticField() const {return fAliceMagneticField;}
  Int_t         GetRunStartTime() const {return fAliceRunStartTime;}
  Int_t         GetRunStopTime() const {return fAliceRunStopTime;}
  Int_t         GetRecoVersion() const {return fAliceReconstructionVersion;}
  Int_t         GetRunQuality() const {return fAliceRunQuality;}
  Float_t       GetBeamEnergy() const {return fAliceBeamEnergy;}
  const char   *GetBeamType() const {return fAliceBeamType.Data();}
  Int_t         GetCalibVersion() const {return fAliceCalibrationVersion;}
  
  Int_t GetDataType() const {return fAliceDataType;}

  Int_t         GetNEvents() const {return fNumEvents;}
  
  AliLHCTag    *GetLHCTag() {return &fLHCTag; } 
  const TClonesArray *GetEventTags() const {return &fEventTag;}
  const TClonesArray *GetDetectorTags() const {return &fDetectorTag;}

 private:
  Int_t    fAliceRunId;                   //the run id
  Float_t  fAliceMagneticField;           //value of the magnetic field
  Int_t    fAliceRunStartTime;            //run start date
  Int_t    fAliceRunStopTime;             //run stop date
  Int_t    fAliceReconstructionVersion;   //reco version
  Bool_t   fAliceRunQuality;              //validation script
  Float_t  fAliceBeamEnergy;              //beam energy cm
  TString  fAliceBeamType;                //run type (pp, AA, pA)
  Int_t    fAliceCalibrationVersion;      //calibration version
  
  Int_t  fAliceDataType;                  //0: simulation -- 1: data
  
  Int_t    fNumEvents;                    //number of events per file
  Int_t    fNumDetectors;                 //number of detector configs per file
  TClonesArray  fEventTag;                //array with all event tags
  TClonesArray  fDetectorTag;             //array with all the detector tags
  
  AliLHCTag   fLHCTag;                    //LHC tag object
  
  ClassDef(AliRunTag,1)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
