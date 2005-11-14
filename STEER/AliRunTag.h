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

#include <stdlib.h>
#include <Riostream.h>

#include "TObject.h"
#include "TClonesArray.h"

#include "AliEventTag.h"
#include "AliLHCTag.h"
#include "AliDetectorTag.h"


//______________________________________________________________________________
class AliRunTag : public TObject
{
 private:
  Int_t    fAliceRunId;                   //the run id
  Float_t  fAliceMagneticField;           //value of the magnetic field
  Int_t    fAliceRunStartTime;            //run start date
  Int_t    fAliceRunStopTime;             //run stop date
  Int_t    fAliceReconstructionVersion;   //reco version
  Bool_t   fAliceRunQuality;              //validation script
  Float_t  fAliceBeamEnergy;              //beam energy cm
  Char_t   fAliceBeamType[5];             //run type (pp, AA, pA)
  Int_t    fAliceCalibrationVersion;      //calibration version
  
  Int_t  fAliceDataType;              //0: simulation -- 1: data
  
  Int_t    fNumEvents;                    //number of events per file
  Int_t    fNumDetectors;                 //number of detector configs per file
  TClonesArray  *fEventTag;               //array with all event tags
  TClonesArray  *fDetectorTag;            //array with all the detector tags
  
  AliLHCTag   fLHCTag;
  
  static TClonesArray *fgEvents;
  static TClonesArray *fgDetectors;
  
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
  void          SetBeamType(char *Ptype) {strcpy(fAliceBeamType,Ptype);}
  void          SetCalibVersion(Int_t Pn) {fAliceCalibrationVersion = Pn;}
  
  void          SetDataType(Int_t i) {fAliceDataType = i;}
  
  void          SetNEvents(Int_t Pn) { fNumEvents = Pn; }
  
  void          SetLHCTag(Float_t Plumin, char *type);
  void          SetDetectorTag(AliDetectorTag *t);
  void          AddEventTag(const AliEventTag &t);
  void          Clear(const char * opt = "");
  
  
  Int_t         GetRunId() {return fAliceRunId;}
  Float_t       GetMagneticField() {return fAliceMagneticField;}
  Int_t         GetRunStartTime() {return fAliceRunStartTime;}
  Int_t         GetRunStopTime() {return fAliceRunStopTime;}
  Int_t         GetRecoVersion() {return fAliceReconstructionVersion;}
  Int_t         GetRunQuality() {return fAliceRunQuality;}
  Float_t       GetBeamEnergy() {return fAliceBeamEnergy;}
  char         *GetBeamType() {return fAliceBeamType;}
  Int_t         GetCalibVersion() {return fAliceCalibrationVersion;}
  
  Int_t GetDataType() {return fAliceDataType;}

  Int_t         GetNEvents() const {return fNumEvents;}
  
  AliLHCTag    *GetLHCTag() { return &fLHCTag; } 
  TClonesArray *GetEventTags() const {return fEventTag;}

  ClassDef(AliRunTag,1)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
