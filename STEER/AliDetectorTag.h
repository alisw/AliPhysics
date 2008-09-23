#ifndef ALIDETECTORTAG_H
#define ALIDETECTORTAG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliDetectorTag
//   This is the class to deal with the tags for the detector level
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include "TObject.h"
#include "TObjArray.h"

//___________________________________________________________________________
class AliDetectorTag : public TObject {
 public:
  AliDetectorTag();
  AliDetectorTag(const AliDetectorTag & t);

  AliDetectorTag &operator=(const AliDetectorTag &rhs);
  virtual ~AliDetectorTag();
  
  //____________________________________________________//
  void SetDetectorMask(UInt_t mask) {fMask = mask; Int2Bin();}
  TObjArray *GetDetectorMask() {return fDetectorArray;}
  UInt_t GetIntDetectorMask();
  void PrintDetectorMask();

  //____________________________________________________//
  Bool_t GetITSSPD() const {return fITSSPD;}
  Bool_t GetITSSDD() const {return fITSSDD;}
  Bool_t GetITSSSD() const {return fITSSSD;}
  Bool_t GetTPC() const {return fTPC;}
  Bool_t GetTRD() const {return fTRD;}
  Bool_t GetTOF() const {return fTOF;}
  Bool_t GetHMPID() const {return fHMPID;}
  Bool_t GetPHOS() const {return fPHOS;}
  Bool_t GetPMD() const {return fPMD;}
  Bool_t GetMUON() const {return fMUON;}
  Bool_t GetFMD() const {return fFMD;}
  Bool_t GetTZERO() const {return fTZERO;}
  Bool_t GetVZERO() const {return fVZERO;}
  Bool_t GetZDC() const {return fZDC;}
  Bool_t GetEMCAL() const {return fEMCAL;}
  
  //____________________________________________________//
 private:
  void Int2Bin();
  void SetDetectorConfiguration();

  void SetITSSPD() {fITSSPD = kTRUE;}
  void SetITSSDD() {fITSSDD = kTRUE;}
  void SetITSSSD() {fITSSSD = kTRUE;}
  void SetTPC() {fTPC = kTRUE;}
  void SetTRD() {fTRD = kTRUE;}
  void SetTOF() {fTOF = kTRUE;}
  void SetHMPID() {fHMPID = kTRUE;}
  void SetPHOS() {fPHOS = kTRUE;}
  void SetPMD() {fPMD = kTRUE;}
  void SetMUON() {fMUON = kTRUE;}
  void SetFMD() {fFMD = kTRUE;}
  void SetTZERO() {fTZERO = kTRUE;}
  void SetVZERO() {fVZERO = kTRUE;}
  void SetZDC() {fZDC = kTRUE;}
  void SetEMCAL() {fEMCAL = kTRUE;}
  
  TObjArray *fDetectorArray; //detectors' names - active
  UInt_t     fMask;          //detector mask
  UInt_t     fDetectors[32]; //detector mask
  Bool_t     fITSSPD;        //ITS-SPD active = 1
  Bool_t     fITSSDD;        //ITS-SDD active = 1
  Bool_t     fITSSSD;        //ITS-SSD active = 1
  Bool_t     fTPC;           //TPC active = 1
  Bool_t     fTRD;           //TRD active = 1
  Bool_t     fTOF;           //TOF active = 1
  Bool_t     fHMPID;         //HMPID active = 1
  Bool_t     fPHOS;          //PHOS active = 1
  Bool_t     fPMD;           //PMD active = 1
  Bool_t     fMUON;          //MUON active = 1
  Bool_t     fFMD;           //FMD active = 1
  Bool_t     fTZERO;         //TZERO active = 1
  Bool_t     fVZERO;         //VZERO active = 1
  Bool_t     fZDC;           //ZDC active = 1
  Bool_t     fEMCAL;         //EMCAL active = 1

  ClassDef(AliDetectorTag,4)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
