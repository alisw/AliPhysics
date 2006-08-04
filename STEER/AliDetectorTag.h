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
//___________________________________________________________________________
class AliDetectorTag : public TObject {
 public:
  AliDetectorTag();
  AliDetectorTag(const AliDetectorTag & t);

  AliDetectorTag &operator=(const AliDetectorTag &rhs);
  virtual ~AliDetectorTag();
  
  //____________________________________________________//
  void SetITS() {fITS = kTRUE;}
  void SetTPC() {fTPC = kTRUE;}
  void SetTRD() {fTRD = kTRUE;}
  void SetTOF() {fTOF = kTRUE;}
  void SetHMPID() {fHMPID = kTRUE;}
  void SetPHOS() {fPHOS = kTRUE;}
  void SetZDC() {fZDC = kTRUE;}
  void SetMUON() {fMUON = kTRUE;}
  void SetPMD() {fPMD = kTRUE;}
  void SetEMCAL() {fEMCAL = kTRUE;}
  void SetVZERO() {fVZERO = kTRUE;}
  void SetTZERO() {fTZERO = kTRUE;}
  
  //____________________________________________________//
  Bool_t GetITS() const {return fITS;}
  Bool_t GetTPC() const {return fTPC;}
  Bool_t GetTRD() const {return fTRD;}
  Bool_t GetTOF() const {return fTOF;}
  Bool_t GetHMPID() const {return fHMPID;}
  Bool_t GetPHOS() const {return fPHOS;}
  Bool_t GetZDC() const {return fZDC;}
  Bool_t GetMUON() const {return fMUON;}
  Bool_t GetPMD() const {return fPMD;}
  Bool_t GetEMCAL() const {return fEMCAL;}
  Bool_t GetVZERO() const {return fVZERO;}
  Bool_t GetTZERO() const {return fTZERO;}
  
  //____________________________________________________//
 private:
  Bool_t   fITS;      //ITS active = 1
  Bool_t   fTPC;      //TPC active = 1
  Bool_t   fTRD;      //TRD active = 1
  Bool_t   fTOF;      //TOF active = 1
  Bool_t   fHMPID;    //HMPID active = 1
  Bool_t   fPHOS;     //PHOS active = 1
  Bool_t   fZDC;      //ZDC active = 1
  Bool_t   fMUON;     //MUON active = 1
  Bool_t   fPMD;      //PMD active = 1
  Bool_t   fEMCAL;    //EMCAL active = 1
  Bool_t   fVZERO;    //VZERO active = 1
  Bool_t   fTZERO;    //TZERO active = 1

  ClassDef(AliDetectorTag,2)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
