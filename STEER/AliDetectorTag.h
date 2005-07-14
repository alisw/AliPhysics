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
#include "TClonesArray.h"
//______________________________________________________________________________
class AliDetectorTag : public TObject
{
 private:
  Bool_t   fITS;		//ITS active = 1
  Bool_t   fTPC;		//TPC active = 1
  Bool_t   fTRD;		//TRD active = 1
  Bool_t   fTOF;		//TOF active = 1
  Bool_t   fHMPID;		//HMPID active = 1
  Bool_t   fPHOS;		//PHOS active = 1
  Bool_t   fZDC;		//ZDC active = 1
  Bool_t   fMUON;		//MUON active = 1
  Bool_t   fABSORBER;		//ABSORBER active = 1
  Bool_t   fPMD;		//PMD active = 1
  Bool_t   fRICH;		//RICH active = 1
  Bool_t   fEMCAL;		//EMCAL active = 1
  Bool_t   fVZERO;		//VZERO active = 1
  Bool_t   fTZERO;		//TZERO active = 1

 public:
  AliDetectorTag();
  AliDetectorTag(AliDetectorTag *t);
  virtual ~AliDetectorTag();
  
  void          SetITS(Int_t n) {fITS = n;}
  void          SetTPC(Int_t n) {fTPC = n;}
  void          SetTRD(Int_t n) {fTRD = n;}
  void          SetTOF(Int_t n) {fTOF = n;}
  void          SetHMPID(Int_t n) {fHMPID = n;}
  void          SetPHOS(Int_t n) {fPHOS = n;}
  void          SetZDC(Int_t n) {fZDC = n;}
  void          SetMUON(Int_t n) {fMUON = n;}
  void          SetABSORBER(Int_t n) {fABSORBER = n;}
  void          SetPMD(Int_t n) {fPMD = n;}
  void          SetRICH(Int_t n) {fRICH = n;}
  void          SetEMCAL(Int_t n) {fEMCAL = n;}
  void          SetVZERO(Int_t n) {fVZERO = n;}
  void          SetTZERO(Int_t n) {fTZERO = n;}
  
  Bool_t        GetITS() {return fITS;}
  Bool_t        GetTPC() {return fTPC;}
  Bool_t        GetTRD() {return fTRD;}
  Bool_t        GetTOF() {return fTOF;}
  Bool_t        GetHMPID() {return fHMPID;}
  Bool_t        GetPHOS() {return fPHOS;}
  Bool_t        GetZDC() {return fZDC;}
  Bool_t        GetMUON() {return fMUON;}
  Bool_t        GetABSORBER() {return fABSORBER;}
  Bool_t        GetPMD() {return fPMD;}
  Bool_t        GetRICH() {return fRICH;}
  Bool_t        GetEMCAL() {return fEMCAL;}
  Bool_t        GetVZERO() {return fVZERO;}
  Bool_t        GetTZERO() {return fTZERO;}
  
  ClassDef(AliDetectorTag,1)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
