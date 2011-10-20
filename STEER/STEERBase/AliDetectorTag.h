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
#include "AliDAQ.h"

//___________________________________________________________________________
class AliDetectorTag : public TObject {
 public:
  AliDetectorTag();
  AliDetectorTag(const AliDetectorTag & t);

  AliDetectorTag &operator=(const AliDetectorTag &rhs);
  virtual ~AliDetectorTag();
  
  void UpdateFromRunTable(AliDetectorTag &detTag);

  //____________________________________________________//
  void SetDetectorMask(UInt_t mask)     {fMaskDAQ = mask; fMaskReco = mask; }
  void SetDetectorMaskDAQ(UInt_t mask)  {fMaskDAQ = mask;}
  void SetDetectorMaskReco(UInt_t mask) {fMaskReco = mask;}
  void SetDetectorValidityRange(UChar_t idet, UShort_t vr) {fDetectorValidityRange[idet] = vr; }
  void SetDetectorStatus(UChar_t idet, TString co) { fDetectorStatus[idet] = co; }
  TObjArray *GetDetectorMask() { return 0; } // {return fDetectorArray;}
  UInt_t GetIntDetectorMask() { return fMaskDAQ; }
  UInt_t GetIntDetectorMaskDAQ() { return fMaskDAQ; }
  UInt_t GetIntDetectorMaskReco() { return fMaskReco; }
  UShort_t GetDetectorValidityRange(UChar_t idet) const { return fDetectorValidityRange[idet]; }
  TString  GetDetectorStatus(UChar_t idet) const { return fDetectorStatus[idet]; }
  const char *GetDetectorMaskDAQ() { return AliDAQ::ListOfTriggeredDetectors(fMaskDAQ); }
  const char *GetDetectorMaskReco() { return AliDAQ::ListOfTriggeredDetectors(fMaskReco); }
  void PrintDetectorMask();

  //____________________________________________________//
  Bool_t GetITSSPD() const {return fMaskDAQ & AliDAQ::kSPD;}
  Bool_t GetITSSDD() const {return fMaskDAQ & AliDAQ::kSSD;}
  Bool_t GetITSSSD() const {return fMaskDAQ & AliDAQ::kSSD;}
  Bool_t GetTPC()    const {return fMaskDAQ & AliDAQ::kTPC;}
  Bool_t GetTRD()    const {return fMaskDAQ & AliDAQ::kTRD;}
  Bool_t GetTOF()    const {return fMaskDAQ & AliDAQ::kTOF;}
  Bool_t GetHMPID()  const {return fMaskDAQ & AliDAQ::kHMPID;}
  Bool_t GetPHOS()   const {return fMaskDAQ & AliDAQ::kPHOS;}
  Bool_t GetPMD()    const {return fMaskDAQ & AliDAQ::kPMD;}
  Bool_t GetMUON()   const {return fMaskDAQ & AliDAQ::kMUON;}
  Bool_t GetFMD()    const {return fMaskDAQ & AliDAQ::kFMD;}
  Bool_t GetTZERO()  const {return fMaskDAQ & AliDAQ::kT0;}
  Bool_t GetVZERO()  const {return fMaskDAQ & AliDAQ::kVZERO;}
  Bool_t GetZDC()    const {return fMaskDAQ & AliDAQ::kZDC;}
  Bool_t GetEMCAL()  const {return fMaskDAQ & AliDAQ::kEMCAL;}
//   #ifdef MFT_UPGRADE
//   Bool_t GetMFT()    const {return fMaskDAQ & AliDAQ::kMFT;}
//   #endif
  Bool_t GetMFT()    const {return fMaskDAQ & AliDAQ::kMFT;}   // AU
  //____________________________________________________//
 private:
  //  void Int2Bin();
  //   void SetDetectorConfiguration();

  void SetITSSPD() {fMaskDAQ |= AliDAQ::kSPD  ;}
  void SetITSSDD() {fMaskDAQ |= AliDAQ::kSDD  ;}
  void SetITSSSD() {fMaskDAQ |= AliDAQ::kSSD  ;}
  void SetTPC()    {fMaskDAQ |= AliDAQ::kTPC  ;}
  void SetTRD()    {fMaskDAQ |= AliDAQ::kTRD  ;}
  void SetTOF()    {fMaskDAQ |= AliDAQ::kTOF  ;}
  void SetHMPID()  {fMaskDAQ |= AliDAQ::kHMPID;}
  void SetPHOS()   {fMaskDAQ |= AliDAQ::kPHOS ;}
  void SetPMD()    {fMaskDAQ |= AliDAQ::kPMD  ;}
  void SetMUON()   {fMaskDAQ |= AliDAQ::kMUON ;}
  void SetFMD()    {fMaskDAQ |= AliDAQ::kFMD  ;}
  void SetTZERO()  {fMaskDAQ |= AliDAQ::kT0   ;}
  void SetVZERO()  {fMaskDAQ |= AliDAQ::kVZERO;}
  void SetZDC()    {fMaskDAQ |= AliDAQ::kZDC  ;}
  void SetEMCAL()  {fMaskDAQ |= AliDAQ::kEMCAL;}
//   #ifdef MFT_UPGRADE
//   void SetMFT()    {fMaskDAQ |= AliDAQ::kMFT;}
//   #endif
  void SetMFT()    {fMaskDAQ |= AliDAQ::kMFT;}   // AU
	
  //   TObjArray *fDetectorArray; //detectors' names - active
  UInt_t     fMaskDAQ;          //detector mask in DAQ
  UInt_t     fMaskReco;         //detector mask in Reco
  //   UInt_t     fDetectors[32]; //detector mask
  //   Bool_t     fITSSPD;        //ITS-SPD active = 1
  //   Bool_t     fITSSDD;        //ITS-SDD active = 1
  //   Bool_t     fITSSSD;        //ITS-SSD active = 1
  //   Bool_t     fTPC;           //TPC active = 1
  //   Bool_t     fTRD;           //TRD active = 1
  //   Bool_t     fTOF;           //TOF active = 1
  //   Bool_t     fHMPID;         //HMPID active = 1
  //   Bool_t     fPHOS;          //PHOS active = 1
  //   Bool_t     fPMD;           //PMD active = 1
  //   Bool_t     fMUON;          //MUON active = 1
  //   Bool_t     fFMD;           //FMD active = 1
  //   Bool_t     fTZERO;         //TZERO active = 1
  //   Bool_t     fVZERO;         //VZERO active = 1
  //   Bool_t     fZDC;           //ZDC active = 1
  //   Bool_t     fEMCAL;         //EMCAL active = 1

  UShort_t   fDetectorValidityRange[AliDAQ::kHLTId];
  TString    fDetectorStatus[AliDAQ::kHLTId];

  ClassDef(AliDetectorTag, 6)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
