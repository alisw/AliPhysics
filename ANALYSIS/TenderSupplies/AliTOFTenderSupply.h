#ifndef ALITOFTENDERSUPPLY_H
#define ALITOFTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  TOF tender, reapply pid on the fly                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <AliTenderSupply.h>
#include <AliLog.h>
#include <AliESDpid.h>

class AliESDpid;
class AliTOFcalib;
class AliTOFT0maker;
class AliESDEevent;
class AliESDtrack;
class AliTOFTenderSupply: public AliTenderSupply {

public:
  AliTOFTenderSupply();
  AliTOFTenderSupply(const char *name, const AliTender *tender=NULL);

  virtual ~AliTOFTenderSupply(){;}

  virtual void              Init();
  virtual void              ProcessEvent();

  // TOF tender methods
  void SetIsMC(Bool_t flag=kFALSE){fIsMC=flag;}
  void SetCorrectExpTimes(Bool_t flag=kTRUE){fCorrectExpTimes=flag;}
  void SetT0DetectorAdjust(Bool_t flag=kFALSE){fT0DetectorAdjust=flag;}
  void SetDebugLevel(Int_t flag=0){fDebugLevel=flag;}
  void SetLHC10dPatch(Bool_t flag=kFALSE){ 
    if (flag == kTRUE) {
      AliInfo(" **** TOF Tender: special setting LHC10d patch is ON");
      AliInfo(" **** TOF Tender: this setting is valid only on LHC10d pass2");
    }
    fLHC10dPatch=flag;
    return;
  }
  void SetAutomaticSettings(Bool_t flag=kTRUE){fAutomaticSettings=flag;}
  void SetForceCorrectTRDBug(Bool_t flag=kTRUE){fForceCorrectTRDBug=flag;}
  void SetUserRecoPass(Int_t flag=0){fUserRecoPass=flag;}
  Int_t GetRecoPass(void){return fRecoPass;}
  void DetectRecoPass();

  /* theoretical expected time: related stuff for LHC10d patch */
  static Float_t GetBetaTh(Float_t m, Float_t p) {return TMath::Sqrt(1. / (1. + m * m / (p * p)));}; // get beta th
  static Float_t GetExpTimeTh(Float_t m, Float_t p, Float_t L) {return L / 2.99792457999999984e-02 / GetBetaTh(m, p);}; // get exp time th
  void RecomputeTExp(AliESDEvent *event) const;
  void RecomputeTExp(AliESDtrack *track) const;
  void FixTRDBug(AliESDEvent *event);
  void FixTRDBug(AliESDtrack *track);
  void InitGeom();
  void FindTRDFix(AliESDtrack *track,Double_t *corr);
  Double_t EstimateLengthInTRD1(AliESDtrack *track);
  Double_t EstimateLengthInTRD2(AliESDtrack *track);
  Double_t EstimateLengthOutTRD(AliESDtrack *track);
  void CorrectDeltaTimes(Double_t pT, Double_t length, Bool_t isTRDout, Double_t *corrections);
  Double_t CorrectExpectedProtonTime(Double_t pT,Double_t length, Bool_t isTRDout);
  Double_t CorrectExpectedKaonTime(Double_t pT,Double_t length, Bool_t isTRDout);
  Double_t CorrectExpectedPionTime(Double_t pT,Double_t length, Bool_t isTRDout);
  Int_t GetOCDBVersion(Int_t runNumber);
  void LoadTOFPIDParams(Int_t runNumber);

private:
  AliESDpid          *fESDpid;         //! ESD pid object

  
  Bool_t fTenderNoAction;    // flag for periods when tender action is not requested/not supported
  Bool_t fIsMC;              // flag for MC data
  Bool_t fCorrectExpTimes;   // flag to apply Expected Time correction 
  Bool_t fCorrectTRDBug;     // flag to fix wrong dE/dx inside TRD
  Bool_t fLHC10dPatch;       // flag to apply special patch for LHC10d (reconstructed with wrong geometry)
  Bool_t fT0DetectorAdjust;  // flag to apply offsets to T0 data (works only on some periods)
  Int_t  fDebugLevel;        // debug purposes 0= no output, 1 Info, 2 lot of info....
  Bool_t fAutomaticSettings; // enable/disable automatic (per run) settings
  Int_t  fRecoPass;          // reconstruction pass: the tender applies different recipes depending on the pass
  Int_t  fUserRecoPass;      // when reco pass is selected by user
  Bool_t fForceCorrectTRDBug; // force TRD bug correction (for some bad MC production...)


  // variables for TOF calibrations and timeZero setup
  AliTOFPIDParams *fTOFPIDParams;   //! TOF PID Params - period depending (OADB loaded)
  AliTOFcalib     *fTOFCalib;       // recalibrate TOF signal with OCDB
  AliTOFT0maker   *fTOFT0maker;     // computation of TOF-T0
  Float_t fT0shift[4];              // T0 detector correction from OCDB
  Float_t fT0IntercalibrationShift; // extra-shift to adjust TOF/TO intercalibration issue in some period

  // variables to parametrize MC
  static Float_t fgT0Aresolution;   // T0 resolution A-Side (MC)
  static Float_t fgT0Cresolution;   // T0 resolution C-Side (MC)

  // variables to steer TRD bug fix
  Bool_t fGeomSet;                 // steer loading GRP entry
  Bool_t fIsEnteringInTRD;
  Bool_t fInTRD;
  Bool_t fIsComingOutTRD;
  Bool_t fOutTRD;
  Float_t fRhoTRDin;                // cm
  Float_t fRhoTRDout;               // cm
  Float_t fStep;                    // cm
  Double_t fMagField;               // magnetic field value [kGauss]
  ULong64_t fCDBkey;

  AliTOFTenderSupply(const AliTOFTenderSupply&c);
  AliTOFTenderSupply& operator= (const AliTOFTenderSupply&c);

  ClassDef(AliTOFTenderSupply, 11);
};


#endif 
