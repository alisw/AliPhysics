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
  void SetTOFres(Float_t res){fTOFres=res;}
  void SetIsMC(Bool_t flag=kFALSE){fIsMC=flag;}
  void SetCorrectExpTimes(Bool_t flag=kTRUE){fCorrectExpTimes=flag;}
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
  virtual void SetTimeZeroType(AliESDpid::EStartTimeType_t tofTimeZeroType) {fTimeZeroType = tofTimeZeroType;}

  /* theoretical expected time: related stuff for LHC10d patch */
  static Float_t GetBetaTh(Float_t m, Float_t p) {return TMath::Sqrt(1. / (1. + m * m / (p * p)));}; // get beta th
  static Float_t GetExpTimeTh(Float_t m, Float_t p, Float_t L) {return L / 2.99792457999999984e-02 / GetBetaTh(m, p);}; // get exp time th
  void RecomputeTExp(AliESDEvent *event) const;
  void RecomputeTExp(AliESDtrack *track) const;

private:
  AliESDpid          *fESDpid;         //! ESD pid object

  Bool_t fIsMC;              // flag for MC data
  Int_t  fTimeZeroType;      // flag to select timeZero type 
  Bool_t fCorrectExpTimes;   // flag to apply Expected Time correction 
  Bool_t fLHC10dPatch;       // flag to apply special patch for LHC10d (reconstructed with wrong geometry)
  Int_t  fDebugLevel;        // debug purposes 0= no output, 1 Info, 2 lot of info....
  Bool_t fAutomaticSettings; // enable/disable automatic (per run) settings

  // variables for TOF calibrations and timeZero setup
  AliTOFcalib     *fTOFCalib;       // recalibrate TOF signal with OCDB
  AliTOFT0maker   *fTOFT0maker;     // computation of TOF-T0
  Float_t fTOFres;                  // TOF resolution
  Float_t fT0shift[4];              // T0 detector correction from OCDB
  Float_t fT0IntercalibrationShift; // extra-shift to adjust TOF/TO intercalibration issue in some period

  // variables to parametrize MC
  static Float_t fgT0Aresolution;   // T0 resolution A-Side (MC)
  static Float_t fgT0Cresolution;   // T0 resolution C-Side (MC)


  AliTOFTenderSupply(const AliTOFTenderSupply&c);
  AliTOFTenderSupply& operator= (const AliTOFTenderSupply&c);

  ClassDef(AliTOFTenderSupply, 6);
};


#endif 
