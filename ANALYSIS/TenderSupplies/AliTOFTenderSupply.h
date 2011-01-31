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
#include "AliESDpid.h"

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

  // TOF method
  void SetTOFres(Float_t res){fTOFres=res;}
  void SetApplyT0(Bool_t flag=kTRUE){fApplyT0=flag;}
  void SetCorrectExpTimes(Bool_t flag=kTRUE){fCorrectExpTimes=flag;}
  void SetLHC10dPatch(Bool_t flag=kFALSE){ 
    if (flag == kTRUE) {
      Print(" **** TOF Tender: special setting LHC10d patch is ON");
      Print(" **** TOF Tender: this setting is valid only on LHC10d pass2");
    }
    fLHC10dPatch=flag;
    return;
  }
  
  virtual void SetTimeZeroType(AliESDpid::EStartTimeType_t tofTimeZeroType) {fTimeZeroType = tofTimeZeroType;}

  /* theoretical expected time related stuff for LHC10d patch */
  static Float_t GetBetaTh(Float_t m, Float_t p) {return TMath::Sqrt(1. / (1. + m * m / (p * p)));}; // get beta th
  static Float_t GetExpTimeTh(Float_t m, Float_t p, Float_t L) {return L / 2.99792457999999984e-02 / GetBetaTh(m, p);}; // get exp time th
  void RecomputeTExp(AliESDEvent *event) const;
  void RecomputeTExp(AliESDtrack *track) const;

private:
  AliESDpid          *fESDpid;         //! ESD pid object

  Bool_t fIsMC;              // flag for MC data
  Bool_t fApplyT0;           // flag to subtract the T0-TOF (deprecated)
  Int_t  fTimeZeroType;      // flag to discriminate the time zero type 
  Bool_t fCorrectExpTimes;   // flag to apply Expected Time correction 
  Bool_t fLHC10dPatch;       // flag to apply special patch for LHC10d (reconstructed with wrong geometry)

  // variables for TOF calibrations
  AliTOFcalib     *fTOFCalib;    //! recalibrate TOF signal with OCDB
  AliTOFT0maker   *fTOFT0maker;     //! TOF maker objects (apply all the correction for T0)

  Float_t fTOFres;                   // TOF resolution
  Float_t fT0shift[4];               // T0 detector correction from OCDB

  AliTOFTenderSupply(const AliTOFTenderSupply&c);
  AliTOFTenderSupply& operator= (const AliTOFTenderSupply&c);

  ClassDef(AliTOFTenderSupply, 3);
};


#endif 
