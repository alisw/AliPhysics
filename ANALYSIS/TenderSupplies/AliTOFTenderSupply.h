#ifndef ALITOFTENDERSUPPLY_H
#define ALITOFTENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  TPC tender, reapply pid on the fly                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <AliTenderSupply.h>
#include "AliESDpid.h"

class AliESDpid;
class AliTOFcalib;
class AliTOFT0maker;

class AliTOFTenderSupply: public AliTenderSupply {

public:
  AliTOFTenderSupply();
  AliTOFTenderSupply(const char *name, const AliTender *tender=NULL);

  virtual ~AliTOFTenderSupply(){;}

  virtual void              Init();
  virtual void              ProcessEvent();

  // TOF method
  void SetTOFres(Float_t res){fTOFres=res;}
  void SetApplyT0(Bool_t flag=kTRUE){fApplyT0=flag;};
  void SetCorrectExpTimes(Bool_t flag=kTRUE){fCorrectExpTimes=flag;};
  
  virtual void SetTimeZeroType(AliESDpid::EStartTimeType_t tofTimeZeroType) {fTimeZeroType = tofTimeZeroType;}

private:
  AliESDpid          *fESDpid;         //! ESD pid object

  Bool_t fIsMC;              // flag for MC data
  Bool_t fApplyT0;           // flag to subtract the T0-TOF (deprecated)
  Int_t  fTimeZeroType;      // flag to discriminate the time zero type 
  Bool_t fCorrectExpTimes;   // flag to apply Expected Time correction 

  // variables for TOF calibrations
  AliTOFcalib     *fTOFCalib;    //! recalibrate TOF signal with OCDB
  AliTOFT0maker   *fTOFT0maker;     //! TOF maker objects (apply all the correction for T0)

  Float_t fTOFres;                   // TOF resolution
  Float_t fT0shift[4];               // T0 detector correction from OCDB

  AliTOFTenderSupply(const AliTOFTenderSupply&c);
  AliTOFTenderSupply& operator= (const AliTOFTenderSupply&c);

  ClassDef(AliTOFTenderSupply, 2);
};


#endif 

