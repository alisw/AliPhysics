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

class AliESDpid;
class AliTOFcalibESD;
class AliTOFT0makerANA;

class AliTOFTenderSupply: public AliTenderSupply {
  
public:
  AliTOFTenderSupply();
  AliTOFTenderSupply(const char *name, const AliTender *tender=NULL);
  
  virtual ~AliTOFTenderSupply(){;}
  
  virtual void              Init();
  virtual void              ProcessEvent();

  // TOF method
  void SetTOFres(Float_t res){fTOFres=res;}
  
private:
  AliESDpid          *fESDpid;         //! ESD pid object
  
  // variables for TOF calibrations
  AliTOFcalibESD     *fTOFesdCalib;    //! recalibrate TOF signal with OCDB
  AliTOFT0makerANA   *fTOFT0maker;     //! TOF maker objects (apply all the correction for T0)

  Float_t fTOFres;                   // TOF resolution
  
  AliTOFTenderSupply(const AliTOFTenderSupply&c);
  AliTOFTenderSupply& operator= (const AliTOFTenderSupply&c);
  
  ClassDef(AliTOFTenderSupply, 1);  // TPC tender task
};


#endif

