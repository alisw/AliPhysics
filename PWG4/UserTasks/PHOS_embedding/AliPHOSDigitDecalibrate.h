#ifndef ALIPHOSDIGITDECALIBRATE_H
#define ALIPHOSDIGITDECALIBRATE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  PHOS tender, apply corrections to PHOS clusters                   //
//  and do track matching                                             //
//  Author : Dmitri Peressounko (RRC KI)                              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include <AliTenderSupply.h>

class TVector3;
class AliPHOSGeometry;
class AliPHOSCalibData ; 
class AliPHOSDigitDecalibrate: public AliTenderSupply {
  
public:
  AliPHOSDigitDecalibrate();
  AliPHOSDigitDecalibrate(const char *name, const AliTender *tender=NULL);
  virtual ~AliPHOSDigitDecalibrate();

  virtual void   Init();
  virtual void   ProcessEvent();
  
  void  SetDecalibration(Int_t mod, TH2F * dec) ;

protected:
  AliPHOSDigitDecalibrate(const AliPHOSDigitDecalibrate&c);
  AliPHOSDigitDecalibrate& operator= (const AliPHOSDigitDecalibrate&c);
private:

  TH2F * hDec[5] ;                           //! Decalibration coeff.
  AliPHOSGeometry   *fPHOSGeo;               //! PHOS geometry
  AliPHOSCalibData *fPHOSCalibData;          //! PHOS calibration object

 
  ClassDef(AliPHOSDigitDecalibrate, 1); // PHOS tender task
};


#endif

