#ifndef ALITOFT0MAKER_H
#define ALITOFT0MAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFT0maker.h,v 1.8 2010/01/19 16:32:20 noferini Exp $ */

///////////////////////////////////////////////
//					     //
//  Manager class for time zero evaluation   //
//  with TOF informations                    //
//					     //
///////////////////////////////////////////////


#include "TObject.h"
#include "TString.h"
#include "AliESDEvent.h"
#include "AliStack.h"

class AliTOFcalibHisto;
class AliTOFT0v1;

class AliTOFT0maker : public TObject {
public:
  
  AliTOFT0maker() ;
  virtual ~AliTOFT0maker() ; // dtor
  AliTOFT0maker(const AliTOFT0maker & t);
  AliTOFT0maker & operator=(const AliTOFT0maker & t);
 
  void SetESDdata(Bool_t val=kTRUE){fESDswitch=val;};

  // return (fCalculated[0]=event time -- fCalculated[1]=sigma event time in ps -- fCalculated[2]=mean event time for each fill) if you can subtruct the event time; return NULL if there is no event time
  Double_t *RemakePID(AliESDEvent *esd,Double_t t0time=0.,Double_t t0sigma=1000.); // t0time and t0sigma in ps

  void      SetTimeResolution(Double_t timeresolution){fTimeResolution=timeresolution;};// TOF timeresolution in [s] e.g. for 120 ps -> 1.2e-10
  Double_t  GetTimeResolution() const {return fTimeResolution;}
  
 private:
  void TakeTimeRawCorrection(AliESDEvent * const esd);
  void RemakeTOFpid(AliESDEvent *esd,Float_t timezero);
  Double_t GetT0Fill(Int_t nrun) const ;

  AliTOFcalibHisto *fCalib; // TOF calibration object pointer

  Bool_t fESDswitch; // if you want take the ESD time instead of the raw + time slewing correction

  Double_t fCalculated[3]; // contains the parameters with the event time
  Double_t fTimeResolution;  // global time resolution used to calculate T0

  Float_t fT0sigma; // T0 resolution
  
  ClassDef(AliTOFT0maker,1);  // Calculate the time zero using TOF detector */
  
};

#endif // ALITOFT0MAKER_H
