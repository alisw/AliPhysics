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

class TH1F;
class AliESDEvent;

class AliTOFcalibHisto;
class AliTOFT0v1;

class AliTOFT0maker : public TObject {
public:
  
  AliTOFT0maker() ;
  virtual ~AliTOFT0maker() ; // dtor
  AliTOFT0maker(const AliTOFT0maker & t);
  AliTOFT0maker & operator=(const AliTOFT0maker & t);
 
  void SetESDdata(Bool_t val=kTRUE){fESDswitch=val;};

  // return (fCalculated[0]=event time -- fCalculated[1]=sigma event time in ps -- fCalculated[2]=mean event time for each fill -- fCalculated[3]=number of tracks at the TOF level) if you can subtruct the event time; return NULL if there is no event time
  Double_t *RemakePID(AliESDEvent *esd,Double_t t0time=0.,Double_t t0sigma=1000.); // t0time and t0sigma in ps

  void      SetTimeResolution(Double_t timeresolution){fTimeResolution=timeresolution;};// TOF timeresolution in [s] e.g. for 120 ps -> 1.2e-10
  Double_t  GetTimeResolution() const {return fTimeResolution;}
  
  void LoadChannelMap(char *filename="$ALICE_ROOT/TOF/enableMap.104892.root"); //load the enable channel map
  void ApplyMask(AliESDEvent * const esd);
  
  void SetNoTOFT0(Bool_t status=kTRUE){fNoTOFT0=status;}; // disable the TOF T0 info
  void SetMaskOffChannel(Bool_t status=kTRUE){fKmask=status;}; // swith for the map off channel
  
 private:
  void TakeTimeRawCorrection(AliESDEvent * const esd);
  void RemakeTOFpid(AliESDEvent *esd,Float_t timezero);
  Double_t GetT0Fill() const;
  
  AliTOFcalibHisto *fCalib; // TOF calibration object pointer
  
  Int_t fnT0; // total number of T0-TOF
  Int_t fiT0; // last T0-TOF used for T0 fill
  Double_t fT0fill[1000];  // array for dynamical t0 fill calculation
  Double_t fT0sigmaTOF[1000]; // array for dynamical t0 fill resolution

  Bool_t fNoTOFT0;   // swithc to avoid T0-TOF is used
  Bool_t fESDswitch; // if you want take the ESD time instead of the raw + time slewing correction
  
  Double_t fCalculated[8]; // contains the parameters with the event time
  Double_t fTimeResolution;  // global time resolution used to calculate T0
  
  Float_t fT0sigma; // T0 resolution
  
  TH1F *fHmapChannel; // histo with the channel map
  Bool_t fKmask; // switch if you want apply a channel filter
  
  static const Int_t fgkNmaxT0step = 500; //number of steps in the t0 fill calculation

  ClassDef(AliTOFT0maker,2);  // Calculate the time zero using TOF detector */
  
};

#endif // ALITOFT0MAKER_H
