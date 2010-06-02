#ifndef ALITOFT0MAKERANA_H
#define ALITOFT0MAKERANA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFT0makerANA.h,v 1.8 2010/01/19 16:32:20 noferini Exp $ */

///////////////////////////////////////////////
//					     //
//  Manager class for time zero evaluation   //
//  with TOF informations                    //
//					     //
///////////////////////////////////////////////


#include "TObject.h"

class TH1F;
class AliESDEvent;

class AliTOFT0v2;
class AliESDpid;

class AliTOFT0makerANA : public TObject {
public:
  
  AliTOFT0makerANA() ;
  AliTOFT0makerANA(AliESDpid * const externalPID);
  virtual ~AliTOFT0makerANA() ; // dtor
  AliTOFT0makerANA(const AliTOFT0makerANA & t);
  AliTOFT0makerANA & operator=(const AliTOFT0makerANA & t);
 
  // return (fCalculated[0]=event time -- fCalculated[1]=sigma event time in ps -- fCalculated[2]=mean event time for each fill -- fCalculated[3]=number of tracks at the TOF level) if you can subtruct the event time -- ...
  Double_t *RemakePID(AliESDEvent *esd,Double_t t0time=0.,Double_t t0sigma=1000.); // t0time and t0sigma in ps

  void      SetTimeResolution(Double_t timeresolution){fTimeResolution=timeresolution;};// TOF timeresolution in [s] e.g. for 120 ps -> 1.2e-10
  Double_t  GetTimeResolution() const {return fTimeResolution;}
  
  void LoadChannelMap(char *filename="$ALICE_ROOT/TOF/enableMap.104892.root"); //load the enable channel map
  void ApplyMask(AliESDEvent * const esd);
  
  void SetNoTOFT0(Bool_t status=kTRUE){fNoTOFT0=status;}; // disable the TOF T0 info
  void SetMaskOffChannel(Bool_t status=kTRUE){fKmask=status;}; // swith for the map off channel
  
 private:
  AliESDpid *fPIDesd; //!ESD pid object
  void TakeTimeRawCorrection(AliESDEvent * const esd);
  void RemakeTOFpid(/*AliESDEvent *esd,*/Float_t timezero);
  Double_t GetT0Fill() const;
  
  Int_t fnT0; //! total number of T0-TOF
  Int_t fiT0; //! last T0-TOF used for T0 fill
  Double_t fT0fill[1000];  //! array for dynamical t0 fill calculation
  Double_t fT0sigmaTOF[1000]; //! array for dynamical t0 fill resolution

  Bool_t fNoTOFT0;   //! switch to avoid T0-TOF is used
  
  Double_t fCalculated[8]; //! contains the parameters with the event time
  Double_t fTimeResolution;  //! global time resolution used to calculate T0
  
  Float_t fT0sigma; //! T0 resolution
  
  TH1F *fHmapChannel; //! histo with the channel map
  Bool_t fKmask; //! switch if you want apply a channel filter
  
  static const Int_t fgkNmaxT0step = 500; //number of steps in the t0 fill calculation

  ClassDef(AliTOFT0makerANA,2);  // Calculate the time zero using TOF detector */
  
};

#endif // ALITOFT0MAKERANA_H
