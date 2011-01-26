#ifndef ALITOFT0MAKERANA_H
#define ALITOFT0MAKERANA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliRsnTOFT0maker.h,v 1.8 2010/01/19 16:32:20 noferini Exp $ */

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
#include "TH1F.h"

class AliTOFcalibHisto;
class AliTOFT0v1;

class AliRsnTOFT0maker : public TObject {
public:

  enum ESettings
  {
    kNone,
    kPass2,
    kPass4,
    kLHC09d10
  };

  ESettings fSettings;
  TString Settings()
  {
    TString out;
    switch (fSettings)
    {
      case kPass2:    out = "pass 2"; break;
      case kPass4:    out = "pass 4"; break;
      case kLHC09d10: out = "LHC09d10"; break;
      default:        out = "none specific"; break;
    }
    return out;
  }
      
    
  AliRsnTOFT0maker() ;
  virtual ~AliRsnTOFT0maker() ; // dtor
  AliRsnTOFT0maker(const AliRsnTOFT0maker & t);
  AliRsnTOFT0maker & operator=(const AliRsnTOFT0maker & t);
 
  void SetESDdata(Bool_t val=kTRUE){fESDswitch=val;};

  // return (fCalculated[0]=event time -- fCalculated[1]=sigma event time in ps -- fCalculated[2]=mean event time for each fill -- fCalculated[3]=number of tracks at the TOF level) if you can subtruct the event time; return NULL if there is no event time
  Double_t *RemakePID(AliESDEvent *esd,Double_t t0time=0.,Double_t t0sigma=1000.); // t0time and t0sigma in ps

  void      SetTimeResolution(Double_t timeresolution){fTimeResolution=timeresolution;};// TOF timeresolution in [s] e.g. for 120 ps -> 1.2e-10
  Double_t  GetTimeResolution() const {return fTimeResolution;}
  
  void LoadChannelMap(char *filename="$ALICE_ROOT/TOF/enableMap.104892.root"); //load the enable channel map
  void ApplyMask(AliESDEvent *esd);
  
  void SetNoTOFT0(Bool_t status=kTRUE){fNoTOFT0=status;}; // disable the TOF T0 info
  void SetMaskOffChannel(Bool_t status=kTRUE){fKmask=status;}; // swith for the map off channel
  
 private:
  void TakeTimeRawCorrection(AliESDEvent * const esd);
  void RemakeTOFpid(AliESDEvent *esd,Float_t timezero);
  Double_t GetT0Fill(Int_t nrun) const ;
  
  AliTOFcalibHisto *fCalib; // TOF calibration object pointer
  
  Bool_t fESDswitch; // if you want take the ESD time instead of the raw + time slewing correction
  
  Double_t fCalculated[4]; // contains the parameters with the event time
  Double_t fTimeResolution;  // global time resolution used to calculate T0
  
  Float_t fT0sigma; // T0 resolution
  
  TH1F *fHmapChannel; // histo with the channel map
  Bool_t fKmask;
  Bool_t fNoTOFT0;
  
  ClassDef(AliRsnTOFT0maker,1);  // Calculate the time zero using TOF detector */
  
};

#endif // ALITOFT0MAKERANA_H
