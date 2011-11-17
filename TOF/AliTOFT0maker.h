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
#include "AliTOFT0v1.h"

class TH1F;
class AliESDEvent;

class AliESDpid;
class AliTOFcalib;

class AliTOFT0maker : public TObject {
public:
  
  AliTOFT0maker() ; // default constructor
  AliTOFT0maker(AliESDpid *externalPID, AliTOFcalib *tofCalib=NULL); // overloaded constructor
  virtual ~AliTOFT0maker() ; // dtor
  
  // return (fCalculated[0]=event time -- fCalculated[1]=sigma event time in ps -- fCalculated[2]=mean event time for each fill -- fCalculated[3]=number of tracks at the TOF level) if you can subtruct the event time; return NULL if there is no event time
  Double_t *ComputeT0TOF(AliESDEvent *esd,Double_t t0time=0.,Double_t t0sigma=1000.); // t0time and t0sigma in ps
  void ApplyT0TOF(AliESDEvent *esd);
  Float_t GetExpectedSigma(Float_t mom, Float_t tof, Float_t mass);
  Double_t  *GetT0p(Float_t p);
  
  void      SetTimeResolution(Double_t timeresolution){fTimeResolution=timeresolution;};// TOF timeresolution in [ps]
  Double_t  GetTimeResolution() const {return fTimeResolution;} // Get TOF Time Resolution
  void SetT0FillWidth(Float_t width){if(width > 50) fT0width = width; else fT0width=150;}; // in ps
  
  void LoadChannelMap(const char *filename="$ALICE_ROOT/TOF/enableMap.104892.root"); //load the enable channel map
  void ApplyMask(AliESDEvent * const esd); // Apply the channel mask
  
  void SetNoTOFT0(Bool_t status=kTRUE){fNoTOFT0=status;}; // disable the TOF T0 info
  void SetMaskOffChannel(Bool_t status=kTRUE){fKmask=status;}; // switch for the map off channel
  
  Float_t  TuneForMC(AliESDEvent *esd); // set enabled channeld, add a smeared t0, add a TOF smearing, return true(smeared) T0 event
  
  void SetT0spread(Float_t t0spread){fT0spreadExt=t0spread;}; // get T0spread
  Float_t GetT0spread() const {return fT0spreadExt;} // get T0spread

  void SetT0fill(Float_t t0fill){fT0fillExt=t0fill;};
  
  void WriteInESD(AliESDEvent *esd);

  Int_t GetPileUpCandidate(){return fT0TOF->GetPileUpCandidate();};

 private:
  
  AliTOFT0maker(const AliTOFT0maker &);
  AliTOFT0maker & operator=(const AliTOFT0maker &);
  
  void SetTOFResponse();

  AliTOFT0v1 *fT0TOF; // T0-TOF
  AliESDpid *fPIDesd; // PID esd
  Bool_t fExternalPIDFlag; // external PID flag
  AliTOFcalib *fTOFcalib; // TOF calibration

  Bool_t fNoTOFT0;   // switch to avoid T0-TOF is used

  Int_t fNmomBins;

  Double_t fCalculated[10]; // contains the parameters with the event time

  Double_t fT0cur[2]; // current T0 and T0 sigma
 
  Double_t fTimeResolution;  // global time resolution used to calculate T0
  
  Float_t fT0sigma; // T0 resolution
  
  TH1F *fHmapChannel; // histo with the channel map
  Bool_t fKmask; // switch if you want apply a channel filter
  
  Float_t fT0width; // T0 FILL width

  Float_t fT0spreadExt;
  Float_t fT0fillExt; // t0spread if set 

  ClassDef(AliTOFT0maker,2);  // Calculate the time zero using TOF detector */
  
};

#endif // ALITOFT0MAKER_H
