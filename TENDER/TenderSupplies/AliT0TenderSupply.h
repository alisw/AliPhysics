#ifndef ALIT0TENDERSUPPLY_H
#define ALIT0TENDERSUPPLY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//      //
//   //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <AliTenderSupply.h>

class AliT0TenderSupply: public AliTenderSupply {
  
 public:
  AliT0TenderSupply();
  AliT0TenderSupply(const char *name, const AliTender *tender=NULL);
  virtual ~AliT0TenderSupply();

  virtual void          Init();
  virtual void          ProcessEvent();
  void SetCorrectMeanTime (Bool_t flag=kFALSE){fCorrectMeanTime=flag;};
  void SetAmplutudeCorrection (Bool_t flag=kFALSE){fCorrectStartTimeOnAmplSatur=flag;};
  void SetPass4LHC11aCorrection (Bool_t flag=kFALSE){fPass4LHC11aCorrection=flag;};
  void SetLHC16jMC (const Double32_t* time, Float_t shift);

 private:
  
  AliT0TenderSupply(const AliT0TenderSupply&c);
  AliT0TenderSupply& operator= (const AliT0TenderSupply&c);


  Bool_t  fCorrectMeanTime; //! mean time shift will be corrected
  Float_t fTimeOffset[4]; //! time offset to be used for fCorrectMeanTime
  Bool_t  fCorrectStartTimeOnAmplSatur; //!  fix start times suffering from saturated amplitude in pmts
  Float_t fAmplitudeThreshold; //! above this value pmt suffer from saturation
  Bool_t fPass4LHC11aCorrection; //! above this value pmt suffer from saturation
  Bool_t fLHC16j;  //  MC production to be fix
  Float_t fFixMeanCFD[24];  // mean CFD in MC
  AliESDEvent *fEvent;   //ESD event
  
  Int_t fnEvent; //event number
  
  ClassDef(AliT0TenderSupply, 3);  // T0 tender supply
};

#endif
