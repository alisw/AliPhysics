#ifndef ALIT0DIGIT_H
#define ALIT0DIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
class TClonesArray;
class TArrayI;

//___________________________________________
class AliT0digit: public TObject {
  ////////////////////////////////////////////////////////////////////////
 public:
  AliT0digit();
  virtual ~AliT0digit();
  Int_t BestTimeRight() {return fBestTimeRight;}
  Int_t BestTimeLeft() {return fBestTimeLeft;}
  Int_t MeanTime() {return fTimeAverage;}
  Int_t TimeDiff() {return fTimeDiff;}
  Int_t SumMult() {return fSumMult;}
  void SetTimeBestRight( Int_t time) {fBestTimeRight = time;}
  void SetTimeBestLeft( Int_t time) {fBestTimeLeft = time;}
  void SetMeanTime(Int_t time) {fTimeAverage=time;}
  void SetDiffTime(Int_t time) {fTimeDiff=time;}
  void SetSumMult(Int_t time) {fSumMult=time;}
  
  virtual void SetTime (TArrayI &o);
  virtual void GetTime (TArrayI &o);
  virtual void SetADC (TArrayI &o);
  virtual void GetADC (TArrayI &o);
  
  virtual void SetTimeAmp (TArrayI &o);
  virtual void GetTimeAmp (TArrayI &o);
  virtual void SetADCAmp (TArrayI &o);
  virtual void GetADCAmp (TArrayI &o);
 private: 

  Int_t fBestTimeRight;        // TOF first particle on the right 
  Int_t fBestTimeLeft;         // TOF first particle on the left
  Int_t fTimeAverage;             // mean time (start signal)
  Int_t fTimeDiff;             // time difference (vertex position)

  TArrayI *fTime;    // array's TDC
  TArrayI *fADC;    // array's ADC
  TArrayI *fTimeAmp;    // array's TDC
  TArrayI *fADCAmp;    // array's ADC
  Int_t fSumMult;   //multiplisity
  AliT0digit( const AliT0digit& );
  AliT0digit& operator=(const AliT0digit&); 
  
  ClassDef(AliT0digit,4)  //Digit (Header) object for set:T0
};

typedef AliT0digit AliSTARTdigit; // for backward compatibility

#endif




