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

  Int_t BestTimeA()   {return fBestTimeA;}
  Int_t BestTimeC()   {return fBestTimeC;}
  Int_t MeanTime()    {return fTimeAverage;}
  Int_t TimeDiff()    {return fTimeDiff;}
  Int_t SumMult()     {return fSumMult;}
  Int_t RefPoint()    {return fRefPoint;}
  void SetTimeBestA( Int_t time) {fBestTimeA = time;}
  void SetTimeBestC( Int_t time) {fBestTimeC = time;}
  void SetMeanTime(Int_t time) {fTimeAverage=time;}
  void SetDiffTime(Int_t time) {fTimeDiff=time;}
  void SetSumMult(Int_t time) {fSumMult=time;}
  void SetRefPoint(Int_t time) {fRefPoint=time;}
  
  virtual void SetTimeCFD (TArrayI &o);
  virtual void GetTimeCFD (TArrayI &o);
  virtual void SetQT0 (TArrayI &o);
  virtual void GetQT0 (TArrayI &o);
  
  virtual void SetTimeLED (TArrayI &o);
  virtual void GetTimeLED (TArrayI &o);
  virtual void SetQT1 (TArrayI &o);
  virtual void GetQT1 (TArrayI &o);

 private: 

  TArrayI *fTimeCFD;    // array's TDC
  TArrayI *fQT0;    // array's ADC
  TArrayI *fTimeLED;    // array's TDC
  TArrayI *fQT1;    // array's ADC
  Int_t fTimeAverage;             // mean time (start signal)
  Int_t fTimeDiff;             // time difference (vertex position)
  Int_t fBestTimeA;        // TOF first particle on the right 
  Int_t fBestTimeC;         // TOF first particle on the left
  Int_t fSumMult;   //multiplisity
  Int_t fRefPoint; // reference point
  AliT0digit( const AliT0digit& );
  AliT0digit& operator=(const AliT0digit&); 
  
  ClassDef(AliT0digit,6)  //Digit (Header) object for set:T0
};


typedef AliT0digit AliSTARTdigit; // for backward compatibility

#endif




