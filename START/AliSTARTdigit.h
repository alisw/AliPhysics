#ifndef ALISTARTDIGIT_H
#define ALISTARTDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
class TArrayI;

//___________________________________________
class AliSTARTdigit: public TObject  {
////////////////////////////////////////////////////////////////////////
 public:
    AliSTARTdigit();
    virtual ~AliSTARTdigit();
    void SetMeanTime(Int_t time) {fTimeAverage=time;}
    Int_t  GetMeanTime() {return fTimeAverage;}
    Int_t  GetBestTimeRight() {return fTimeBestRight ;}
    Int_t  GetBestTimeLeft() {return fTimeBestLeft ;}
    Int_t  GetSumADCRight() {return fSumADCRight ;}
    void SetTimeBestRight( Int_t time) {fTimeBestRight = time;}
    void SetTimeBestLeft( Int_t time) {fTimeBestLeft = time;}
    void SetSumADCRight( Int_t ADC) {fSumADCRight = ADC;}
    //    void SetProcessKoef( Float_t pp) {fProcessKoef = pp;}
    virtual void SetTime (TArrayI &o);
    virtual void GetTime (TArrayI &o);
    virtual void SetADC (TArrayI &o);
    virtual void GetADC (TArrayI &o);
     virtual const char* GetName() const {return "START_D";}
  private: 
    //    Float_t fProcessKoef;  // for pp fProcessKoef=1 ; for Pb-Pb - 0.001
    Int_t fTimeAverage;     // Average time
    Int_t fTimeBestRight;   //TOF first particle on the right
    Int_t fTimeBestLeft;    //TOF first particle on the left
    Int_t fSumADCRight;    // multiplicity on the right side
    TArrayI *fTime;    // array's TDC
    TArrayI *fADC;    // array's ADC

    ClassDef(AliSTARTdigit,1)  //Digit (Header) object for set:START
};


#endif



