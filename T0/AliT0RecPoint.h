#ifndef ALIT0RECPOINT_H
#define ALIT0RECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>


//___________________________________________
class AliT0RecPoint: public TObject  {
////////////////////////////////////////////////////////////////////////
 public:
    AliT0RecPoint();
    AliT0RecPoint(const AliT0RecPoint &o);
    AliT0RecPoint& operator= (const AliT0RecPoint &)  { return *this;}
    virtual ~AliT0RecPoint() {}

     Int_t  GetMeanTime()   const {return fTimeAverage;}
     Int_t  GetOnlineMean() const {return fTimeOnlineMean;}
     Int_t  GetBestTimeA()  const {return fTimeBestA ;}
     Int_t  GetBestTimeC()  const {return fTimeBestC ;}
     Int_t GetMultC()       const {return fMultC;}
     Int_t GetMultA()       const {return fMultA;}
     Float_t  GetVertex()   const {return fVertexPosition;}


    void SetMeanTime(Int_t time) {fTimeAverage=time;}
    void SetOnlineMean(Int_t time) {fTimeOnlineMean=time;}
    void SetTimeBestA( Int_t time) {fTimeBestA = time;}
    void SetTimeBestC( Int_t time) {fTimeBestC = time;}
    void SetVertex( Float_t vertex) {fVertexPosition= vertex;}
    void SetMultC(Int_t mult) {fMultC = mult;}
    void SetMultA(Int_t mult) {fMultA = mult;}

    void SetTime (Int_t ipmt, Float_t time) { fTime[ipmt] = time;}
    Float_t GetTime (Int_t ipmt)const { return fTime[ipmt];}
    void SetAmp (Int_t ipmt, Float_t adc) { fADC[ipmt] = adc;}
    Float_t GetAmp (Int_t ipmt) const{ return fADC[ipmt];}
    void SetAmpLED (Int_t ipmt, Float_t adc) { fADCLED[ipmt] = adc;}
    Float_t AmpLED (Int_t ipmt) const{ return fADCLED[ipmt];}

  private: 
    Int_t fTimeAverage;     // Average time
    Int_t fTimeOnlineMean; // online mean signal
    Float_t fVertexPosition;     // Diffrence time between C and A
    Int_t fTimeBestA;   //TOF first particle on the A
    Int_t fTimeBestC;    //TOF first particle on the C
    Int_t fMultC; // multiplicity on the 
    Int_t fMultA; // multiplicity on the 
 
    Float_t fTime[24];    // array's TDC
    Float_t fADC[24];    // array's amplitude
    Float_t fADCLED[24];    // array's LED amplitude


    ClassDef(AliT0RecPoint,4)  //Digit (Header) object for set:T0
};

typedef AliT0RecPoint AliSTARTRecPoint; // for backward compatibility

#endif



