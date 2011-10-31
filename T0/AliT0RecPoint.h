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
    AliT0RecPoint& operator= (const AliT0RecPoint &) { return *this;}
    virtual ~AliT0RecPoint() {}

    Double32_t  GetMeanTime()   const {return fTimeAverage;}
    Int_t  GetOnlineMean() const {return fTimeOnlineMean;}
    Double32_t  GetBestTimeA()  const {return fTimeBestA ;}
    Double32_t  GetBestTimeC()  const {return fTimeBestC ;}
    Float_t GetMultC()       const {return fMultC;}
    Double32_t  Get1stTimeA()  const {return fTime1stA ;}
    Double32_t  Get1stTimeC()  const {return fTime1stC ;}
    Float_t GetMultA()       const {return fMultA;}
    Double32_t  GetVertex()   const {return fVertexPosition;}


    void SetMeanTime(Double32_t time) {fTimeAverage=time;}
    void SetOnlineMean(Int_t time) {fTimeOnlineMean=time;}
    void SetTimeBestA( Double32_t time) {fTimeBestA = time;}
    void SetTimeBestC( Double32_t time) {fTimeBestC = time;}
    void SetTime1stA( Double32_t time) {fTime1stA = time;}
    void SetTime1stC( Double32_t time) {fTime1stC = time;}
    void SetVertex( Double32_t vertex) {fVertexPosition= vertex;}
    void SetMultC(Float_t mult) {fMultC = mult;}
    void SetMultA(Float_t mult) {fMultA = mult;}

    void    SetTime (Int_t ipmt, Double32_t time) { fTime[ipmt] = time;}
    Double32_t GetTime (Int_t ipmt)const { return fTime[ipmt];}
    void    SetAmp (Int_t ipmt, Double32_t adc) { fADC[ipmt] = adc;}
    Double32_t GetAmp (Int_t ipmt) const{ return fADC[ipmt];}
    void    SetAmpLED (Int_t ipmt, Double32_t adc) { fADCLED[ipmt] = adc;}
    Double32_t AmpLED (Int_t ipmt) const{ return fADCLED[ipmt];}
    
    void    SetT0clock (Double32_t time) { fT0clock = time;}
    Double32_t GetT0clock () const{ return fT0clock;}

    Bool_t GetT0Trig(Int_t i) {return (fT0trig&(1<<i)) != 0;}
    Int_t GetT0Trig() {return fT0trig;}
    void   SetT0Trig(Bool_t *tr );
    void PrintTriggerSignals(Int_t trig);

    Float_t GetTimeFull(Int_t ch, Int_t hit) {return fTimeFull[ch][hit];}
    Float_t GetOrA(Int_t hit) {return fOrA[hit];}
    Float_t GetOrC(Int_t hit) {return fOrC[hit];}
    Float_t GetTVDC(Int_t hit) {return fTVDC[hit];}

    void SetTimeFull(Int_t ch, Int_t hit, Float_t time) {fTimeFull[ch][hit] = time;}
    void SetOrA (Int_t hit, Float_t time) { fOrA[hit] = time ;}
    void SetOrC (Int_t hit, Float_t time) { fOrC[hit] = time;}
    void SetTVDC(Int_t hit, Float_t time) { fTVDC[hit] = time;}

  private: 
    Double32_t fTimeAverage;     // Average time with best particles
    Int_t   fTimeOnlineMean; // online mean signal
    Double32_t fVertexPosition;     // Diffrence time between C and A
    Double32_t fTimeBestA;   //TOF first particle on the A
    Double32_t fTimeBestC;    //TOF first particle on the C
    Float_t fMultC; // multiplicity on the 
    Float_t fMultA; // multiplicity on the 
    Double32_t fT0clock; // T0 with best particles
    Int_t   fT0trig;    // T0 trigger signals
 
    Double32_t fTime[24];    // array's TDC
    Double32_t fADC[24];    // array's amplitude
    Double32_t fADCLED[24];    // array's LED amplitude

 
    Float_t fTimeFull[24][5];    // array's TDC no-correction; centred around 0
    Float_t fOrA[5];  //hardware OrA centred around 0
    Float_t fOrC[5];  //hardware OrC centred around 0
    Float_t fTVDC[5]; //hardware TVDC centred around 0
    Bool_t fPileup;
    Bool_t fSattelite;
    Double32_t fTime1stA;   //TOF first particle on the A
    Double32_t fTime1stC;    //TOF first particle on the C

    ClassDef(AliT0RecPoint,8)  // RecPoints (Header) object for set:T0
};

typedef AliT0RecPoint AliSTARTRecPoint; // for backward compatibility

#endif



