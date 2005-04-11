#ifndef ALISTARTRECPOINT_H
#define ALISTARTRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
class TArrayI;

//___________________________________________
class AliSTARTRecPoint: public TObject  {
////////////////////////////////////////////////////////////////////////
 public:
    AliSTARTRecPoint();
    virtual ~AliSTARTRecPoint() {}
    Int_t  GetMeanTime() {return fTimeAverage;}
    Int_t  GetBestTimeRight() {return fTimeBestRight ;}
    Int_t  GetBestTimeLeft() {return fTimeBestLeft ;}
    Int_t GetMultC() {return fMultC;}
    Int_t GetMultA() {return fMultA;}
    Int_t GetMult() {return fMult;}
    Float_t  GetVertex() {return fVertexPosition;}

    void SetMeanTime(Int_t time) {fTimeAverage=time;}
    void SetTimeBestRight( Int_t time) {fTimeBestRight = time;}
    void SetTimeBestLeft( Int_t time) {fTimeBestLeft = time;}
    void SetVertex( Float_t vertex) {fVertexPosition= vertex;}
    void SetMultC(Int_t mult) {fMultC = mult;}
    void SetMultA(Int_t mult) {fMultA = mult;}
    void SetMult(Int_t mult) {fMult = mult;}

  private: 
    //    Float_t fProcessKoef;  // for pp fProcessKoef=1 ; for Pb-Pb - 0.001
    Int_t fTimeAverage;     // Average time
    Float_t fVertexPosition;     // Diffrence time between left and right
    Int_t fTimeBestRight;   //TOF first particle on the right
    Int_t fTimeBestLeft;    //TOF first particle on the left
    Int_t fMultC; // multiplicity on the 
    Int_t fMultA; // multiplicity on the 
    Int_t fMult; // multiplicity A && C 
 
    ClassDef(AliSTARTRecPoint,2)  //Digit (Header) object for set:START
};


#endif



