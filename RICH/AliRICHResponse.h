#ifndef ALIRICHRESPONSE_H
#define ALIRICHRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


//----------------------------------------------
//
// Chamber response virtual base class
//
#include <TObject.h>
class AliRICHSegmentation;


class AliRICHResponse :
public TObject {
 public:
    //
    // Configuration methods
    //
    // Number of sigmas over which cluster didintegration is performed
    virtual void    SetSigmaIntegration(Float_t p1)           =0;
    virtual Float_t SigmaIntegration()                        =0;
    // charge slope in ADC/e
    virtual void    SetChargeSlope(Float_t p1)                =0;
    virtual Float_t ChargeSlope()                             =0;
    // sigma of the charge spread function
    virtual void    SetChargeSpread(Float_t p1, Float_t p2)   =0;
    virtual Float_t ChargeSpreadX()                           =0;
    virtual Float_t ChargeSpreadY()                           =0;
    // Adc-count saturation value
    virtual void    SetMaxAdc(Float_t p1)                     =0;
    virtual Float_t MaxAdc()                                  =0;
    // anode cathode Pitch
    virtual void    SetPitch(Float_t)                         =0;
    virtual Float_t Pitch()                                   =0;
    // alpha feedback
    virtual void    SetAlphaFeedback(Float_t)                 =0;
    virtual Float_t AlphaFeedback()                           =0;
    // ionisation enrgy
    virtual void    SetEIonisation(Float_t)                   =0;
    virtual Float_t EIonisation()                             =0;
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss)                       =0;
    virtual Float_t IntPH()                                    =0;
    // Charge disintegration
    virtual Float_t IntXY(AliRICHSegmentation *)                 =0;
    virtual Int_t   FeedBackPhotons(Float_t *source, Float_t qtot) =0;
    //
    // Mathieson parameters
    virtual void   SetSqrtKx3(Float_t p1)                        =0;
    virtual void   SetKx2(Float_t p1)                            =0;
    virtual void   SetKx4(Float_t p1)                            =0;
    virtual void   SetSqrtKy3(Float_t p1)                        =0;
    virtual void   SetKy2(Float_t p1)                            =0;
    virtual void   SetKy4(Float_t p1)                            =0;
    ClassDef(AliRICHResponse,1)
};

#endif
