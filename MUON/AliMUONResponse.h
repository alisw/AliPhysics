#ifndef ALIMUONRESPONSE_H
#define ALIMUONRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include <TObject.h>

class TF1;
class AliSegmentation;
class AliMUONTransientDigit;

class AliMUONResponse : public TObject 
{
 public:
    AliMUONResponse();
    virtual ~AliMUONResponse();
 
    //
    // Configuration methods
    //
    // Set number of sigmas over which cluster disintegration is performed
    virtual void    SetSigmaIntegration(Float_t p1)           =0;
    // Get number of sigmas over which cluster disintegration is performed
    virtual Float_t SigmaIntegration() const                  =0;
    // Set single electron pulse height (ADCcounts/e)
    virtual void    SetChargeSlope(Float_t p1)                =0;
    // Get single electron pulse height (ADCcounts/e)
    virtual Float_t ChargeSlope() const                       =0;
    // Set sigmas of the charge spread function
    virtual void    SetChargeSpread(Float_t p1, Float_t p2)   =0;
    // Get sigma_X of the charge spread function
    virtual Float_t ChargeSpreadX() const                     =0;
    // Get sigma_Y of the charge spread function
    virtual Float_t ChargeSpreadY() const                     =0;
    // Set maximum Adc-count value
    virtual void    SetMaxAdc(Int_t p1)                       =0;
    // Set saturation value
    virtual void    SetSaturation(Int_t p1)                   =0;
    // Set zero suppression threshold
    virtual void    SetZeroSuppression(Int_t val)             =0;
    // Get maximum Adc-count value
    virtual Int_t MaxAdc() const                              =0;
    // Get saturation value
    virtual Int_t Saturation() const                          =0;
    // Get maximum zero suppression threshold
    virtual Int_t ZeroSuppression() const                     =0;
    // Set anode cathode Pitch
    virtual void    SetPitch(Float_t)                         =0;
    // Get anode cathode Pitch
    virtual Float_t Pitch() const                             =0;
    // Set the charge correlation
    virtual void SetChargeCorrel(Float_t correl)              =0;
    // Get the charge correlation
    virtual Float_t ChargeCorrel() const                      =0;
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t eloss)                      =0;
    // Charge disintegration 
    virtual Float_t IntXY(AliSegmentation *)                  =0;
    // Noise, zero-suppression, adc saturation
    //virtual Int_t DigitResponse(Int_t digit)                =0;
    virtual Int_t DigitResponse(Int_t digit, 
                                AliMUONTransientDigit* where) =0;
    // 
    ClassDef(AliMUONResponse,1) // Chamber response virtual base class 
};
#endif







