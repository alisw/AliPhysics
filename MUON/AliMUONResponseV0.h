#ifndef ALIMUONRESPONSEV0_H
#define ALIMUONRESPONSEV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMUONResponse.h"

class AliMUONResponseV0 : 
public AliMUONResponse {
 public:
    AliMUONResponseV0(){}
    virtual ~AliMUONResponseV0(){}
    //
    // Configuration methods
    //
    // Set number of sigmas over which cluster didintegration is performed
    virtual void    SetSigmaIntegration(Float_t p1) {fSigmaIntegration=p1;}
    // Get number of sigmas over which cluster didintegration is performed   
    virtual Float_t SigmaIntegration() {return fSigmaIntegration;}    
    // Set single electron pulse height (ADCcounts/e)
    virtual void    SetChargeSlope(Float_t p1) {fChargeSlope=p1;}
    // Get Set single electron pulse height (ADCcounts/e)
    virtual Float_t ChargeSlope()      {return fChargeSlope;}
    // Set sigmas of the charge spread function
    virtual void    SetChargeSpread(Float_t p1, Float_t p2)
	{fChargeSpreadX=p1; fChargeSpreadY=p2;}
    // Get sigma_X of the charge spread function
    virtual Float_t ChargeSpreadX()    {return fChargeSpreadX;}
    // Get sigma_Y of the charge spread function
    virtual Float_t ChargeSpreadY()    {return fChargeSpreadY;}        
    // Set maximum Adc-count value
    virtual void    SetMaxAdc(Int_t p1) {fMaxAdc=p1;}
    // Set zero suppression threshold
    virtual void    SetZeroSuppression(Int_t p1) {fZeroSuppression=p1;}
    // Get maximum Adc-count value   
    virtual Int_t   MaxAdc()           {return fMaxAdc;}
    // Get zero suppression threshold
    virtual Int_t   ZeroSuppression() {return fZeroSuppression;}
    // Set anode cathode Pitch
    virtual Float_t Pitch()            {return fPitch;}
    // Get anode cathode Pitch
    virtual void    SetPitch(Float_t p1) {fPitch=p1;};
    // Set Mathieson parameters
    // Mathieson \sqrt{Kx3}
    virtual void    SetSqrtKx3(Float_t p1) {fSqrtKx3=p1;};
    // Mathieson Kx2
    virtual void    SetKx2(Float_t p1) {fKx2=p1;};
    // Mathieson Kx4
    virtual void    SetKx4(Float_t p1) {fKx4=p1;};
    // Mathieson \sqrt{Ky3}
    virtual void    SetSqrtKy3(Float_t p1) {fSqrtKy3=p1;};
    // Mathieson Ky2
    virtual void    SetKy2(Float_t p1) {fKy2=p1;};
    // Mathieson Ky4
    virtual void    SetKy4(Float_t p1) {fKy4=p1;};
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t  IntPH(Float_t eloss);
    // Charge disintegration
    virtual Float_t  IntXY(AliMUONSegmentation * segmentation);
    // Noise, zero-suppression, adc saturation
    virtual Int_t DigitResponse(Int_t digit);

    ClassDef(AliMUONResponseV0,1) // Implementation of Mathieson response
 protected:
    Float_t fChargeSlope;              // Slope of the charge distribution
    Float_t fChargeSpreadX;            // Width of the charge distribution in x
    Float_t fChargeSpreadY;            // Width of the charge distribution in y
    Float_t fSigmaIntegration;         // Number of sigma's used for charge distribution
    Int_t   fMaxAdc;                   // Maximum ADC channel
    Int_t   fZeroSuppression;          // Zero suppression threshold
    Float_t fSqrtKx3;                  // Mathieson \Sqrt{Kx3)
    Float_t fKx2;                      // Mathieson Kx2
    Float_t fKx4;                      // Mathieson Kx4    
    Float_t fSqrtKy3;                  // Mathieson \Sqrt{Kx3)
    Float_t fKy2;                      // Mathieson Ky2
    Float_t fKy4;                      // Mathieson Ky4
    Float_t fPitch;                    // anode-cathode pitch
};
#endif











