#ifndef ALIMUONRESPONSEV0_H
#define ALIMUONRESPONSEV0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup sim
/// \class AliMUONResponseV0
/// \brief Implementation of Mathieson response

#include "AliMUONResponse.h"
#include "AliMUONMathieson.h"

class AliMUONResponseV0 : public AliMUONResponse
{
 public:
  AliMUONResponseV0();
    virtual ~AliMUONResponseV0();
    //
    // Configuration methods
    //
    // Set number of sigmas over which cluster didintegration is performed
    virtual void    SetSigmaIntegration(Float_t p1) {fSigmaIntegration=p1;}
    // Get number of sigmas over which cluster didintegration is performed   
    virtual Float_t SigmaIntegration() const {return fSigmaIntegration;}    
    // Set single electron pulse height (ADCcounts/e)
    virtual void    SetChargeSlope(Float_t p1) {fChargeSlope=p1;}
    // Get Set single electron pulse height (ADCcounts/e)
    virtual Float_t ChargeSlope() const     {return fChargeSlope;}
    // Set sigmas of the charge spread function
    virtual void    SetChargeSpread(Float_t p1, Float_t p2)
	{fChargeSpreadX=p1; fChargeSpreadY=p2;}
    // Get sigma_X of the charge spread function
    virtual Float_t ChargeSpreadX() const    {return fChargeSpreadX;}
    // Get sigma_Y of the charge spread function
    virtual Float_t ChargeSpreadY() const    {return fChargeSpreadY;}        
    // Set maximum Adc-count value
    virtual void    SetMaxAdc(Int_t p1) {fMaxAdc=p1;}
    // Set saturation value
    virtual void    SetSaturation(Int_t p1) {fSaturation=p1;}
    // Set zero suppression threshold
    virtual void    SetZeroSuppression(Int_t p1) {fZeroSuppression=p1;}
    // Get maximum Adc-count value   
    virtual Int_t   MaxAdc() const          {return fMaxAdc;}
    // Get saturation value   
    virtual Int_t   Saturation() const      {return fSaturation;}

    // Get zero suppression threshold
    virtual Int_t   ZeroSuppression() const {return fZeroSuppression;}
    // Set the charge correlation
    virtual void SetChargeCorrel(Float_t correl){fChargeCorrel = correl;}
    // Get the charge correlation
    virtual Float_t ChargeCorrel() const {return fChargeCorrel;}


    // Set anode cathode Pitch
    virtual Float_t Pitch() const           {return fMathieson->Pitch();}
    // Get anode cathode Pitch
    virtual void    SetPitch(Float_t p1)    {fMathieson->SetPitch(p1);};

    // Set Mathieson parameters
    // Mathieson \sqrt{Kx3} and derived Kx2 and Kx4 
    // passing pointer to class Mathieson for backward compatibility
    virtual void    SetSqrtKx3AndDeriveKx2Kx4(Float_t SqrtKx3);
    // Mathieson \sqrt{Kx3}
    virtual void    SetSqrtKx3(Float_t p1) {fMathieson->SetSqrtKx3(p1);};
    // Mathieson Kx2
    virtual void    SetKx2(Float_t p1)     {fMathieson->SetKx2(p1);};
    // Mathieson Kx4
    virtual void    SetKx4(Float_t p1)     {fMathieson->SetKx4(p1);};
    // Mathieson \sqrt{Ky3} and derived Ky2 and Ky4
    virtual void SetSqrtKy3AndDeriveKy2Ky4(Float_t SqrtKy3);
    // Mathieson \sqrt{Ky3}
    virtual void    SetSqrtKy3(Float_t p1) {fMathieson->SetSqrtKy3(p1);};
    // Mathieson Ky2
    virtual void    SetKy2(Float_t p1)     {fMathieson->SetKy2(p1);};
    // Mathieson Ky4
      virtual void SetKy4(Float_t p1)     {fMathieson->SetKy4(p1);};
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t  IntPH(Float_t eloss) const;
    // Charge disintegration
    virtual Float_t  IntXY(Int_t idDE, 
			   AliMUONGeometrySegmentation* segmentation) const;

    virtual Float_t GetAnod(Float_t x) const;
    
    virtual void DisIntegrate(const AliMUONHit& hit, TList& digits);
    
    virtual void Print(Option_t* opt="") const;
     
 protected:
    Float_t fChargeSlope;              ///< Slope of the charge distribution
    Float_t fChargeSpreadX;            ///< Width of the charge distribution in x
    Float_t fChargeSpreadY;            ///< Width of the charge distribution in y
    Float_t fSigmaIntegration;         ///< Number of sigma's used for charge distribution
    Int_t   fMaxAdc;                   ///< Maximum ADC channel
    Int_t   fSaturation;               ///< Pad saturation in ADC channel
    Int_t   fZeroSuppression;          ///< Zero suppression threshold
    Float_t fChargeCorrel;             ///< \brief amplitude of charge correlation on 2 cathods
                                       ///  is RMS of ln(q1/q2)
    AliMUONMathieson* fMathieson;      ///< pointer to mathieson fct
    Float_t fChargeThreshold;          ///< Charges below this threshold are = 0  

  private:
    AliMUONResponseV0(const AliMUONResponseV0& rhs);
    AliMUONResponseV0& operator = (const AliMUONResponseV0& rhs);

   
    ClassDef(AliMUONResponseV0,2) // Implementation of detector response
};

#endif











