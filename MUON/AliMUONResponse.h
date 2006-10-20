#ifndef ALIMUONRESPONSE_H
#define ALIMUONRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup sim
/// \class AliMUONResponse
/// \brief Chamber response base class

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliMUONDigit;
class AliMUONGeometrySegmentation;
class AliMUONHit;
class TF1;
class TList;

class AliMUONResponse : public TObject 
{
 public:
    AliMUONResponse();
    virtual ~AliMUONResponse();
 
    //
    // Configuration methods
    //
    // Set number of sigmas over which cluster disintegration is performed
    virtual void    SetSigmaIntegration(Float_t)           {return;}
    // Get number of sigmas over which cluster disintegration is performed
    virtual Float_t SigmaIntegration() const                  {return 1.;}
    // Set single electron pulse height (ADCcounts/e)
    virtual void    SetChargeSlope(Float_t )                {return;}
    // Get single electron pulse height (ADCcounts/e)
    virtual Float_t ChargeSlope() const                       {return 1.;}
    // Set sigmas of the charge spread function
    virtual void    SetChargeSpread(Float_t , Float_t )   {return;}
    // Get sigma_X of the charge spread function
    virtual Float_t ChargeSpreadX() const                     {return 1.;}
    // Get sigma_Y of the charge spread function
    virtual Float_t ChargeSpreadY() const                     {return 1.;}
    // Set maximum Adc-count value
    virtual void    SetMaxAdc(Int_t )                       {return;}
    // Set saturation value
    virtual void    SetSaturation(Int_t )                   {return;}
    // Set zero suppression threshold
    virtual void    SetZeroSuppression(Int_t )             {return;}
    // Get maximum Adc-count value
    virtual Int_t MaxAdc() const                              {return kTRUE;}
    // Get saturation value
    virtual Int_t Saturation() const                          {return kTRUE;}
    // Get maximum zero suppression threshold
    virtual Int_t ZeroSuppression() const                     {return kTRUE;}
    // Set anode cathode Pitch
    virtual void    SetPitch(Float_t)                         {return;}
    // Get anode cathode Pitch
    virtual Float_t Pitch() const                             {return 1.;}
    // Set the charge correlation
    virtual void SetChargeCorrel(Float_t)                     {return;}
    // Get the charge correlation
    virtual Float_t ChargeCorrel() const                      {return 1.;}
    //  
    // Chamber response methods
    // Pulse height from scored quantity (eloss)
    virtual Float_t IntPH(Float_t) const                      {return 1.;}
    // Charge disintegration 
    virtual Float_t IntXY(Int_t, AliMUONGeometrySegmentation*) const {return 1.;}
    
    /// Go from one hit to several digits, applying charge spreading.
    virtual void DisIntegrate(const AliMUONHit& hit, TList& digits);

    // 
    ClassDef(AliMUONResponse,1) // Chamber response virtual base class 
};
#endif







