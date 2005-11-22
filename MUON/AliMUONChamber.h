#ifndef ALIMUONCHAMBER_H
#define ALIMUONCHAMBER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONChamber
/// \brief MUON tracking chamber class

#include <TObject.h>
#include <TObjArray.h>

#include "AliMUONResponse.h"
#include "AliMUONGeometrySegmentation.h"

class AliMUONClusterFinderVS;
class AliMUON;
class AliMUONHit;


class AliMUONChamber : public TObject
{
 public:
    AliMUONChamber();
    AliMUONChamber(Int_t id);
    virtual ~AliMUONChamber();
    
//
// Get chamber Id
  virtual Int_t   GetId() const {return fId;}
//

//  
// Set response model
  virtual void    SetResponseModel(AliMUONResponse* thisResponse) {fResponse=thisResponse;}
//  
//  Get pointer to response model
  virtual AliMUONResponse* &ResponseModel(){return fResponse;}

//
// Member function forwarding to the segmentation and response models
//
// Calculate pulse height from energy loss  
  virtual Float_t IntPH(Float_t eloss) {return fResponse->IntPH(eloss);}

// Initialisation of charge fluctuation for given hit
  virtual void    ChargeCorrelationInit();

// Configuration forwarding
//
// Define signal distribution region
// by number of sigmas of the distribution function
  virtual void   SetSigmaIntegration(Float_t p1)
      {fResponse->SetSigmaIntegration(p1);}
// Set the single electron pulse-height (ADCchan/e)  
  virtual void   SetChargeSlope(Float_t p1)              {fResponse->SetChargeSlope(p1);}
// Set width of charge distribution function  
  virtual void   SetChargeSpread(Float_t p1, Float_t p2) {fResponse->SetChargeSpread(p1,p2);}
// Set maximum ADC count value
  virtual void   SetMaxAdc(Int_t p1)                   {fResponse->SetMaxAdc(p1);}
//  
// Cluster formation method (charge disintegration)
  virtual void   DisIntegration(AliMUONHit* hit,
				Int_t& x, Float_t newclust[6][500]);
//
  virtual void SetChargeCorrel(Float_t correl) {fResponse->SetChargeCorrel(correl);}

 protected:
  AliMUONChamber(const AliMUONChamber & rChamber);
  // assignment operator  
  AliMUONChamber& operator =(const AliMUONChamber& rhs);

  Int_t   fId;   // chamber number
  Float_t fCurrentCorrel; //! charge correlation for current hit.

  AliMUONResponse        *fResponse;        // pointer to response
  AliMUON                *fMUON;            // pointer to MUON
  ClassDef(AliMUONChamber,3) // Muon tracking chamber class
};

#endif
