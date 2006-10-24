/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// -----------------------------
// Class AliMUONResponseFactory
// -----------------------------
// Factory for muon response
// Class separated from AliMUONFactoryV4

/* $Id$ */

#include "AliMUONResponseFactory.h"
#include "AliRun.h"
#include "AliLog.h"

#include "AliMpPlaneType.h"

#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONChamber.h"
#include "AliMUONResponseV0.h"
#include "AliMUONResponseTrigger.h"
#include "AliMUONResponseTriggerV1.h"

ClassImp(AliMUONResponseFactory)

//__________________________________________________________________________
  AliMUONResponseFactory::AliMUONResponseFactory(const char* name)
    : TNamed(name, ""),
      fMUON(0),
      fResponse0(0)
{
/// Standard constructor
  
  AliDebug(1,Form("ctor this=%p",this));
}

//__________________________________________________________________________
  AliMUONResponseFactory::AliMUONResponseFactory()
    : TNamed(),
      fMUON(0),
      fResponse0(0)
{
/// Default constructor

  AliDebug(1,Form("default (empty) ctor this=%p",this));
}

//__________________________________________________________________________

AliMUONResponseFactory::~AliMUONResponseFactory()
{
/// Destructor
	AliDebug(1,Form("dtor this=%p",this));
}
          
//__________________________________________________________________________
void AliMUONResponseFactory::BuildCommon() 
{
/// Construct the default response.

  // Default response: 5 mm of gas
  fResponse0 = new AliMUONResponseV0;
  fResponse0->SetSqrtKx3AndDeriveKx2Kx4(0.7131); // sqrt(0.5085)
  fResponse0->SetSqrtKy3AndDeriveKy2Ky4(0.7642); // sqrt(0.5840)
  fResponse0->SetPitch(AliMUONConstants::Pitch()); // anode-cathode distance
  fResponse0->SetSigmaIntegration(10.);
  fResponse0->SetChargeSlope(10);
  fResponse0->SetChargeSpread(0.18, 0.18);
  fResponse0->SetMaxAdc(4096);
  fResponse0->SetSaturation(3000);
  fResponse0->SetZeroSuppression(6);
}       
        
//__________________________________________________________________________
void AliMUONResponseFactory::BuildStation1() 
{
/// Configuration for Chamber TC1/2  (Station 1) ----------           

  // Response for 4 mm of gas (station 1)
  // automatic consistency with width of sensitive medium in CreateGeometry ????
  AliMUONResponseV0* responseSt1 = new AliMUONResponseV0;
  // Mathieson parameters from L.Kharmandarian's thesis, page 190
  responseSt1->SetSqrtKx3AndDeriveKx2Kx4(0.7000); // sqrt(0.4900)
  responseSt1->SetSqrtKy3AndDeriveKy2Ky4(0.7550); // sqrt(0.5700)
  responseSt1->SetPitch(AliMUONConstants::PitchSt1()); // anode-cathode distance
  responseSt1->SetSigmaIntegration(10.);
  // ChargeSlope larger to compensate for the smaller anode-cathode distance
  // and keep the same most probable ADC channel for mip's
  responseSt1->SetChargeSlope(62.5); 
  // assumed proportionality to anode-cathode distance for ChargeSpread
  responseSt1->SetChargeSpread(0.144, 0.144);
  responseSt1->SetMaxAdc(4096);
  responseSt1->SetSaturation(3000);
  responseSt1->SetZeroSuppression(6);

   for (Int_t chamber = 0; chamber < 2; chamber++) {
    fMUON->SetResponseModel(chamber, responseSt1); // special response      
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
  }
}

//__________________________________________________________________________
void AliMUONResponseFactory::BuildStation2() 
{
/// Configuration for Chamber TC3/4 (Station 2) -----------

  for (Int_t chamber = 2; chamber < 4; chamber++) {
    fMUON->SetResponseModel(chamber, fResponse0); // normal response        
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
  }
}       
        
//__________________________________________________________________________
void AliMUONResponseFactory::BuildStation3() 
{
/// Configuration for Chamber TC5/6  (Station 3) ----------          

  for (Int_t chamber = 4; chamber < 6; chamber++) {
    fMUON->SetResponseModel(chamber, fResponse0); // normal response        
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
  }
}       
        
//__________________________________________________________________________
void AliMUONResponseFactory::BuildStation4() 
{
/// Configuration for Chamber TC7/8  (Station 4) ----------          

  for (Int_t chamber = 6; chamber < 8; chamber++) {
    fMUON->SetResponseModel(chamber, fResponse0); // normal response        
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
  }
}       
        
//__________________________________________________________________________
void AliMUONResponseFactory::BuildStation5() 
{
/// Configuration for Chamber TC9/10  (Station 5) ---------           

  for (Int_t chamber = 8; chamber < 10; chamber++) {
    fMUON->SetResponseModel(chamber, fResponse0); // normal response        
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
  }
}       
        
//__________________________________________________________________________
void AliMUONResponseFactory::BuildStation6() 
{
/// Configuration for Trigger Chambers   (Station 6,7) ---------           

    Bool_t resTrigV1 = fMUON->GetTriggerResponseV1();    

    for (Int_t chamber = 10; chamber < 14; chamber++) {
	
	if (!resTrigV1) 
	    fMUON->SetResponseModel(chamber, new AliMUONResponseTrigger);
    	else 
	    fMUON->SetResponseModel(chamber, new AliMUONResponseTriggerV1);
	
	fMUON->Chamber(chamber).SetChargeCorrel(0); // same charge on both cathodes
  }
}       

//__________________________________________________________________________
void AliMUONResponseFactory::Build(AliMUON* where) 
{
/// Construct MUON responses

  fMUON = where;

  // Set default parameters
  fMUON->SetIshunt(0);
  fMUON->SetMaxStepGas(0.1);
  fMUON->SetMaxStepAlu(0.1);

  // Build stations
  BuildCommon();
  BuildStation1();
  BuildStation2();
  BuildStation3();
  BuildStation4();
  BuildStation5();
  BuildStation6();
}

//__________________________________________________________________________
void AliMUONResponseFactory::BuildStation(AliMUON* where, Int_t stationNumber) 
{
/// Construct MUON responses for given station

  fMUON = where;
  if (!fResponse0) BuildCommon(); 
    
  switch (stationNumber) {    
    case 1:  BuildStation1(); break;
    case 2:  BuildStation2(); break;
    case 3:  BuildStation3(); break;
    case 4:  BuildStation4(); break;
    case 5:  BuildStation5(); break;
    case 6:  BuildStation6(); break;
    
    default: AliFatal("Wrong station number");
  }  
}         
