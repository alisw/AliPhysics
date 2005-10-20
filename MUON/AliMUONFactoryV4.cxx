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

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response //
////////////////////////////////////////////////////////////

/* $Id$ */

#include "AliMUONFactoryV4.h"
#include "AliRun.h"
#include "AliLog.h"

#include "AliMpPlaneType.h"

#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONTriggerConstants.h"
#include "AliMUONChamber.h"
#include "AliMUONResponseV0.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONVGeometryDEIndexing.h"
#include "AliMUONSegmentationManager.h"
#include "AliMUONSt12QuadrantSegmentation.h"
#include "AliMUONSt345SlatSegmentationV2.h"
#include "AliMUONTriggerSegmentationV2.h"
#include "AliMUONResponseTrigger.h"

ClassImp(AliMUONFactoryV4)

//__________________________________________________________________________
  AliMUONFactoryV4::AliMUONFactoryV4(const char* name)
    : TNamed(name, ""),
      fMUON(0),
      fResponse0(0),
      fDESegmentations(0)
{  
  AliDebug(1,Form("ctor this=%p",this));
  fDESegmentations = new TObjArray();
  fDESegmentations->SetOwner(kTRUE);
}

//__________________________________________________________________________
  AliMUONFactoryV4::AliMUONFactoryV4()
    : TNamed(),
      fMUON(0),
      fResponse0(0),
      fDESegmentations(0)
{
  AliDebug(1,Form("default (empty) ctor this=%p",this));
// Default constructor
}

//__________________________________________________________________________
AliMUONFactoryV4::AliMUONFactoryV4(const AliMUONFactoryV4& rhs)
 : TNamed(rhs)
{
  // Protected copy constructor

  AliFatal("Not implemented.");
}

//__________________________________________________________________________

AliMUONFactoryV4::~AliMUONFactoryV4()
{
// Destructor
	AliDebug(1,Form("dtor this=%p",this));
  delete fDESegmentations;
}

//__________________________________________________________________________
AliMUONFactoryV4&  AliMUONFactoryV4::operator=(const AliMUONFactoryV4& rhs)
{
  // Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    
          
//__________________________________________________________________________
Bool_t AliMUONFactoryV4::IsGeometryDefined(Int_t ichamber)
{
// Return true, if det elements for the chamber with the given ichamber Id
// are defined in geometry (the geometry builder for this chamber was activated)

  if ( ! fMUON ||
       ! fMUON->Chamber(ichamber).GetGeometry() ||
       ! fMUON->Chamber(ichamber).GetGeometry()->GetDEIndexing() ||
       ! fMUON->Chamber(ichamber).GetGeometry()->GetDEIndexing()->GetNofDetElements() )
       
    return kFALSE;
  
  return kTRUE;
}  

//__________________________________________________________________________
void AliMUONFactoryV4::BuildCommon() 
{
  //
  // Construct the default response.
  //

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
void AliMUONFactoryV4::BuildStation1() 
{
  //--------------------------------------------------------
  // Configuration for Chamber TC1/2  (Station 1) ----------           
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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

  // Quadrant segmentations:
  AliMUONSt12QuadrantSegmentation* bendSt1
    = new AliMUONSt12QuadrantSegmentation(kStation1, kBendingPlane);
  AliMUONSt12QuadrantSegmentation* nonbendSt1
    = new AliMUONSt12QuadrantSegmentation(kStation1, kNonBendingPlane);
  
  // Add in the array (for safe deleting)  
  fDESegmentations->Add(bendSt1);  
  fDESegmentations->Add(nonbendSt1);  

  AliMUONGeometrySegmentation* segmentation[2];

  for (Int_t chamber = 0; chamber < 2; chamber++) {

    segmentation[0] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
    segmentation[1] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
        
    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    //--------------------------------------------------------
    // Configuration for Chamber TC1/2  (Station 1) ----------           


    // cathode 0
    segmentation[0]->Add(id0,      bendSt1);
    segmentation[0]->Add(id0 +  3, nonbendSt1);
    segmentation[0]->Add(id0 +  2, bendSt1);
    segmentation[0]->Add(id0 +  1, nonbendSt1); 
    fMUON->SetSegmentationModel(chamber, 1, segmentation[0]);   

    // cathode 1
    segmentation[1]->Add(id0,      nonbendSt1);
    segmentation[1]->Add(id0 +  3, bendSt1);
    segmentation[1]->Add(id0 +  2, nonbendSt1);
    segmentation[1]->Add(id0 +  1, bendSt1);
    fMUON->SetSegmentationModel(chamber, 2, segmentation[1]);
        
    fMUON->SetResponseModel(chamber, responseSt1); // special response      
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
        
  }
}

//__________________________________________________________________________
void AliMUONFactoryV4::BuildStation2() 
{
  //
  //--------------------------------------------------------
  // Configuration for Chamber TC3/4 (Station 2) -----------
  ///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  // Quadrant segmentations:
  AliMUONSt12QuadrantSegmentation* bendSt2
    = new AliMUONSt12QuadrantSegmentation(kStation2, kBendingPlane);
  AliMUONSt12QuadrantSegmentation* nonbendSt2
    = new AliMUONSt12QuadrantSegmentation(kStation2, kNonBendingPlane);

  // Add in the array (for safe deleting)  
  fDESegmentations->Add(bendSt2);  
  fDESegmentations->Add(nonbendSt2);  

  AliMUONGeometrySegmentation* segmentation[2];

  for (Int_t chamber = 2; chamber < 4; chamber++) {

    segmentation[0] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
    segmentation[1] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
        
    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    //--------------------------------------------------------
    // Configuration for Chamber TC3/4  (Station 2) ----------           


    // cathode 0
    segmentation[0]->Add(id0,      bendSt2);
    segmentation[0]->Add(id0 +  3, nonbendSt2);
    segmentation[0]->Add(id0 +  2, bendSt2);
    segmentation[0]->Add(id0 +  1, nonbendSt2); 
    fMUON->SetSegmentationModel(chamber, 1, segmentation[0]);   

    // cathode 1
    segmentation[1]->Add(id0,      nonbendSt2);
    segmentation[1]->Add(id0 +  3, bendSt2);
    segmentation[1]->Add(id0 +  2, nonbendSt2);
    segmentation[1]->Add(id0 +  1, bendSt2);
    fMUON->SetSegmentationModel(chamber, 2, segmentation[1]);
        
    fMUON->SetResponseModel(chamber, fResponse0); // normal response        
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
        
  }
}       
        
//_____________________________________________________________________________
void
AliMUONFactoryV4::BuildChamber345(Int_t firstDetElemId, Int_t lastDetElemId)
{
  // Build a single chamber for stations 345.
  // The first and lastDetElemId must correspond to the same chamber.
	
  Int_t ichamber = firstDetElemId/100 - 1;
  Int_t test = lastDetElemId/100-1;
  
  if ( test != ichamber )
	{
		AliFatal(Form("DetElemIds %d and %d not part of the same chamber !",
									firstDetElemId,lastDetElemId));
	}
	
  const Int_t kNPLANES = 2;
  const AliMpPlaneType kptypes[kNPLANES] = { kBendingPlane, kNonBendingPlane };
  
  AliMUONChamber& chamber = fMUON->Chamber(ichamber);
	
  for ( Int_t iplane = 0; iplane < kNPLANES; ++iplane )
	{
		AliMUONGeometrySegmentation* segmentation = 
		new AliMUONGeometrySegmentation(chamber.GetGeometry());
		
		for ( Int_t d = firstDetElemId; d <= lastDetElemId; ++d ) 
		{
			if ( !AliMUONSegmentationManager::IsValidDetElemId(d) )
	    {
	      AliWarning(Form("You are requesting an invalid detElemId = %d, I am skipping it",d));
	      continue;
	    }
			
			AliMUONVGeometryDESegmentation* slatSeg = 
	    new AliMUONSt345SlatSegmentationV2(d,kptypes[iplane]);
	    
			fDESegmentations->Add(slatSeg);
			
			segmentation->Add(d,slatSeg);
		}
		
		fMUON->SetSegmentationModel(ichamber,iplane+1,segmentation);
	}
	
  fMUON->SetResponseModel(ichamber,fResponse0);
	
  chamber.SetChargeCorrel(0.11); // 11% charge spread  
}

//_____________________________________________________________________________
void
AliMUONFactoryV4::BuildChamberTrigger(Int_t firstDetElemId, Int_t lastDetElemId)
{
  // Build a single chamber for trigger stations.
  // The first and lastDetElemId must correspond to the same chamber.
	
  Int_t ichamber = firstDetElemId/100 - 1;
  Int_t test = lastDetElemId/100-1;
  
  if ( test != ichamber )
	{
		AliFatal(Form("DetElemIds %d and %d not part of the same chamber !",
									firstDetElemId,lastDetElemId));
	}
	
  const Int_t kNPLANES = 2;
  const AliMpPlaneType kptypes[kNPLANES] = { kBendingPlane, kNonBendingPlane };
  
  AliMUONChamber& chamber = fMUON->Chamber(ichamber);
	
  for ( Int_t iplane = 0; iplane < kNPLANES; ++iplane )
	{
		AliMUONGeometrySegmentation* segmentation = 
		new AliMUONGeometrySegmentation(chamber.GetGeometry());
		
		for ( Int_t d = firstDetElemId; d <= lastDetElemId; ++d ) 
		{
			if ( !AliMUONSegmentationManager::IsValidDetElemId(d) )
	    {
	      AliWarning(Form("You are requesting an invalid detElemId = %d, I am skipping it",d));
	      continue;
	    }
			
			AliMUONVGeometryDESegmentation* slatSeg = 
	    new AliMUONTriggerSegmentationV2(d,kptypes[iplane]);
			      
			fDESegmentations->Add(slatSeg);
			
			segmentation->Add(d,slatSeg);
		}
		
		fMUON->SetSegmentationModel(ichamber,iplane+1,segmentation);
	}
	
  fMUON->SetResponseModel(ichamber,new AliMUONResponseTrigger);
	
  chamber.SetChargeCorrel(0); // same charge on both cathodes
}

//__________________________________________________________________________
void AliMUONFactoryV4::BuildStation3() 
{
  BuildChamber345(500,517);
  BuildChamber345(600,617);
}

//__________________________________________________________________________
void AliMUONFactoryV4::BuildStation4() 
{
  BuildChamber345(700,725);
  BuildChamber345(800,825);
}

//__________________________________________________________________________
void AliMUONFactoryV4::BuildStation5() 
{       
  BuildChamber345(900,925);
  BuildChamber345(1000,1025);
}

//__________________________________________________________________________
void AliMUONFactoryV4::BuildStation6() 
{ 
  BuildChamberTrigger(1100,1117);
  BuildChamberTrigger(1200,1217);
  BuildChamberTrigger(1300,1317);
  BuildChamberTrigger(1400,1417);
}

//__________________________________________________________________________
void AliMUONFactoryV4::Build(AliMUON* where, const char* what) 
{
  //
  // Construct MUON from chambers, segmentation and responses
  //

  fMUON = where;
  char tmp[20];
  strcpy(tmp, what);

  if (strcmp(tmp, "default")==0) {
    // Set default parameters
    fMUON->SetIshunt(0);
    fMUON->SetMaxStepGas(0.1);
    fMUON->SetMaxStepAlu(0.1);

    // Build stations
    BuildCommon();
    if (IsGeometryDefined(0))  BuildStation1();
    if (IsGeometryDefined(2))  BuildStation2();
    if (IsGeometryDefined(4))  BuildStation3();
    if (IsGeometryDefined(6))  BuildStation4();
    if (IsGeometryDefined(8))  BuildStation5();
    if (IsGeometryDefined(10)) BuildStation6();
  } 
  else
    AliDebug(0,"Non default version of MUON selected. You have to construct yourself the MUON elements !!");
}

//__________________________________________________________________________
void AliMUONFactoryV4::BuildStation(AliMUON* where, Int_t stationNumber) 
{
  //
  // Construct MUON from chambers, segmentation and responses
  //

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
