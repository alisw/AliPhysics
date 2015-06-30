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

//=============================================================================================
//
//      Class describing geometry of one MFT half-disk support
//
//      Contact author: raphael.tieulent@cern.ch
//
//=============================================================================================

#include "TGeoManager.h"
#include "TGeoTube.h"
#include "AliMFTSupport.h"

ClassImp(AliMFTSupport)

//=============================================================================================

AliMFTSupport::AliMFTSupport():
TNamed(),fZin(0.), fZout(0.), fRmin(0.), fRmax(0.),fIsBottom(kFALSE){
  
  // default constructor
  
}
//=============================================================================================

AliMFTSupport::AliMFTSupport(Double_t zIn, Double_t zOut, Double_t rMin, Double_t rMax, Bool_t isBottom):
TNamed(),fZin(zIn), fZout(zOut), fRmin(rMin), fRmax(rMax),fIsBottom(isBottom){

  // default constructor
  
}


//=============================================================================================

AliMFTSupport::~AliMFTSupport() {
  
}


//=============================================================================================
TGeoVolume *  AliMFTSupport::CreateVolume(){
  // Create Shapes
  Double_t phiMin =0., phiMax=0.;
  
  phiMin = 180.*fIsBottom;
  phiMax = 180.*(fIsBottom+1.);
  
  // Create Shapes
  TGeoTubeSeg *support = new TGeoTubeSeg(fRmin, fRmax, (fZin - fZout)/2., phiMin, phiMax);
  
  // Get Mediums
  TGeoMedium *medBe  = gGeoManager->GetMedium("MFT_Be");
  
  // Create Volumes
  // Chip Volume
  TGeoVolume *supportVol = new TGeoVolume(Form("%s_Support",GetName()), support, medBe);
  supportVol->SetVisibility(kTRUE);
  supportVol->SetLineColor(kGray);
  
  
  return supportVol;
}
