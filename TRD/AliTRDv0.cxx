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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 0 -- fast simulator                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h> 

#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h> 
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliRun.h"
#include "AliMC.h"

#include "AliTRDgeometry.h"
#include "AliTRDhit.h"
#include "AliTRDv0.h"

ClassImp(AliTRDv0)
  
//_____________________________________________________________________________
AliTRDv0::AliTRDv0()
  :AliTRD() 
  ,fHitsOn(0)
{
  //
  // AliTRDv0 default constructor
  //

}

//_____________________________________________________________________________
AliTRDv0::AliTRDv0(const char *name, const char *title) 
  :AliTRD(name,title) 
  ,fHitsOn(0)
{
  //
  // Standard constructor for Transition Radiation Detector version 0
  //

}

//_____________________________________________________________________________
AliTRDv0::~AliTRDv0()
{
  //
  // AliTRDv0 destructor
  //

}

//_____________________________________________________________________________
void AliTRDv0::CreateGeometry()
{
  //
  // Create the GEANT geometry for the Transition Radiation Detector - Version 0
  // This version covers the full azimuth. 
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* frame = gAlice->GetModule("FRAME");
  if (!frame) {
    AliError("TRD needs FRAME to be present\n");
    return;
  }

  // Define the chambers
  AliTRD::CreateGeometry();

}

//_____________________________________________________________________________
void AliTRDv0::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector
  //

  AliTRD::CreateMaterials();

}

//_____________________________________________________________________________
void AliTRDv0::Init() 
{
  //
  // Initialize Transition Radiation Detector after geometry is built
  //

  AliTRD::Init();

  AliDebug(1,"          Fast simulator\n\n");
  AliDebug(1,"++++++++++++++++++++++++++++++++++++++++++++++");
  
}

//_____________________________________________________________________________
void AliTRDv0::StepManager()
{
  //
  // Procedure called at every step in the TRD
  // Fast simulator. If switched on, a hit is produced when a track
  // crosses the border between amplification region and pad plane.
  //

  Int_t   pla = 0; 
  Int_t   cha = 0;
  Int_t   sec = 0; 

  Float_t hits[3];
  Int_t   det;

  TLorentzVector p;

  // Use pad plane as sensitive volume
  TString  cIdSens = "L";
  TString  cIdCurrent;
  Char_t   cIdChamber[3];
           cIdChamber[2] = 0;

  const Int_t kNplan = AliTRDgeometry::Nplan();

  // Writing out hits enabled?
  if (!(fHitsOn)) {
    return;
  }

  // Use only charged tracks and count them only once per volume
  if (gMC->TrackCharge()    && 
      gMC->IsTrackEntering()) {
    
    // Check on sensitive volume
    cIdCurrent = gMC->CurrentVolName();
    if (cIdSens == cIdCurrent[1]) {

      gMC->TrackPosition(p);
      for (Int_t i = 0; i < 3; i++) {
        hits[i] = p[i];
      }

      // The sector number (0 - 17)
      // The numbering goes clockwise and starts at y = 0
      Float_t phi = kRaddeg*TMath::ATan2(hits[0],hits[1]);
      if (phi < 90.0) {
        phi = phi + 270.0;
      }
      else {
        phi = phi -  90.0;
      }
      sec = ((Int_t) (phi / 20.0));

      // The plane and chamber number
      cIdChamber[0]   = cIdCurrent[2];
      cIdChamber[1]   = cIdCurrent[3];
      Int_t idChamber = atoi(cIdChamber);
      cha = ((Int_t) idChamber / kNplan);
      pla = ((Int_t) idChamber % kNplan);
      det = fGeometry->GetDetector(pla,cha,sec);

      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),det,hits,0,kTRUE);       

    }

  }  

}
