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

/*
$Log$
Revision 1.14  2000/02/28 19:10:26  cblume
Include the new TRD classes

Revision 1.13.4.1  2000/02/28 18:01:53  cblume
Change to new hit version and introduce geometry class

Revision 1.13  1999/11/05 22:50:28  fca
Do not use Atan, removed from ROOT too

Revision 1.12  1999/11/02 16:35:56  fca
New version of TRD introduced

Revision 1.11  1999/11/01 20:41:51  fca
Added protections against using the wrong version of FRAME

Revision 1.10  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 0 -- fast simulator                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDfullClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>

#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"
  
#include "AliTRDv0.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDv0)
 
//_____________________________________________________________________________
AliTRDv0::AliTRDv0(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 0
  //

  fIdSens     = 0;
  fHitsOn     = 0;

  fIdChamber1 = 0;
  fIdChamber2 = 0;
  fIdChamber3 = 0;

}

//_____________________________________________________________________________
void AliTRDv0::CreateGeometry()
{
  //
  // Create the GEANT geometry for the Transition Radiation Detector - Version 0
  // This version covers the full azimuth. 
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* FRAME = gAlice->GetModule("FRAME");
  if (!FRAME) return;

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

  // Identifier of the sensitive volume (amplification region)
  fIdSens     = gMC->VolId("UL06");

  // Identifier of the TRD-driftchambers
  fIdChamber1 = gMC->VolId("UCIO");
  fIdChamber2 = gMC->VolId("UCIM");
  fIdChamber3 = gMC->VolId("UCII");

  printf("          Fast simulator\n\n");
  for (Int_t i = 0; i < 80; i++) printf("*");
  printf("\n");
  
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
  Int_t   iIdSens, icSens; 
  Int_t   iIdChamber, icChamber;

  Float_t hits[4];
  Int_t   det[1];

  TLorentzVector p;
  TClonesArray  &lhits = *fHits;

  // Writing out hits enabled?
  if (!(fHitsOn)) return;

  // Use only charged tracks and count them only once per volume
  if (gMC->TrackCharge()    && 
      gMC->IsTrackExiting()) {
    
    // Check on sensitive volume
    iIdSens = gMC->CurrentVolID(icSens);
    if (iIdSens == fIdSens) { 

      gMC->TrackPosition(p);
      for (Int_t i = 0; i < 3; i++) hits[i] = p[i];
      // No charge created
      hits[3] = 0;

      // The sector number (0 - 17)
      // The numbering goes clockwise and starts at y = 0
      Float_t phi = kRaddeg*TMath::ATan2(hits[0],hits[1]);
      if (phi < 90.) 
        phi = phi + 270.;
      else
        phi = phi -  90.;
      sec = ((Int_t) (phi / 20));

      // The chamber number 
      //   0: outer left
      //   1: middle left
      //   2: inner
      //   3: middle right
      //   4: outer right
      iIdChamber = gMC->CurrentVolOffID(1,icChamber);
      if      (iIdChamber == fIdChamber1)
        cha = (hits[2] < 0 ? 0 : 4);
      else if (iIdChamber == fIdChamber2)       
        cha = (hits[2] < 0 ? 1 : 3);
      else if (iIdChamber == fIdChamber3)       
        cha = 2;

      // The plane number (0 - 5)
      pla = icChamber - TMath::Nint((Float_t) (icChamber / 7)) * 6 - 1;

      det[0] = fGeometry->GetDetector(pla,cha,sec);
      new(lhits[fNhits++]) AliTRDhit(fIshunt
                                    ,gAlice->CurrentTrack()
                                    ,det
                                    ,hits);

    }

  }  

}
