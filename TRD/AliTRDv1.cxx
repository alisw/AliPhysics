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
Revision 1.10  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 1 -- coarse simulation             //
//  This version has two detector arms, leaving the space in front of the    //
//  HMPID and PHOS empty                                                     //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDv1Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>

#include "AliTRDv1.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"
 
ClassImp(AliTRDv1)

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for the Transition Radiation Detector version 1
  //

  fIdSens     = 0;
  fHitsOn     = 0;

  fIdSpace1   = 0;
  fIdSpace2   = 0;
  fIdSpace3   = 0;

  fIdChamber1 = 0;
  fIdChamber2 = 0;
  fIdChamber3 = 0;

}
 
//_____________________________________________________________________________
void AliTRDv1::CreateGeometry()
{
  //
  // Create the GEANT geometry for the Transition Radiation Detector - Version 1
  // This version covers only part of the azimuth. 
  //
  // Author:  Christoph Blume (C.Blume@gsi.de) 20/07/99 
  //

  Float_t xpos, ypos, zpos;

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* FRAME = gAlice->GetModule("FRAME");
  if (!FRAME) return;

  // Define the chambers
  AliTRD::CreateGeometry();

  // Position the the TRD-sectors only in one TRD-volume in the spaceframe
  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gspos("TRD ",1,"BTR1",xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDv1::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector version 1
  //

  AliTRD::CreateMaterials();

}

//_____________________________________________________________________________
void AliTRDv1::Init() 
{
  //
  // Initialise the Transition Radiation Detector after the geometry is built
  //

  printf("**************************************"
	 "  TRD  "
	 "**************************************\n");
  printf("\n     Version 1 of TRD initialing, "
	 "with openings for PHOS and RICH\n\n");

  AliTRD::Init();

  //
  // Check that FRAME is there otherwise we have no place where to
  // put TRD
  AliModule* FRAME=gAlice->GetModule("FRAME");
  if(!FRAME) {
    Error("Ctor","TRD needs FRAME to be present\n");
    exit(1);
  } else 
    if(FRAME->IsVersion()!=0) {
      Error("Ctor","FRAME version 0 needed with this version of TRD\n");
      exit(1);
    }

  for (Int_t i = 0; i < 80; i++) printf("*");
  printf("\n");
  
  // Identifier of the sensitive volume (amplification region)
  fIdSens     = gMC->VolId("UL06");

  // Identifier of the TRD-spaceframe volumina
  fIdSpace1   = gMC->VolId("B028");
  fIdSpace2   = gMC->VolId("B029");
  fIdSpace3   = gMC->VolId("B030");

  // Identifier of the TRD-driftchambers
  fIdChamber1 = gMC->VolId("UCIO");
  fIdChamber2 = gMC->VolId("UCIM");
  fIdChamber3 = gMC->VolId("UCII");

  printf("**************************************"
	 "  TRD  "
	 "**************************************\n");
}

//_____________________________________________________________________________
void AliTRDv1::StepManager() 
{
  //
  // Procedure called at every step in the TRD
  // Fast simulator. If switched on, a hit is produced when a track
  // crosses the border between amplification region and pad plane.
  //

  Int_t   vol[3]; 
  Int_t   iIdSens, icSens; 
  Int_t   iIdSpace, icSpace;
  Int_t   iIdChamber, icChamber;

  Int_t   secMap1[10] = {  3,  7,  8,  9, 10, 11,  2,  1, 18, 17 };
  Int_t   secMap2[ 5] = { 16, 15, 14, 13, 12 };
  Int_t   secMap3[ 3] = {  5,  6,  4 };

  Float_t hits[4];

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

      iIdSpace   = gMC->CurrentVolOffID(4,icSpace  );
      iIdChamber = gMC->CurrentVolOffID(1,icChamber);

      // The sector number
      if      (iIdSpace == fIdSpace1) 
        vol[0] = secMap1[icSpace-1];
      else if (iIdSpace == fIdSpace2) 
        vol[0] = secMap2[icSpace-1];
      else if (iIdSpace == fIdSpace3) 
        vol[0] = secMap3[icSpace-1];

      // The chamber number 
      //   1: outer left
      //   2: middle left
      //   3: inner
      //   4: middle right
      //   5: outer right
      if      (iIdChamber == fIdChamber1)
        vol[1] = (hits[2] < 0 ? 1 : 5);
      else if (iIdChamber == fIdChamber2)       
        vol[1] = (hits[2] < 0 ? 2 : 4);
      else if (iIdChamber == fIdChamber3)       
        vol[1] = 3;

      // The plane number
      vol[2] = icChamber - TMath::Nint((Float_t) (icChamber / 7)) * 6;

      new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);

    }

  }  

}
