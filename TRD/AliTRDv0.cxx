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
Revision 1.21  2002/02/20 14:01:40  hristov
Compare a TString with a string, otherwise the conversion cannot be done on Sun

Revision 1.20  2002/02/13 16:58:37  cblume
Bug fix reported by Jiri. Make atoi input zero terminated in StepManager()

Revision 1.19  2002/02/11 14:25:27  cblume
Geometry update, compressed hit structure

Revision 1.18  2000/11/30 17:38:08  cblume
Changes to get in line with new STEER and EVGEN

Revision 1.17  2000/11/01 14:53:21  cblume
Merge with TRD-develop

Revision 1.14.2.3  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.14.2.2  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.14.2.1  2000/09/18 13:48:18  cblume
Adapt to new AliTRDhit

Revision 1.16  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.15  2000/06/07 16:25:37  cblume
Try to remove compiler warnings on Sun and HP

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

#include <stdlib.h> 

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h> 
#include <TLorentzVector.h>

#include "AliRun.h"
#include "AliConst.h"
  
#include "AliTRDv0.h"
#include "AliTRDhit.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDv0)
  
//_____________________________________________________________________________
AliTRDv0::AliTRDv0():AliTRD() 
{
  //
  // AliTRDv0 default constructor
  //

  fHitsOn     = 0;

}

//_____________________________________________________________________________
AliTRDv0::AliTRDv0(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 0
  //

  fHitsOn     = 0;

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
  if (!frame) return;

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
  if (!(fHitsOn)) return;

  // Use only charged tracks and count them only once per volume
  if (gMC->TrackCharge()    && 
      gMC->IsTrackEntering()) {
    
    // Check on sensitive volume
    cIdCurrent = gMC->CurrentVolName();
    if (cIdSens == cIdCurrent[1]) {

      gMC->TrackPosition(p);
      for (Int_t i = 0; i < 3; i++) hits[i] = p[i];

      // The sector number (0 - 17)
      // The numbering goes clockwise and starts at y = 0
      Float_t phi = kRaddeg*TMath::ATan2(hits[0],hits[1]);
      if (phi < 90.) 
        phi = phi + 270.;
      else
        phi = phi -  90.;
      sec = ((Int_t) (phi / 20));

      // The plane and chamber number
      cIdChamber[0] = cIdCurrent[2];
      cIdChamber[1] = cIdCurrent[3];
      Int_t idChamber = atoi(cIdChamber);
      cha = ((Int_t) idChamber / kNplan);
      pla = ((Int_t) idChamber % kNplan);
      det = fGeometry->GetDetector(pla,cha,sec);

      AddHit(gAlice->CurrentTrack(),det,hits,0,kTRUE);       

    }

  }  

}
