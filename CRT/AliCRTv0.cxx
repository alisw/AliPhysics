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
Revision 1.2  2002/07/09 08:45:35  hristov
Old style include files needed on HP (aCC)

Revision 1.1  2002/06/16 17:08:19  hristov
First version of CRT


*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// ALICE Cosmic Ray Trigger                                                  //
//                                                                           //
//  This class contains the functions for version 0 of the ALICE Cosmic Ray  //
//  Trigger.                                                                 //
//
//   Authors:
//
//   Arturo Fernandez <afernand@fcfm.buap.mx>
//   Enrique Gamez    <egamez@fcfm.buap.mx>
//
//   Universidad Autonoma de Puebla
//
//
//Begin_Html
/*
<img src="picts/AliCRTv0Class.gif">
</pre>
<br clear=left>
<p>The responsible person for this module is
<a href="mailto:egamez@fcfm.buap.mx">Enrique Gamez</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <iostream.h>

#include <TMath.h>
#include <TGeometry.h>
#include <TTUBE.h>
#include <TNode.h>
#include <TLorentzVector.h>

#include "AliCRTv0.h"
#include "AliCRTConstants.h"

#include "AliRun.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliConst.h"
#include "AliPDG.h"

ClassImp(AliCRTv0)
 
//_____________________________________________________________________________
AliCRTv0::AliCRTv0() : AliCRT()
{
  //
  // Default constructor for CRT
  //
  fMucur = 0;
}
 
//_____________________________________________________________________________
AliCRTv0::AliCRTv0(const char *name, const char *title)
  : AliCRT(name,title)
{
  //
  // Standard constructor for CRT
  //
  //Begin_Html
  /*
    <img src="picts/AliCRTv0.gif">
  */
  //End_Html
}

//_____________________________________________________________________________
void AliCRTv0::BuildGeometry()
{

}

//_____________________________________________________________________________
void AliCRTv0::CreateGeometry()
{
  //
  // Create geometry for the CRT array
  //

  Int_t  idrotm[2499];    // The rotation matrix.

  Int_t * idtmed = fIdtmed->GetArray() - 1099 ;

  //
  // Molasse
  //

  // Exactly above the hall
  Float_t tspar[5];
  tspar[0] = 1170.;
  tspar[1] = 1170. + 375.;
  tspar[2] = (1900.+1150.)/2.+100.;
  tspar[3] = 0.;
  tspar[4] = 180.;
  gMC->Gsvolu("CMO1", "TUBS", idtmed[1103], tspar, 5);
  gMC->Gspos("CMO1", 1, "ALIC", 0., 500., 1900.-tspar[2]+400., 0, "MANY");

  Float_t tbox[3];
  tbox[0] = 1250.;
  tbox[1] = (4420. - 1670.)/2.;
  tbox[2] = (1900.+1150.)/2. + 200.;
  gMC->Gsvolu("CM12", "BOX", idtmed[1103], tbox, 3);
  gMC->Gspos("CM12", 1, "ALIC", 0., 4420. -tbox[1], 1900.-tbox[2]+400., 0, "MANY");

  AliMatrix(idrotm[2003], 0., 0., 90., 0., 90., 90.);
  // Along the PM25
  Float_t tube[3];
  tube[0] = 455. + 100.;
  tube[1] = 555. + 375.;
  tube[2] = (5150. - 1166.)/2.;
  gMC->Gsvolu("CMO2", "TUBE", idtmed[1103], tube, 3);
  gMC->Gspos("CMO2", 1, "ALIC", -2100., 4420.-tube[2], 0., idrotm[2003], "MANY");


  // Along the PGC2
  tube[0] = 650.;
  tube[1] = 2987.7;
  tube[2] = (5150. - 690.)/2.;
  gMC->Gsvolu("CMO3", "TUBE", idtmed[1103], tube, 3);
  gMC->Gspos("CMO3", 1, "ALIC", 375., 4420.-tube[2], 1900.+2987.7, idrotm[2003], "MANY");
  // Behind the PGC2 up to the end of the M. volume.
  tbox[0] = 12073.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (12073. - 1900.+2987.7+650.)/2.;
  gMC->Gsvolu("CMO7", "BOX", idtmed[1103], tbox, 3);
  gMC->Gspos("CMO7", 1, "ALIC", 0., 4420.-tbox[1], 1900.+2987.7+650.+tbox[2], 0, "MANY");

  // Along the PX24 , upper part.
  tube[0] = 1250.;
  tube[1] = 2300;
  tube[2] = 2575. - 1300. + 95.;
  gMC->Gsvolu("CMO4", "TUBE", idtmed[1103], tube, 3);
  gMC->Gspos("CMO4", 1, "ALIC", 0., 404.+1300.+tube[2], -2300., idrotm[2003], "MANY");

  // Along the PX24 , lower part
  tspar[0] = 1250.;
  tspar[1] = 2300;
  tspar[2] = 1300.;
  tspar[3] = kRaddeg*TMath::ASin(1070./1150.);
  tspar[4] = 360. - tspar[3];
  gMC->Gsvolu("CMO5", "TUBS", idtmed[1103], tspar, 5);
  gMC->Gspos("CMO5", 1, "ALIC", 0., 404., -2300., idrotm[2003], "MANY");
  // behind the PX24
  tbox[0] = 12073.;
  tbox[1] = 2575. + 95.;
  tbox[2] = 8523./2.;
  gMC->Gsvolu("CMO6", "BOX", idtmed[1103], tbox, 3);
  gMC->Gspos("CMO6", 1, "ALIC", 0., 4420.-tbox[1], -3550.-tbox[2], 0, "MANY");


  // On the right side of th hall
  tbox[0] = (12073. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (8437.7+650.)/2.;
  gMC->Gsvolu("CMO8", "BOX", idtmed[1103], tbox, 3);
  gMC->Gspos("CMO8", 1, "ALIC", 1250.+tbox[0], 4420.-tbox[1], -3550.+tbox[2], 0, "MANY");

  // on the left side of the hall, behind 
  tbox[0] = (12073. - 2755.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (8437.7+650.)/2.;
  gMC->Gsvolu("CMO9", "BOX", idtmed[1103], tbox, 3);
  gMC->Gspos("CMO9", 1, "ALIC", -2755.-tbox[0], 4420.-tbox[1], -3550.+tbox[2], 0, "MANY");


  // Molasse betwen the PX24 & PM25 on the left side.
  tbox[0] = (2755. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (3550. - 555.)/2.;
  gMC->Gsvolu("CM10", "BOX", idtmed[1103], tbox, 3);
  gMC->Gspos("CM10", 1, "ALIC", -1250.-tbox[0], 4420.-tbox[1], -tbox[2]-555., 0, "MANY");


  // Molasse betwen the PGC2 & PM25 on the left side.
  tbox[0] = (2755. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (1900.+2987.7 - 555. + 650.)/2.;
  gMC->Gsvolu("CM11", "BOX", idtmed[1103], tbox, 3);
  gMC->Gspos("CM11", 1, "ALIC", -1250.-tbox[0], 4420.-tbox[1], 555.+tbox[2], 0, "MANY");

  //
  // Scintillators

  Float_t box[3];
  box[0] = AliCRTConstants::fgCageLenght/2.; // Half Length of the box along the X axis, cm.
  box[1] = AliCRTConstants::fgCageHeight/2.; // Half Length of the box along the Y axis, cm.
  box[2] = AliCRTConstants::fgCageWidth/2.;  // Half Length of the box along the Z axis, cm.


  // Define the Scintillators. as a big box.
  Float_t scint[3];
  scint[0] = AliCRTConstants::fgActiveAreaLenght/2.;       // Half Length in X
  scint[1] = AliCRTConstants::fgActiveAreaHeight/2.;       // Half Length in Y
  scint[2] = AliCRTConstants::fgActiveAreaWidth/2.;        // Half Length in Z
  gMC->Gsvolu("CRT1", "BOX ", idtmed[1102], scint, 3);  // Scintillators
  // Divide the modules in 2 planes.
  gMC->Gsdvn("CRT2", "CRT1", 2, 2);
  // Now divide each plane in 8 palettes
  gMC->Gsdvn("CRT3", "CRT1", 8, 3);


  //
  // Define the coordinates where the draw will begin.
  //

  //
  // -- X axis.
  // we'll start dawing from the center.
  Float_t initX = 0.;

  //
  // -- Y axis
  Float_t gapY   = 30.;        // 30 cms. above the barrel.
  // For the height we staimate the from the center of the ceiling,
  // if were a cilinder, must be about 280cm.
  Float_t barrel = 790.; // Barrel radius.
  Float_t height  = barrel + gapY - 30.;
  Float_t initY = height;

  //
  // -- Z axis.
  // we'll start dawing from the center.

  //
  // Put 4 modules on the top of the magnet
  Int_t step = 4;
  for ( Int_t i = 1 ; i <= 4 ; i++ ) {
    gMC->Gspos("CRT1", i, "ALIC", initX, initY, (i-step)*box[2], 0, "ONLY");
    step--;
  }

  // Modules on the barrel sides.
  // Because the openenig angle for each face is 22.5, and if we want to
  //    put the modules right in the middle
  Float_t xtragap = 10.;
  Float_t initXside = (height+xtragap)*TMath::Sin(2*22.5*kDegrad); //rigth side
  Float_t initYside = (height+xtragap)*TMath::Cos(2*22.5*kDegrad);

  // Put 4 modules on the left side of the magnet
  // The rotation matrix parameters, for the left side.
  AliMatrix(idrotm[232], 90., 315., 90., 45., 0., 337.5);
  Int_t stepl = 4;
  for ( Int_t i = 1 ; i <= 4 ; i++ ) {
    gMC->Gspos("CRT1", i+4, "ALIC", initXside, initYside, (i-stepl)*box[2],
	       idrotm[232], "ONLY");
    stepl--;
  }

  // Put 4 modules on the right side of the magnet
  // The rotation matrix parameters for the right side.
  AliMatrix(idrotm[231], 90., 45., 90., 315., 180., 202.5);
  Int_t stepr = 4;
  for ( Int_t i = 1 ; i <= 4 ; i++ ) {
    gMC->Gspos("CRT1", i+8, "ALIC", -initXside, initYside, (i-stepr)*box[2],
	       idrotm[231], "ONLY");
    stepr--;
  }

}

//_____________________________________________________________________________
void AliCRTv0::CreateMaterials()
{
  //
  //--
  //

  // Use the standard materials.
  AliCRT::CreateMaterials();
}


//_____________________________________________________________________________
void AliCRTv0::DrawDetector()
{

}

//_____________________________________________________________________________
void AliCRTv0::DrawModule()
{
  //
  // Draw a shaded view of the L3 magnet
  //
   cout << "AliCRTv0::DrawModule() : Drawing the module" << endl;

   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);

   gMC->Gsatt("ALIC","seen",0);
   gMC->Gsatt("L3MO","seen",1); // L3 Magnet
   gMC->Gsatt("CRT1","seen",1); // Scintillators

   // Draw the molasse volumes
   gMC->Gsatt("CMO1","seen",0); // Exactly above the HALL
   gMC->Gsatt("CMO2","seen",0); // Molasse, along the PM25
   gMC->Gsatt("CMO3","seen",0); // molasse along the PGC2
   gMC->Gsatt("CMO4","seen",0); // Molasse, behind the PX24 upper part
   gMC->Gsatt("CMO5","seen",0); // molasse behind px24, lower part
   gMC->Gsatt("CMO6","seen",0); // behind the PX24
   gMC->Gsatt("CMO7","seen",0); // behind the PGC2
   gMC->Gsatt("CMO8","seen",0); // on the right side.
   gMC->Gsatt("CMO9","seen",0); // on the left side.
   gMC->Gsatt("CM10","seen",0); // betwen PX24 & PM25.
   gMC->Gsatt("CM11","seen",0); // betwen PGC2 & PM25.
   gMC->Gsatt("CM12","seen",0); // box above the hall.

   gMC->Gdopt("hide", "on");
   gMC->Gdopt("edge","off");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox("ALIC", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 10, 9.5, .009, .009);
   gMC->Gdhead(1111, "View of CRT(ACORDE)");
   gMC->Gdman(18, 4, "MAN");


}

//_____________________________________________________________________________
void AliCRTv0::Init()
{
  //
  // Initialise L3 magnet after it has been built
  Int_t i;
  //
  if(fDebug) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" CRTv0_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the CRTv0 initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }

}

//_____________________________________________________________________________
void AliCRTv0::StepManager()
{
  //
  // Called for every step in the CRT Detector
  //
  Float_t hits[12];
  Int_t   vol[5];

  // Check if this is the last step of the track in the current volume
  Bool_t laststepvol = gMC->IsTrackEntering();
  // Obtain the medium
  TLorentzVector xyz;
  gMC->TrackPosition(xyz);
  TLorentzVector pxyz;
  gMC->TrackMomentum(pxyz);

  if ( laststepvol && (strcmp(gMC->CurrentVolName(),"CRT1") == 0) ) {
    if (  gMC->TrackCharge() != 0 || gMC->TrackPid() == kGamma ) {
      Float_t vert[3];

      hits[0] = fMucur++;

      if ( (gMC->TrackPid() != kMuonPlus) && (gMC->TrackPid() != kMuonMinus)) {
	hits[1] = -(Float_t)gMC->TrackPid();
      } else {
	hits[1] = (Float_t)gMC->TrackPid();
      }

      TLorentzVector xyz;
      gMC->TrackPosition(xyz);
      TLorentzVector pxyz;
      gMC->TrackMomentum(pxyz);

      hits[2]  = xyz[0]; // X pit
      hits[3]  = xyz[1]; // Y pit
      hits[4]  = xyz[2]; // Z pit
      hits[5]  = pxyz[0]; // pxug
      hits[6] = pxyz[1]; // pyug
      hits[7] = pxyz[2]; // pzug

      hits[8] = gMC->GetMedium(); // layer
      hits[9] = vert[0]; // xver
      hits[10] = vert[1]; // yver
      hits[11] = vert[2]; // zver
    }
  }

  // Store the hit.
  AddHit(gAlice->CurrentTrack(),vol, hits);
}
