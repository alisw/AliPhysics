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
Revision 1.1  2002/06/16 17:08:19  hristov
First version of CRT


*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Alice COsmic Ray Trigger                                                  //
//                                                                           //
//  This class contains the functions for version 0 of the Cosmic Rays ALICE //
//  detector.                                                                //
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

  //
  //-- Create the hall
  CreateHall();

  Int_t  idrotm[2499];    // The rotation matrix.

  // idtmed[1099->1198] equivalent to fIdtmed[0->100]
  Int_t * idtmed = fIdtmed->GetArray() - 1099 ;

  // In order to generate the more correctly the modules (more datails)
  // we will create a box (air box) as a (sub)mother volume.
  // Inside this box we'll put the scintillator tiles, the PMTs the frame
  // and, maybe, some other things.

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
  //gMC->Gsdvn("CRT2", "CRT1", 2, 2);
  // Now divide each plane in 8 palettes
  //gMC->Gsdvn("CRT3", "CRT1", 8, 3);


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
    gMC->Gspos("CRT1", i, "CRTA", initX, initY, (i-step)*box[2], 0, "ONLY");
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
    gMC->Gspos("CRT1", i+4, "CRTA", initXside, initYside, (i-stepl)*box[2],
	       idrotm[232], "ONLY");
    stepl--;
  }

  // Put 4 modules on the right side of the magnet
  // The rotation matrix parameters for the right side.
  AliMatrix(idrotm[231], 90., 45., 90., 315., 180., 202.5);
  Int_t stepr = 4;
  for ( Int_t i = 1 ; i <= 4 ; i++ ) {
    gMC->Gspos("CRT1", i+8, "CRTA", -initXside, initYside, (i-stepr)*box[2],
	       idrotm[231], "ONLY");
    stepr--;
  }

}

//_____________________________________________________________________________
void AliCRTv0::CreateHall()
{

  Float_t r2;
  Float_t phid, phim, pbox[3], h, r, tspar[5];
  Float_t w1, dh, am, bm, dl,cm, hm, dr, dx, xl;
  Int_t idrotm[1999];
  Float_t trdpar[4], trapar[11];
  Float_t phi;
  
  Int_t *idtmed = fIdtmed->GetArray()-1899;

  // Because the BODY,  was "filled with vaccum", we're goin to superimpose
  // a volume with molasse, then using the boolenaq operations with volumes
  // we are going the replace, with air, the part closed by the hall walls

  // For the moment, just use the simple mother volume ...
  Float_t mola[3];
  mola[0] = 2000.;
  mola[1] = 2000.;
  mola[2] = 3000.;
  gMC->Gsvolu("CRTA", "BOX", idtmed[1905], mola, 3);
  gMC->Gspos("CRTA", 1, "ALIC", 0., 0., 0., 0, "MANY");

  // Now make a cilinder with the size of the hall (raughtly)
  // and fill it with air.
  Float_t hall[3];
  hall[0] = 0.;    // Inner radius
  hall[1] = 1170.; // Outer radius
  hall[2] = 2625. + 200.; // Hall lenght
  gMC->Gsvolu("CRTB", "TUBE", idtmed[1914], hall, 3);
  // perform the substraction SMOL - SALL, asigning SMOL as MANY and 
  // SALL as ONLY.
  gMC->Gspos("CRTB", 1, "CRTA", 0., 500., -2*375., 0, "ONLY");

  //     RB24/26 TUNNEL FLOOR 
  
  r   = 220.;
  h   = 140.;
  phi = TMath::ACos(h / r);
  xl  = r * TMath::Sin(phi);
  dr  = 100.;
  dh  = dr * TMath::Cos(phi);
  dl  = dr * TMath::Sin(phi);

  AliMatrix(idrotm[1900], 90., 0., 0., 0., 90., 90.);
  AliMatrix(idrotm[1901], 270., 0., 90., 90., 0., 0.);
  //     END WALL 
  
  pbox[0] = 1200.;
  pbox[1] = 1300.;
  pbox[2] = 60.;
  gMC->Gsvolu("CRTC", "BOX ", idtmed[1956], pbox, 3);
  gMC->Gspos("CRTC", 1, "CRTA", 0., 404., 1960, 0, "ONLY");
  
  //     hall floor, the interior part (left)
  
  phid      = 16.197;
  trdpar[0] = 700.;
  trdpar[1] = TMath::Tan(phid * kDegrad) * 190. + 700.;
  trdpar[2] = 550.;
  trdpar[3] = 95.;
  gMC->Gsvolu("CRTD", "TRD1", idtmed[1956], trdpar, 4);
  gMC->Gspos("CRTD", 1, "CRTA", 0., -801., 1350., idrotm[1900], "ONLY");

  //     hall floor, the outside part
  
  phid      = 16.197;
  trdpar[0] = 700.;
  trdpar[1] = TMath::Tan(phid * kDegrad) * 190. + 700.;
  trdpar[2] = 1325.;
  trdpar[3] = 95.;
  gMC->Gsvolu("CRTE", "TRD1", idtmed[1956], trdpar, 4);
  gMC->Gspos("CRTE", 2, "CRTA", 0., -801., -2125., idrotm[1900], "ONLY");
  
  //     hall side walls 

  // Interior walls  
  trapar[0] = 550.;
  trapar[1] = 0.;
  trapar[2] = 0.;
  trapar[3] = 1273.78/2;
  trapar[4] = 235.;
  trapar[5] = 50.;
  trapar[6] = TMath::ATan((trapar[4] - trapar[5]) / 2. / trapar[3]) * kRaddeg;
  trapar[7] = trapar[3];
  trapar[8] = trapar[4];
  trapar[9] = trapar[5];
  trapar[10] = trapar[6];
  dx = trapar[4] * 1.5 + 700. - trapar[5] * .5;
  gMC->Gsvolu("CRTF", "TRAP", idtmed[1956], trapar, 11);// interior wall
  gMC->Gspos("CRTF", 1, "CRTA", dx, -896+trapar[3],  1350., 0, "ONLY");
  gMC->Gspos("CRTF", 2, "CRTA",-dx, -896+trapar[3],  1350., idrotm[1901], "ONLY");

  // Exterior walls
  Float_t trapare[11];
  trapare[0] = 275.;
  for ( Int_t i = 1 ; i <= 10 ; i++ ) {
    trapare[i] = trapar[i];
  }
  gMC->Gsvolu("CRTG", "TRAP", idtmed[1956], trapare, 11);// exterior wall
  gMC->Gspos("CRTG", 1, "CRTA", dx, -896+trapar[3],  -1075., 0, "ONLY");
  gMC->Gspos("CRTG", 2, "CRTA",-dx, -896+trapar[3],  -1075., idrotm[1901], "ONLY");

  pbox[0] = 50.;
  pbox[1] = (500. - (trapar[3] * 2. - 896.)) / 2.;
  pbox[2] = 1625.;
  gMC->Gsvolu("CRTH", "BOX ", idtmed[1956], pbox, 3);
  gMC->Gspos("CRTH", 1, "CRTA",  1120., 500-pbox[1], 275., 0, "ONLY");
  gMC->Gspos("CRTH", 2, "CRTA", -1120., 500-pbox[1], 275., 0, "ONLY");

  //     slanted wall close to L3 magnet 
  
  phim = 45.;
  hm   = 790.;
  //rm   = hm / TMath::Cos(phim / 2. * kDegrad);
  am   = hm * TMath::Tan(phim / 2. * kDegrad);
  bm   = (hm + 76.) / hm * am;
  cm   = bm * 2. / TMath::Sqrt(2.);
  trapar[0] = 800.;
  trapar[1] = 0.;
  trapar[2] = 0.;
  trapar[3] = (1273.78 - cm) / 2.;
  trapar[4] = 235. - cm * TMath::Tan(phid * kDegrad) / 2.;
  trapar[5] = 50.;
  trapar[6] = TMath::ATan((trapar[4] - trapar[5]) / 2. / trapar[3]) * kRaddeg;
  trapar[7] = trapar[3];
  trapar[8] = trapar[4];
  trapar[9] = trapar[5];
  trapar[10] = trapar[6];
  w1 = trapar[4];
  dx = cm*TMath::Tan(phid * kDegrad) + 700. + trapar[4] * 1.5 - trapar[5] * .5;
  gMC->Gsvolu("CRTI", "TRAP", idtmed[1956], trapar, 11);
  r2 = cm - 896. + trapar[3];
  gMC->Gspos("CRTI", 1, "CRTA", dx, r2, 0., 0, "ONLY");
  gMC->Gspos("CRTI", 2, "CRTA",-dx, r2, 0., idrotm[1901], "ONLY");
  trapar[3]  = cm / 2.;
  trapar[4]  = w1 + cm / 2.;
  trapar[5]  = w1;
  trapar[6]  = TMath::ATan(.5) * kRaddeg;
  trapar[7]  = trapar[3];
  trapar[8]  = trapar[4];
  trapar[9]  = trapar[5];
  trapar[10] = trapar[6];
  dx = 1170. - trapar[4] * .5 - trapar[5] * .5;
  gMC->Gsvolu("CRTJ", "TRAP", idtmed[1956], trapar, 11);
  r2 = trapar[3] - 896.;
  gMC->Gspos("CRTJ", 1, "CRTA", dx, r2, 0., 0, "ONLY");
  gMC->Gspos("CRTJ", 2, "CRTA",-dx, r2, 0., idrotm[1901], "ONLY");

  tspar[0] = 1070.;
  tspar[1] = 1170.;
  tspar[2] = pbox[2];
  tspar[3] = 0.;
  tspar[4] = 180.;
  gMC->Gsvolu("CRTK", "TUBS", idtmed[1956], tspar, 5);
  gMC->Gspos("CRTK", 1, "CRTA", 0., 500., 300., 0, "ONLY");
  trdpar[0] = 1170 - trapar[4] * 2.;
  trdpar[1] = trdpar[0] + TMath::Tan(phim * kDegrad) * 76.;
  trdpar[2] = 800.;
  trdpar[3] = 38.;
  gMC->Gsvolu("CRTL", "TRD1", idtmed[1956], trdpar, 4);
  gMC->Gspos("CRTL", 1, "CRTA", 0., -858., 0., idrotm[1900], "ONLY");



  // Define the setion tube of the PX24, at the same level of hall
  // rotate the tubes around X, Z'=Y, Y'=-Z
  AliMatrix(idrotm[2001], 0., 0., 90., 0., 90., 90.);
  Float_t pxi[5];
  pxi[0] = 1150.;              // inside radius
  pxi[1] = 1250.;              // outside radius
  pxi[2] = 1300.;              // half lenght in Z
  pxi[3] = kRaddeg*TMath::ASin(tspar[0]/pxi[0]);//starting angle of the segment
  pxi[4] = 360.-pxi[3];               // ending angle of the segment
  gMC->Gsvolu("CRTM", "TUBS", idtmed[1956], pxi, 5);
  gMC->Gspos("CRTM", 1, "CRTA", 0., 404., -2300., idrotm[2001], "MANY");

  // Define the setion tube of the PX24, above the hall
  Float_t pxa[3];
  pxa[0] = pxi[0];
  pxa[1] = pxi[1];
  pxa[2] = 2550. - pxi[2]; // Half lenght
  gMC->Gsvolu("CRTN", "TUBE", idtmed[1956], pxa, 3);
  gMC->Gspos("CRTN", 1, "CRTA", 0.,pxi[2]+404+pxa[2], -2300., idrotm[2001], "MANY");
  // Fill this section with air.
  Float_t pxb[3];
  pxb[0] = 0.;
  pxb[1] = pxa[0];
  pxb[2] = pxa[2];
  gMC->Gsvolu("CRTO", "TUBE", idtmed[1914], pxb, 3);
  gMC->Gspos("CRTO", 1, "CRTA", 0., pxi[2]+404+pxa[2], -2300., idrotm[2001], "ONLY");


  // PM25 Acces shaft.
  Float_t pma[3];
  pma[0] = 910./2.;// Inner radius 
  pma[1] = pma[0] + 100.; // Outer Radius
  pma[2] = 5100./2.; // Half lenght 
  gMC->Gsvolu("CRTP", "TUBE", idtmed[1956], pma, 3);
  gMC->Gspos("CRTP", 1, "CRTA", -2100., 1654., 0., idrotm[2001], "ONLY");
  // Fill it with air.
  Float_t pmb[3];
  pmb[0] = 0.;
  pmb[1] = pma[0];
  pmb[2] = pma[2];
  gMC->Gsvolu("CRTQ", "TUBE", idtmed[1914], pmb, 3);
  gMC->Gspos("CRTQ", 1, "CRTA", -2100., 1654., 0., idrotm[2001], "ONLY");


  // PGC2 Acces shaft.
  Float_t pgc[3];
  pgc[0] = 1200./2.;// Inner Radius 
  pgc[1] = pgc[0] + 100.; // outer Radius
  pgc[2] = 5100./2.; // Half lenght 
  gMC->Gsvolu("CRTR", "TUBE", idtmed[1956], pgc, 3);
  gMC->Gspos("CRTR", 1, "CRTA", 375., 1654., 4850., idrotm[2001], "ONLY");
  // Fill it with air.
  Float_t pgd[3];
  pgd[0] = 0.;
  pgd[1] = pgc[0];
  pgd[2] = pgc[2];
  gMC->Gsvolu("CRTS", "TUBE", idtmed[1914], pgd, 3);
  gMC->Gspos("CRTS", 1, "CRTA", 375., 1654., 4850., idrotm[2001], "ONLY");

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

   gMC->Gsatt("ALIC","seen",0); // Mother volume, pit ceiling
   gMC->Gsatt("L3MO","seen",1); // L3 Magnet
   gMC->Gsatt("CRT1","seen",1); // Scintillators (air) box.

   // Draw the volumes for all the hall.
   gMC->Gsatt("CRTA","seen",1);
   gMC->Gsatt("CRTB","seen",1);
   gMC->Gsatt("CRTC","seen",1);
   gMC->Gsatt("CRTD","seen",1);
   gMC->Gsatt("CRTE","seen",1);
   gMC->Gsatt("CRTF","seen",1);
   gMC->Gsatt("CRTG","seen",1);
   gMC->Gsatt("CRTH","seen",1);
   gMC->Gsatt("CRTI","seen",1);
   gMC->Gsatt("CRTJ","seen",1);
   gMC->Gsatt("CRTK","seen",1);
   gMC->Gsatt("CRTL","seen",1);
   gMC->Gsatt("CRTM","seen",1);
   gMC->Gsatt("CRTN","seen",1);
   gMC->Gsatt("CRTO","seen",1);
   gMC->Gsatt("CRTP","seen",1);
   gMC->Gsatt("CRTQ","seen",1);
   gMC->Gsatt("CRTR","seen",1);
   gMC->Gsatt("CRTS","seen",1);


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
