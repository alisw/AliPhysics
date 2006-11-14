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
// ALICE Cosmic Ray Trigger                                                  //
//                                                                           //
//  This class contains the functions for version 0 of the ALICE Cosmic Ray  //
//  Trigger. This vesion is suposed to work as standalone module             //
//                                                                           //
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
<img src="picts/AliCRTv1Class.gif">
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

#include "AliCRTv1.h"

#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TVirtualMC.h>

#include "AliRun.h"
#include "AliConst.h"

#include "AliCRThit.h"
#include "AliCRTModule.h"
#include "AliCRTConstants.h"
#include "AliMC.h"
#include "AliLog.h"

ClassImp(AliCRTv1)
 
//_____________________________________________________________________________
AliCRTv1::AliCRTv1()
  : AliCRT()
{
  //
  // Default constructor
  //
  fIshunt = 0;
  fHits = 0;
}
 
//_____________________________________________________________________________
AliCRTv1::AliCRTv1(const char *name, const char *title)
  : AliCRT(name, title)
{
  //
  // Standard constructor
  //
  //Begin_Html
  /*
    <img src="picts/AliCRTv1.gif">
  */
  //End_Html
  fIshunt = 1; // All hits are associated with primary particles  

  fHits =  new TClonesArray("AliCRThit",400);
  gAlice->GetMCApp()->AddHitList(fHits);

  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}

//_____________________________________________________________________________
AliCRTv1::~AliCRTv1()
{
  //
  // Default destructor
  //
}

//_____________________________________________________________________________
void AliCRTv1::CreateMaterials()
{
  //
  // Create Materials.
  // Use the parent class definition of the materials
  //
  AliCRT::CreateMaterials();
}

//_____________________________________________________________________________
void AliCRTv1::CreateGeometry()
{
  //
  // Create geometry for the CRT array
  //

  Int_t  idrotm[2499];    // The rotation matrix.
  Int_t* idtmed = fIdtmed->GetArray() - 1099 ;
  AliCRTConstants* crtConstants = AliCRTConstants::Instance();

  // Create the mother volume, the one which will contain all the material
  // above the hall.
  Float_t pbox[3];
  pbox[0] = AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  //pbox[0] = 12073;
  pbox[1] = crtConstants->Depth();
  pbox[2] = pbox[0];
  gMC->Gsvolu("CRT", "BOX", idtmed[1114], pbox, 3);
  gMC->Gspos("CRT", 1, "ALIC", 0, 0, 0, 0, "ONLY");

  // Shafts.
  this->CreateShafts();

  // Molasse.
  this->CreateMolasse();

  // This volume can be seen as the volume which ACORDE will ocupate
  // above the upper face of the L3 magnet. Inside this volume the detectors
  // aboce the magnet will be, then there will be two copies of this volume,
  // one for each side.
  Float_t box[3];
  //box[0] = 2*crtConstants->MagMinRadius()*TMath::Sin(kDegrad*22.5);
  box[0] = crtConstants->MagMinRadius()*TMath::Sin(kDegrad*22.5);
  box[1] = crtConstants->MagMaxRadius() - crtConstants->MagMinRadius();
  box[2] = crtConstants->MagnetLenght()/2;
  gMC->Gsvolu("CRT1", "BOX", idtmed[1134], box, 3);

  // Check if the AliCRTModule instance have been set, otherwise
  // use the default values
  if ( !fModule ) {
    Info("CreateGeometry", "Using default dimensions");
    fModule = new AliCRTModule("CRTmod", "Default module dimensions");
  }

  // The full module volume.
  // This volume will be ocupied by all the material of the module
  // the scintillators, the aluminium frame, etc.
  box[0] = fModule->FrameLength()/2;
  box[1] = fModule->FrameThickness()/2;
  box[2] = fModule->FrameWidth()/2;
  gMC->Gsvolu("CRT2", "BOX", idtmed[1114], box, 3);

  // The scintillators
  box[0] = crtConstants->SinglePaletteLenght()/4;
  box[1] = crtConstants->SinglePaletteHeight();
  box[2] = crtConstants->SinglePaletteWidth()/2;
  gMC->Gsvolu("CRT3", "BOX", idtmed[1112], box, 3);
  gMC->Gspos("CRT3", 1, "CRT2", 0, 2, 0, 0, "ONLY");

  // The metallic frame
  box[0] = fModule->FrameLength()/2;
  box[1] = fModule->FrameThickness()/2;
  box[2] = 2;
  gMC->Gsvolu("CRT4", "BOX", idtmed[1108], box, 3);
  gMC->Gspos("CRT4", 1, "CRT2", 0, 0,  13 - box[2], 0, "MANY");
  gMC->Gspos("CRT4", 2, "CRT2", 0, 0, -13 + box[2], 0, "MANY");

  box[0] = 2;
  box[1] = fModule->FrameThickness()/2;
  box[2] = fModule->FrameWidth()/2;
  gMC->Gsvolu("CRT5", "BOX", idtmed[1108], box, 3);
  gMC->Gspos("CRT5", 1, "CRT2",  140 - box[0], 0, 0, 0, "MANY");
  gMC->Gspos("CRT5", 2, "CRT2", -140 + box[0], 0, 0, 0, "MANY");

  // The support bars
  box[0] = 2;
  box[1] = fModule->FrameThickness()/2;
  box[2] = 500;
  gMC->Gsvolu("CRT6", "BOX", idtmed[1108], box, 3);

  // Now put into the volume CR11 all the above volumes.
  // 20 scintillation modules
  // 4 support bars
  Int_t copyNumber = 0;
  for ( Int_t k = 0; k < fModule->NumberOfRows(); k++ ) {
    Float_t zCoordinate = k*fModule->ZGap() - 450;
    gMC->Gspos("CRT2",++copyNumber,"CRT1",-150, 15, zCoordinate, 0, "MANY");
    gMC->Gspos("CRT2",++copyNumber,"CRT1",150, 15, zCoordinate, 0, "MANY");

  }

  // Put the support bars
  gMC->Gspos("CRT6", 1, "CRT1",  -75, 5, 0, 0, "ONLY");
  gMC->Gspos("CRT6", 2, "CRT1", -225, 5, 0, 0, "ONLY");
  gMC->Gspos("CRT6", 3, "CRT1",   75, 5, 0, 0, "ONLY");
  gMC->Gspos("CRT6", 4, "CRT1",  225, 5, 0, 0, "ONLY");

  // Now put a copy of CR11 on the 3 upper faces of the magnet
  // In the right side side of the magnet
  AliMatrix(idrotm[231], 90, 45, 90, 135, 0, 0);
  // In the left side side of the magnet
  AliMatrix(idrotm[232], 90, 315, 90, 45, 0, 0);

  Float_t x = crtConstants->MagMaxRadius();
  gMC->Gspos("CRT1", 1, "ALIC", 0, x, 0, 0, "MANY");
  gMC->Gspos("CRT1", 2, "ALIC", -x*TMath::Sin(kDegrad*45), x*TMath::Cos(kDegrad*45), 0, idrotm[231], "MANY");
  gMC->Gspos("CRT1", 3, "ALIC",  x*TMath::Sin(kDegrad*45), x*TMath::Cos(kDegrad*45), 0, idrotm[232], "MANY");

}

//_____________________________________________________________________________
void AliCRTv1::CreateMolasse()
{
  //
  //
  //
  Int_t  idrotm[2499];    // The rotation matrix.
  Int_t* idtmed = fIdtmed->GetArray() - 1099 ;

  Float_t px24radius = 2300/2;
  Float_t px24X = 0;
  //Float_t px24Y = ;
  Float_t px24Z = 2300;

  Float_t pm25radius = 910/2;
  Float_t pm25X = 2100;
  //Float_t pm25Y = ;
  Float_t pm25Z = 0;

  Float_t pgc2radius = 1100/2;
  Float_t pgc2X = -375;
  //Float_t pgc2Y = ;
  Float_t pgc2Z = -(1900 + 2987.7);

  Float_t concreteWidth = 100; // Standard width of the hall walls.


  // Create a local mother volume.
  Float_t pbox[3];
  pbox[0] = AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[1] = AliCRTConstants::Instance()->Depth()/2;
  pbox[2] = pbox[0];
  gMC->Gsvolu("CMO1", "BOX", idtmed[1114], pbox, 3);

  // Now put the molasse exactly above the hall. OK
  // Above the ceiling
  Float_t ptubs[5];
  ptubs[0] = 1170;
  ptubs[1] = 2100 - pm25radius;
  ptubs[2] = 1900/2 + px24radius;
  ptubs[3] = 0;
  ptubs[4] = 180;
  gMC->Gsvolu("CMO2", "TUBS", idtmed[1123], ptubs, 5);
  gMC->Gspos("CMO2", 1, "CMO1", 0, 500-AliCRTConstants::Instance()->Depth()/2, ptubs[2]-1900, 0, "MANY");

  // Molasse around the RB24/26 Wall. OK
  ptubs[0] = 220 + 1600;
  ptubs[1] = AliCRTConstants::Instance()->Depth() - ptubs[0];
  ptubs[2] = 2987.7/2 - 1100/4 - concreteWidth/2;
  ptubs[3] = 0;
  ptubs[4] = 180;
  gMC->Gsvolu("CMO3", "TUBS", idtmed[1123], ptubs, 5);
  gMC->Gspos("CMO3", 1, "CMO1", 70, 40-AliCRTConstants::Instance()->Depth()/2, -1900 - ptubs[2], 0, "MANY");

  // A big block above the RB24/26 wall. OK
  pbox[0] = AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[1] = (AliCRTConstants::Instance()->Depth() - 220 - 1600)/2;
  pbox[2] = 2987.7/2 - 1100/4 - concreteWidth/2;
  gMC->Gsvolu("CMO4", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CMO4", 1, "CMO1", 0, AliCRTConstants::Instance()->Depth()/2 - pbox[1], -1900 - pbox[2], 0, "MANY");
  // Small blocks below the volume CMO4 on both sides of the wall RB24/26. OK
  pbox[0] = (AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - ptubs[0])/2;
  pbox[1] = AliCRTConstants::Instance()->Depth()/2 - pbox[1];
  gMC->Gsvolu("CM17", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM17", 1, "CMO1", AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - pbox[0], -AliCRTConstants::Instance()->Depth()/2 + pbox[1], -1900 - pbox[2], 0, "MANY");
  gMC->Gspos("CM17", 2, "CMO1", -AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)+ pbox[0], -AliCRTConstants::Instance()->Depth()/2 + pbox[1], -1900 - pbox[2], 0, "MANY");

  // And a big block of molasse above the hall up to the surface. OK
  pbox[0] = pm25X - pm25radius;
  pbox[1] = (AliCRTConstants::Instance()->Depth()-500-1170)/2;
  pbox[2] = (1900 + 1150)/2;
  gMC->Gsvolu("CMO5", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CMO5", 1, "CMO1", 0,AliCRTConstants::Instance()->Depth()/2-pbox[1], pbox[2]-1900, 0, "MANY");
  // Small blocks of molasse betwen the blocks CMO2, CMO5 and PM25. Ok
  pbox[0] = (pm25X - pm25radius - 1170)/2;
  pbox[1] = 1000;
  gMC->Gsvolu("CM16", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM16", 1, "CMO1", 1170 + pbox[0], -AliCRTConstants::Instance()->Depth()/2+pbox[1], pbox[2] - 1900, 0, "MANY");

  // Molasse around the shafts.
  AliMatrix(idrotm[2003], 0, 0, 90, 0, 90, 90);
  // Around the PX24, the open section. OK
  ptubs[0] = px24radius + concreteWidth;
  ptubs[1] = ptubs[0] + 1000;
  ptubs[2] = (2300 - (5150 - AliCRTConstants::Instance()->Depth()))/2;
  ptubs[3] = 180 + kRaddeg*TMath::ASin(1070/ptubs[0]);
  ptubs[4] = 180 -  kRaddeg*TMath::ASin(1070/ptubs[0]);
  gMC->Gsvolu("CMO6", "TUBS", idtmed[1123], ptubs, 5);
  gMC->Gspos("CMO6", 1, "CMO1", px24X, ptubs[2] - AliCRTConstants::Instance()->Depth()/2, px24Z, idrotm[2003], "MANY");

  // Around the PX24, the closed section. OK
  Float_t ptube[3];
  ptube[0] = px24radius + concreteWidth;
  ptube[1] = ptube[0] + 1000;
  ptube[2] = (5150 - 2300)/2;
  gMC->Gsvolu("CMO7", "TUBE", idtmed[1123], ptube, 3);
  gMC->Gspos("CMO7", 1, "CMO1", px24X, AliCRTConstants::Instance()->Depth()/2 - ptube[2], px24Z, idrotm[2003], "MANY");

  // Around PM25. OK
  ptube[0] = pm25radius + concreteWidth;
  ptube[1] = ptube[0] + 400;
  ptube[2] = AliCRTConstants::Instance()->Depth()/2;
  gMC->Gsvolu("CMO8", "TUBE", idtmed[1123], ptube, 3);
  gMC->Gspos("CMO8", 1, "CMO1", pm25X, 0, pm25Z, idrotm[2003], "MANY");
  // On both sides of the PM25 along the HALL.
  pbox[0] = (2100 + pm25radius - 1170)/2;
  pbox[1] = AliCRTConstants::Instance()->Depth()/2;
  pbox[2] = (3*px24radius - pm25radius)/2;
  gMC->Gsvolu("CM18", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM18", 1, "CMO1", 2100, 0, pbox[2] + pm25radius, 0, "MANY");

  pbox[2] = (1900 - pm25radius)/2;
  gMC->Gsvolu("CM19", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM19", 1, "CMO1", 2100, 0, -pbox[2] - pm25radius, 0, "MANY");

  // Around the PGC2. OK
  ptube[0] = pgc2radius + concreteWidth;
  ptube[1] = 2987.7 - 740;
  ptube[2] = AliCRTConstants::Instance()->Depth()/2;
  gMC->Gsvolu("CMO9", "TUBE", idtmed[1123], ptube, 3);
  gMC->Gspos("CMO9", 1, "CMO1", pgc2X, 0, pgc2Z, idrotm[2003], "MANY");

  // On both sides of the PGC2.OK
  pbox[0] = (AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - 1100 - 375)/2;
  pbox[1] = AliCRTConstants::Instance()->Depth()/2;
  pbox[2] = pgc2radius + concreteWidth;
  gMC->Gsvolu("CM10", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM10", 1, "CMO1", AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - pbox[0], 0, pgc2Z, 0, "MANY");
  gMC->Gspos("CM10", 2, "CMO1", -AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) + pbox[0], 0, pgc2Z, 0, "MANY");

  // big block of molasse behind the PX24. OK
  pbox[0] = AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[1] = AliCRTConstants::Instance()->Depth()/2;
  pbox[2] = (pbox[0] - (2300 + 1150 + 100))/2;
  gMC->Gsvolu("CM12", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM12", 1, "CMO1", px24X, 0, px24Z + px24radius + concreteWidth + pbox[2], 0, "MANY");

  // big block of molasse in the opposite side of the PM25. OK
  pbox[0] = (AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - 1150)/2;
  pbox[1] = AliCRTConstants::Instance()->Depth()/2;
  pbox[2] = (1900 + 2300 + 1150)/2;
  gMC->Gsvolu("CM13", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM13", 1, "CMO1", -1150 - pbox[0], 0, pbox[2] - 1900, 0, "MANY");

  // big block of molasse behind the PM25. OK
  pbox[0] = (AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - (2100 + 910/2 + 100))/2;
  pbox[1] = AliCRTConstants::Instance()->Depth()/2;
  pbox[2] = (1900 + 2300 + 1150)/2;
  gMC->Gsvolu("CM14", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM14", 1, "CMO1", pm25X + pm25radius + concreteWidth + pbox[0], 0, pbox[2] - 1900, 0, "MANY");

  // big block of molasse behind the PGC2. OK
  pbox[0] = AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[1] = AliCRTConstants::Instance()->Depth()/2;
  pbox[2] = (pbox[0] - (2987.7 + 1900 + 1100/2 + 100))/2;
  gMC->Gsvolu("CM15", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM15", 1, "CMO1", 0, 0, -pbox[0] + pbox[2], 0, "MANY");

  gMC->Gspos("CMO1",1,"CRT",0,AliCRTConstants::Instance()->Depth()/2,0,0,"MANY");

}

//_____________________________________________________________________________
void AliCRTv1::CreateShafts()
{
  //
  //
  //
  Int_t  idrotm[2499];    // The rotation matrix.
  Int_t* idtmed = fIdtmed->GetArray() - 1099 ;

  //
  // Acces shafts
  //
  AliMatrix(idrotm[2001], 0, 0, 90, 0, 90, 90);


  // Create a bing cilinder to hold the main structures in the shaft.
  //   All the structures relative to the shaft will be put into
  // this volume.
  //   This shaft is composed by an open tube down in the hall, and
  // a cilinder avobe the level of the ceiling.
  Float_t ptube[3];
  ptube[0] = 0;    // inner radius
  ptube[1] = 1250; // outer radius
  ptube[2] = 5150/2; // Half lenght in Z
  gMC->Gsvolu("CSF1", "TUBE", idtmed[1114], ptube, 3);

  Float_t ptubs[5];
  // The open section of the PX24
  ptubs[0] = 1150; // Inner radius
  ptubs[1] = 1250; // Outer radius
  ptubs[2] = 1300; // Half length
  ptubs[3] = 180 + kRaddeg*TMath::ASin(1070/ptubs[0]); // starting angle
  ptubs[4] = 180 -  kRaddeg*TMath::ASin(1070/ptubs[0]);
  gMC->Gsvolu("CSF2", "TUBS", idtmed[1116], ptubs, 5);
  gMC->Gspos("CSF2", 1, "CSF1", 0, 0, -ptube[2] + ptubs[2], 0, "MANY");

  // The other part of the shaft.
  ptube[0] = ptubs[0]; // Inner radius
  ptube[1] = ptubs[1]; // Outer radius
  ptube[2] = 5150/2 - ptubs[2]; // Half lenght
  gMC->Gsvolu("CSF3", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF3", 1, "CSF1", 0, 0, 5150/2 - ptube[2], 0, "MANY");

  Float_t pbox[3];
  // Concrete walls along the shaft (next to the elevator.)
  pbox[0] = 480/2;  // Half length in X
  pbox[1] = 120/2;  // Half length in Y
  pbox[2] = 5150/2; // Half length in Z
  gMC->Gsvolu("CSW1", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW1", 1, "CSF1", 820+pbox[0],  150+pbox[1], 0, 0, "MANY");
  gMC->Gspos("CSW1", 2, "CSF1", 820+pbox[0], -300-pbox[1], 0, 0, "MANY");

  //
  pbox[0] = 120/2;  // Half length in X
  pbox[1] = 750/2;  // Half length in Y
  pbox[2] = 5150/2; // Half length in Z
  gMC->Gsvolu("CSW2", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW2", 1, "CSF1", 820-60, 150+pbox[1], 0, 0, "MANY");

  //
  pbox[0] = 120/2;  // Half length in X
  pbox[1] = 600/2;  // Half lenght in Y
  pbox[2] = 5150/2; // Half length in Z
  gMC->Gsvolu("CSW3", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW3", 1, "CSF1", 820-60, -300-pbox[1], 0, 0, "MANY");

  // Material below the counting rooms.
  pbox[0] = 400/2;
  pbox[1] = 2300/2;
  pbox[2] = 300/2;
  gMC->Gsvolu("CSW4", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW4",1,"CSF1",2300/2-pbox[0],0,3000-5150/2-pbox[2], 0, "MANY");

  // Shielding plug.
  pbox[0] = 1400/2;
  pbox[1] = 2300/2;
  pbox[2] = 170/2;
  gMC->Gsvolu("CSW5", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW5", 1, "CSF1", 0, 0, 3000-5150/2-130, 0, "MANY");

  // The end of the support for the shielding plug.
  pbox[0] = 170/2;
  pbox[1] = 2300/2;
  pbox[2] = 300/2;
  gMC->Gsvolu("CSW6", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW6",1,"CSF1",-1400/2-pbox[0],0,3000-5150/2-pbox[2],0,"MANY");

  // ...
  pbox[0] = 100/2;
  pbox[1] = 2300/2;
  pbox[2] = 450/2;
  gMC->Gsvolu("CSW7", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW7",1,"CSF1",-1400/2-170-pbox[0],0,3000-5150/2+pbox[2],0,"MANY");

  // Material close to the pipe.
  pbox[0] = 300/2;
  pbox[1] = 2300/2;
  pbox[2] = 170/2;
  gMC->Gsvolu("CSW8", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW8",1,"CSF1",-2300/2+pbox[0],0,2500-5150/2,0,"MANY");

  // Now put the shaft into the mother volume.
  gMC->Gspos("CSF1", 1, "CRT", 0, AliCRTConstants::Instance()->Depth() - 5150/2, 2300, idrotm[2001], "MANY");

  // PM25 Access Shaft
  ptube[0] = 910/2;
  ptube[1] = ptube[0] + 100;
  ptube[2] = (5150 - 1166)/2;
  gMC->Gsvolu("CSF4", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF4", 1, "CRT", 2100, AliCRTConstants::Instance()->Depth()-ptube[2], 0, idrotm[2001], "MANY");

  // PGC2 Access Shaft
  ptube[0] = 1100/2;
  ptube[1] = ptube[0] + 100;
  ptube[2] = (5150 - 690)/2;
  gMC->Gsvolu("CSF5", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF5", 1, "CRT", -375, AliCRTConstants::Instance()->Depth()-ptube[2], -1900 - 2987.7, idrotm[2001], "MANY");

}

//_____________________________________________________________________________
void AliCRTv1::DrawDetector() const
{
  //
  // Draw a shaded view of the L3 magnet
  //
  Info("DrawDetector", "Drawing CRT module");

  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("ALIC", "seen", 0);
  /*
  gMC->Gsatt("L3MO","seen",0); // L3 Magnet, Mother
  gMC->Gsatt("L3CO","seen",1); // Coils
  gMC->Gsatt("L3C1","seen",1); // Coils
  gMC->Gsatt("L3YO","seen",1); // Yoke
  gMC->Gsatt("L3DO","seen",0); // return Yoke (DOOR)
  gMC->Gsatt("L3FR","seen",1); // DOOR
  gMC->Gsatt("L3IR","seen",0); // Inner layer
  gMC->Gsatt("L3O1","seen",1); // Door opening
  gMC->Gsatt("L3O2","seen",1); // Door opening
  */
  gMC->Gsatt("CRT", "seen",0); // CRT mother volume.

  gMC->Gsatt("CMO1","seen",0); // Molasse.

  gMC->Gsatt("CSF1","seen",0); // PX24 access shaft.
  gMC->Gsatt("CSF2", "seen", 1); // PX24 open section
  gMC->Gsatt("CSF3", "seen", 1); // PX24, upper part.
  gMC->Gsatt("CSW1", "seen", 1);
  gMC->Gsatt("CSW2", "seen", 1);
  gMC->Gsatt("CSW3", "seen", 1);
  gMC->Gsatt("CSW4", "seen", 1);
  gMC->Gsatt("CSW5", "seen", 1);
  gMC->Gsatt("CSW6", "seen", 1);
  gMC->Gsatt("CSW7", "seen", 1);
  gMC->Gsatt("CSW8", "seen", 1);

  gMC->Gsatt("CSF4","seen",1); // PM25 access shaft.
  gMC->Gsatt("CSF5","seen",1); // PGC2 access shaft.

  gMC->Gsatt("CRT",  "seen", 0); // CRT Mother volume.
  gMC->Gsatt("CRT1", "seen", 0); // ?
  gMC->Gsatt("CRT2", "seen", 0); // Module air box
  gMC->Gsatt("CRT3", "seen", 1); // Scintillators
  gMC->Gsatt("CRT3", "colo", 2); // Scintillators
  gMC->Gsatt("CRT4", "seen", 1); // Aluminium frame (long bars)
  gMC->Gsatt("CRT4", "colo", 3); //
  gMC->Gsatt("CRT5", "seen", 1); // Aluminium frame (short bars)
  gMC->Gsatt("CRT5", "colo", 3); //
  gMC->Gsatt("CRT6", "seen", 1); // Module support
  gMC->Gsatt("CRT6", "colo", 3); //

  gMC->Gdopt("hide", "on");
  gMC->Gdopt("edge","off");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox("ALIC", 0, 3000, -3000, 3000, -6000, 6000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 70, 30, 0, 10, 9.5, .001, .001);
  gMC->Gdhead(1111, "View of CRT(ACORDE)");
  gMC->Gdman(18, 4, "MAN");

}

//_____________________________________________________________________________
void AliCRTv1::Init()
{
  //
  // Initialise L3 magnet after it has been built
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" CRTv1_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the CRTv1 initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }

}

//____________________________________________________________________________
void AliCRTv1::StepManager()
{
  //
  // Called for every step in the Cosmic Ray Trigger
  //
  static Int_t   vol[1];
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;

  static Float_t hits[14];
  static Float_t eloss;

  if ( gMC->TrackPid() != kMuonMinus ) return;

  // Only charged tracks
  if ( !(gMC->TrackCharge()) ) return;

  if (gMC->IsNewTrack()) {
    // Reset the deposited energy
    eloss = 0;
  }

  // Add th energy loss in each step.
  eloss += gMC->Edep();

  if ( ( (strcmp(gMC->CurrentVolName(),"CRT4") == 0) || // Magnet
	 (strcmp(gMC->CurrentVolName(),"CRT5") == 0) || // CRT
	 (strcmp(gMC->CurrentVolName(),"CRT6") == 0) || // Magnet Doors
	 (strcmp(gMC->CurrentVolName(),"CSF2") == 0) || // PX24
	 (strcmp(gMC->CurrentVolName(),"CSF3") == 0) || // PM25
	 (strcmp(gMC->CurrentVolName(),"CSF4") == 0) )  // PGC2
       && gMC->IsTrackEntering() ) {

  /*
  if ( (strcmp(gMC->CurrentVolName(),"CRT3") == 0)
       && gMC->IsTrackEntering() ) {
  */
    // Get current particle id(ipart),track position (pos) and momentum (mom)
    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);
    ipart = gMC->TrackPid();

    ipart = gMC->TrackPid();
    hits[0]  = (Float_t)ipart; //                 (fId)
    
    hits[1] = pos[0]; // X coordinate (fX)
    hits[2] = pos[1]; // Y coordinate (fY)
    hits[3] = pos[2]; // Z coordinate (fZ)
    hits[4] = mom[0]; // Px           (fpxug)
    hits[5] = mom[1]; // Py           (fpyug)
    hits[6] = mom[2]; // Pz           (fpzug)
    hits[7] = eloss;              // Energy loss

    Info("StepManager", "X=%f", pos[0]);

    // Tag the volumes
    if      ( (strcmp(gMC->CurrentVolName(),"CRT4")==0) ) vol[0] = 1; // Magnet
    else if ( (strcmp(gMC->CurrentVolName(),"CRT5")==0) ) vol[0] = 2; // CRT
    else if ( (strcmp(gMC->CurrentVolName(),"CRT6")==0) ) vol[0] = 3; // Doors
    else if ( (strcmp(gMC->CurrentVolName(),"CSF2")==0) ) vol[0] = 4; // PX24
    else if ( (strcmp(gMC->CurrentVolName(),"CSF3")==0) ) vol[0] = 5; // PM25
    else if ( (strcmp(gMC->CurrentVolName(),"CSF4")==0) ) vol[0] = 6; // PGC2
    else                                                  vol[0] = -1;// ?
    //vol[0]  = gMC->GetMedium();  //layer(flay)
    Info("StepManager", "Adding hit");
    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol, hits);
    Info("StepManager", "Hit added");
    // Reset the deposited energy only when you reach the Magnet
    if ( (strcmp(gMC->CurrentVolName(),"CRT4")==0) ) eloss = 0;

  } else {
    return;
  }

}

//_____________________________________________________________________________
void AliCRTv1::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a CRT hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliCRThit(fIshunt,track,vol,hits);
}
