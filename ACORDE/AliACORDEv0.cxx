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
//  Trigger. This version will be used to simulation comic rays in alice with//
//  all the detectors. It include geometry and hits (posicion and momentum)  //
//  									     //	
//   Author: Enrique Gamez                                                   //
//                                                                           //
//                  Send comments to:                                        //
//      Arturo Fernandez <afernand@fcfm.buap.mx>                             //
//      Eleazar Cuautle  <ecuautle@nucleares.unam.mx>                        //
//									     //
//	Last update: Nov. 17th. 2009					     //
//	Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch		     //
//	FCFM-BUAP, Puebla, Pue. Mexico					     //
//									     //
///////////////////////////////////////////////////////////////////////////////


#include "AliACORDEv0.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TVirtualMC.h>
#include <TPDGCode.h>


#include "AliRun.h"
#include "AliConst.h"
#include "AliACORDEhit.h"
#include "AliACORDEConstants.h"
#include "AliMC.h"
#include "AliLog.h"

ClassImp(AliACORDEv0)
 
//_____________________________________________________________________________
AliACORDEv0::AliACORDEv0()
  : AliACORDE()
{
  //
  // Default constructor
  fIshunt = 0;
  fHits = 0;
  //
} 
//_____________________________________________________________________________
AliACORDEv0::AliACORDEv0(const char *name, const char *title)
  : AliACORDE(name, title)
{
  //
  // Standard constructor
  //
  fIshunt = 1; // All hits are associated with primary particles 
  fHits =  new TClonesArray("AliACORDEhit",400);
  gAlice->GetMCApp()->AddHitList(fHits);
}
//_____________________________________________________________________________
AliACORDEv0::~AliACORDEv0()
{
  //
  // Default destructor
  //
}

//_____________________________________________________________________________
void AliACORDEv0::CreateGeometry()
{
  CreateAcorde();
  if (GetCreateCavern()) CreateCavern();
}

void AliACORDEv0::CreateCavern()
{
  Int_t* idtmed = fIdtmed->GetArray() - 1099 ;
    // Create the mother volume, the one which will contain all the material
  // above the hall.
  Float_t pbox[3];
  pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  //pbox[0] = 12073;
  pbox[1] = AliACORDEConstants::Instance()->Depth();
  pbox[2] = pbox[0];
  gMC->Gsvolu("ACORDE", "BOX", idtmed[1114], pbox, 3);
  gMC->Gspos("ACORDE", 1, "ALIC", 0, 0, 0, 0, "ONLY");
  CreateShafts();
  CreateMolasse();
}

void AliACORDEv0::CreateShafts()

{

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
  gMC->Gspos("CSF1", 1, "ACORDE", 0, AliACORDEConstants::Instance()->Depth() - 5150/2, 2300, idrotm[2001], "MANY");

  // PM25 Access Shaft
  ptube[0] = 910/2;
  ptube[1] = ptube[0] + 100;
  ptube[2] = (5150 - 1166)/2;
  gMC->Gsvolu("CSF4", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF4", 1, "ACORDE", 2100, AliACORDEConstants::Instance()->Depth()-ptube[2], 0, idrotm[2001], "MANY");

  // PGC2 Access Shaft
  ptube[0] = 1100/2;
  ptube[1] = ptube[0] + 100;
  ptube[2] = (5150 - 690)/2;
  gMC->Gsvolu("CSF5", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF5", 1, "ACORDE", -375, AliACORDEConstants::Instance()->Depth()-ptube[2], -1900 - 2987.7, idrotm[2001], "MANY");

}


void AliACORDEv0::CreateMolasse()

{

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
  pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
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
  gMC->Gspos("CMO2", 1, "CMO1", 0, 500-AliACORDEConstants::Instance()->Depth()/2, ptubs[2]-1900, 0, "MANY");

  // Molasse around the RB24/26 Wall. OK
  ptubs[0] = 220 + 1600;
  ptubs[1] = AliACORDEConstants::Instance()->Depth() - ptubs[0];
  ptubs[2] = 2987.7/2 - 1100/4 - concreteWidth/2;
  ptubs[3] = 0;
  ptubs[4] = 180;
  gMC->Gsvolu("CMO3", "TUBS", idtmed[1123], ptubs, 5);
  gMC->Gspos("CMO3", 1, "CMO1", 70, 40-AliACORDEConstants::Instance()->Depth()/2, -1900 - ptubs[2], 0, "MANY");

  // A big block above the RB24/26 wall. OK
  pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[1] = (AliACORDEConstants::Instance()->Depth() - 220 - 1600)/2;
  pbox[2] = 2987.7/2 - 1100/4 - concreteWidth/2;
  gMC->Gsvolu("CMO4", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CMO4", 1, "CMO1", 0, AliACORDEConstants::Instance()->Depth()/2 - pbox[1], -1900 - pbox[2], 0, "MANY");
  // Small blocks below the volume CMO4 on both sides of the wall RB24/26. OK
  pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) -
ptubs[0])/2;
  pbox[1] = AliACORDEConstants::Instance()->Depth()/2 - pbox[1];
  gMC->Gsvolu("CM17", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM17", 1, "CMO1", AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - pbox[0], -AliACORDEConstants::Instance()->Depth()/2 + pbox[1], -1900 - pbox[2], 0, "MANY");
  gMC->Gspos("CM17", 2, "CMO1", -AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad)+ pbox[0], -AliACORDEConstants::Instance()->Depth()/2 + pbox[1], -1900 - pbox[2], 0, "MANY");

  // And a big block of molasse above the hall up to the surface. OK
  pbox[0] = pm25X - pm25radius;
  pbox[1] = (AliACORDEConstants::Instance()->Depth()-500-1170)/2;
  pbox[2] = (1900 + 1150)/2;
  gMC->Gsvolu("CMO5", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CMO5", 1, "CMO1", 0,AliACORDEConstants::Instance()->Depth()/2-pbox[1], pbox[2]-1900, 0, "MANY");
  // Small blocks of molasse betwen the blocks CMO2, CMO5 and PM25. Ok
  pbox[0] = (pm25X - pm25radius - 1170)/2;
  pbox[1] = 1000;
  gMC->Gsvolu("CM16", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM16", 1, "CMO1", 1170 + pbox[0], -AliACORDEConstants::Instance()->Depth()/2+pbox[1], pbox[2] - 1900, 0, "MANY");

  // Molasse around the shafts.
  AliMatrix(idrotm[2003], 0, 0, 90, 0, 90, 90);
  // Around the PX24, the open section. OK
  ptubs[0] = px24radius + concreteWidth;
  ptubs[1] = ptubs[0] + 1000;
  ptubs[2] = (2300 - (5150 - AliACORDEConstants::Instance()->Depth()))/2;
  ptubs[3] = 180 + kRaddeg*TMath::ASin(1070/ptubs[0]);
  ptubs[4] = 180 -  kRaddeg*TMath::ASin(1070/ptubs[0]);
  gMC->Gsvolu("CMO6", "TUBS", idtmed[1123], ptubs, 5);
  gMC->Gspos("CMO6", 1, "CMO1", px24X, ptubs[2] - AliACORDEConstants::Instance()->Depth()/2, px24Z, idrotm[2003], "MANY");
  // Around the PX24, the closed section. OK
  Float_t ptube[3];
  ptube[0] = px24radius + concreteWidth;
  ptube[1] = ptube[0] + 1000;
  ptube[2] = (5150 - 2300)/2;
  gMC->Gsvolu("CMO7", "TUBE", idtmed[1123], ptube, 3);
  gMC->Gspos("CMO7", 1, "CMO1", px24X, AliACORDEConstants::Instance()->Depth()/2 - ptube[2], px24Z, idrotm[2003], "MANY");

  // Around PM25. OK
  ptube[0] = pm25radius + concreteWidth;
  ptube[1] = ptube[0] + 400;
  ptube[2] = AliACORDEConstants::Instance()->Depth()/2;
  gMC->Gsvolu("CMO8", "TUBE", idtmed[1123], ptube, 3);
  gMC->Gspos("CMO8", 1, "CMO1", pm25X, 0, pm25Z, idrotm[2003], "MANY");
  // On both sides of the PM25 along the HALL.
  pbox[0] = (2100 + pm25radius - 1170)/2;
  pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
  pbox[2] = (3*px24radius - pm25radius)/2;
  gMC->Gsvolu("CM18", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM18", 1, "CMO1", 2100, 0, pbox[2] + pm25radius, 0, "MANY");

  pbox[2] = (1900 - pm25radius)/2;
  gMC->Gsvolu("CM19", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM19", 1, "CMO1", 2100, 0, -pbox[2] - pm25radius, 0, "MANY");

  // Around the PGC2. OK
  ptube[0] = pgc2radius + concreteWidth;
  ptube[1] = 2987.7 - 740;
  ptube[2] = AliACORDEConstants::Instance()->Depth()/2;
  gMC->Gsvolu("CMO9", "TUBE", idtmed[1123], ptube, 3);
  gMC->Gspos("CMO9", 1, "CMO1", pgc2X, 0, pgc2Z, idrotm[2003], "MANY");

  // On both sides of the PGC2.OK
  pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) -
1100 - 375)/2;
  pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
  pbox[2] = pgc2radius + concreteWidth;
  gMC->Gsvolu("CM10", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM10", 1, "CMO1", AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) - pbox[0], 0, pgc2Z, 0, "MANY");
  gMC->Gspos("CM10", 2, "CMO1", -AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) + pbox[0], 0, pgc2Z, 0, "MANY");

  // big block of molasse behind the PX24. OK
  pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
  pbox[2] = (pbox[0] - (2300 + 1150 + 100))/2;
  gMC->Gsvolu("CM12", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM12", 1, "CMO1", px24X, 0, px24Z + px24radius + concreteWidth + pbox[2], 0, "MANY");

  // big block of molasse in the opposite side of the PM25. OK
  pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) -
1150)/2;
  pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
  pbox[2] = (1900 + 2300 + 1150)/2;
  gMC->Gsvolu("CM13", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM13", 1, "CMO1", -1150 - pbox[0], 0, pbox[2] - 1900, 0, "MANY");

  // big block of molasse behind the PM25. OK
  pbox[0] = (AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad) -
(2100 + 910/2 + 100))/2;
  pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
  pbox[2] = (1900 + 2300 + 1150)/2;
  gMC->Gsvolu("CM14", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM14", 1, "CMO1", pm25X + pm25radius + concreteWidth + pbox[0], 0, pbox[2] - 1900, 0, "MANY");

  // big block of molasse behind the PGC2. OK
  pbox[0] = AliACORDEConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[1] = AliACORDEConstants::Instance()->Depth()/2;
  pbox[2] = (pbox[0] - (2987.7 + 1900 + 1100/2 + 100))/2;
  gMC->Gsvolu("CM15", "BOX", idtmed[1123], pbox, 3);
  gMC->Gspos("CM15", 1, "CMO1", 0, 0, -pbox[0] + pbox[2], 0, "MANY");

  gMC->Gspos("CMO1",1,"ACORDE",0,AliACORDEConstants::Instance()->Depth()/2,0,0,"MANY");

}

void AliACORDEv0::CreateAcorde()
{
  //
  // Create geometry for the ACORDE array
  // done in two main steps
  //  1.- definition of the modules
  //  2.- placement of the modules
  //
  Int_t  idrotm[2499];    // The rotation matrix.
  Int_t* idtmed = fIdtmed->GetArray() - 1099;
  AliACORDEConstants* constants = AliACORDEConstants::Instance();
  Float_t box[3];
  Float_t placed_at;
  Float_t placed_at2;
  Float_t small = 0.05; // to separate slightly some volumes
                        // by half a mm so that they do not overlap


  // 1.- Definition of a module
  // *  ACORDE1 => volume filled with air, representing a module
  //               it contains all other volumes defining the module
  //               there are 60 copies of it
  // *  ACORDE2 => volume defining one scintillator pad
  //               there are 2 copies of it per module
  // *  ACORDE3-6 => volumes representing the Al walls of box
  //               surrounding the plastic
  //               3: long wall, 2 copies (front, back)
  //               4: end caps, 2 copies (left, right)
  //               5: long stripe to model the profile 
  //                  4 copies (upper front and back, lower)
  //               6: short stripe to model the profile
  //                  4 copies (upper left, right; lower)

  // The full module volume.
  // This volume will be ocupied by all the material of the module
  // the scintillators, the aluminium frame, etc.
  box[0] = constants->ModuleLength()/2;
  box[1] = constants->ModuleHeight()/2;
  box[2] = constants->ModuleWidth()/2;
  gMC->Gsvolu("ACORDE1", "BOX", idtmed[1114], box, 3);

  // The scintillators
  box[0] = constants->PlasticLength()/2;
  box[1] = constants->PlasticHeight()/2;
  box[2] = constants->PlasticWidth()/2;
  gMC->Gsvolu("ACORDE2", "BOX", idtmed[1112], box, 3);

  // it is important to keep this order for easy assignment of 
  // a volume to a physical module:
  placed_at = box[1]+constants->ProfileThickness()
    - constants->ModuleHeight()/2+small;
  gMC->Gspos("ACORDE2", 1, "ACORDE1", 0, placed_at, 0, 0, "MANY");
  placed_at = placed_at + 2.0*box[1]+small;
  gMC->Gspos("ACORDE2", 2, "ACORDE1", 0, placed_at, 0, 0, "MANY");


  // The metallic frame: long walls of box
  // back,front,left,right, defined looking
  // from the + z diraction into alice; i.e.
  // back ==> z<0, front ==> z>0
  // left ==> x<0, right ==> x>0
  // up ==> increasing y, down ==> decreasing y
  box[0] = constants->ModuleLength()/2;
  box[1] = constants->ModuleHeight()/2;
  box[2] = constants->ProfileThickness()/2.0; 
  gMC->Gsvolu("ACORDE3", "BOX", idtmed[1108], box, 3);
  // front wall
  placed_at = constants->ModuleWidth()/2-constants->ProfileThickness()/2.0;
  gMC->Gspos("ACORDE3", 1, "ACORDE1", 0, 0, placed_at, 0, "MANY");
  // back wall
  gMC->Gspos("ACORDE3", 2, "ACORDE1", 0, 0, -placed_at , 0, "MANY");

  // The metallic frame: end caps
  box[0] = constants->ProfileThickness()/2.0;
  box[1] = constants->ModuleHeight()/2;
  box[2] = constants->ModuleWidth()/2;
  gMC->Gsvolu("ACORDE4", "BOX", idtmed[1108], box, 3);
  // right cap
  placed_at = constants->ModuleLength()/2-constants->ProfileThickness()/2.0;
  gMC->Gspos("ACORDE4", 1, "ACORDE1", placed_at, 0, 0, 0, "MANY");
  // left cap
  gMC->Gspos("ACORDE4", 2, "ACORDE1", -placed_at, 0, 0, 0, "MANY");

  // The metallic frame: the profile, long stripes
  box[0] = constants->ModuleLength()/2.0;
  box[1] = constants->ProfileThickness()/2;
  box[2] = constants->ProfileWidth()/2;
  gMC->Gsvolu("ACORDE5", "BOX", idtmed[1108], box, 3);
  // upper front
  placed_at = constants->ModuleHeight()/2-box[1];
  placed_at2 = constants->ModuleWidth()/2-
    constants->ProfileThickness()-box[2];
  gMC->Gspos("ACORDE5", 1, "ACORDE1",0,placed_at,placed_at2, 0, "MANY");
  // upper back
  gMC->Gspos("ACORDE5", 2, "ACORDE1",0,placed_at,-placed_at2, 0, "MANY");
  // lower front
  gMC->Gspos("ACORDE5", 3, "ACORDE1",0,-placed_at,placed_at2, 0, "MANY");
  // lower back
  gMC->Gspos("ACORDE5", 4, "ACORDE1",0,-placed_at,-placed_at2, 0, "MANY");

  // The metallic frame: the profile, long stripes
  box[0] = constants->ProfileWidth()/2.0;
  box[1] = constants->ProfileThickness()/2;
  box[2] = constants->ModuleWidth()/2-constants->ProfileWidth();
  gMC->Gsvolu("ACORDE6", "BOX", idtmed[1108], box, 3);
  // upper right
  placed_at = constants->ModuleHeight()/2-box[1];
  placed_at2 = constants->ModuleLength()/2-
    constants->ProfileThickness()-box[0];
  gMC->Gspos("ACORDE6", 1, "ACORDE1",placed_at2,placed_at,0, 0, "MANY");
  // upper left
  gMC->Gspos("ACORDE6", 2, "ACORDE1",-placed_at2,placed_at,0, 0, "MANY");
  // lower right
  gMC->Gspos("ACORDE6", 3, "ACORDE1",placed_at2,-placed_at,0, 0, "MANY");
  // lower left
  gMC->Gspos("ACORDE6", 4, "ACORDE1",-placed_at2,-placed_at,0, 0, "MANY");

  // End of MODULE definition

  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  // 2.- placement of the module
  // Now put all of them in the right position in 
  // master volume ALIC

  // rotation matrices (see Geant manual for conventions)
  // for columns 4 and 5
  AliMatrix(idrotm[231], 90, 45, 90, 135, 0, 0);
  // for columns 0 and 1
  AliMatrix(idrotm[232], 90, 315, 90, 45, 0, 0);

  // place each one of the 6 columns in turn
  // for the first and the last column the position
  // of the two last modules depends on the value 
  // of the fITSGeometry variable

  // it is important to keep this order because
  // the copy number defines the module!

  // first column, except first and last  modules
  for (Int_t copy = 2; copy < 10; copy++)
    gMC->Gspos("ACORDE1",copy,"ALIC",
	       constants->OldModulePositionX(copy-1),
	       constants->OldModulePositionY(copy-1),
	       constants->OldModulePositionZ(copy-1),
	       idrotm[232], "MANY");
  // second column
  for (Int_t copy = 11; copy < 21; copy++)
    gMC->Gspos("ACORDE1",copy,"ALIC",
	       constants->OldModulePositionX(copy-1),
	       constants->OldModulePositionY(copy-1),
	       constants->OldModulePositionZ(copy-1),
	       idrotm[232], "MANY");
  // third and fourth columns
  for (Int_t copy = 21; copy < 41; copy++)
    gMC->Gspos("ACORDE1",copy,"ALIC",
	       constants->OldModulePositionX(copy-1),
	       constants->OldModulePositionY(copy-1),
	       constants->OldModulePositionZ(copy-1),
	       0, "MANY");
  // fifth column
  for (Int_t copy = 41; copy < 51; copy++)
    gMC->Gspos("ACORDE1",copy,"ALIC",
	       constants->OldModulePositionX(copy-1),
	       constants->OldModulePositionY(copy-1),
	       constants->OldModulePositionZ(copy-1),
	       idrotm[231], "MANY");
  // last column, except first and last  modules
  for (Int_t copy = 52; copy < 60; copy++)
    gMC->Gspos("ACORDE1",copy,"ALIC",
	       constants->OldModulePositionX(copy-1),
	       constants->OldModulePositionY(copy-1),
	       constants->OldModulePositionZ(copy-1),
	       idrotm[231], "MANY");
  // the last four modules
  if (Get4CentralModulesGeometry()) {
    gMC->Gspos("ACORDE1",1,"ALIC",
	       constants->OldExtraModulePositionX(),
	       constants->OldExtraModulePositionY(),
	       constants->OldExtraModulePositionZ(0),
	       0, "MANY");  
    gMC->Gspos("ACORDE1",10,"ALIC",
	       constants->OldExtraModulePositionX(),
	       constants->OldExtraModulePositionY(),
	       constants->OldExtraModulePositionZ(1),
	       0, "MANY");  
    gMC->Gspos("ACORDE1",51,"ALIC",
	       constants->OldExtraModulePositionX(),
	       constants->OldExtraModulePositionY(),
	       constants->OldExtraModulePositionZ(2),
	       0, "MANY");  
    gMC->Gspos("ACORDE1",60,"ALIC",
	       constants->OldExtraModulePositionX(),
	       constants->OldExtraModulePositionY(),
	       constants->OldExtraModulePositionZ(3),
	       0, "MANY");  
  } else {
    gMC->Gspos("ACORDE1",1,"ALIC",
	       constants->OldModulePositionX(0),
	       constants->OldModulePositionY(0),
	       constants->OldModulePositionZ(0),
	       idrotm[232], "MANY");
    gMC->Gspos("ACORDE1",10,"ALIC",
	       constants->OldModulePositionX(9),
	       constants->OldModulePositionY(9),
	       constants->OldModulePositionZ(9),
	       idrotm[232], "MANY");
    gMC->Gspos("ACORDE1",51,"ALIC",
	       constants->OldModulePositionX(50),
	       constants->OldModulePositionY(50),
	       constants->OldModulePositionZ(50),
	       idrotm[231], "MANY");
    gMC->Gspos("ACORDE1",60,"ALIC",
	       constants->OldModulePositionX(59),
	       constants->OldModulePositionY(59),
	       constants->OldModulePositionZ(59),
	       idrotm[231], "MANY");
  } // end if (fITSGeometry)

}

//____________________________________________________________________________

void AliACORDEv0::Init()
{
  // Initialise L3 magnet after it has been built
  Int_t i;
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" ACORDEv0_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    // Here the ACORDEv initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
 // AliACORDE::Init();  
}
//____________________________________________________________________________
void AliACORDEv0::StepManager()
{
  //
  // Called for every step in the Cosmic Ray Trigger
  //


  // volume: 
  //  [0] = module number 1-60 (1==>(0-0), 60 (5-9)
  //  [1] = Plastic number: 0 (down) to 1 (up)
  static Int_t   vol[2]; 
  //
  // hit
  // [0] = PID
  // [1-3] = x, y, z 
  // [4] = time 
  // [5-7] = px, py, pz
  // [8] = energy 
  // [9] = energy loss
  // [10] = length of track through plastic
  static Float_t hits[11];

  // local static variables
  static Float_t eloss;
  static Float_t step;
  // scintillator volume
  static Int_t idScint = gMC->VolId("ACORDE2");

  // local variables
  Int_t copy;
  TLorentzVector pos;
  TLorentzVector mom;

  // only charged tracks
  if ( !gMC->TrackCharge() || !gMC->IsTrackAlive() ) return;

  // only in sensitive material
  if (gMC->CurrentVolID(copy) == idScint) {
    step  += gMC->TrackStep();
    eloss += gMC->Edep();
    // set all hit variables except eloss which is resetted
    // set volume variables
    if (gMC->IsTrackEntering()) {
      eloss = 0.0;
      step = 0.0;
      gMC->TrackPosition(pos);
      gMC->TrackMomentum(mom);
      // hit
      // [0] = PID
      // [1-3] = x, y, z 
      // [4] = time 
      // [5-7] = px, py, pz
      // [8] = energy 
      // [9] = energy loss
      hits[0]  = (Float_t ) gMC->TrackPid(); 
      hits[1] = pos[0]; 
      hits[2] = pos[1]; 
      hits[3] = pos[2]; 
      hits[4] = gMC->TrackTime();
      hits[5] = mom[0]; 
      hits[6] = mom[1]; 
      hits[7] = mom[2]; 
      hits[8] = gMC->Etot();
      // volume: 
      //  [0] = module number 1-60 (1==>(0-0), 60 (5-9)
      //  [1] = Plastic number: 0 (down) to 1 (up)
      Int_t copyPlastic; // plastic: down=1, up=2
      Int_t copyModule; // module: 1-60
      gMC->CurrentVolID(copyPlastic);
      gMC->CurrentVolOffID(1, copyModule);
      // module
      vol[0] = copyModule;
      // plastic: 0 = down, 1 = up
      vol[1] = copyPlastic;
    } // end if gMC->IsTrackEntering()

    // set hit[9] = total energy loss and book hit
    if( gMC->IsTrackExiting() || 
	gMC->IsTrackStop() || 
	gMC->IsTrackDisappeared()){
      hits[9] = eloss;
      hits[10] = step;
      eloss = 0.0;
      step = 0.0;
      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol, hits);
    }
  } // end if in scintillator

}

//_____________________________________________________________________________
void AliACORDEv0::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a ACORDE hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliACORDEhit(fIshunt,track,vol,hits);
}

