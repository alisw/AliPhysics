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
Revision 1.10  2002/10/29 17:20:37  hristov
Corrections for subscript out of range (Alpha)

Revision 1.9  2002/10/23 06:47:56  alibrary
Introducing Riostream.h

Revision 1.8  2002/10/14 14:55:34  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.4.2.4  2002/10/10 14:40:31  hristov
Updating VirtualMC to v3-09-02

Revision 1.7  2002/10/07 11:13:25  gamez
Access shafts added

Revision 1.6  2002/07/26 06:21:12  gamez
CRT3 volume taken as sensitive volume

Revision 1.5  2002/07/25 12:52:34  morsch
AddHit call only if hit has been defined.

Revision 1.4  2002/07/12 12:57:29  gamez
Division of CRT1 corrected

Revision 1.3.2.1  2002/07/12 12:32:50  gamez
Division of CRT1 corrected

Revision 1.3  2002/07/10 15:57:04  gamez
CreateHall() removed, and new Molasse volumes

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
//  Trigger. This version will be used to simulation comic rays in alice     //
//  with all the detectors.                                                  //
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

#include <Riostream.h>

#include <TGeometry.h>
#include <TBRIK.h>
#include <TNode.h>
#include <TLorentzVector.h>

#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h"
#include "AliPDG.h"

#include "AliCRTv0.h"
#include "AliCRTConstants.h"

ClassImp(AliCRTv0)
 
//_____________________________________________________________________________
AliCRTv0::AliCRTv0() : AliCRT()
{
  //
  // Default constructor for CRT v0
  //
}
 
//_____________________________________________________________________________
AliCRTv0::AliCRTv0(const char *name, const char *title)
  : AliCRT(name,title)
{
  //
  // Standard constructor for CRT v0
  //
  //Begin_Html
  /*
    <img src="picts/AliCRTv0.gif">
  */
  //End_Html
}

//_____________________________________________________________________________
AliCRTv0::AliCRTv0(const AliCRTv0& crt)
{
  //
  // Copy ctor.
  //
  crt.Copy(*this);
}

//_____________________________________________________________________________
AliCRTv0& AliCRTv0::operator= (const AliCRTv0& crt)
{
  //
  // Asingment operator.
  //
  crt.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliCRTv0::BuildGeometry()
{
  //
  // Create the ROOT TNode geometry for the CRT
  //

  TNode *node, *top;

  const Int_t kColorCRT = kRed;

  // Find the top node alice.
  top = gAlice->GetGeometry()->GetNode("alice");

  new TBRIK("S_CRT_A", "CRT box", "void", 
	    AliCRTConstants::fgActiveAreaLenght/2., 
	    AliCRTConstants::fgActiveAreaHeight/2., 
	    AliCRTConstants::fgActiveAreaWidth/2.);

  
  new TRotMatrix("Left", "Left", 90., 315., 90., 45., 0., 337.5);
  new TRotMatrix("Right", "Right", 90., 45., 90., 315., 180., 202.5);
  new TRotMatrix("Up", "Up", 90., 0., 90., 90., 0., 90.);
  top->cd();

  //
  // Put 4 modules on the top of the magnet
  Float_t box = AliCRTConstants::fgCageWidth/2.;
  top->cd();
  node = new TNode("upper1", "upper1", "S_CRT_A", 0., 790.,  3.*box, "Up");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper2", "upper2", "S_CRT_A", 0., 790.,    box, "Up");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper3", "upper3", "S_CRT_A", 0., 790., -1.*box, "Up");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper4", "upper4", "S_CRT_A", 0., 790., -3.*box, "Up");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);


  // Modules on the left side.
  Float_t xtragap = 10.;
  Float_t initXside = (790.+xtragap)*TMath::Sin(2*22.5*kDegrad); //rigth side
  Float_t initYside = (790.+xtragap)*TMath::Cos(2*22.5*kDegrad);
  top->cd();
  node = new TNode("upper5", "upper5", "S_CRT_A", initXside, initYside,  3.*box, "Left");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper6", "upper6", "S_CRT_A", initXside, initYside,    box, "Left");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper7", "upper7", "S_CRT_A", initXside, initYside, -1.*box, "Left");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper8", "upper8", "S_CRT_A", initXside, initYside, -3.*box, "Left");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);


  // Modules on the right side.
  top->cd();
  node = new TNode("upper9", "upper9", "S_CRT_A", -initXside, initYside,  3.*box, "Right");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper10", "upper10", "S_CRT_A", -initXside, initYside,    box, "Right");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper11","upper11", "S_CRT_A", -initXside, initYside, -1.*box, "Right");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);

  top->cd();
  node = new TNode("upper12","upper12", "S_CRT_A", -initXside, initYside, -3.*box, "Right");
  node->SetLineColor(kColorCRT);
  fNodes->Add(node);


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
  CreateMolasse();

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
  gMC->Gsvolu("CRT1", "BOX ", idtmed[1112], scint, 3);     // Scintillators

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

  // Divide the modules in 2 planes.
  //gMC->Gsdvn("CRT2", "CRT1", 2, 2);
  // Now divide each plane in 8 palettes
  //gMC->Gsdvn("CRT3", "CRT2", 8, 3);

}

//_____________________________________________________________________________
void AliCRTv0::CreateMolasse()
{
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
  gMC->Gsvolu("CMO1", "TUBS", idtmed[1123], tspar, 5);
  gMC->Gspos("CMO1", 1, "ALIC", 0., 500., 1900.-tspar[2]+400., 0, "MANY");

  Float_t tbox[3];
  tbox[0] = 1250.;
  tbox[1] = (4420. - 1670.)/2.;
  tbox[2] = (1900.+1150.)/2. + 200.;
  gMC->Gsvolu("CM12", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM12", 1, "ALIC", 0., 4420. -tbox[1], 1900.-tbox[2]+400., 0, "MANY");

  AliMatrix(idrotm[2003], 0., 0., 90., 0., 90., 90.);
  // Along the PM25
  Float_t tube[3];
  tube[0] = 455. + 100.;
  tube[1] = 555. + 375.;
  tube[2] = (5150. - 1166.)/2.;
  gMC->Gsvolu("CMO2", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO2", 1, "ALIC", -2100., 4420.-tube[2], 0., idrotm[2003], "MANY");


  // Along the PGC2
  tube[0] = 650.;
  tube[1] = 2987.7;
  tube[2] = (5150. - 690.)/2.;
  gMC->Gsvolu("CMO3", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO3", 1, "ALIC", 375., 4420.-tube[2], 1900.+2987.7, idrotm[2003], "MANY");
  // Behind the PGC2 up to the end of the M. volume.
  tbox[0] = 12073.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (12073. - 1900.-2987.7-650.)/2.;
  gMC->Gsvolu("CMO7", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO7", 1, "ALIC", 0., 4420.-tbox[1], 1900.+2987.7+650.+tbox[2], 0, "MANY");

  // Along the PX24 , upper part.
  tube[0] = 1250.;
  tube[1] = 2300;
  tube[2] = 2575. - 1300. + 95.;
  gMC->Gsvolu("CMO4", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO4", 1, "ALIC", 0., 404.+1300.+tube[2], -2300., idrotm[2003], "MANY");

  // Along the PX24 , lower part
  tspar[0] = 1250.;
  tspar[1] = 2300;
  tspar[2] = 1300.;
  tspar[3] = kRaddeg*TMath::ASin(1070./1150.);
  tspar[4] = 360. - tspar[3];
  gMC->Gsvolu("CMO5", "TUBS", idtmed[1123], tspar, 5);
  gMC->Gspos("CMO5", 1, "ALIC", 0., 404., -2300., idrotm[2003], "MANY");
  // behind the PX24
  tbox[0] = 12073.;
  tbox[1] = 2575. + 95.;
  tbox[2] = 8523./2.;
  gMC->Gsvolu("CMO6", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO6", 1, "ALIC", 0., 4420.-tbox[1], -3550.-tbox[2], 0, "MANY");


  // On the right side of th hall
  tbox[0] = (12073. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (8437.7+650.)/2.;
  gMC->Gsvolu("CMO8", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO8", 1, "ALIC", 1250.+tbox[0], 4420.-tbox[1], -3550.+tbox[2], 0, "MANY");

  // on the left side of the hall, behind 
  tbox[0] = (12073. - 2755.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (8437.7+650.)/2.;
  gMC->Gsvolu("CMO9", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO9", 1, "ALIC", -2755.-tbox[0], 4420.-tbox[1], -3550.+tbox[2], 0, "MANY");


  // Molasse betwen the PX24 & PM25 on the left side.
  tbox[0] = (2755. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (3550. - 555.)/2.;
  gMC->Gsvolu("CM10", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM10", 1, "ALIC", -1250.-tbox[0], 4420.-tbox[1], -tbox[2]-555., 0, "MANY");


  // Molasse betwen the PGC2 & PM25 on the left side.
  tbox[0] = (2755. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (1900.+2987.7 - 555. + 650.)/2.;
  gMC->Gsvolu("CM11", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM11", 1, "ALIC", -1250.-tbox[0], 4420.-tbox[1], 555.+tbox[2], 0, "MANY");


}

//_____________________________________________________________________________
void AliCRTv0::CreateShafts()
{
  //
  //
  //
  Int_t  idrotm[2499];    // The rotation matrix.

  Int_t * idtmed = fIdtmed->GetArray() - 1099 ;

  // HAll ceiling
  Float_t ptubs[5];
  ptubs[0] = 1070.;
  ptubs[1] = 1170.;
  ptubs[2] = 1900.;
  ptubs[3] = 0.;
  ptubs[4] = 180.;
  gMC->Gsvolu("CHC1", "TUBS", idtmed[1116], ptubs, 5);
  gMC->Gspos("CHC1", 1, "ALIC", 0., 500., 0., 0, "ONLY");


  //
  // Acces shafts
  //
  AliMatrix(idrotm[2001], 0., 0., 90., 0., 90., 90.);
  
  // PX24
  ptubs[0] = 1150.;
  ptubs[1] = 1250.;
  ptubs[2] = 1300.;
  ptubs[3] = kRaddeg*TMath::ASin(1070./ptubs[0]);
  ptubs[4] = 360 - ptubs[3];
  gMC->Gsvolu("CSF1", "TUBS", idtmed[1116], ptubs, 5);
  gMC->Gspos("CSF1", 1, "ALIC", 0., 404., -2300., idrotm[2001], "MANY");

  Float_t ptube[3];
  ptube[0] = ptubs[0];
  ptube[1] = ptubs[1];
  ptube[2] = 2575. - ptubs[2] + 95.;
  gMC->Gsvolu("CSF2", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF2", 1, "ALIC", 0., 404.+ptubs[2]+ptube[2], -2300., idrotm[2001], "MANY");
  
  // Concrete walls along the shaft
  Float_t pbox[3];
  pbox[0] = 585./2.;
  pbox[1] = 2575. + 95.;
  pbox[2] = 20.;
  gMC->Gsvolu("CSW1", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW1", 1, "ALIC", -290-pbox[0], 404.-1300.+pbox[1], -3450.+210.*2, 0, "MANY");
  
  //
  pbox[0] = 750./2.;
  pbox[1] = 2575. + 95.;
  pbox[2] = 20.;
  gMC->Gsvolu("CSW3", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW3", 1, "ALIC", 420.-290.+pbox[0], 404.-1300.+pbox[1], -3450.+210.*2, 0, "MANY");
  
  //
  pbox[0] = 60.;
  pbox[1] = 2575. + 95.;
  pbox[2] = 210.;
  gMC->Gsvolu("CSW2", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW2", 1, "ALIC", -290-pbox[0], 404.-1300.+pbox[1], -3450.+pbox[2], 0, "MANY");
  gMC->Gspos("CSW2", 2, "ALIC", 420.-290.+pbox[0], 404.-1300.+pbox[1], -3450.+pbox[2], 0, "MANY");
  
  
  // 
  pbox[0] = 1000.;
  pbox[1] = 80.;
  pbox[2] = 200.;
  gMC->Gsvolu("CSP1", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP1", 1, "ALIC", 0., 2600.-700., -1150-pbox[2], 0, "MANY");
  
  //
  pbox[0] = 340.8;
  pbox[1] = 300./2.;
  pbox[2] = 460./2.;
  gMC->Gsvolu("CSP2", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP2", 1, "ALIC", 0., 2950.-700., -3450+pbox[2], 0, "MANY");
  
  //
  pbox[0] = 600.;
  pbox[1] = 150.;
  pbox[2] = 75.;
  gMC->Gsvolu("CSP3", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP3", 1, "ALIC", 0., 2950.-700., -1150.-210.-pbox[2], 0, "MANY");
  
  //
  pbox[0] = 600.;
  pbox[1] = 250.;
  pbox[2] = 38.;
  gMC->Gsvolu("CSP4", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP4", 1, "ALIC", 0., 2950.-700.+155.+pbox[1], -1150.-210.-pbox[2], 0, "MANY");
  
  
  // Shielding plug
  pbox[0] = 850.;
  pbox[1] = 90.;
  pbox[2] = 720.;
  gMC->Gsvolu("CSP5", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP5", 1, "ALIC", 0., 2950.-700., -3450.+460.+pbox[2], 0, "MANY");
  
  //
  pbox[0] = 80.;
  pbox[1] = 150.;
  pbox[2] = 720.;
  gMC->Gsvolu("CSP6", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP6", 1, "ALIC", 1150.-600., 2950.-700., -3450.+460.+pbox[2], 0, "MANY");
  gMC->Gspos("CSP6", 2, "ALIC", -1150.+600., 2950.-700., -3450.+460.+pbox[2], 0, "MANY");
  
  
  //
  pbox[0] = 130.;
  pbox[1] = 60.;
  pbox[2] = 750.;
  gMC->Gsvolu("CSP7", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP7", 1, "ALIC", 850.+pbox[0], 2950.-700.+100., -3450.+460.+pbox[2], 0, "MANY");
  gMC->Gspos("CSP7", 2, "ALIC", -850.-pbox[0], 2950.-700.+100., -3450.+460.+pbox[2], 0, "MANY");
  
  
  // PM25 Acces Shaft
  ptube[0] = 910./2.;
  ptube[1] = ptube[0] + 100.;
  ptube[2] = (5150. - 1166.)/2.;
  gMC->Gsvolu("CSF3", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF3", 1, "ALIC", -2100., AliCRTConstants::fgDepth-ptube[2], 0., idrotm[2001], "MANY");
  
  // PGC2 Access Shaft
  ptube[0] = 1100./2.;
  ptube[1] = ptube[0] + 100.;
  ptube[2] = (5150. - 690.)/2.;
  gMC->Gsvolu("CSF4", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF4", 1, "ALIC", 375., AliCRTConstants::fgDepth-ptube[2], 1900.+2987.7, idrotm[2001], "MANY");

}

//_____________________________________________________________________________

void AliCRTv0::CreateMaterials()
{
  // Use the standard materials.
  AliCRT::CreateMaterials();  
}


//_____________________________________________________________________________
void AliCRTv0::DrawDetector()
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
  // Called for every step in the Cosmic Ray Trigger
  //
  static Int_t   vol[5];
  Int_t		 copy;
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;

  static Float_t hits[14];
  Int_t tracknumber = gAlice->CurrentTrack();

  static Float_t eloss;
  static Float_t tlength;
  Float_t theta;
  Float_t phi;

  if ( !gMC->IsTrackAlive() ) return;

  if (gMC->IsNewTrack()) {
    // Reset the deposited energy
    eloss = 0.;
  }
  
  eloss += gMC->Edep(); // Store the energy loss along the trajectory.
  tlength += gMC->TrackStep();

  if (gMC->IsTrackEntering() && (strcmp(gMC->CurrentVolName(),"CM12") == 0) ) {

  // Get current particle id (ipart), track position (pos) and momentum (mom)
    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);
    ipart = gMC->TrackPid();

    Double_t tc   = mom[0]*mom[0]+mom[1]*mom[1];
    Double_t pt   = TMath::Sqrt(tc);
    theta   = Float_t(TMath::ATan2(pt,Double_t(mom[2])))*kRaddeg;
    phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;


    vol[0]    = gMC->CurrentVolOffID(1, vol[1]);
    vol[2]    = gMC->CurrentVolID(copy);
    vol[3]    = copy;
    
    hits[0]  = 0.f; //                 (fnmou)
    hits[1]  = (Float_t)ipart; //      (fId)

    hits[2]  = pos[0]; // X coordinate (fX)
    hits[3]  = pos[1]; // Y coordinate (fY)
    hits[4]  = pos[2]; // Z coordinate (fZ)
    hits[5]  = mom[0]; // Px           (fpxug)
    hits[6]  = mom[1]; // Py           (fpyug)
    hits[7]  = mom[2]; // Pz           (fpzug)
    
    hits[8]  = gMC->GetMedium();//layer(flay)
    hits[9]  = theta;   // arrival angle
    hits[10] = phi;     // 
    hits[11] = eloss;   // Energy loss
    hits[12] = tlength; // Trajectory lenght
    hits[13] = (Float_t)tracknumber;

    AddHit(gAlice->CurrentTrack(),vol, hits);

  }

}

