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

#include <TGeometry.h>
#include <TBRIK.h>
#include <TNode.h>
#include <TVirtualMC.h>

#include "AliRun.h"
#include "AliConst.h"

#include "AliCRTv0.h"
#include "AliCRTConstants.h"

ClassImp(AliCRTv0)
 
//_____________________________________________________________________________
AliCRTv0::AliCRTv0()
  : AliCRT()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliCRTv0::AliCRTv0(const char *name, const char *title)
  : AliCRT(name, title)
{
  //
  // Standard constructor
  //
  //Begin_Html
  /*
    <img src="picts/AliCRTv0.gif">
  */
  //End_Html
  SetMarkerColor(kRed);
  SetMarkerStyle(kRed);
  SetMarkerSize(0.4);
}

//_____________________________________________________________________________
AliCRTv0::AliCRTv0(const AliCRTv0& crt)
  : AliCRT(crt)
{
  //
  // Copy constructor
  //
  crt.Copy(*this);
}

//_____________________________________________________________________________
AliCRTv0::~AliCRTv0()
{
  //
  // Default destructor
  //
}

//_____________________________________________________________________________
AliCRTv0& AliCRTv0::operator=(const AliCRTv0& crt)
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
  Int_t* idtmed = fIdtmed->GetArray() - 1099 ;

  // Create a mother volume.
  Float_t pbox[3];
  //pbox[0] = AliCRTConstants::fgDepth*TMath::Tan(67.5*kDegrad);
  pbox[0] = 12073.;
  pbox[1] = AliCRTConstants::fgDepth;
  pbox[2] = pbox[0];
  gMC->Gsvolu("CRT", "BOX", idtmed[1114], pbox, 3);
  gMC->Gspos("CRT", 1, "ALIC", 0, 0, 0, 0, "ONLY");

  // Create a big volume with air barrel above the magnet
  Float_t barrel[10];
  Float_t magnetSides = 3.;
  Float_t planesPerpendicularToZ = 2.;
  barrel[0] = 22.5;
  barrel[1] = 45*magnetSides;
  barrel[2] = magnetSides;
  barrel[3] = planesPerpendicularToZ;
  barrel[4] = -700.;
  barrel[5] = AliCRTConstants::fgMagMinRadius;
  barrel[6] = AliCRTConstants::fgMagMinRadius + 2.; // 2 cm width
  barrel[7] = -barrel[4];
  barrel[8] = barrel[5];
  barrel[9] = barrel[6];
  gMC->Gsvolu("CRT4", "PGON", idtmed[1112], barrel, 10);
  gMC->Gspos("CRT4", 1 , "CRT", 0., -30., 0., 0, "MANY");

  //
  Float_t box[3];
  box[0] = AliCRTConstants::fgSinglePaletteLenght/4;
  box[1] = AliCRTConstants::fgSinglePaletteHeight/2;
  box[2] = AliCRTConstants::fgSinglePaletteWidth/2;
  gMC->Gsvolu("CRT6", "BOX", idtmed[1112], box, 3);

  // In the right side side of the magnet
  AliMatrix(idrotm[231], 90., 45., 90., 315., 180., 202.5);

  // In the left side side of the magnet
  AliMatrix(idrotm[232], 90., 315., 90., 315., 0.0000040, 263.0707092);

  // Now put them into the volume created above
  // First above the magnet.
  const Float_t away = (2.*barrel[5]*TMath::Sin(kDegrad*22.5))/4.;
  const Int_t nModules = 10;
  for (Int_t i = 0; i < nModules; i++) {
    Float_t zCoordinate = i*100 - 450;
    // In the lef side
    gMC->Gspos("CRT6", i, "CRT4", -away, barrel[5]+1., zCoordinate, 0, "MANY");
    // In the rigth side
    gMC->Gspos("CRT6",i+10,"CRT4", away, barrel[5]+1., zCoordinate, 0, "MANY");

    // The most away part (left side)
    gMC->Gspos("CRT6", i+20, "CRT4", 3*away, barrel[5]+31 - away, zCoordinate, idrotm[232], "MANY");
    // The inner part (left side)
    gMC->Gspos("CRT6", i+30, "CRT4", 4*away, barrel[5]+31 - 2*away, zCoordinate, idrotm[232], "MANY");

    // The most away part (rigth side)
    gMC->Gspos("CRT6", i+40, "CRT4", -3*away, barrel[5]+31 - away, zCoordinate, idrotm[231], "MANY");
    // The inner part (rigth side)
    gMC->Gspos("CRT6", i+50, "CRT4", -4*away, barrel[5]+31 - 2*away, zCoordinate, idrotm[231], "MANY");
  }

  // Now the magnet doors
  magnetSides = 8.;
  barrel[1] = 45*magnetSides;
  barrel[2] = magnetSides;
  barrel[4] = 700.;
  barrel[5] = 0;
  barrel[6] = 790;
  barrel[7] = barrel[4] + 2.;
  barrel[8] = barrel[5];
  barrel[9] = barrel[6];
  gMC->Gsvolu("CRT5", "PGON", idtmed[1111], barrel, 10);
  gMC->Gspos("CRT5", 1, "CRT", 0., -30., 0., 0, "ONLY");

  AliMatrix(idrotm[300], 90., 0., 90., 90., 180., 0.);
  gMC->Gspos("CRT5", 2, "CRT", 0., -30., 0., idrotm[300], "ONLY");

}

//_____________________________________________________________________________
void AliCRTv0::DrawDetector()
{
  //
  // Draw a shaded view of the L3 magnet
  //

  Info("DrawDetector", "Drawing the module");

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
