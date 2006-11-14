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

#include "AliCRTv0.h"

#include <TGeometry.h>
#include <TBRIK.h>
#include <TNode.h>
#include <TVirtualMC.h>

#include "AliRun.h"
#include "AliConst.h"

#include "AliCRTConstants.h"
#include "AliCRTModule.h"

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
  //PH  SetMarkerColor(kRed);
  //PH  SetMarkerStyle(kRed);
  //PH  SetMarkerSize(0.4);
}

//_____________________________________________________________________________
AliCRTv0::~AliCRTv0()
{
  //
  // Default destructor
  //
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

  AliCRTConstants* crtConstants = AliCRTConstants::Instance();

  new TBRIK("S_CRT_A", "CRT box", "void", 
	    crtConstants->ActiveAreaLenght()/2., 
	    crtConstants->ActiveAreaHeight()/2., 
	    crtConstants->ActiveAreaWidth()/2.);

  
  new TRotMatrix("Left", "Left", 90., 315., 90., 45., 0., 337.5);
  new TRotMatrix("Right", "Right", 90., 45., 90., 315., 180., 202.5);
  new TRotMatrix("Up", "Up", 90., 0., 90., 90., 0., 90.);
  top->cd();

  //
  // Put 4 modules on the top of the magnet
  Float_t box = crtConstants->CageWidth()/2.;
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
  Int_t* idtmed = fIdtmed->GetArray() - 1099;
  AliCRTConstants* crtConstants = AliCRTConstants::Instance();

  // Create the mother volume.
  // This volume can be seen as the volume which ACORDE will ocupate
  // above the upper face of the L3 magnet. Inside this volume the detectors
  // aboce the magnet will be, then there will be two copies of this volume,
  // one for each side.
  Float_t box[3];
  //box[0] = 2*crtConstants->MagMinRadius()*TMath::Sin(kDegrad*22.5);
  box[0] = crtConstants->MagMinRadius()*TMath::Sin(kDegrad*22.5);
  box[1] = crtConstants->MagMaxRadius() - crtConstants->MagMinRadius();
  box[2] = crtConstants->MagnetLenght()/2;
  gMC->Gsvolu("CRT1", "BOX", idtmed[1112], box, 3);

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

  Float_t x = crtConstants->MagMinRadius()+10;
  gMC->Gspos("CRT1", 1, "ALIC", 0, x, 0, 0, "MANY");
  gMC->Gspos("CRT1", 2, "ALIC", -x*TMath::Sin(kDegrad*45), x*TMath::Cos(kDegrad*45), 0, idrotm[231], "MANY");
  gMC->Gspos("CRT1", 3, "ALIC",  x*TMath::Sin(kDegrad*45), x*TMath::Cos(kDegrad*45), 0, idrotm[232], "MANY");

}

//_____________________________________________________________________________
void AliCRTv0::DrawDetector() const
{
  //
  // Draw a shaded view of the L3 magnet
  //

  Info("DrawDetector", "Drawing the module");

  gMC->Gsatt("*", "seen", -1);

  gMC->Gsatt("ALIC","seen",0);

  gMC->Gsatt("L3MO","seen",0); // L3 Magnet, Mother
  gMC->Gsatt("L3CO","seen",1); // Coils
  gMC->Gsatt("L3C1","seen",1); // Coils
  gMC->Gsatt("L3YO","seen",1); // Yoke
  gMC->Gsatt("L3DO","seen",0); // return Yoke (DOOR)
  gMC->Gsatt("L3FR","seen",1); // DOOR
  gMC->Gsatt("L3IR","seen",0); // Inner layer
  gMC->Gsatt("L3O1","seen",1); // Door opening
  gMC->Gsatt("L3O2","seen",1); // Door opening

  gMC->Gsatt("CRT1", "seen", 0); // CRT Mother
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
  //gMC->Gdraw("alic", 40, 30, 0, 10, 9.5, .009, .009);
  gMC->Gdraw("alic", 30, 40, 0, -30, -60, .09, .09);
  gMC->Gdhead(1111, "View of CRT(ACORDE)");
  gMC->Gdman(18, 4, "MAN");
}
