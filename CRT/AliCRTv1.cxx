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
#include "AliCRTConstants.h"
#include "AliMC.h"

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
  fIshunt =  1; // All hits are associated with primary particles  

  fHits =  new TClonesArray("AliCRThit",400);
  gAlice->GetMCApp()->AddHitList(fHits);

  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}

//_____________________________________________________________________________
AliCRTv1::AliCRTv1(const AliCRTv1& crt)
  : AliCRT(crt)
{
  //
  // Copy ctor.
  //
  crt.Copy(*this);
}

//_____________________________________________________________________________
AliCRTv1::~AliCRTv1()
{
  //
  // Default destructor
  //
}

//_____________________________________________________________________________
AliCRTv1& AliCRTv1::operator=(const AliCRTv1& crt)
{
  //
  // Asingment operator
  //
  crt.Copy(*this);
  return *this;
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

  // Shafts.
  this->CreateShafts();

  // Molasse.
  this->CreateMolasse();

  AliCRTConstants* crtConstants = AliCRTConstants::Instance();

  // Create a big volume with air barrel above the magnet
  Float_t barrel[10];
  Float_t magnetSides = 3.;
  Float_t planesPerpendicularToZ = 2.;
  barrel[0] = 22.5;
  barrel[1] = 45*magnetSides;
  barrel[2] = magnetSides;
  barrel[3] = planesPerpendicularToZ;
  barrel[4] = -700.;
  barrel[5] = crtConstants->MagMinRadius();
  barrel[6] = crtConstants->MagMinRadius() + 2.; // 2 cm width
  barrel[7] = -barrel[4];
  barrel[8] = barrel[5];
  barrel[9] = barrel[6];
  gMC->Gsvolu("CRT4", "PGON", idtmed[1112], barrel, 10);
  gMC->Gspos("CRT4", 1 , "CRT", 0., -30., 0., 0, "ONLY");

  //
  Float_t box[3];
  box[0] = crtConstants->SinglePaletteLenght()/4;
  box[1] = crtConstants->SinglePaletteHeight()/2;
  box[2] = crtConstants->SinglePaletteWidth()/2;
  gMC->Gsvolu("CRT6", "BOX", idtmed[1113], box, 3);

  // In the right side side of the magnet
  AliMatrix(idrotm[231], 90., 45., 90., 315., 180., 202.5);

  // In the left side side of the magnet
  //AliMatrix(idrotm[232], 90., 315., 90., 315., 0.0000040, 263.0707092);
  AliMatrix(idrotm[232], 90, 315, 90, 315, 0, 263);

  // Now put them into the volume created above
  // First above the magnet.
  Float_t away = (2.*barrel[5]*TMath::Sin(kDegrad*22.5))/4.;
  Int_t nModules = 10;
  for (Int_t i = 0; i < nModules; i++) {
    Float_t zCoordinate = i*100 - 450;
    // In the lef side
    gMC->Gspos("CRT6", i, "CRT4", -away, barrel[5]+1., zCoordinate, 0, "ONLY");
    // In the rigth side
    gMC->Gspos("CRT6",i+10,"CRT4", away, barrel[5]+1., zCoordinate, 0, "ONLY");

    // The most away part (left side)
    gMC->Gspos("CRT6", i+20, "CRT4", 3*away, barrel[5]+21 - away, zCoordinate, idrotm[232], "ONLY");
    // The inner part (left side)
    gMC->Gspos("CRT6", i+30, "CRT4", 4*away, barrel[5]+21 - 2*away, zCoordinate, idrotm[232], "ONLY");

    // The most away part (rigth side)
    gMC->Gspos("CRT6", i+40, "CRT4", -3*away, barrel[5]+21. - away, zCoordinate, idrotm[231], "ONLY");
    // The inner part (rigth side)
    gMC->Gspos("CRT6", i+50, "CRT4", -4*away, barrel[5]+21 - 2*away, zCoordinate, idrotm[231], "ONLY");
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
void AliCRTv1::CreateMolasse()
{
  //
  //
  //
  Int_t  idrotm[2499];    // The rotation matrix.
  Int_t* idtmed = fIdtmed->GetArray() - 1099 ;

  // Exactly above the hall
  Float_t tspar[5];
  tspar[0] = 1170;
  tspar[1] = 1170 + 375;
  tspar[2] = (1900 + 1150)/2 + 100;
  tspar[3] = 0;
  tspar[4] = 180;
  gMC->Gsvolu("CMO1", "TUBS", idtmed[1123], tspar, 5);
  gMC->Gspos("CMO1", 1, "CRT", 0, 500., 1900 - tspar[2] + 400, 0, "MANY");

  Float_t tbox[3];
  tbox[0] = 1250;
  tbox[1] = (4420 - 1670)/2;
  tbox[2] = (1900 + 1150)/2 + 200;
  gMC->Gsvolu("CM12", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM12",1,"CRT",0, 4420 - tbox[1], 1900 - tbox[2] + 400, 0,"MANY");

  AliMatrix(idrotm[2003], 0., 0., 90., 0., 90., 90.);
  // Along the PM25
  Float_t tube[3];
  tube[0] = 455 + 100;
  tube[1] = 555 + 375;
  tube[2] = (5150 - 1166)/2;
  gMC->Gsvolu("CMO2", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO2", 1, "CRT", -2100, 4420 - tube[2], 0, idrotm[2003], "MANY");

  // Along the PGC2
  tube[0] = 650;
  tube[1] = 2987.7;
  tube[2] = (5150 - 690)/2;
  gMC->Gsvolu("CMO3", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO3",1,"CRT", 375, 4420 - tube[2], 1900 + 2987.7,idrotm[2003],"MANY");

  // Behind the PGC2 up to the end of the M. volume.
  tbox[0] = 12073;
  tbox[1] = 2575 + 95;
  tbox[2] = (12073 - 1900 - 2987.7 - 650)/2.;
  gMC->Gsvolu("CMO7", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO7", 1, "CRT", 0, 4420 - tbox[1], 1900 + 2987.7 + 650 + tbox[2], 0, "MANY");

  // Along the PX24 , upper part.
  tube[0] = 1250;
  tube[1] = 2300;
  tube[2] = 2575 - 1300 + 95;
  gMC->Gsvolu("CMO4", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO4", 1, "CRT", 0, 404 + 1300 + tube[2], -2300, idrotm[2003], "MANY");

  // Along the PX24 , lower part
  tspar[0] = 1250;
  tspar[1] = 2300;
  tspar[2] = 1300;
  tspar[3] = kRaddeg*TMath::ASin(1070./1150.);
  tspar[4] = 360 - tspar[3];
  gMC->Gsvolu("CMO5", "TUBS", idtmed[1123], tspar, 5);
  gMC->Gspos("CMO5", 1, "CRT", 0., 404, -2300, idrotm[2003], "MANY");
  // behind the PX24
  tbox[0] = 12073;
  tbox[1] = 2575 + 95;
  tbox[2] = 8523/2;
  gMC->Gsvolu("CMO6", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO6", 1, "CRT", 0., 4420 - tbox[1], -3550 - tbox[2], 0, "MANY");

  // On the right side of th hall
  tbox[0] = (12073 - 1250)/2;
  tbox[1] = 2575 + 95;
  tbox[2] = (8437.7+650)/2;
  gMC->Gsvolu("CMO8", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO8", 1, "CRT", 1250 + tbox[0], 4420 - tbox[1], -3550 + tbox[2], 0, "MANY");

  // on the left side of the hall, behind 
  tbox[0] = (12073 - 2755)/2;
  tbox[1] = 2575 + 95;
  tbox[2] = (8437.7 + 650)/2.;
  gMC->Gsvolu("CMO9", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO9", 1, "CRT", -2755 - tbox[0], 4420 - tbox[1], -3550 + tbox[2], 0, "MANY");

  // Molasse betwen the PX24 & PM25 on the left side.
  tbox[0] = (2755 - 1250)/2;
  tbox[1] = 2575 + 95;
  tbox[2] = (3550 - 555)/2;
  gMC->Gsvolu("CM10", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM10", 1, "CRT", -1250 - tbox[0], 4420 - tbox[1], -tbox[2] - 555, 0, "MANY");

  // Molasse betwen the PGC2 & PM25 on the left side.
  tbox[0] = (2755 - 1250)/2;
  tbox[1] = 2575 + 95;
  tbox[2] = (1900 + 2987.7 - 555 + 650)/2;
  gMC->Gsvolu("CM11", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM11", 1, "CRT", -1250 - tbox[0], 4420 - tbox[1], 555 + tbox[2], 0, "MANY");

}

//_____________________________________________________________________________
void AliCRTv1::CreateShafts()
{
  //
  //
  //
  Int_t  idrotm[2499];    // The rotation matrix.
  Int_t* idtmed = fIdtmed->GetArray() - 1099 ;

  // Create a mother volume.
  Float_t pbox[3];
  //pbox[0] = AliCRTConstants::Instance()->Depth()*TMath::Tan(67.5*kDegrad);
  pbox[0] = 12073.;
  pbox[1] = AliCRTConstants::Instance()->Depth();
  pbox[2] = pbox[0];
  gMC->Gsvolu("CRT", "BOX", idtmed[1114], pbox, 3);
  gMC->Gspos("CRT", 1, "ALIC", 0, 0, 0, 0, "ONLY");

  // HAll ceiling
  Float_t ptubs[5];
  ptubs[0] = 1070;
  ptubs[1] = 1170;
  ptubs[2] = 1900;
  ptubs[3] = 0;
  ptubs[4] = 180;
  gMC->Gsvolu("CHC1", "TUBS", idtmed[1116], ptubs, 5);
  gMC->Gspos("CHC1", 1, "CRT", 0, 500, 0, 0, "ONLY");

  //
  // Acces shafts
  //
  AliMatrix(idrotm[2001], 0., 0., 90., 0., 90., 90.);

  // PX24
  ptubs[0] = 1150;
  ptubs[1] = 1250;
  ptubs[2] = 1300;
  ptubs[3] = kRaddeg*TMath::ASin(1070/ptubs[0]);
  ptubs[4] = 360 - ptubs[3];
  gMC->Gsvolu("CSF1", "TUBS", idtmed[1116], ptubs, 5);
  gMC->Gspos("CSF1", 1, "CRT", 0., 404, -2300, idrotm[2001], "MANY");

  Float_t ptube[3];
  ptube[0] = ptubs[0];
  ptube[1] = ptubs[1];
  ptube[2] = 2575 - ptubs[2] + 95;
  gMC->Gsvolu("CSF2", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF2", 1, "CRT", 0, 404 + ptubs[2] + ptube[2], -2300, idrotm[2001], "MANY");

  // Concrete walls along the shaft
  pbox[0] = 585/2;
  pbox[1] = 2575 + 95;
  pbox[2] = 20;
  gMC->Gsvolu("CSW1", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW1", 1, "CRT", -290 - pbox[0], 404 - 1300 + pbox[1], -3450 + 210*2, 0, "MANY");

  //
  pbox[0] = 750/2;
  pbox[1] = 2575 + 95;
  pbox[2] = 20;
  gMC->Gsvolu("CSW3", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW3", 1, "CRT", 420 - 290 +pbox[0], 404 - 1300 + pbox[1], -3450 + 210*2, 0, "MANY");

  //
  pbox[0] = 60;
  pbox[1] = 2575 + 95;
  pbox[2] = 210;
  gMC->Gsvolu("CSW2", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW2", 1, "CRT", -290 - pbox[0], 404 - 1300 + pbox[1], -3450 + pbox[2], 0, "MANY");
  gMC->Gspos("CSW2", 2, "CRT", 420 - 290 + pbox[0], 404 - 1300 + pbox[1], -3450 + pbox[2], 0, "MANY");

  // 
  pbox[0] = 1000;
  pbox[1] = 80;
  pbox[2] = 200;
  gMC->Gsvolu("CSP1", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP1", 1, "CRT", 0, 2600 - 700, -1150 - pbox[2], 0, "MANY");

  //
  pbox[0] = 340.8;
  pbox[1] = 300/2.;
  pbox[2] = 460/2.;
  gMC->Gsvolu("CSP2", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP2", 1, "CRT", 0, 2950.-700., -3450+pbox[2], 0, "MANY");

  //
  pbox[0] = 600;
  pbox[1] = 150;
  pbox[2] = 75;
  gMC->Gsvolu("CSP3", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP3", 1, "CRT", 0, 2950.-700., -1150.-210.-pbox[2], 0, "MANY");

  //
  pbox[0] = 600;
  pbox[1] = 250;
  pbox[2] = 38;
  gMC->Gsvolu("CSP4", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP4", 1, "CRT", 0, 2950 - 700 + 155+pbox[1], -1150 - 210 - pbox[2], 0, "MANY");

  // Shielding plug
  pbox[0] = 850;
  pbox[1] = 90;
  pbox[2] = 720;
  gMC->Gsvolu("CSP5", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP5", 1, "CRT", 0, 2950 - 700, -3450 + 460 + pbox[2], 0,"MANY");

  //
  pbox[0] = 80;
  pbox[1] = 150;
  pbox[2] = 720;
  gMC->Gsvolu("CSP6", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP6", 1, "CRT", 1150 - 600 , 2950 - 700, -3450 + 460 + pbox[2], 0, "MANY");
  gMC->Gspos("CSP6", 2, "CRT", -1150 + 600, 2950 - 700, -3450 + 460 + pbox[2], 0, "MANY");

  //
  pbox[0] = 130;
  pbox[1] = 60;
  pbox[2] = 750;
  gMC->Gsvolu("CSP7", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP7", 1, "CRT", 850 + pbox[0], 2950 - 700 + 100, -3450 + 460 + pbox[2], 0, "MANY");
  gMC->Gspos("CSP7", 2, "CRT", -850 - pbox[0], 2950 - 700+ 100, -3450 + 460 + pbox[2], 0, "MANY");

  // PM25 Acces Shaft
  ptube[0] = 910/2;
  ptube[1] = ptube[0] + 100;
  ptube[2] = (5150 - 1166)/2;
  gMC->Gsvolu("CSF3", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF3", 1, "CRT", -2100, AliCRTConstants::Instance()->Depth()-ptube[2], 0, idrotm[2001], "MANY");

  // PGC2 Access Shaft
  ptube[0] = 1100/2;
  ptube[1] = ptube[0] + 100;
  ptube[2] = (5150 - 690)/2;
  gMC->Gsvolu("CSF4", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF4", 1, "CRT", 375, AliCRTConstants::Instance()->Depth()-ptube[2], 1900 + 2987.7, idrotm[2001], "MANY");

}

//_____________________________________________________________________________
void AliCRTv1::DrawDetector() const
{
  //
  // Draw a shaded view of the L3 magnet
  //
  //cout << "AliCRTv1::DrawModule() : Drawing the module" << endl;
  
  
  Int_t able = 1;
  Int_t enable = 0;
  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("alic", "seen", 0);
  
  gMC->Gsatt("ALIC","seen",enable);
  gMC->Gsatt("CRT", "seen",enable);
  gMC->Gsatt("L3MO","seen",  able); // L3 Magnet
  //gMC->Gsatt("CRT1","seen",  able); // Scintillators
  gMC->Gsatt("CRT4","seen",  able); // Scintillators barrel
  
  // Draw the molasse volumes
  gMC->Gsatt("CMO1","seen",enable); // Exactly above the HALL
  gMC->Gsatt("CMO2","seen",enable); // Molasse, along the PM25
  gMC->Gsatt("CMO3","seen",enable); // molasse along the PGC2
  gMC->Gsatt("CMO4","seen",enable); // Molasse, behind the PX24 upper part
  gMC->Gsatt("CMO5","seen",enable); // molasse behind px24, lower part
  gMC->Gsatt("CMO6","seen",enable); // behind the PX24
  gMC->Gsatt("CMO7","seen",enable); // behind the PGC2
  gMC->Gsatt("CMO8","seen",enable); // on the right side.
  gMC->Gsatt("CMO9","seen",enable); // on the left side.
  gMC->Gsatt("CM10","seen",enable); // betwen PX24 & PM25.
  gMC->Gsatt("CM11","seen",enable); // betwen PGC2 & PM25.
  gMC->Gsatt("CM12","seen",enable); // box above the hall.
  
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
void AliCRTv1::Init()
{
  //
  // Initialise L3 magnet after it has been built
  Int_t i;
  //
  if(fDebug) {
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

    // Tag the volumes
    if      ( (strcmp(gMC->CurrentVolName(),"CRT4")==0) ) vol[0] = 1; // Magnet
    else if ( (strcmp(gMC->CurrentVolName(),"CRT5")==0) ) vol[0] = 2; // CRT
    else if ( (strcmp(gMC->CurrentVolName(),"CRT6")==0) ) vol[0] = 3; // Doors
    else if ( (strcmp(gMC->CurrentVolName(),"CSF2")==0) ) vol[0] = 4; // PX24
    else if ( (strcmp(gMC->CurrentVolName(),"CSF3")==0) ) vol[0] = 5; // PM25
    else if ( (strcmp(gMC->CurrentVolName(),"CSF4")==0) ) vol[0] = 6; // PGC2
    else                                                  vol[0] = -1;// ?
    //vol[0]  = gMC->GetMedium();  //layer(flay)

    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol, hits);

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

//_____________________________________________________________________________
void AliCRTv1::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
  AliDetector::ResetHits();
}

//_____________________________________________________________________________
void AliCRTv1::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  AliDetector::ResetDigits();
} 

//____________________________________________________________________________
void AliCRTv1::FinishEvent()
{
  //
  //
}
