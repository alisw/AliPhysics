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


#include <TBRIK.h>
#include <TGeometry.h>
#include <TLorentzVector.h>
#include <TNode.h>
#include <TPDGCode.h>

#include "AliCRTConstants.h"
#include "AliCRTv1.h"
#include "AliConst.h"
#include "AliMagF.h"
#include "AliRun.h"

ClassImp(AliCRTv1)
 
//_____________________________________________________________________________
AliCRTv1::AliCRTv1() : AliCRTv0()
{
  //
  // Default constructor for CRT
  //
  fCRTStatus = kTRUE;
  fRICHStatus = kFALSE;
  fTPCStatus = kFALSE;
  fMagnetStatus = kTRUE;

  fCRTModule = kFALSE;
}
 
//_____________________________________________________________________________
AliCRTv1::AliCRTv1(const char *name, const char *title)
  : AliCRTv0(name,title)
{
  //
  // Standard constructor for CRT
  //
  //Begin_Html
  /*
    <img src="picts/AliCRTv1.gif">
  */
  //End_Html
  fCRTStatus = kTRUE;
  fCRTModule = kFALSE;

  fRICHStatus = kFALSE;
  fTPCStatus = kFALSE;
  fMagnetStatus = kFALSE;
}

//_____________________________________________________________________________
AliCRTv1::AliCRTv1(const AliCRTv1& crt)
{
  //
  // Copy ctor.
  //
  crt.Copy(*this);
}

//_____________________________________________________________________________
AliCRTv1& AliCRTv1::operator= (const AliCRTv1& crt)
{
  //
  // Asingment operator
  //
  crt.Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliCRTv1::CreateGeometry()
{
  //
  // Create geometry for the CRT array
  //

  Int_t  idrotm[2499];    // The rotation matrix.

  Int_t * idtmed = fIdtmed->GetArray() - 1099 ;

  //
  // Shafts.
  this->CreateShafts();

  //
  // Molasse.
  this->CreateMolasse();


  //
  // Scintillators

  Float_t box[3];
  box[0] = AliCRTConstants::fgCageLenght/2.; // Half Length of the box along the X axis, cm.
  box[1] = AliCRTConstants::fgCageHeight/2.; // Half Length of the box along the Y axis, cm.
  box[2] = AliCRTConstants::fgCageWidth/2.;  // Half Length of the box along the Z axis, cm.

  //
  // Create a big voluem with air barrel above the magnet
  Float_t barrel[10];
  Float_t magnetSides = 3.;
  Float_t planesPerpendicularToZ = 2.;
  Float_t rMin = 790.;
  Float_t rMax = rMin + 20.; // 20 cm width
  barrel[0] = 22.5;
  barrel[1] = 45*magnetSides;
  barrel[2] = magnetSides;
  barrel[3] = planesPerpendicularToZ;
  barrel[4] = -600.;
  barrel[5] = rMin;
  barrel[6] = rMax;
  barrel[7] = 600.;
  barrel[8] = rMin;
  barrel[9] = rMax;
  gMC->Gsvolu("CRT4", "PGON", idtmed[1114], barrel, 10);
  gMC->Gspos("CRT4", 1 , "CRT", 0., -30., 0., 0, "ONLY");
  

  // Create  the current sicuiitllator arry
  // Define the Scintillators. as a big box.
  Float_t scint[3];
  scint[0] = AliCRTConstants::fgActiveAreaLenght/2.;   // Half Length in X
  scint[1] = AliCRTConstants::fgActiveAreaHeight/2.;   // Half Length in Y
  scint[2] = AliCRTConstants::fgActiveAreaWidth/2.;    // Half Length in Z
  gMC->Gsvolu("CRT1", "BOX ", idtmed[1112], scint, 3); // Scintillators
  //
  // -- X axis.
  // we'll start dawing from the center.
  Float_t initX = 0.;
  
  //
  // -- Y axis
  Float_t gapY   = 30.;        // 30 cms. above the barrel.
  // For the height we staimate the from the center of the ceiling,
  // if were a cilinder, must be about 280cm.
  Float_t barrelc = 790.; // Barrel radius.
  Float_t height  = barrelc + gapY - 30.;
  Float_t initY = height;
  
  //
  // -- Z axis.
  // we'll start dawing from the center.
  
  //
  // Put 4 modules on the top of the magnet
  Int_t step = 4;
  for ( Int_t i = 1 ; i <= 4 ; i++ ) {
    gMC->Gspos("CRT1", i, "CRT", initX, initY, (i-step)*box[2], 0, "ONLY");
    step--;
  }
  
  // Modules on the barrel sides.
  // Because the openenig angle for each face is 22.5, and if we want to
  //    put the modules right in the middle
  Float_t xtragap = 10.;
  Float_t initXside = (height+xtragap)*TMath::Sin(2*22.5*kDegrad);//rigthside
  Float_t initYside = (height+xtragap)*TMath::Cos(2*22.5*kDegrad);
  
  // Put 4 modules on the left side of the magnet
  // The rotation matrix parameters, for the left side.
  AliMatrix(idrotm[232], 90., 315., 90., 45., 0., 337.5);
  Int_t stepl = 4;
  for ( Int_t i = 1 ; i <= 4 ; i++ ) {
    gMC->Gspos("CRT1", i+4, "CRT", initXside, initYside, (i-stepl)*box[2],
	       idrotm[232], "ONLY");
    stepl--;
  }
  
  // Put 4 modules on the right side of the magnet
  // The rotation matrix parameters for the right side.
  AliMatrix(idrotm[231], 90., 45., 90., 315., 180., 202.5);
  Int_t stepr = 4;
  for ( Int_t i = 1 ; i <= 4 ; i++ ) {
    gMC->Gspos("CRT1", i+8, "CRT", -initXside, initYside, (i-stepr)*box[2],
	       idrotm[231], "ONLY");
    stepr--;
  }
  
  this->CreateMagnetGeometry();
  this->CreateRICHGeometry();
  this->CreateTPCGeometry();
  
}

//_____________________________________________________________________________
void AliCRTv1::CreateMagnetGeometry()
{
  
  cout<<"\n\n\tYou are requiring the CRT with the Magnet Activated!\n\n";
  
  Int_t  idrotm[2499];    // The rotation matrix.

  Int_t * idtmed = fIdtmed->GetArray() - 1099 ;

  // Disable the CRT StepManager method.
  fCRTStatus = kFALSE;

  Float_t barrel[10];
  Float_t magnetSides = 3.;
  Float_t planesPerpendicularToZ = 2.;
  //Float_t rMin = 790.;
  //Float_t rMax = rMin + 20.; // 20 cm width

  // MAgnet
  // Create the upper faces of the magnet.
  barrel[0] = 22.5;
  barrel[1] = 360.;
  barrel[2] = 8.;
  barrel[3] = 2.;
  barrel[4] = -600.;
  barrel[5] = 580.;
  barrel[6] = 790.;
  barrel[7] = 600.;
  barrel[8] = 580.;
  barrel[9] = 790.;
  gMC->Gsvolu("C3MO", "PGON", idtmed[1114], barrel, 10);
  gMC->Gspos("C3MO", 1, "CRT", 0., -30., 0., 0, "ONLY");

  // Define coils 
  
  barrel[5] = 585.;
  barrel[6] = 690.;
  barrel[8] = 585.;
  barrel[9] = 690.;
  gMC->Gsvolu("C3CO", "PGON", idtmed[1108], barrel, 10); //Aluminium
  gMC->Gspos("C3CO", 1, "C3MO", 0., 0., 0., 0, "ONLY");
  
  barrel[5] = 580.;
  barrel[6] = 585.;
  barrel[8] = 580.;
  barrel[9] = 585.;
  gMC->Gsvolu("C3C1", "PGON", idtmed[1128], barrel, 10);// Aluminium
  gMC->Gspos("C3C1", 1, "C3MO", 0., 0., 0., 0, "ONLY");

  // Define yoke 
  
  barrel[5] = 690.;
  barrel[6] = 790.;
  barrel[8] = 690.;
  barrel[9] = 790.;
  gMC->Gsvolu("C3YO", "PGON", idtmed[1109], barrel, 10); // Iron
  gMC->Gspos("C3YO", 1, "C3MO", 0., 0., 0., 0, "ONLY");


  // Now create one inside the magnet as L3C1
  // voulme for tracking.
  barrel[0] = 22.5;
  barrel[1] = 45*magnetSides;
  barrel[2] = magnetSides;
  barrel[3] = planesPerpendicularToZ;
  barrel[4] = -600.;
  barrel[5] = 575.;
  barrel[6] = 580.;
  barrel[7] = 600.;
  barrel[8] = 575.;
  barrel[9] = 580.;
  gMC->Gsvolu("C3CI", "PGON", idtmed[1134], barrel, 10);
  gMC->Gspos("C3CI", 1 , "CRT", 0., -30., 0., 0, "ONLY");

  // And a detector layer in the door 10 cm thick
  // Volume for tracking.
  barrel[0] = 22.5;
  barrel[1] = 360.;
  barrel[2] = 8.;
  barrel[3] = 2.;
  barrel[4] = 590.;
  barrel[5] = 0.;
  barrel[6] = 580.;
  barrel[7] = 600.;
  barrel[8] = barrel[5];
  barrel[9] = barrel[6];
  gMC->Gsvolu("C3C2", "PGON", idtmed[1154], barrel, 10); // Air
  gMC->Gspos("C3C2", 1, "CRT",  0., -30., 0., 0, "ONLY");
  AliMatrix(idrotm[1010], 90., 0., 90., 90., 180., 0.);
  gMC->Gspos("C3C2", 2, "CRT",  0., -30., 0., idrotm[1010], "ONLY");



  barrel[4] = 600.;
  barrel[5] = 0.;
  barrel[6] = 790.;
  barrel[7] = 700.;
  barrel[8] = barrel[5];
  barrel[9] = barrel[6];
  gMC->Gsvolu("C3DO", "PGON", idtmed[1174], barrel, 10); // Air
  gMC->Gspos("C3DO", 1, "CRT", 0., -30., 0., 0, "ONLY");
  AliMatrix(idrotm[1010], 90., 0., 90., 90., 180., 0.);
  gMC->Gspos("C3DO", 2, "CRT", 0., -30., 0., idrotm[1010], "ONLY");

  barrel[4] = 610.;
  barrel[5] = 0.;
  barrel[6] = 790.;
  barrel[7] = 700.;
  barrel[8] = barrel[5];
  barrel[9] = barrel[6];
  gMC->Gsvolu("C3FR", "PGON", idtmed[1149], barrel, 10); // Iron
  gMC->Gspos("C3FR", 1, "C3DO", 0., 0., 0., 0, "ONLY");
  // INNER LAYER 
  
  barrel[4] = 600.;
  barrel[7] = 610.;
  gMC->Gsvolu("C3IR", "PGON", idtmed[1149], barrel, 10); //Iron
  gMC->Gspos("C3IR", 1, "C3DO", 0., 0., 0., 0, "ONLY");

}

//_____________________________________________________________________________
void AliCRTv1::CreateTPCGeometry()
{
  cout<<"\n\n\tYou are requiring the CRT with the TPC Activated!\n\n";
  Int_t * idtmed = fIdtmed->GetArray() - 1099 ;

  // Disable the CRT StepManager method.
  fCRTStatus = kFALSE;
  // Disable the MAgnet
  fMagnetStatus = kFALSE;
  // Disable th RICH
  fRICHStatus = kFALSE;

  // TPC
  // Tpc SAndwich 1 - Al
  // TSA1
  Float_t tube[5];
  tube[0]=274.8124;
  tube[1]=278.;
  tube[2]=252.1;
  tube[3] = 0.;
  tube[4] = 180.;
  gMC->Gsvolu("CSA1","TUBS",idtmed[1154],tube,5);
  // TSA1->TOCV (0.,0.,3.) ->TOIN (0.,0.,0.)->TPC (0.,0.,0.)->ALIC(0.,0.,0.)
  gMC->Gspos("CSA1 ",1,"CRT",0.,0.,0.,0,"ONLY");

}

//_____________________________________________________________________________
void AliCRTv1::CreateRICHGeometry()
{

  cout<<"\n\n\tYou are requiring the CRT with the RICH Activated!\n\n";

  Int_t  idrotm[2499];    // The rotation matrix.

  Int_t * idtmed = fIdtmed->GetArray() - 1099 ;

  // Disable the CRT StepManager method.
  fCRTStatus = kFALSE;
  // Disable the MAgnet
  fMagnetStatus = kFALSE;


  // now create  volume to simulate the HMPID volume. CSI
  Float_t csi_length = 160*.8 + 2.6;
  Float_t csi_width = 144*.84 + 2*2.6;
  Float_t tbox[3];
  tbox[0] = csi_width/2;
  tbox[1] = 11.5;
  tbox[2] = csi_length/2;
  gMC->Gsvolu("CRIC ", "BOX ", idtmed[1174], tbox, 3);

  Double_t dOffset = 490+1.267 - 8/2;  // distance from center of mother volume ALIC to methane
  
  Double_t dAlpha = 19.5; // angle between centers of chambers - y-z plane
  Double_t dAlphaRad = dAlpha*kDegrad;
  
  Double_t dBeta = 20.;   // angle between center of chambers - y-x plane
  Double_t dBetaRad = dBeta*kDegrad;
   
  Double_t dRotAngle = 60.;     // the whole RICH is to be rotated in x-y plane + means clockwise rotation 
  Double_t dRotAngleRad = dRotAngle*kDegrad;
   
    
   TRotMatrix *pRotMatrix; // tmp pointer
   
   TVector3 vector(0,dOffset,0); // Position of chamber 2 without rotation

   // Chamber 0  standalone (no other chambers in this row) 
   AliMatrix(idrotm[1000],90, -dRotAngle+360,90-dAlpha, 90-dRotAngle, dAlpha, -90+300);
   pRotMatrix=new TRotMatrix("rot993","rot993",90,-dRotAngle, 90-dAlpha,90-dRotAngle,dAlpha, -90);
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateX(dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("CRIC",1,"CRT",vector.X(),vector.Y(),vector.Z(),idrotm[1000], "ONLY");
   
   // Chamber 1   
   AliMatrix(idrotm[1001],90,-dBeta-dRotAngle,90,90-dBeta-dRotAngle, 0,0);

   pRotMatrix=new TRotMatrix("rot994","rot994",90,-dBeta-dRotAngle,90,90-dBeta-dRotAngle,0,0);  
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(-dBetaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("CRIC",2,"CRT",vector.X(),vector.Y(),vector.Z(),idrotm[1001], "ONLY");           
   
   // Chamber 2   the top one with no Alpha-Beta rotation
   AliMatrix(idrotm[1002],90,-dRotAngle,90,90-dRotAngle,0,0);

   pRotMatrix=new TRotMatrix("rot995","rot995",90,-dRotAngle,90,90-dRotAngle,0,0);
   
   vector.SetXYZ(0,dOffset,0);
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("CRIC",3,"CRT",vector.X(),vector.Y(),vector.Z(),idrotm[1002], "ONLY");           
   
   // Chamber 3
   AliMatrix(idrotm[1003],90,dBeta-dRotAngle,90.,90+dBeta-dRotAngle,0,0);
   pRotMatrix=new TRotMatrix("rot996","rot996", 90,dBeta-dRotAngle,90.,90+dBeta-dRotAngle,0,0);
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(dBetaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("CRIC",4,"CRT",vector.X(),vector.Y(),vector.Z(),idrotm[1003], "ONLY");

   // Chamber 4
   AliMatrix(idrotm[1004],90,360-dBeta-dRotAngle,108.2,90-dBeta-dRotAngle,18.2,90-dBeta-60);
   pRotMatrix=new TRotMatrix("rot997","rot997",90,360-dBeta-dRotAngle,108.2,90-dBeta-dRotAngle,18.2,90-dBeta);
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(-dBetaRad); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("CRIC",5,"CRT",vector.X(),vector.Y(),vector.Z(),idrotm[1004], "ONLY");

   // Chamber 5
   AliMatrix(idrotm[1005],90,-dRotAngle+360,90+dAlpha,90-dRotAngle,dAlpha,90-60);     

   pRotMatrix=new TRotMatrix("rot998","rot998",90,-dRotAngle,90+dAlpha,90-dRotAngle,dAlpha,90);     
   
   vector.SetXYZ(0,dOffset,0); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("CRIC",6,"CRT",vector.X(),vector.Y(),vector.Z(),idrotm[1005], "ONLY");           
   
   // Chamber 6           
   AliMatrix(idrotm[1006],90,dBeta-dRotAngle+360,108.2,90+dBeta-dRotAngle,18.2,90+dBeta-60);

   pRotMatrix=new TRotMatrix("rot999","rot999",90,dBeta-dRotAngle,108.2,90+dBeta-dRotAngle,18.2,90+dBeta);    
   
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(dBetaRad); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("CRIC",7,"CRT",vector.X(),vector.Y(),vector.Z(),idrotm[1006], "ONLY");
   
}

//_____________________________________________________________________________
void AliCRTv1::CreateMolasse()
{
  //
  //
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
  gMC->Gsvolu("CMO1", "TUBS", idtmed[1123], tspar, 5);
  gMC->Gspos("CMO1", 1, "CRT", 0., 500., 1900.-tspar[2]+400., 0, "MANY");

  Float_t tbox[3];
  tbox[0] = 1250.;
  tbox[1] = (4420. - 1670.)/2.;
  tbox[2] = (1900.+1150.)/2. + 200.;
  gMC->Gsvolu("CM12", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM12", 1, "CRT", 0., 4420. -tbox[1], 1900.-tbox[2]+400., 0, "MANY");

  AliMatrix(idrotm[2003], 0., 0., 90., 0., 90., 90.);
  // Along the PM25
  Float_t tube[3];
  tube[0] = 455. + 100.;
  tube[1] = 555. + 375.;
  tube[2] = (5150. - 1166.)/2.;
  gMC->Gsvolu("CMO2", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO2", 1, "CRT", -2100., 4420.-tube[2], 0., idrotm[2003], "MANY");


  // Along the PGC2
  tube[0] = 650.;
  tube[1] = 2987.7;
  tube[2] = (5150. - 690.)/2.;
  gMC->Gsvolu("CMO3", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO3", 1, "CRT", 375., 4420.-tube[2], 1900.+2987.7, idrotm[2003], "MANY");
  // Behind the PGC2 up to the end of the M. volume.
  tbox[0] = 12073.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (12073. - 1900.-2987.7-650.)/2.;
  gMC->Gsvolu("CMO7", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO7", 1, "CRT", 0., 4420.-tbox[1], 1900.+2987.7+650.+tbox[2], 0, "MANY");

  // Along the PX24 , upper part.
  tube[0] = 1250.;
  tube[1] = 2300;
  tube[2] = 2575. - 1300. + 95.;
  gMC->Gsvolu("CMO4", "TUBE", idtmed[1123], tube, 3);
  gMC->Gspos("CMO4", 1, "CRT", 0., 404.+1300.+tube[2], -2300., idrotm[2003], "MANY");

  // Along the PX24 , lower part
  tspar[0] = 1250.;
  tspar[1] = 2300;
  tspar[2] = 1300.;
  tspar[3] = kRaddeg*TMath::ASin(1070./1150.);
  tspar[4] = 360. - tspar[3];
  gMC->Gsvolu("CMO5", "TUBS", idtmed[1123], tspar, 5);
  gMC->Gspos("CMO5", 1, "CRT", 0., 404., -2300., idrotm[2003], "MANY");
  // behind the PX24
  tbox[0] = 12073.;
  tbox[1] = 2575. + 95.;
  tbox[2] = 8523./2.;
  gMC->Gsvolu("CMO6", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO6", 1, "CRT", 0., 4420.-tbox[1], -3550.-tbox[2], 0, "MANY");


  // On the right side of th hall
  tbox[0] = (12073. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (8437.7+650.)/2.;
  gMC->Gsvolu("CMO8", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO8", 1, "CRT", 1250.+tbox[0], 4420.-tbox[1], -3550.+tbox[2], 0, "MANY");

  // on the left side of the hall, behind 
  tbox[0] = (12073. - 2755.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (8437.7+650.)/2.;
  gMC->Gsvolu("CMO9", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CMO9", 1, "CRT", -2755.-tbox[0], 4420.-tbox[1], -3550.+tbox[2], 0, "MANY");


  // Molasse betwen the PX24 & PM25 on the left side.
  tbox[0] = (2755. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (3550. - 555.)/2.;
  gMC->Gsvolu("CM10", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM10", 1, "CRT", -1250.-tbox[0], 4420.-tbox[1], -tbox[2]-555., 0, "MANY");


  // Molasse betwen the PGC2 & PM25 on the left side.
  tbox[0] = (2755. - 1250.)/2.;
  tbox[1] = 2575. + 95.;
  tbox[2] = (1900.+2987.7 - 555. + 650.)/2.;
  gMC->Gsvolu("CM11", "BOX", idtmed[1123], tbox, 3);
  gMC->Gspos("CM11", 1, "CRT", -1250.-tbox[0], 4420.-tbox[1], 555.+tbox[2], 0, "MANY");


}

//_____________________________________________________________________________
void AliCRTv1::CreateShafts()
{
  //
  //
  //
  Int_t  idrotm[2499];    // The rotation matrix.

  Int_t * idtmed = fIdtmed->GetArray() - 1099 ;

  // Create a mother volume.
  Float_t pbox[3];
  //pbox[0] = AliCRTConstants::fgDepth*TMath::Tan(67.5*kDegrad);
  pbox[0] = 12073.;
  pbox[1] = AliCRTConstants::fgDepth;
  pbox[2] = pbox[0];
  gMC->Gsvolu("CRT", "BOX", idtmed[1114], pbox, 3);
  gMC->Gspos("CRT", 1, "ALIC", 0., 0., 0., 0, "ONLY");

  // HAll ceiling
  Float_t ptubs[5];
  ptubs[0] = 1070.;
  ptubs[1] = 1170.;
  ptubs[2] = 1900.;
  ptubs[3] = 0.;
  ptubs[4] = 180.;
  gMC->Gsvolu("CHC1", "TUBS", idtmed[1116], ptubs, 5);
  gMC->Gspos("CHC1", 1, "CRT", 0., 500., 0., 0, "ONLY");


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
  gMC->Gspos("CSF1", 1, "CRT", 0., 404., -2300., idrotm[2001], "MANY");

  Float_t ptube[3];
  ptube[0] = ptubs[0];
  ptube[1] = ptubs[1];
  ptube[2] = 2575. - ptubs[2] + 95.;
  gMC->Gsvolu("CSF2", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF2", 1, "CRT", 0., 404.+ptubs[2]+ptube[2], -2300., idrotm[2001], "MANY");
  
  // Concrete walls along the shaft
  pbox[0] = 585./2.;
  pbox[1] = 2575. + 95.;
  pbox[2] = 20.;
  gMC->Gsvolu("CSW1", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW1", 1, "CRT", -290-pbox[0], 404.-1300.+pbox[1], -3450.+210.*2, 0, "MANY");
  
  //
  pbox[0] = 750./2.;
  pbox[1] = 2575. + 95.;
  pbox[2] = 20.;
  gMC->Gsvolu("CSW3", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW3", 1, "CRT", 420.-290.+pbox[0], 404.-1300.+pbox[1], -3450.+210.*2, 0, "MANY");
  
  //
  pbox[0] = 60.;
  pbox[1] = 2575. + 95.;
  pbox[2] = 210.;
  gMC->Gsvolu("CSW2", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSW2", 1, "CRT", -290-pbox[0], 404.-1300.+pbox[1], -3450.+pbox[2], 0, "MANY");
  gMC->Gspos("CSW2", 2, "CRT", 420.-290.+pbox[0], 404.-1300.+pbox[1], -3450.+pbox[2], 0, "MANY");
  
  
  // 
  pbox[0] = 1000.;
  pbox[1] = 80.;
  pbox[2] = 200.;
  gMC->Gsvolu("CSP1", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP1", 1, "CRT", 0., 2600.-700., -1150-pbox[2], 0, "MANY");
  
  //
  pbox[0] = 340.8;
  pbox[1] = 300./2.;
  pbox[2] = 460./2.;
  gMC->Gsvolu("CSP2", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP2", 1, "CRT", 0., 2950.-700., -3450+pbox[2], 0, "MANY");
  
  //
  pbox[0] = 600.;
  pbox[1] = 150.;
  pbox[2] = 75.;
  gMC->Gsvolu("CSP3", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP3", 1, "CRT", 0., 2950.-700., -1150.-210.-pbox[2], 0, "MANY");
  
  //
  pbox[0] = 600.;
  pbox[1] = 250.;
  pbox[2] = 38.;
  gMC->Gsvolu("CSP4", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP4", 1, "CRT", 0., 2950.-700.+155.+pbox[1], -1150.-210.-pbox[2], 0, "MANY");
  
  
  // Shielding plug
  pbox[0] = 850.;
  pbox[1] = 90.;
  pbox[2] = 720.;
  gMC->Gsvolu("CSP5", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP5", 1, "CRT", 0., 2950.-700., -3450.+460.+pbox[2], 0, "MANY");
  
  //
  pbox[0] = 80.;
  pbox[1] = 150.;
  pbox[2] = 720.;
  gMC->Gsvolu("CSP6", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP6", 1, "CRT", 1150.-600., 2950.-700., -3450.+460.+pbox[2], 0, "MANY");
  gMC->Gspos("CSP6", 2, "CRT", -1150.+600., 2950.-700., -3450.+460.+pbox[2], 0, "MANY");
  
  
  //
  pbox[0] = 130.;
  pbox[1] = 60.;
  pbox[2] = 750.;
  gMC->Gsvolu("CSP7", "BOX", idtmed[1116], pbox, 3);
  gMC->Gspos("CSP7", 1, "CRT", 850.+pbox[0], 2950.-700.+100., -3450.+460.+pbox[2], 0, "MANY");
  gMC->Gspos("CSP7", 2, "CRT", -850.-pbox[0], 2950.-700.+100., -3450.+460.+pbox[2], 0, "MANY");
  
  
  // PM25 Acces Shaft
  ptube[0] = 910./2.;
  ptube[1] = ptube[0] + 100.;
  ptube[2] = (5150. - 1166.)/2.;
  gMC->Gsvolu("CSF3", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF3", 1, "CRT", -2100., AliCRTConstants::fgDepth-ptube[2], 0., idrotm[2001], "MANY");
  
  // PGC2 Access Shaft
  ptube[0] = 1100./2.;
  ptube[1] = ptube[0] + 100.;
  ptube[2] = (5150. - 690.)/2.;
  gMC->Gsvolu("CSF4", "TUBE", idtmed[1116], ptube, 3);
  gMC->Gspos("CSF4", 1, "CRT", 375., AliCRTConstants::fgDepth-ptube[2], 1900.+2987.7, idrotm[2001], "MANY");

}

//_____________________________________________________________________________
void AliCRTv1::DrawDetector()
{
  //
  // Draw a shaded view of the L3 magnet
  //
  cout << "AliCRTv1::DrawModule() : Drawing the module" << endl;
  
  
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
  static Int_t   vol[5];
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;

  static Float_t hits[14];
  static Float_t eloss;
  static Float_t elossMag;

  if ( !gMC->IsTrackAlive() ) return;

  if (gMC->IsNewTrack()) {
    // Reset the deposited energy
    eloss = 0.;
    elossMag = 0.;
  }

  // Add th energy loss in each step.
  eloss += gMC->Edep();

  gMC->TrackPosition(pos);

  //
  // CRT
  //

  if ( gMC->IsTrackEntering() && (strcmp(gMC->CurrentVolName(),"CRT4") == 0)
       &&(gMC->TrackPid() == kMuonMinus || gMC->TrackPid() == kMuonPlus) ) {
    
    // Get current particle id(ipart),track position (pos) and momentum (mom)
    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);
    ipart = gMC->TrackPid();
    
    vol[0] = 1;
    vol[1] = 0;
    vol[2] = 0;
    vol[3] = 0;
    vol[4] = 0;

    ipart = gMC->TrackPid();
    hits[0]  = (Float_t)ipart; //                 (fId)
    
    hits[1]  = pos[0]; // X coordinate (fX)
    hits[2]  = pos[1]; // Y coordinate (fY)
    hits[3]  = pos[2]; // Z coordinate (fZ)
    hits[4]  = mom[0]; // Px           (fpxug)
    hits[5]  = mom[1]; // Py           (fpyug)
    hits[6]  = mom[2]; // Pz           (fpzug)
    
    hits[7]  = gMC->GetMedium();  //layer(flay)
    hits[8] = eloss;              // Energy loss
    
    hits[9] = 1; // CRT mother activated.
    hits[10] = 0;
    hits[11] = 0;
    hits[12] = 0;
    hits[13] = 0;

    //hits[9] = gAlice->GetCurrentTrackNumber();
    
    AddHit(gAlice->GetCurrentTrackNumber(),vol, hits);
    
    eloss = 0.;

  } else if (gMC->IsTrackEntering()&&(strcmp(gMC->CurrentVolName(),"CRT1")==0)
	     &&(gMC->TrackPid()==kMuonMinus || gMC->TrackPid()==kMuonPlus)) {
    
    vol[0] = 0;
    vol[1] = 1;
    vol[2] = 0;
    vol[3] = 0;
    vol[4] = 0;

    hits[9] = 0; // CRT mother activated.
    hits[10] = 1;
    hits[11] = 0;
    hits[12] = 0;
    hits[13] = 0;

    //hits[10] = 1;
    
    //AddHit(gAlice->GetCurrentTrackNumber(),vol, hits);
    
    //eloss = 0.;


  } else if (gMC->IsTrackEntering()&&(strcmp(gMC->CurrentVolName(),"C3CI")==0)
      &&(gMC->TrackPid()==kMuonMinus || gMC->TrackPid()==kMuonPlus)) {

    //
    // Inside the magnet, upper part.
    //
  
    // Get current particle id(ipart),track position (pos) and momentum (mom)

    vol[0] = 0;
    vol[1] = 0;
    vol[2] = 1;
    vol[3] = 0;
    vol[4] = 0;

    hits[9] = 0; // CRT mother activated.
    hits[10] = 0;
    hits[11] = 1;
    hits[12] = 0;
    hits[13] = 0;
      
    AddHit(gAlice->GetCurrentTrackNumber(),vol, hits);
    
    //eloss = 0.;

  } else if ( gMC->IsTrackEntering()&&(strcmp(gMC->CurrentVolName(),"CRIC")==0)
       && (gMC->TrackPid()==kMuonMinus || gMC->TrackPid()==kMuonPlus) ) {

    //
    // HMPID
    //
    
    // Get current particle id(ipart),track position (pos) and momentum (mom)

    vol[0] = 0;
    vol[1] = 0;
    vol[2] = 0;
    vol[3] = 1;
    vol[4] = 0;

    hits[9] = 0;
    hits[10] = 0;
    hits[11] = 0;
    hits[12] = 1;
    hits[13] = 0;
    
    AddHit(gAlice->GetCurrentTrackNumber(),vol, hits);
    
    //eloss = 0.;


  } else if (gMC->IsTrackEntering()&&(strcmp(gMC->CurrentVolName(),"CSA1")==0)
	     &&(gMC->TrackPid()==kMuonMinus || gMC->TrackPid()==kMuonPlus)) {

    //
    // TPC
    //
    
    // Get current particle id(ipart),track position (pos) and momentum (mom)
    
    vol[0] = 0;
    vol[1] = 0;
    vol[2] = 0;
    vol[3] = 0;
    vol[4] = 1;

    hits[9] = 0;
    hits[10] = 0;
    hits[11] = 0;
    hits[12] = 0;
    hits[13] = 1;

    
    AddHit(gAlice->GetCurrentTrackNumber(),vol, hits);
    
    //eloss = 0.;

  } else {
    return;
  }


}
