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
//  Beam pipe class. Test version                                            //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPIPEClass.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TVirtualMC.h>

#include "AliMagF.h"
#include "AliPIPEv1.h"
#include "AliRun.h"
 
ClassImp(AliPIPEv1)
 
//_____________________________________________________________________________
AliPIPEv1::AliPIPEv1()
{
  //
  // Default constructor for beam pipe
  //
}
 
//_____________________________________________________________________________
AliPIPEv1::AliPIPEv1(const char *name, const char *title)
       : AliPIPE(name,title)
{
  //
  // Standard constructor for beam pipe
  //
}
 
//_____________________________________________________________________________
void AliPIPEv1::CreateGeometry()
{
  //
  // Create Beam Pipe geometry
  //
  //Begin_Html
  /*
    <img src="picts/AliPIPE.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliPIPETree.gif">
  */
  //End_Html

  Float_t tpar[3], dzmo, zpos, absorberDistance, absorberEnd;
  Float_t r2, dr;

  const Double_t kZFlange = 150;
  
  Int_t *idtmed = fIdtmed->GetArray()-1999;
  
  
  absorberDistance   = 90.;  // DEFINES DRIFT LENGTH 
  //z_nose  = 102.;
  //z_cone  = 285.;
  //theta1  = 24.;  // 1. angle defining the front absorber 
  //theta2  = 5.;   // 2. angle defining the front absorbe 
  //acc_max = 9.;   // ANGLE POLAIRE MAXIMUM 
  //acc_min = 2.;   // ANGLE POLAIRE MINIMUM DE DETECTION 
  absorberEnd   = 503.;
  //d_steel = 1.;   // THICKNESS OF STEEL SUPPORT 
  //d_poly  = 7.5;
  //d_pb    = 2.5;
  //abs_cc  = 315.; // DEFINES LENGTH OF CARBON 
  //abs_c   = 358.;
  //abs_s   = 150.; // DEFINES W-SHIELD LENGTH 
  //abs_n   = 80.;  // START OF NOSE 
  //r_abs   = 4.;
  //r_pb    = .1;
  //epsilon = .01;
  //theta_r = 3.;
  //d_rear  = 35.;
  //theta_open = .75;
  
  //z_l3      = 700.;
  //zmag_in   = 725.;
  //zmag_out  = 1225.;
  //zfil_in   = 1471.;
  //zfil_out  = 1591.;
  //zcon_in   = 1900.;
  //zcon_out  = 2e3;
  //zcone_e   = 859.0875;
  //spec_l    = 1800.;
  //zplug_in  = 1780.;
  //zplug_out = 1900.;
  
  //     Chamber position 
  //      CZ1=515.5 
  //cz1 = 511.;
  //cz2 = 686.;
  //cz3 = 971.;
  //cz4 = 1245.;
  //cz5 = 1445.;
  //cz6 = 1610.;
  //cz7 = 1710.;
  
  
  //     the mother of all beam pipes 
  
  tpar[0] = 0.;
  tpar[1] = 3.;
  tpar[2] = (absorberDistance + 700.) / 2.;
  dzmo = tpar[2] - absorberDistance;
  gMC->Gsvolu("QQMO", "TUBE", idtmed[2015], tpar, 3);
  gMC->Gspos("QQMO", 1, "ALIC", 0., 0., -dzmo, 0, "ONLY");
  
  //       BEAM PIPE IN DRIFT SPACE 
  
  //     -30-kZFlange 
  tpar[0] = 0.;
  tpar[1] = 3.;
  tpar[2] = 30;
  gMC->Gsvolu("QDT1", "TUBE", idtmed[2015], tpar, 3);
  
  tpar[0] = 2.9;
  gMC->Gsvolu("QTB1", "TUBE", idtmed[2004], tpar, 3);
  gMC->Gspos("QTB1", 1, "QDT1", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("QDT1", 1, "QQMO", 0., 0., dzmo, 0, "ONLY");
  
  
  //     30-90 
  tpar[0] = 0.;
  tpar[1] = 3.;
  tpar[2] = 30.;
  gMC->Gsvolu("QDT2", "TUBE", idtmed[2015], tpar, 3);
  
  tpar[0] = 2.9;
  gMC->Gsvolu("QTB2", "TUBE", idtmed[2004], tpar, 3);
  gMC->Gspos("QTB2", 1, "QDT2", 0., 0., 0.,   0, "ONLY");
  gMC->Gspos("QDT2", 1, "QQMO", 0., 0., dzmo, 0, "ONLY");
  
  //       beam pipe outside absorber on the left side 
  
  
  
  //     -30 - kZFlange 
  tpar[0] = 0.;
  tpar[1] = 3.;
  tpar[2] = (kZFlange - 30)/2;
  gMC->Gsvolu("QDT5", "TUBE", idtmed[2015], tpar, 3);
  
  tpar[0] = 2.9;
  zpos    = -30. - tpar[2] + dzmo;
  gMC->Gsvolu("QTB5", "TUBE", idtmed[2004], tpar, 3);
  gMC->Gspos("QTB5", 1, "QDT5", 0., 0., 0.,   0, "ONLY");
  gMC->Gspos("QDT5", 1, "QQMO", 0., 0., zpos, 0, "ONLY");
  
  //     STRAIGHT STEEL PIECE 
  
  zpos    = -kZFlange;
  r2      = 2.9;
  dr      = .015;
  tpar[0] = 0.;
  tpar[1] = r2 + dr;
  tpar[2] = (zpos + 700.) / 2.;
  gMC->Gsvolu("QDT7", "TUBE", idtmed[2015], tpar, 3);
  tpar[0] = r2;
  gMC->Gsvolu("QTB7", "TUBE", idtmed[2018], tpar, 3);
  gMC->Gspos("QTB7", 1, "QDT7", 0., 0., 0.,   0, "ONLY");
  zpos = zpos - tpar[2] + dzmo;
  gMC->Gspos("QDT7", 1, "QQMO", 0., 0., zpos, 0, "ONLY");
  
  //     flange dn 63 
  
  tpar[0] = 3.;
  tpar[1] = 5.7;
  tpar[2] = 2.;
  gMC->Gsvolu("QN63", "TUBE", idtmed[2018], tpar, 3);
  zpos = tpar[2] - kZFlange;
  gMC->Gspos("QN63", 1, "ALIC", 0., 0., zpos, 0, "ONLY");
  
  
  //     Replace Absorber or Shield by Beam-Pipe 
  //     in case they are not selected 
  
  if (gAlice->GetModule("ABSO") == 0) {
    
    gMC->Gspos("QN63", 2, "ALIC", 0., 0., kZFlange, 0, "ONLY");
    r2      = 2.9;
    dr      = .1;
    tpar[0] = 0.;
    tpar[1] = r2 + dr;
    tpar[2] = (kZFlange - absorberDistance) / 2.;
    gMC->Gsvolu("QDT8", "TUBE", idtmed[2015], tpar, 3);
    tpar[0] = r2;
    gMC->Gsvolu("QTB8", "TUBE", idtmed[2004], tpar, 3);
    gMC->Gspos("QTB8", 1, "QDT8", 0., 0., 0., 0, "ONLY");
    zpos    = absorberDistance + tpar[2];
    gMC->Gspos("QDT8", 1, "ALIC", 0., 0., zpos, 0, "ONLY");
    dr      = .015;
    tpar[0] = 0.;
    tpar[1] = r2 + dr;
    tpar[2] = (absorberEnd - kZFlange) / 2.;
    gMC->Gsvolu("QDTS", "TUBE", idtmed[2015], tpar, 3);
    tpar[0] = r2;
    gMC->Gsvolu("QTBS", "TUBE", idtmed[2018], tpar, 3);
    gMC->Gspos("QTBS", 1, "QDTS", 0., 0., 0., 0, "ONLY");
    zpos = tpar[2] + kZFlange;
    gMC->Gspos("QDTS", 1, "ALIC", 0., 0., zpos, 0, "ONLY");
  }
  if (gAlice->GetModule("SHIL") == 0) {
    r2      = 2.9;
    dr      = .015;
    tpar[0] = 0.;
    tpar[1] = r2 + dr;
    tpar[2] = (700. - absorberEnd) / 2.;
    gMC->Gsvolu("QDT9", "TUBE", idtmed[2015], tpar, 3);
    tpar[0] = r2;
    gMC->Gsvolu("QTB9", "TUBE", idtmed[2018], tpar, 3);
    gMC->Gspos("QTB9", 1, "QDT9", 0., 0., 0., 0, "ONLY");
    zpos = absorberEnd + tpar[2];
    gMC->Gspos("QDT9", 1, "ALIC", 0., 0., zpos, 0, "ONLY");
  }
}

//_____________________________________________________________________________
void AliPIPEv1::CreateMaterials()
{
  //
  // Create materials for beam pipe
  //

  Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  
  //     STEEL 
  
  
  // --- Define the various materials for GEANT --- 
  AliMaterial(5, "BERILLIUM$", 9.01, 4., 1.848, 35.3, 36.7);
  AliMaterial(16, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  
  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  epsil  = .001;  // Tracking precision, 
  stemax = -1.;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  
  //    Air 
  
  AliMedium(15, "AIR_L3_US", 15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Beryllium 
  
  AliMedium(5, "BE_L3_US", 5, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Vacuum 
  
  AliMedium(16, "VA_L3_US", 16, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Steel 
  
  AliMedium(19, "ST_L3_US", 19, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
}

