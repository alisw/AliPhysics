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
Revision 1.6.2.1  1999/12/03 16:38:51  fca
Correct overlap in magnet

Revision 1.6  1999/09/29 09:24:30  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Magnetic Dipole version 1                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliDIPOv2Class.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:andreas.morsch@cern.ch">Andreas Morsch</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliDIPOv2.h"
#include "AliRun.h"
#include "AliConst.h"
 
ClassImp(AliDIPOv2)
 
//_____________________________________________________________________________
AliDIPOv2::AliDIPOv2()
{
  //
  // Default constructor for magnetic dipole version 2
  //
}
 
//_____________________________________________________________________________
AliDIPOv2::AliDIPOv2(const char *name, const char *title)
  : AliDIPO(name,title)
{
  //
  // Standard constructor for the magnetic dipole version 2    
   SetMarkerColor(7);
   SetMarkerStyle(2);
   SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliDIPOv2::CreateGeometry()
{
  //
  // Creation of the geometry of the magnetic DIPOLE version 2
  //
  //Begin_Html
  /*
    <img src="picts/AliDIPOv2Tree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliDIPOv2.gif">
  */
  //End_Html

  Float_t cpar[5], tpar[3], ypar[4];
  Float_t dz, dx, dy;
  Int_t idrotm[1899];
  Float_t acc_max, the1, phi1, the2, phi2, the3, phi3;

  Int_t *idtmed = fIdtmed->GetArray()-1799;
  
  //abs_d   = 90.;  // DEFINES DRIFT LENGTH 
  //z_nose  = 102.;
  //z_cone  = 285.;
  //theta1  = 24.;  // 1. angle defining the front absorber 
  //theta2  = 5.;   // 2. angle defining the front absorbe 
  acc_max = 9.;   // ANGLE POLAIRE MAXIMUM 
  //acc_min = 2.;   // ANGLE POLAIRE MINIMUM DE DETECTION 
  //abs_l   = 503.;
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
  
  
  //       DIPOLE MAGNET 
  
  tpar[0] = 300.;
  tpar[1] = 415.;
  tpar[2] = 250.;
  gMC->Gsvolu("DDIP", "BOX ", idtmed[1814], tpar, 3);
  
  //       COILS 
  
  // air - m.f. 
  cpar[0] = 210.;
  cpar[1] = 263.;
  cpar[2] = 83.2/2.;
  cpar[3] = 120.;
  cpar[4] = 240.;
  gMC->Gsvolu("DC1 ", "TUBS", idtmed[1813], cpar, 5);
  cpar[3] = -60.;
  cpar[4] = 60.;
  gMC->Gsvolu("DC2 ", "TUBS", idtmed[1813], cpar, 5);
  // ... define regions for higher cuts 
  cpar[0] += 10.;
  cpar[1] += -10.;
  cpar[2] += -10.;
  cpar[3]  = 120.;
  cpar[4]  = 240.;
  gMC->Gsvolu("DC3 ", "TUBS", idtmed[1833], cpar, 5);
  gMC->Gspos("DC3 ", 1, "DC1 ", 0., 0., 0., 0, "ONLY");
  cpar[3] = -60.;
  cpar[4] = 60.;
  gMC->Gsvolu("DC4 ", "TUBS", idtmed[1833], cpar, 5);
  gMC->Gspos("DC4 ", 1, "DC2 ", 0., 0., 0., 0, "ONLY");
  // ... 
  dz = 83.2/2. - 250.;
  gMC->Gspos("DC1 ", 1, "DDIP", 0., 0.,  dz, 0, "ONLY");
  gMC->Gspos("DC1 ", 2, "DDIP", 0., 0., -dz, 0, "ONLY");
  gMC->Gspos("DC2 ", 1, "DDIP", 0., 0.,  dz, 0, "ONLY");
  gMC->Gspos("DC2 ", 2, "DDIP", 0., 0., -dz, 0, "ONLY");
  the1 = 180.;
  phi1 = 0.;
  the2 = 90.;
  phi2 = 150.;
  the3 = 90.;
  phi3 = 60.;
  AliMatrix(idrotm[1800], the1, phi1, the2, phi2, the3, phi3);
  phi2 = 30.;
  the3 = -90.;
  phi3 = -60.;
  AliMatrix(idrotm[1801], the1, phi1, the2, phi2, the3, phi3);
  the1 = 0.;
  phi1 = 0.;
  the2 = 90.;
  phi2 = 150.;
  the3 = 90.;
  phi3 = 60.;
  AliMatrix(idrotm[1802], the1, phi1, the2, phi2, the3, phi3);
  phi2 = 30.;
  the3 = -90.;
  phi3 = -60.;
  AliMatrix(idrotm[1803], the1, phi1, the2, phi2, the3, phi3);
  cpar[0] = 25.;
  cpar[1] = 108.2;
  cpar[2] = 26.5;
  cpar[3] = 270.;
  cpar[4] = 360.;
  gMC->Gsvolu("DC11", "TUBS", idtmed[1813], cpar, 5);
  // ... higher cuts 
  cpar[0] += 10.;
  cpar[1] += -10.;
  cpar[2] += -10.;
  gMC->Gsvolu("DC21", "TUBS", idtmed[1833], cpar, 5);
  gMC->Gspos("DC21", 1, "DC11", 0., 0., 0., 0, "ONLY");
  // ... 
  dx = TMath::Sin(30*kDegrad) * -236.5;
  dy = TMath::Cos(30*kDegrad) * -236.5;
  dz = cpar[1] + 10. - 250.;
  gMC->Gspos("DC11", 1, "DDIP",  dx, dy,  dz, idrotm[1800], "ONLY");
  gMC->Gspos("DC11", 2, "DDIP",  dx, dy, -dz, idrotm[1802], "ONLY");
  gMC->Gspos("DC11", 3, "DDIP", -dx, dy,  dz, idrotm[1801], "ONLY");
  gMC->Gspos("DC11", 4, "DDIP", -dx, dy, -dz, idrotm[1803], "ONLY");
  cpar[0] = 25.;
  cpar[1] = 25.+83.2;
  cpar[2] = 53./2.;
  cpar[3] = 0.;
  cpar[4] = 90.;
  gMC->Gsvolu("DC12", "TUBS", idtmed[1813], cpar, 5);
  // ... higher cuts 
  cpar[0] += 10.;
  cpar[1] += -10.;
  cpar[2] += -10.;
  gMC->Gsvolu("DC22", "TUBS", idtmed[1833], cpar, 5);
  gMC->Gspos("DC22", 1, "DC12", 0., 0., 0., 0, "ONLY");
  // ... 
  dx = TMath::Sin(30*kDegrad) * -236.5;
  dy = TMath::Cos(30*kDegrad) * 236.5;
  dz = cpar[1] + 10. - 250.;
  gMC->Gspos("DC12", 1, "DDIP",  dx, dy,  dz, idrotm[1801], "ONLY");
  gMC->Gspos("DC12", 2, "DDIP",  dx, dy, -dz, idrotm[1803], "ONLY");
  gMC->Gspos("DC12", 3, "DDIP", -dx, dy,  dz, idrotm[1800], "ONLY");
  gMC->Gspos("DC12", 4, "DDIP", -dx, dy, -dz, idrotm[1802], "ONLY");
  the1 = 90.;
  phi1 = 60.;
  the2 = 90.;
  phi2 = 150.;
  the3 = 0.;
  phi3 = 0.;
  AliMatrix(idrotm[1804], the1, phi1, the2, phi2, the3, phi3);
  the1 = 90.;
  phi1 = 120.;
  the2 = 90.;
  phi2 = 210.;
  the3 = 0.;
  phi3 = 0.;
  AliMatrix(idrotm[1805], the1, phi1, the2, phi2, the3, phi3);
  tpar[0] = 53./2.;
  tpar[1] = 83.2/2.;
  tpar[2] = 283.6/2.;
  gMC->Gsvolu("DL1 ", "BOX ", idtmed[1813], tpar, 3);
  // ... higher cuts 
  tpar[0] -= 10.;
  tpar[1] -= 10.;
  gMC->Gsvolu("DL2 ", "BOX ", idtmed[1833], tpar, 3);
  gMC->Gspos("DL2 ", 1, "DL1 ", 0., 0., 0., 0, "ONLY");
  // ... 
  dx = -60.5;
  dy = -238.;
  dz = 0.;
  gMC->Gspos("DL1 ", 1, "DDIP", dx,  dy, dz, idrotm[1804], "ONLY");
  gMC->Gspos("DL1 ", 2, "DDIP", dx, -dy, dz, idrotm[1805], "ONLY");
  gMC->Gspos("DL1 ", 3, "DDIP",-dx,  dy, dz, idrotm[1805], "ONLY");
  gMC->Gspos("DL1 ", 4, "DDIP",-dx, -dy, dz, idrotm[1804], "ONLY");
  
  //   YOKE 
  
  ypar[1] = 275.8;
  ypar[2] = 5.;
  ypar[3] = 156.8;
  ypar[0] = ypar[1] - ypar[3] * 2. * TMath::Tan(acc_max * kDegrad);
  gMC->Gsvolu("DY1 ", "TRD1", idtmed[1809], ypar, 4);
  // iron - 
  dy = 283.5;
  dy = 279;
  gMC->Gspos("DY1 ", 1, "DDIP", 0.,  dy, 0., 0, "ONLY");
  gMC->Gspos("DY1 ", 2, "DDIP", 0., -dy, 0., 0, "ONLY");
  ypar[2] = 60.;
  gMC->Gsvolu("DY2 ", "TRD1", idtmed[1829], ypar, 4);
  // iron - 
  dy = ypar[2] + 284.;
  gMC->Gspos("DY2 ", 1, "DDIP", 0.,  dy, 0., 0, "ONLY");
  gMC->Gspos("DY2 ", 2, "DDIP", 0., -dy, 0., 0, "ONLY");
  the1 = 99.;
  phi1 = 0.;
  the2 = 90.;
  phi2 = 90.;
  the3 = 9.;
  phi3 = 0.;
  AliMatrix(idrotm[1806], the1, phi1, the2, phi2, the3, phi3);
  the1 = 261.;
  phi1 = 0.;
  the3 = 171.;
  AliMatrix(idrotm[1807], the1, phi1, the2, phi2, the3, phi3);
  tpar[0] = 60.;
  tpar[1] = 283.;
  tpar[2] = 156.8;
  gMC->Gsvolu("DYL ", "BOX ", idtmed[1814], tpar, 3);
  tpar[0] = 5.;
  tpar[1] = 73.;
  gMC->Gsvolu("DY3 ", "BOX ", idtmed[1809], tpar, 3);
  dx = tpar[0] - 60.;
  dy = tpar[1] + 137.;
  gMC->Gspos("DY3 ", 1, "DYL ", dx,  dy, 0., 0, "ONLY");
  gMC->Gspos("DY3 ", 2, "DYL ", dx, -dy, 0., 0, "ONLY");
  tpar[0] = 55.;
  gMC->Gsvolu("DY4 ", "BOX ", idtmed[1829], tpar, 3);
  dx = dx + 5. + tpar[0];
  gMC->Gspos("DY4 ", 1, "DYL ", dx,  dy, 0., 0, "ONLY");
  gMC->Gspos("DY4 ", 2, "DYL ", dx, -dy, 0., 0, "ONLY");
  tpar[0] = 37.7;
  tpar[1] = 137.;
  gMC->Gsvolu("DY5 ", "BOX ", idtmed[1829], tpar, 3);
  dx = 60. - tpar[0];
  gMC->Gspos("DY5 ", 1, "DYL ", dx, 0., 0., 0, "ONLY");
  tpar[0] = 5.;
  gMC->Gsvolu("DY6 ", "BOX ", idtmed[1809], tpar, 3);
  dx = dx - 37.7 - tpar[0];
  gMC->Gspos("DY6 ", 1, "DYL ", dx, 0., 0., 0, "ONLY");
  tpar[0] = 17.3;
  tpar[1] = 5.;
  gMC->Gsvolu("DY7 ", "BOX ", idtmed[1809], tpar, 3);
  dx = tpar[0] - 60.;
  dy = tpar[1] + 127.;
  gMC->Gspos("DY7 ", 1, "DYL ", dx,  dy, 0., 0, "ONLY");
  gMC->Gspos("DY7 ", 2, "DYL ", dx, -dy, 0., 0, "ONLY");
  
  dx = ypar[0] + ypar[3] * TMath::Tan(acc_max * kDegrad) - 60.;
  gMC->Gspos("DYL ", 1, "DDIP", dx, 0., 0., idrotm[1806], "ONLY");
  gMC->Gspos("DYL ", 2, "DDIP",-dx, 0., 0., idrotm[1807], "ONLY");
  gMC->Gspos("DDIP", 1, "ALIC", 0., 0., 725.+250., 0, "MANY");
  gMC->Gsatt("DDIP", "SEEN", 0);
  gMC->Gsatt("DC21", "SEEN", 0);
  gMC->Gsatt("DC22", "SEEN", 0);
  gMC->Gsatt("DC3 ", "SEEN", 0);
  gMC->Gsatt("DC4 ", "SEEN", 0);
}

//_____________________________________________________________________________
void AliDIPOv2::DrawModule()
{
  //
  // Draw a shaded view of the muon absorber
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("DDIP","seen",0);
  gMC->Gsatt("DC1 ","seen",1);
  gMC->Gsatt("DC2 ","seen",1);
  gMC->Gsatt("DC3 ","seen",1);
  gMC->Gsatt("DC4 ","seen",1);
  gMC->Gsatt("DC11","seen",1);
  gMC->Gsatt("DC21","seen",1);
  gMC->Gsatt("DC12","seen",1);
  gMC->Gsatt("DC22","seen",1);
  gMC->Gsatt("DL1 ","seen",1);
  gMC->Gsatt("DL2 ","seen",1);
  gMC->Gsatt("DY1 ","seen",1);
  gMC->Gsatt("DY2 ","seen",1);
  gMC->Gsatt("DYL ","seen",1);
  gMC->Gsatt("DY3 ","seen",1);
  gMC->Gsatt("DY4 ","seen",1);
  gMC->Gsatt("DY5 ","seen",1);
  gMC->Gsatt("DY6 ","seen",1);
  gMC->Gsatt("DY7 ","seen",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox(".");
  gMC->DefaultRange();
  gMC->Gdraw("alic", 30, 30, 0, 17, 13.5, .019, .019);
  gMC->Gdhead(1111, "Magnetic Dipole Version 2");
  gMC->Gdman(16, 4, "MAN");
}

//_____________________________________________________________________________
void AliDIPOv2::CreateMaterials()
{
  //
  // Create Materials for Magnetic Dipole version 2
  //
  
  Int_t ISXFLD   = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  Float_t acoil[3]  = { 26.98,1.01,16. };
  Float_t zcoil[3]  = { 13.,1.,8. };
  Float_t wcoil[3]  = { .66,.226,.114 };
  
  Float_t epsil, stmin, deemax, tmaxfd, stemax;
  
  // --- Define the various materials for GEANT --- 
  //     Aluminum 
  AliMaterial(9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(49, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  
  //     Iron 
  AliMaterial(10, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  
  //     Air 
  AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
  AliMaterial(35, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
  AliMaterial(55, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
  
  //     Vacuum 
  AliMaterial(16, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(36, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(56, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  
  //     stainless Steel 
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(39, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(59, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  
  //     Coil 
  AliMixture(14, "Al$", acoil, zcoil, 2.122, 3, wcoil);
  AliMixture(34, "Al$", acoil, zcoil, 2.122, 3, wcoil);
  AliMixture(54, "Al$", acoil, zcoil, 2.122, 3, wcoil);
  
  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  epsil  = .001;  // Tracking precision, 
  stemax = -1.;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 
  
  //    Aluminum 
  AliMedium(9, "ALU_C0          ",  9, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(29, "ALU_C1          ", 29, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(49, "ALU_C2          ", 49, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Iron 
  AliMedium(10, "FE_C0           ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(30, "FE_C1           ", 30, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(50, "FE_C2           ", 50, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Air 
  AliMedium(15, "AIR_C0          ", 15, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_C1          ", 35, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(55, "AIR_C2          ", 55, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Vacuum 
  AliMedium(16, "VA_C0           ", 16, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(36, "VA_C1           ", 36, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(56, "VA_C2           ", 56, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Steel 
  AliMedium(19, "ST_C0           ", 19, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(39, "ST_C1           ", 39, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(59, "ST_C3           ", 59, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Coil 
  AliMedium(14, "Coil_C1         ", 14, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(34, "Coil_C2         ", 34, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(54, "Coil_C3         ", 54, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}

