///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Magnetic Dipole version 1                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliDIPOv1Class.gif">
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

#include "AliDIPOv1.h"
#include "AliRun.h"
 
ClassImp(AliDIPOv1)
 
//_____________________________________________________________________________
AliDIPOv1::AliDIPOv1()
{
  //
  // Default constructor for the magnetic dipole version 1
  //
}
 
//_____________________________________________________________________________
AliDIPOv1::AliDIPOv1(const char *name, const char *title)
       : AliDIPO(name,title)
{
  //
  // Standard constructor for magnetic dipole version 1
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliDIPOv1::CreateGeometry()
{
  //
  // Creation of the geometry of the magnetic DIPOLE version 1
  //
  //Begin_Html
  /*
    <img src="picts/AliDIPOv1Tree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliDIPOv1.gif">
  */
  //End_Html

  AliMC* pMC = AliMC::GetMC();
  
  Float_t par[5];
  
  Int_t *idtmed = fIdtmed->GetArray()-1799;
  
  //abs_d   = 90.;  // DEFINES DRIFT LENGTH 
  //z_nose  = 102.;
  //z_cone  = 285.;
  //theta1  = 24.;  // 1. angle defining the front absorber 
  //theta2  = 5.;   // 2. angle defining the front absorbe 
  //acc_max = 9.;   // ANGLE POLAIRE MAXIMUM 
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
  par[0] = 0.;
  par[1] = 280.;
  par[2] = 250.;
  pMC->Gsvolu("DDIP", "TUBE", idtmed[1801], par, 3);
  
  //       COIL 
  par[0] = 250.;
  par[1] = 125.;
  par[2] = 165.;
  par[3] = 204.;
  par[4] = 244.;
  
  pMC->Gsvolu("DIPC", "CONE", idtmed[1810], par, 5);
  pMC->Gspos("DIPC", 1, "DDIP", 0., 0., 0., 0, "ONLY");
  par[0] = 250.;
  par[1] = 115.;
  par[2] = 125.;
  par[3] = 194.;
  par[4] = 204.;
  pMC->Gsvolu("DIIC", "CONE", idtmed[1807], par, 5);
  pMC->Gspos("DIIC", 1, "DDIP", 0., 0., 0., 0, "ONLY");
  
  //       YOKE 
  par[0] = 250.;
  par[1] = 165.;
  par[2] = 195.;
  par[3] = 244.;
  par[4] = 274.;
  
  pMC->Gsvolu("DIPY", "CONE", idtmed[1834], par, 5);
  pMC->Gspos("DIPY", 1, "DDIP", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("DDIP", 1, "ALIC", 0., 0., 725.+250, 0, "ONLY");
}

//_____________________________________________________________________________
void AliDIPOv1::DrawModule()
{
  //
  // Draw a shaded view of the muon absorber
  //

  AliMC* pMC = AliMC::GetMC();
  
  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  pMC->Gsatt("DDIP","seen",1);
  pMC->Gsatt("DIPC","seen",1);
  pMC->Gsatt("DIIC","seen",1);
  pMC->Gsatt("DIPY","seen",1);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox(".");
  pMC->DefaultRange();
  pMC->Gdraw("alic", 30, 30, 0, 17, 13.5, .019, .019);
  pMC->Gdhead(1111, "Magnetic Dipole Version 1");
  pMC->Gdman(16, 4, "MAN");
}

//_____________________________________________________________________________
void AliDIPOv1::CreateMaterials()
{
  //
  // Create Materials for Dipole Magnet version 1
  //
  
  Int_t ISXFLD   = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  Float_t epsil, stmin, deemax, tmaxfd, stemax;
  
  //     STEEL 
  
  
  // --- Define the various materials for GEANT --- 
  AliMaterial(9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
  AliMaterial(10, "IRON$     ", 55.85, 26., 7.87, 0, 17.1);
  AliMaterial(16, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMixture(24, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  
  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  epsil  = .001; // Tracking precision, 
  stemax = -1.;  // Maximum displacement for multiple scat 
  tmaxfd = -20.; // Maximum angle due to field deflection 
  deemax = -.3;  // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 
  
  //    Air 

  AliMedium(1, "AIR_DI_US         ", 15, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(2, "AIR_DI_US         ", 15, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(3, "AIR_L3_US         ", 15, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Aluminum 
  
  AliMedium(8, "ALU_DI_US         ", 9, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(11, "ALU_DI_SH         ", 9, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Iron 
  
  AliMedium(31, "FE_NF_US          ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(32, "FE_DI_US          ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(33, "FE_L3_US          ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(34, "FE_NF_SH          ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "FE_DI_SH          ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(36, "FE_L3_SH          ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Vacuum 
  
  AliMedium(37, "VA_NF_US          ", 16, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(38, "VA_DI_US          ", 16, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(39, "VA_L3_US          ", 16, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Steel 
  
  AliMedium(75, "ST_L3_US          ", 24, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}

