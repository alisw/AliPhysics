///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Muon ABSOrber                                                            //
//  This class contains the description of the muon absorber geometry        //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliABSOClass.gif">
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

#include "AliABSO.h"
#include "AliRun.h"
#include "AliConst.h"
 
ClassImp(AliABSO)
 
//_____________________________________________________________________________
AliABSO::AliABSO()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliABSO::AliABSO(const char *name, const char *title)
       : AliModule(name,title)
{
  //
  // Standard constructor
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliABSO::CreateGeometry()
{
  //
  // Creation of the geometry of the muon absorber
  //
  //Begin_Html
  /*
    <img src="gif/AliABSOTree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="gif/AliABSO.gif">
  */
  //End_Html

  AliMC* pMC = AliMC::GetMC();
  
  Int_t *idtmed = gAlice->Idtmed();
  
  Float_t d_pb, cpar[5], dpar[12], tpar[3], zpos,
    cpar1[5], cpar2[5], cpar3[5], cpar4[5], cpar5[12], 
    cpar7[5], cpar8[5], cpar9[5], abs_c, abs_d,
    abs_l, cpar10[5], r_abs;
  Float_t theta1, theta2, abs_cc, d_rear, dz, zr,
    z_cone, d_poly, z_nose, theta_open;
  Float_t acc_min, acc_max, par[50], d_steel, z_w,
    theta_r, epsilon;
  //
  abs_d   = 90.;    // DEFINES DRIFT LENGTH 
  z_nose  = 102.;
  z_cone  = 285.;
  theta1  = 24.;    // 1. angle defining the front absorber 
  theta2  = 5.;     // 2. angle defining the front absorbe 
  acc_max = 9.;     // ANGLE POLAIRE MAXIMUM 
  acc_min = 2.;     // ANGLE POLAIRE MINIMUM DE DETECTION 
  abs_l   = 503.;
  d_steel = 1.;     // THICKNESS OF STEEL SUPPORT 
  d_poly  = 7.5;
  d_pb    = 2.5;
  abs_cc  = 315.;   // DEFINES LENGTH OF CARBON 
  abs_c   = 358.;
  //abs_s   = 150.;   // DEFINES W-SHIELD LENGTH 
  //abs_n   = 80.;    // START OF NOSE 
  r_abs   = 4.;
  //r_pb    = .1;
  epsilon = .01;
  theta_r = 3.;
  d_rear  = 35.;
  theta_open = .75;
  //
  //z_l3     = 700.;
  //zmag_in  = 725.;
  //zmag_out = 1225.;
  //zfil_in  = 1471.;
  //zfil_out = 1591.;
  //zcon_in  = 1900.;
  //zcon_out = 2e3;
  //zcone_e  = 859.0875;
  //spec_l   = 1800.;
  //zplug_in = 1780.;
  //zplug_out= 1900.;
  //
  //     Chamber position 
  //      CZ1=515.5 
  //cz1 = 511.;
  //cz2 = 686.;
  //cz3 = 971.;
  //cz4 = 1245.;
  //cz5 = 1445.;
  //cz6 = 1610.;
  //cz7 = 1710.;
  //
  // --- Outer shape of front absorber 
  par[0] = 0.;
  par[1] = 360.;
  par[2] = 4.;
  par[3] = abs_d;
  par[4] = 0.;
  par[5] = abs_d * TMath::Tan(theta1 * kDegrad);
  par[6] = z_nose;
  par[7] = 0.;
  par[8] = par[6] * TMath::Tan(theta1 * kDegrad);
  par[9] = z_cone;
  par[10] = 0.;
  par[11] = par[8] + (par[9] - par[6]) * TMath::Tan(theta2 * kDegrad);
  par[12] = abs_l;
  par[13] = 0.;
  par[14] = par[11] + (par[12] - par[9]) * TMath::Tan(acc_max * kDegrad);
  pMC->Gsvolu("ABSM", "PCON", idtmed[1605], par, 15);
  //
  // --- Now define all elements of the absorber 
  //
  //       TUNGSTEN NOSE SEGMENT BETWEEN Z=90 AND 112 CM 
  //       SHAPED ALONG A 24 DEG LINE 
  //
  cpar1[0] = (z_nose - abs_d) / 2.;
  cpar1[1] = abs_d * TMath::Tan(acc_max * kDegrad) + d_steel;
  cpar1[2] = abs_d * TMath::Tan(theta1 * kDegrad);
  cpar1[3] = z_nose * TMath::Tan(acc_max * kDegrad) + d_steel;
  cpar1[4] = z_nose * TMath::Tan(theta1 * kDegrad);
  pMC->Gsvolu("ANOS", "CONE", idtmed[1611], cpar1, 5);
  //
  dz = cpar1[0] + abs_d;
  pMC->Gspos("ANOS", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //       IRON  SUPPORT STRUCTURE 
  //
  cpar2[0] = (abs_l - abs_d) / 2.;
  cpar2[1] = abs_d * TMath::Tan(acc_max * kDegrad);
  cpar2[2] = cpar2[1] + d_steel;
  cpar2[3] = abs_l * TMath::Tan(acc_max * kDegrad);
  cpar2[4] = cpar2[3] + d_steel;
  pMC->Gsvolu("ASST", "CONE", idtmed[1658], cpar2, 5);
  dz = cpar2[0] + abs_d;
  pMC->Gspos("ASST", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //       PB FRONT SHIELD INNER SEGMENT, ALSO POLYETHYLENE WAS 
  //       CONSIDERED FOR THIS REGION 
  //
  cpar3[0] = (z_cone - z_nose) / 2.;
  cpar3[1] = cpar1[3];
  cpar3[2] = cpar1[3] + d_poly;
  cpar3[3] = z_cone * TMath::Tan(acc_max * kDegrad) + d_steel;
  cpar3[4] = cpar3[3] + d_poly;
  pMC->Gsvolu("AWFS", "CONE", idtmed[1652], cpar3, 5);
  dz = cpar3[0] + z_nose;
  pMC->Gspos("AWFS", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //       PB OUTER SURFACE 
  //
  cpar5[0] = 0.;
  cpar5[1] = 360.;
  cpar5[2] = 3.;
  cpar5[3] = z_nose;
  cpar5[4] = z_nose * TMath::Tan(acc_max * kDegrad) + d_steel + d_poly;
  cpar5[5] = z_nose * TMath::Tan(theta1 * kDegrad);
  cpar5[6] = z_cone;
  cpar5[7] = z_cone * TMath::Tan(acc_max * kDegrad) + d_steel + d_poly;
  cpar5[8] = cpar5[7] + d_pb;
  cpar5[9] = abs_l;
  cpar5[10] = abs_l * TMath::Tan(acc_max * kDegrad) + d_steel + d_poly;
  cpar5[11] = cpar5[10] + d_pb;
  pMC->Gsvolu("APBS", "PCON", idtmed[1612], cpar5, 12);
  dz = 0.;
  pMC->Gspos("APBS", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //     POLYETHYLEN LAYER 
  //
  cpar4[0] = (abs_l - z_cone) / 2.;
  cpar4[1] = z_cone * TMath::Tan(acc_max * kDegrad) + d_steel;
  cpar4[2] = cpar4[1] + d_poly;
  cpar4[3] = abs_l * TMath::Tan(acc_max * kDegrad) + d_steel;
  cpar4[4] = cpar4[3] + d_poly;
  pMC->Gsvolu("APOL", "CONE", idtmed[1657], cpar4, 5);
  dz = cpar4[0] + z_cone;
  pMC->Gspos("APOL", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //     LEAD INNER SHIELD (inner radius const up to z=abs_c) 
  //
  z_w      = r_abs / TMath::Tan(acc_min * kDegrad);
  cpar8[0] = (abs_c - z_w) / 2.;
  cpar8[1] = r_abs;
  cpar8[2] = r_abs + epsilon;
  cpar8[3] = r_abs;
  cpar8[4] = abs_c * TMath::Tan(acc_min * kDegrad);
  pMC->Gsvolu("AWI1", "CONE", idtmed[1652], cpar8, 5);
  dz = cpar8[0] + z_w;
  pMC->Gspos("AWI1", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //     TUNGSTEN OPENING CONE UP TO THE END 
  //
  cpar8[0] = (abs_l - abs_c) / 2.;
  cpar8[1] = r_abs;
  cpar8[2] = abs_c * TMath::Tan(acc_min * kDegrad);
  cpar8[3] = cpar8[1] + cpar8[0] * 2. * TMath::Tan(theta_open * kDegrad);
  cpar8[4] = abs_l * TMath::Tan(acc_min * kDegrad);
  pMC->Gsvolu("AWI2", "CONE", idtmed[1651], cpar8, 5);
  dz = cpar8[0] + abs_c;
  pMC->Gspos("AWI2", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //     CONCRETE CONE 
  //
  cpar7[0] = (abs_l - d_rear - abs_cc) / 2.;
  cpar7[1] = abs_cc * TMath::Tan(acc_min * kDegrad);
  cpar7[2] = abs_cc * TMath::Tan(acc_max * kDegrad);
  cpar7[3] = (abs_l - d_rear) * TMath::Tan(acc_min * kDegrad);
  cpar7[4] = (abs_l - d_rear) * TMath::Tan(acc_max * kDegrad);
  pMC->Gsvolu("ACON", "CONE", idtmed[1656], cpar7, 5);
  dz = cpar7[0] + abs_cc;
  pMC->Gspos("ACON", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //     REAR SHIELD 
  //
  zr = abs_l - d_rear;
  cpar9[0] = 2.5;
  cpar9[1] = zr * TMath::Tan(theta_r * kDegrad);
  cpar9[2] = zr * TMath::Tan(acc_max * kDegrad);
  cpar9[3] = cpar9[1] + TMath::Tan(theta_r * kDegrad) * 5.;
  cpar9[4] = cpar9[2] + TMath::Tan(acc_max * kDegrad) * 5.;
  pMC->Gsvolu("ARE1", "CONE", idtmed[1652], cpar9, 5);
  dz  = cpar9[0] + zr;
  zr += 5.;
  pMC->Gspos("ARE1", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  cpar9[1] = zr * TMath::Tan(theta_r * kDegrad);
  cpar9[2] = zr * TMath::Tan(acc_max * kDegrad);
  cpar9[3] = cpar9[1] + TMath::Tan(theta_r * kDegrad) * 5.;
  cpar9[4] = cpar9[2] + TMath::Tan(acc_max * kDegrad) * 5.;
  pMC->Gsvolu("ARE2", "CONE", idtmed[1657], cpar9, 5);
  dz  = cpar9[0] + zr;
  zr += 5.;
  pMC->Gspos("ARE2", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  cpar9[1] = zr * TMath::Tan(theta_r * kDegrad);
  cpar9[2] = zr * TMath::Tan(acc_max * kDegrad);
  cpar9[3] = cpar9[1] + TMath::Tan(theta_r * kDegrad) * 5.;
  cpar9[4] = cpar9[2] + TMath::Tan(acc_max * kDegrad) * 5.;
  pMC->Gsvolu("ARE3", "CONE", idtmed[1652], cpar9, 5);
  dz  = cpar9[0] + zr;
  zr += 5.;
  pMC->Gspos("ARE3", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  cpar9[1] = zr * TMath::Tan(theta_r * kDegrad);
  cpar9[2] = zr * TMath::Tan(acc_max * kDegrad);
  cpar9[3] = cpar9[1] + TMath::Tan(theta_r * kDegrad) * 5.;
  cpar9[4] = cpar9[2] + TMath::Tan(acc_max * kDegrad) * 5.;
  pMC->Gsvolu("ARE4", "CONE", idtmed[1657], cpar9, 5);
  dz  = cpar9[0] + zr;
  zr += 5.;
  pMC->Gspos("ARE4", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  cpar9[1] = zr * TMath::Tan(theta_r * kDegrad);
  cpar9[2] = zr * TMath::Tan(acc_max * kDegrad);
  cpar9[3] = cpar9[1] + TMath::Tan(theta_r * kDegrad) * 5.;
  cpar9[4] = cpar9[2] + TMath::Tan(acc_max * kDegrad) * 5.;
  pMC->Gsvolu("ARE5", "CONE", idtmed[1652], cpar9, 5);
  dz  = cpar9[0] + zr;
  zr += 5.;
  pMC->Gspos("ARE5", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  cpar9[1] = zr * TMath::Tan(theta_r * kDegrad);
  cpar9[2] = zr * TMath::Tan(acc_max * kDegrad);
  cpar9[3] = cpar9[1] + TMath::Tan(theta_r * kDegrad) * 5.;
  cpar9[4] = cpar9[2] + TMath::Tan(acc_max * kDegrad) * 5.;
  pMC->Gsvolu("ARE6", "CONE", idtmed[1657], cpar9, 5);
  dz  = cpar9[0] + zr;
  zr += 5.;
  pMC->Gspos("ARE6", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  cpar9[1] = zr * TMath::Tan(theta_r * kDegrad);
  cpar9[2] = zr * TMath::Tan(acc_max * kDegrad);
  cpar9[3] = cpar9[1] + TMath::Tan(theta_r * kDegrad) * 5.;
  cpar9[4] = cpar9[2] + TMath::Tan(acc_max * kDegrad) * 5.;
  pMC->Gsvolu("ARE7", "CONE", idtmed[1612], cpar9, 5);
  dz  = cpar9[0] + zr;
  zr += 5.;
  pMC->Gspos("ARE7", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //     TUNGSTEN REAR SHIELD INNER PART 
  //
  zr = abs_l - d_rear;
  cpar10[0] = d_rear / 2.;
  cpar10[1] = zr * TMath::Tan(acc_min * kDegrad);
  cpar10[2] = zr * TMath::Tan(theta_r * kDegrad);
  cpar10[3] = cpar10[1] + d_rear * TMath::Tan(acc_min * kDegrad);
  cpar10[4] = cpar10[2] + d_rear * TMath::Tan(theta_r * kDegrad);
  pMC->Gsvolu("ARIN", "CONE", idtmed[1611], cpar10, 5);
  dz = cpar10[0] + zr;
  pMC->Gspos("ARIN", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //     ELEMENTS OF THE BEAM PIPE TO BE POSITIONED INTO THE ABSORBER 
  //
  //     MOTHER VOLUME 1. SEGMENT 
  //
  tpar[0] = 0.;
  tpar[1] = r_abs;
  tpar[2] = (abs_c - abs_d) / 2.;
  pMC->Gsvolu("AATU", "TUBE", idtmed[1655], tpar, 3);
  //
  tpar[1] = r_abs - .8;
  tpar[0] = tpar[1] - .2;
  tpar[2] = (abs_c - abs_d) / 2.;
  pMC->Gsvolu("ATUB", "TUBE", idtmed[1649], tpar, 3);
  dz = 0.;
  pMC->Gspos("ATUB", 1, "AATU", 0., 0., dz, 0, "ONLY");
  //
  dz = (abs_c - abs_d) / 2. + abs_d;
  pMC->Gspos("AATU", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  //     MOTHER VOLUME 2. SEGMENT 
  //
  cpar[0] = (abs_l - abs_c) / 2.;
  cpar[1] = 0.;
  cpar[2] = r_abs;
  cpar[3] = 0.;
  cpar[4] = cpar[2] + cpar[0] * 2. * TMath::Tan(theta_open * kDegrad);
  pMC->Gsvolu("AAT1", "CONE", idtmed[1655], cpar, 5);
  //
  cpar[0]  = (abs_l - abs_c) / 2.;
  cpar[2] += -.8;
  cpar[1]  = cpar[2] - .2;
  cpar[4] += -.8;
  cpar[3]  = cpar[4] - .2;
  pMC->Gsvolu("ATU1", "CONE", idtmed[1649], cpar, 5);
  dz = 0.;
  pMC->Gspos("ATU1", 1, "AAT1", 0., 0., dz, 0, "ONLY");
  //
  dz = (abs_l - abs_c) / 2. + abs_c;
  pMC->Gspos("AAT1", 1, "ABSM", 0., 0., dz, 0, "ONLY");
  //
  pMC->Gspos("ABSM", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  //
  //       absorber support structure 
  //
  //       attention this element is positioned into ALIC 
  //
  dpar[0] = 0.;
  dpar[1] = 360.;
  dpar[2] = 3.;
  dpar[3] = abs_l;
  dpar[4] = abs_l * TMath::Tan(acc_max * kDegrad);
  dpar[5] = dpar[4] + 4. / TMath::Cos(acc_max * kDegrad);
  dpar[6] = 600.;
  dpar[7] = TMath::Tan(acc_max * kDegrad) * 600;
  dpar[8] = dpar[7] + 4. / TMath::Cos(acc_max * kDegrad);
  dpar[9] = 670.;
  dpar[10] = 159.;
  dpar[11] = 163.5;
  pMC->Gsvolu("ASUP", "PCON", idtmed[1618], dpar, 12);
  dz = 0.;
  pMC->Gspos("ASUP", 1, "ALIC", 0., 0., dz, 0, "ONLY");
  //
  //     Flange at the entrance of the absorber 
  //
  tpar[0] = 3.;
  tpar[1] = 5.7;
  tpar[2] = 2.;
  pMC->Gsvolu("AF63", "TUBE", idtmed[1618], tpar, 3);
  zpos = abs_d + tpar[2];
  pMC->Gspos("AF63", 1, "ABSM", 0., 0., zpos, 0, "ONLY");
}

//_____________________________________________________________________________
void AliABSO::DrawModule()
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
  pMC->Gsatt("ABSM","seen",1);
  pMC->Gsatt("ANOS","seen",1);
  pMC->Gsatt("ASST","seen",1);
  pMC->Gsatt("AWFS","seen",1);
  pMC->Gsatt("APBS","seen",1);
  pMC->Gsatt("APOL","seen",1);
  pMC->Gsatt("AWI1","seen",1);
  pMC->Gsatt("AWI2","seen",1);
  pMC->Gsatt("ACON","seen",1);
  pMC->Gsatt("ARE1","seen",1);
  pMC->Gsatt("ARE2","seen",1);
  pMC->Gsatt("ARE3","seen",1);
  pMC->Gsatt("ARE4","seen",1);
  pMC->Gsatt("ARE5","seen",1);
  pMC->Gsatt("ARE6","seen",1);
  pMC->Gsatt("ARE7","seen",1);
  pMC->Gsatt("ARIN","seen",1);
  pMC->Gsatt("AATU","seen",1);
  pMC->Gsatt("ATUB","seen",1);
  pMC->Gsatt("AAT1","seen",1);
  pMC->Gsatt("ATU1","seen",1);
  pMC->Gsatt("ASUP","seen",1);
  pMC->Gsatt("AF63","seen",1);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 21.5, 15, .04, .04);
  pMC->Gdhead(1111, "Muon Absorber");
  pMC->Gdman(16, 6, "MAN");
  pMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliABSO::CreateMaterials()
{
  //
  // Define materials for muon absorber
  //
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  Float_t apoly[2]  = { 12.01,1. };
  Float_t zpoly[2]  = { 6.,1. };
  Float_t wpoly[2]  = { .33,.67 };
  Float_t aconc[10] = { 1.,12.01,15.994,22.99,24.305,26.98,
			28.086,39.1,40.08,55.85 };
  Float_t zconc[10] = { 1.,6.,8.,11.,12.,13.,14.,19.,20.,26. };
  Float_t wconc[10] = { .01,.001,.529107,.016,.002,.033872,
			.337021,.013,.044,.014 };
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  //
  //     Carbon 
  AliMaterial(6, "CARBON$   ", 12.01, 6., 2.265, 18.8, 49.9);
  AliMaterial(26, "CARBON$   ", 12.01, 6., 2.265, 18.8, 49.9);
  AliMaterial(46, "CARBON$   ", 12.01, 6., 2.265, 18.8, 49.9);
  //
  //     Aluminum 
  AliMaterial(9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(49, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  //
  //     Iron 
  AliMaterial(10, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  //
  //     Tungsten 
  AliMaterial(12, "TUNGSTEN$ ", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(32, "TUNGSTEN$ ", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(52, "TUNGSTEN$ ", 183.85, 74., 19.3, .35, 10.3);
  //
  //     Lead 
  AliMaterial(13, "LEAD$     ", 207.19, 82., 11.35, .56, 18.5);
  AliMaterial(33, "LEAD$     ", 207.19, 82., 11.35, .56, 18.5);
  AliMaterial(53, "LEAD$     ", 207.19, 82., 11.35, .56, 18.5);
  //
  //     Air 
  AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500.);
  AliMaterial(35, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500.);
  AliMaterial(55, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500.);
  //
  //     Vacuum 
  AliMaterial(16, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(36, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(56, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  //
  //     Concrete 
  AliMixture(17, "CONCRETE$", aconc, zconc, 2.35, 10, wconc);
  AliMixture(37, "CONCRETE$", aconc, zconc, 2.35, 10, wconc);
  AliMixture(57, "CONCRETE$", aconc, zconc, 2.35, 10, wconc);
  //
  //     Poly CH2 
  AliMixture(18, "POLYETHYLEN$", apoly, zpoly, .95, -2, wpoly);
  //
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  //
  AliMixture(38, "POLYETHYLEN$", apoly, zpoly, .95, 2, wpoly);
  AliMixture(58, "POLYETHYLEN$", apoly, zpoly, .95, 2, wpoly);
  //
  //     stainless Steel 
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(39, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(59, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  //
  // **************** 
  //     Defines tracking media parameters. 
  //
  epsil  = .001;  // Tracking precision, 
  stemax = -1.;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 
  //
  //    Carbon 
  AliMedium(1606, "C_C0             ", 6, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1626, "C_C1             ", 26, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1646, "C_C2             ", 46, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Aluminum 
  AliMedium(1609, "ALU_C0          ", 9, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1629, "ALU_C1          ", 29, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1649, "ALU_C2          ", 49, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Iron 
  AliMedium(1610, "FE_C0           ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1630, "FE_C1           ", 30, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1650, "FE_C2           ", 50, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Tungsten 
  AliMedium(1612, "W_C0            ", 12, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1632, "W_C1            ", 32, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1652, "W_C2            ", 52, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Lead 
  AliMedium(1613, "PB_C0           ", 13, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1633, "PB_C1           ", 33, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1653, "PB_C2           ", 53, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Air 
  AliMedium(1615, "AIR_C0          ", 15, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1635, "AIR_C1          ", 35, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1655, "AIR_C2          ", 55, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Vacuum 
  AliMedium(1616, "VA_C0           ", 16, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1636, "VA_C1           ", 36, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1656, "VA_C2           ", 56, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Concrete 
  AliMedium(1617, "CC_C0           ", 17, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1637, "CC_C1           ", 37, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1657, "CC_C2           ", 57, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Polyethilene 
  AliMedium(1618, "CH2_C0 B        ", 18, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1638, "CH2_C1          ", 38, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1658, "CH2_C2          ", 58, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Steel 
  AliMedium(1619, "ST_C0           ", 19, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1639, "ST_C1           ", 39, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1659, "ST_C3           ", 59, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
void AliABSO::Init()
{
  //
  // Initialisation of the muon absorber after it has been built
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" ABSO_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}
 
