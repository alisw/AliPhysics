///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Muon Shield Class                                                        //
//  This class contains a description of the muon shield                     //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliSHILClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliSHIL.h"
#include <TNode.h>
#include <TTUBE.h>
#include <TBRIK.h>
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"

ClassImp(AliSHIL)
 
//_____________________________________________________________________________
AliSHIL::AliSHIL()
{
  //
  // Default constructor for muon shield
  //
}
 
//_____________________________________________________________________________
AliSHIL::AliSHIL(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor for muon shield
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliSHIL::BuildGeometry()
{
  //
  // Root TNode geometry only for active detectors
  //
}
 
//_____________________________________________________________________________
void AliSHIL::CreateGeometry()
{
  //
  // Build muon shield geometry
  //
  // Origin N.Van Eijndhoven
  //
  //Begin_Html
  /*
    <img src="gif/AliSHIL.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="gif/AliSHILTree.gif">
  */
  //End_Html

  Float_t r_pb, cpar[5], parm[12], tpar[3], zmag_out,
    abs_c, abs_l, r_abs, pcpar[12];
  Float_t dr, dz, 
    zs, theta_open;
  Float_t cz1, cz2;
  Float_t dzcoch1, dzcoch2, acc_min, dzs, zcone_e, zmag_in,
    z_close, z_l3, shl1, shl4;
    
  AliMC* pMC = AliMC::GetMC();
  
  Int_t *idtmed = gAlice->Idtmed();
  
  //abs_d   = 90.;  // DEFINES DRIFT LENGTH 
  //z_nose  = 102.;
  //z_cone  = 285.;
  //theta1  = 24.;  // 1. angle defining the front absorber 
  //theta2  = 5.;   // 2. angle defining the front absorbe 
  //acc_max = 9.;   // ANGLE POLAIRE MAXIMUM 
  acc_min = 2.;   // ANGLE POLAIRE MINIMUM DE DETECTION 
  abs_l   = 503.;
  //d_steel = 1.;   // THICKNESS OF STEEL SUPPORT 
  //d_poly  = 7.5;
  //d_pb    = 2.5;
  //abs_cc  = 315.; // DEFINES LENGTH OF CARBON 
  abs_c   = 358.;
  //abs_s   = 150.; // DEFINES W-SHIELD LENGTH 
  //abs_n   = 80.;  // START OF NOSE 
  r_abs   = 4.;
  r_pb    = .1;
  //epsilon = .01;
  //theta_r = 3.;
  //d_rear  = 35.;
  theta_open = .75;
  
  z_l3      = 700.;
  zmag_in   = 725.;
  zmag_out  = 1225.;
  //zfil_in   = 1471.;
  //zfil_out  = 1591.;
  //zcon_in   = 1900.;
  //zcon_out  = 2e3;
  zcone_e   = 859.0875;
  //spec_l    = 1800.;
  //zplug_in  = 1780.;
  //zplug_out = 1900.;
  
  //     Chamber position 
  //      CZ1=515.5 
  cz1 = 511.;
  cz2 = 686.;
  //cz3 = 971.;
  //cz4 = 1245.;
  //cz5 = 1445.;
  //cz6 = 1610.;
  //cz7 = 1710.;
  
  
  //  the mother of all shields 
  
  parm[0] = 0.;
  parm[1] = 360.;
  parm[2] = 3.;
  parm[3] = abs_l;
  parm[4] = 0.;
  parm[5] = abs_l * TMath::Tan(acc_min * kDegrad);
  parm[6] = zcone_e;
  parm[7] = 0.;
  parm[8] = 30.;
  parm[9] = 1900.;
  parm[10] = 0.;
  parm[11] = 30.;
  pMC->Gsvolu("YMOT", "PCON", idtmed[1715], parm, 12);
  
  //       beam shield ouside absorber inside l3 
  
  shl1    = (z_l3 - abs_l) / 2.;
  cpar[0] = shl1;
  cpar[1] = 0.;
  cpar[2] = abs_l * TMath::Tan(acc_min * kDegrad);
  cpar[3] = 0.;
  cpar[4] = z_l3 * TMath::Tan(acc_min * kDegrad);
  pMC->Gsvolu("YMS1", "CONE", idtmed[1759], cpar, 5);
  
  //       OUTSIDE STEEL TUBE 
  // Pb/W 
  cpar[1] = cpar[2] - 4.;
  cpar[3] = cpar[4] - 4.;
  pMC->Gsvolu("YSH1", "CONE", idtmed[1718], cpar, 5);
  pMC->Gspos("YSH1", 1, "YMS1", 0., 0., 0., 0, "ONLY");
  
  //      COCHES FOR CHAMBERS 
  dzcoch1 = cz1 + 4. - (abs_l + shl1);
  dzcoch2 = cz2 + 4. - (abs_l + shl1);
  
  //       AIR 
  cpar[0] = 7.5;
  dr      = cpar[0] * 2. * TMath::Tan(acc_min * kDegrad);
  cpar[2] = (cz1 + 4. - cpar[0]) * TMath::Tan(acc_min * kDegrad);
  cpar[1] = cpar[2] - 3.;
  cpar[3] = cpar[1] + dr;
  cpar[4] = cpar[2] + dr;
  pMC->Gsvolu("YAC1", "CONE", idtmed[1714], cpar, 5);
  pMC->Gspos("YAC1", 1, "YMS1", 0., 0., dzcoch1, 0, "ONLY");
  cpar[2] = (cz2 + 4. - cpar[0]) * TMath::Tan(acc_min * kDegrad);
  cpar[1] = cpar[2] - 3.;
  cpar[3] = cpar[1] + dr;
  cpar[4] = cpar[2] + dr;
  pMC->Gsvolu("YAC2", "CONE", idtmed[1714], cpar, 5);
  pMC->Gspos("YAC2", 1, "YMS1", 0., 0., dzcoch2, 0, "ONLY");
  
  //       STEEL 
  cpar[0] = (cz1 + 4. + 20. - abs_l) / 2.;
  dr      = cpar[0] * 2. * TMath::Tan(acc_min * kDegrad);
  //      CPAR(3)=ABS_L*TAN(kDegrad*ACC_MIN)-3. 
  cpar[2] = abs_l * TMath::Tan(acc_min * kDegrad) - 4.;
  cpar[1] = cpar[2] - 3.;
  cpar[3] = cpar[1] + dr;
  cpar[4] = cpar[2] + dr;
  dzcoch1 = abs_l + cpar[0] - (abs_l + shl1);
  pMC->Gsvolu("YSC1", "CONE", idtmed[1718], cpar, 5);
  pMC->Gspos("YSC1", 1, "YMS1", 0., 0., dzcoch1, 0, "ONLY");
  
  cpar[0] = 16.;
  dr      = cpar[0] * 2. * TMath::Tan(acc_min * kDegrad);
  //      CPAR(3)=(CZ2+4.-CPAR(1))*TAN(kDegrad*ACC_MIN)-3. 
  cpar[2] = (cz2 + 4. - cpar[0]) * TMath::Tan(acc_min * kDegrad) - 4.;
  cpar[1] = cpar[2] - 3.;
  cpar[3] = cpar[1] + dr;
  cpar[4] = cpar[2] + dr;
  pMC->Gsvolu("YSC2", "CONE", idtmed[1718], cpar, 5);
  pMC->Gspos("YSC2", 1, "YMS1", -4., 0., dzcoch2, 0, "ONLY");
  
  //       ... beam pipe 
  cpar[0] = shl1;
  cpar[1] = 0.;
  cpar[2] = r_abs + (abs_l - abs_c) * TMath::Tan(theta_open * kDegrad);
  cpar[3] = 0.;
  cpar[4] = cpar[2] + shl1 * 2. * TMath::Tan(theta_open * kDegrad);
  pMC->Gsvolu("YMB1", "CONE", idtmed[1755], cpar, 5);
  
  cpar[0]  = shl1;
  cpar[2] += -.8;
  cpar[1]  = cpar[2] - .2;
  cpar[4] += -.8;
  cpar[3]  = cpar[4] - .2;
  pMC->Gsvolu("YBS1", "CONE", idtmed[1749], cpar, 5);
  pMC->Gspos("YBS1", 1, "YMB1", 0., 0., 0., 0, "ONLY");
  
  
  pMC->Gspos("YMB1", 1, "YMS1", 0., 0., 0., 0, "ONLY");
  dz = shl1 + abs_l;
  pMC->Gspos("YMS1", 1, "YMOT", 0., 0., dz, 0, "ONLY");
  
  //       BEAM SHIELD OUTSIDE L3 
  
  //     L3->DIPOLE 
  cpar[0] = (zmag_in - z_l3) / 2.;
  cpar[1] = 0.;
  cpar[2] = z_l3 * TMath::Tan(acc_min * kDegrad);
  cpar[3] = 0.;
  cpar[4] = zmag_in * TMath::Tan(acc_min * kDegrad);
  pMC->Gsvolu("YMS2", "CONE", idtmed[1759], cpar, 5);
  
  //       OUTSIDE STEEL TUBE 
  
  // Pb/W 
  cpar[1] = cpar[2] - 4.;
  cpar[3] = cpar[4] - 4.;
  pMC->Gsvolu("YSH2", "CONE", idtmed[1718], cpar, 5);
  pMC->Gspos("YSH2", 1, "YMS2", 0., 0., 0., 0, "ONLY");
  
  //       ... beam pipe 

  zs      = z_l3 - abs_c;
  cpar[1] = 0.;
  cpar[2] = r_abs + zs * TMath::Tan(theta_open * kDegrad);
  cpar[3] = 0.;
  cpar[4] = cpar[2] + cpar[0] * 2. * TMath::Tan(theta_open * kDegrad);
  pMC->Gsvolu("YMB2", "CONE", idtmed[1755], cpar, 5);
  
  cpar[2] += -.8;
  cpar[1]  = cpar[2] - .2;
  cpar[4] += -.8;
  cpar[3]  = cpar[4] - .2;
  pMC->Gsvolu("YBS2", "CONE", idtmed[1749], cpar, 5);
  pMC->Gspos("YBS2", 1, "YMB2", 0., 0., 0., 0, "ONLY");
  
  
  pMC->Gspos("YMB2", 1, "YMS2", 0., 0., 0., 0, "ONLY");
  dz = cpar[0] + z_l3;
  pMC->Gspos("YMS2", 1, "YMOT", 0., 0., dz, 0, "ONLY");
  pcpar[0] = 0.;
  pcpar[1] = 360.;
  pcpar[2] = 3.;
  pcpar[3] = 0.;
  pcpar[4] = 0.;
  pcpar[5] = zmag_in * TMath::Tan(acc_min * kDegrad);
  pcpar[6] = zcone_e - zmag_in;
  pcpar[7] = 0.;
  pcpar[8] = 30.;
  pcpar[9] = zmag_out - zmag_in;
  pcpar[10] = 0.;
  pcpar[11] = 30.;
  pMC->Gsvolu("YMS3", "PCON", idtmed[1759], pcpar, 12);
  
  //       OUTSIDE STEEL TUBE 
  
  // Pb/W 
  pcpar[4]  = pcpar[5] - 4.;
  pcpar[7]  = pcpar[8] - 4.;
  pcpar[10] = pcpar[11] - 4.;
  pMC->Gsvolu("YSH3", "PCON", idtmed[1718], pcpar, 12);
  pMC->Gspos("YSH3", 1, "YMS3", 0., 0., 0., 0, "MANY");
  
  //       ... beam pipe up to closing cone 
  
  zs      = zmag_in - abs_c;
  z_close = 804.;
  
  cpar[0] = (z_close - zmag_in) / 2.;
  cpar[1] = 0.;
  cpar[2] = r_abs + zs * TMath::Tan(theta_open * kDegrad);
  cpar[3] = 0.;
  cpar[4] = cpar[2] + cpar[0] * 2. * TMath::Tan(theta_open * kDegrad);
  pMC->Gsvolu("YMB3", "CONE", idtmed[1755], cpar, 5);
  
  cpar[2] += -.8;
  cpar[1]  = cpar[2] - .2;
  cpar[4] += -.8;
  cpar[3]  = cpar[4] - .2;
  pMC->Gsvolu("YBS3", "CONE", idtmed[1749], cpar, 5);
  pMC->Gspos("YBS3", 1, "YMB3", 0., 0., 0.,      0, "ONLY");
  pMC->Gspos("YMB3", 1, "YMS3", 0., 0., cpar[0], 0, "ONLY");
  
  //       .closing cone 
  
  dzs = cpar[0] * 2.;
  zs  = z_close - abs_c;
  
  cpar[0] = 25.;
  cpar[1] = 0.;
  cpar[2] = r_abs + zs * TMath::Tan(theta_open * kDegrad);
  cpar[3] = 0.;
  cpar[4] = r_abs;
  pMC->Gsvolu("YMB5", "CONE", idtmed[1755], cpar, 5);
  
  cpar[2] += -.8;
  cpar[1]  = cpar[2] - .2;
  cpar[4] += -.8;
  cpar[3]  = cpar[4] - .2;
  pMC->Gsvolu("YBS5", "CONE", idtmed[1749], cpar, 5);
  pMC->Gspos("YBS5", 1, "YMB5", 0., 0., 0., 0, "ONLY");
  
  dzs += cpar[0];
  pMC->Gspos("YMB5", 1, "YMS3", 0., 0., dzs, 0, "ONLY");
  dzs += cpar[0];
  
  //       OUTSIDE PB-TUBE 
  
  tpar[0] = 30. - r_pb - 4.;
  tpar[1] = 30. - r_pb;
  tpar[2] = (zmag_out - z_close - 50.) / 2.;
  dzs    += tpar[2];
  pMC->Gsvolu("YNW1", "TUBE", idtmed[1752], tpar, 3);
  pMC->Gspos("YNW1", 1, "YMS3", 0., 0., dzs, 0, "MANY");
  
  //       constant beam pipe up to end of magnet 
  
  tpar[0] = 0.;
  tpar[1] = r_abs;
  tpar[2] = (zmag_out - z_close - 50.) / 2.;
  pMC->Gsvolu("YMB6", "TUBE", idtmed[1755], tpar, 3);
  
  tpar[1] = r_abs - .8;
  tpar[0] = tpar[1] - .2;
  pMC->Gsvolu("YBS6", "TUBE", idtmed[1749], tpar, 3);
  pMC->Gspos("YBS6", 1, "YMB6", 0., 0., 0.,  0, "ONLY");
  pMC->Gspos("YMB6", 1, "YMS3", 0., 0., dzs, 0, "ONLY");
  
  dz = zmag_in;
  pMC->Gspos("YMS3", 1, "YMOT", 0., 0., dz,  0, "ONLY");
  
  //       DIPOLE-> 
  
  tpar[0] = 0.;
  tpar[1] = 30.;
  tpar[2] = (1900. - zmag_out) / 2.;
  shl4 = tpar[2];
  pMC->Gsvolu("YMS4", "TUBE", idtmed[1759], tpar, 3);
  //      CALL GSVOLU('YMS4','TUBE',IDTMED(1752),TPAR,3,IL3) ! W 
  
  //       OUTSIDE STEEL TUBE 
  
  // Pb/W 
  tpar[0] = tpar[1] - 4.;
  pMC->Gsvolu("YSH4", "TUBE", idtmed[1718], tpar, 3);
  pMC->Gspos("YSH4", 1, "YMS4", 0., 0., 0., 0, "MANY");
  
  //       OUTSIDE PB-TUBE 
  
  tpar[0] = 30. - r_pb - 4.;
  tpar[1] = 30. - r_pb;
  pMC->Gsvolu("YNW2", "TUBE", idtmed[1752], tpar, 3);
  pMC->Gspos("YNW2", 1, "YMS4", 0., 0., 0., 0, "MANY");
  
  //       ... beam pipe 
  
  tpar[0] = 0.;
  tpar[1] = r_abs;
  tpar[2] = shl4;
  pMC->Gsvolu("YMB4", "TUBE", idtmed[1755], tpar, 3);
  
  tpar[2] = shl4;
  tpar[1] = r_abs - .8;
  tpar[0] = tpar[1] - .2;
  pMC->Gsvolu("YBS4", "TUBE", idtmed[1749], tpar, 3);
  pMC->Gspos("YBS4", 1, "YMB4", 0., 0., 0., 0, "ONLY");
  
  
  dz = zmag_out + shl4;
  pMC->Gspos("YMB4", 1, "YMS4", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("YMS4", 1, "YMOT", 0., 0., dz, 0, "ONLY");
  pMC->Gspos("YMOT", 1, "ALIC", 0., 0., 0., 0, "ONLY");
}

//_____________________________________________________________________________
void AliSHIL::CreateMaterials()
{
  //
  // Defines materials for the muon shield
  //
  // Origin N. Van Eijndhoven
  //
  
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  Float_t apbw[2]   = { 207.2,183.85 };
  Float_t zpbw[2]   = { 82.,74. };
  Float_t wpbw[2]   = { .5,.5 };
  
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  
  //     STEEL 
  
  
  //     LEAD/TUNGSTEN MIXTURE 
  
  
  // --- Define the various materials for GEANT --- 
  //     Aluminum 
  AliMaterial(9,  "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(49, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  
  //     Iron 
  AliMaterial(10, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON$     ", 55.85, 26., 7.87, 1.76, 17.1);
  
  //     Tungsten 
  AliMaterial(12, "TUNGSTEN$ ", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(32, "TUNGSTEN$ ", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(52, "TUNGSTEN$ ", 183.85, 74., 19.3, .35, 10.3);
  
  //     Lead 
  AliMaterial(13, "LEAD$     ", 207.19, 82., 11.35, .56, 18.5);
  AliMaterial(33, "LEAD$     ", 207.19, 82., 11.35, .56, 18.5);
  AliMaterial(53, "LEAD$     ", 207.19, 82., 11.35, .56, 18.5);
  
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
  
  //     Lead/Tungsten 
  AliMixture(20, "LEAD/TUNGSTEN$", apbw, zpbw, 15.325, 2, wpbw);
  AliMixture(40, "LEAD/TUNGSTEN$", apbw, zpbw, 15.325, 2, wpbw);
  AliMixture(60, "LEAD/TUNGSTEN$", apbw, zpbw, 15.325, 2, wpbw);
  
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
  AliMedium(1709, "ALU_C0          ", 9, 0,  ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1729, "ALU_C1          ", 29, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1749, "ALU_C2          ", 49, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Iron 
  AliMedium(1710, "FE_C0           ", 10, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1730, "FE_C1           ", 30, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1750, "FE_C2           ", 50, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Tungsten 
  AliMedium(1712, "W_C0            ", 12, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1732, "W_C1            ", 32, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1752, "W_C2            ", 52, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Lead 
  AliMedium(1713, "PB_C0           ", 13, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1733, "PB_C1           ", 33, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1753, "PB_C2           ", 53, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Air 
  AliMedium(1715, "AIR_C0          ", 15, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1735, "AIR_C1          ", 35, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1755, "AIR_C2          ", 55, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Vacuum 
  AliMedium(1716, "VA_C0           ", 16, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1736, "VA_C1           ", 36, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1756, "VA_C2           ", 56, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Steel 
  AliMedium(1719, "ST_C0           ", 19, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1739, "ST_C1           ", 39, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1759, "ST_C3           ", 59, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Lead/Tungsten 
  AliMedium(1720, "PB/W0           ", 20, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1740, "PB/W1           ", 40, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1760, "PB/W3           ", 60, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
void AliSHIL::DrawDetector () 
{
  //
  // Draw a shaded view of the muon shield
  //

  AliMC* pMC = AliMC::GetMC();
  
  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  pMC->Gsatt("YMOT","seen",1);
  pMC->Gsatt("YMS1","seen",1);
  pMC->Gsatt("YSH1","seen",1);
  pMC->Gsatt("YAC1","seen",1);
  pMC->Gsatt("YAC2","seen",1);
  pMC->Gsatt("YSC1","seen",1);
  pMC->Gsatt("YSC2","seen",1);
  pMC->Gsatt("YMB1","seen",1);
  pMC->Gsatt("YBS1","seen",1);
  pMC->Gsatt("YMS2","seen",1);
  pMC->Gsatt("YSH2","seen",1);
  pMC->Gsatt("YMB2","seen",1);
  pMC->Gsatt("YBS2","seen",1);
  pMC->Gsatt("YMS3","seen",1);
  pMC->Gsatt("YSH3","seen",1);
  pMC->Gsatt("YMB3","seen",1);
  pMC->Gsatt("YBS3","seen",1);
  pMC->Gsatt("YMB5","seen",1);
  pMC->Gsatt("YBS5","seen",1);
  pMC->Gsatt("YNW1","seen",1);
  pMC->Gsatt("YMB6","seen",1);
  pMC->Gsatt("YBS6","seen",1);
  pMC->Gsatt("YMS4","seen",1);
  pMC->Gsatt("YSH4","seen",1);
  pMC->Gsatt("YNW2","seen",1);
  pMC->Gsatt("YMB4","seen",1);
  pMC->Gsatt("YBS4","seen",1);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 30, 30, 0, 26.5, 18, .03, .03);
  pMC->Gdhead(1111, "Muon Shield");
  pMC->Gdman(16, 6, "MAN");
}

//_____________________________________________________________________________
void AliSHIL::Init()
{
  //
  // Initialise the muon shield after it has been built
  //
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" SHIL_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the ABSO initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

 
//_____________________________________________________________________________
void AliSHIL::StepManager()
{
  //
  // Called at every step in the muon shield
  //
}
