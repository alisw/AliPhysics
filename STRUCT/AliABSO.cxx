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
//  Muon ABSOrber                                                            //
//  This class contains the description of the muon absorber geometry        //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliABSOClass.gif">
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
#include "AliMagF.h"
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
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliABSO::CreateGeometry()
{
  //
  // Creation of the geometry of the muon absorber
  //
}

//_____________________________________________________________________________
void AliABSO::DrawModule() const
{
  //
  // Draw a shaded view of the muon absorber
  //
}


//_____________________________________________________________________________
void AliABSO::CreateMaterials()
{
  //
  // Define materials for muon absorber
  //
  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
//
// Air
//
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  Float_t dAir1 = 1.20479E-10;
//
// Polyethylene
//
  Float_t apoly[2]  = { 12.01,1. };
  Float_t zpoly[2]  = { 6.,1. };
  Float_t wpoly[2]  = { .33,.67 };
//
// Concrete
//
  Float_t aconc[10] = { 1.,12.01,15.994,22.99,24.305,26.98,28.086,39.1,40.08,55.85 };
  Float_t zconc[10] = { 1.,6.,8.,11.,12.,13.,14.,19.,20.,26. };
  Float_t wconc[10] = { .01,.001,.529107,.016,.002,.033872, .337021,.013,.044,.014 };
//
// Steel
//  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
//
// Ni-Cu-W alloy
//
  Float_t aniwcu[3] ={58.6934, 183.84, 63.546};
  Float_t zniwcu[3] ={28., 74., 29};
  Float_t wniwcu[3] ={0.015,0.95,0.035};
//
// Poly Concrete
//                      H     Li     F       C      Al     Si      Ca      Pb     O
  Float_t aPolyCc[9] = {1. ,  6.941, 18.998, 12.01, 26.98, 28.086, 40.078, 207.2, 15.999};
  Float_t zPolyCc[9] = {1. ,  3.   ,  9.   ,  6.  , 13.  , 14.   , 20.   ,  82. ,  8.   };
  Float_t wPolyCc[9] = {4.9,  1.2  ,  1.3  ,  1.1 ,  0.15,  0.02 ,  0.06 ,   0.7,  1.1  };
  Float_t wtot=0;
  Int_t   i=0;

  for (i=0; i<9; i++) wtot+=wPolyCc[i];
  for (i=0; i<9; i++) wPolyCc[i]/=wtot;  

//
// Insulation powder
//                    Si         O       Ti     Al
  Float_t ains[4] ={28.0855, 15.9994, 47.867,  26.982};
  Float_t zins[4] ={14.,      8.    , 22.   ,  13.   };
  Float_t wins[4] ={ 0.3019,  0.4887,  0.1914,  0.018};
  
//
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  //
  //     Carbon 
  AliMaterial( 6, "CARBON0$   ", 12.01, 6., 1.75, 24.4, 49.9);
  AliMaterial(26, "CARBON1$   ", 12.01, 6., 1.75, 24.4, 49.9);
  AliMaterial(46, "CARBON2$   ", 12.01, 6., 1.75, 24.4, 49.9);
  //
  //     Magnesium
  AliMaterial( 7, "MAGNESIUM$ ", 24.31, 12., 1.74, 25.3, 46.0);
  //
  //     Aluminum 
  AliMaterial(9,  "ALUMINIUM0$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "ALUMINIUM1$", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(49, "ALUMINIUM2$", 26.98, 13., 2.7, 8.9, 37.2);
  //
  //     Iron 
  AliMaterial(10, "IRON0$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON1$     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON2$     ", 55.85, 26., 7.87, 1.76, 17.1);
  //
  //     Copper
  AliMaterial(11, "COPPER0$   ", 63.55, 29., 8.96, 1.43, 15.1);
  AliMaterial(31, "COPPER1$   ", 63.55, 29., 8.96, 1.43, 15.1);
  AliMaterial(51, "COPPER2$   ", 63.55, 29., 8.96, 1.43, 15.1);
  //
  //     Tungsten 
  AliMaterial(12, "TUNGSTEN0$ ", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(32, "TUNGSTEN1$ ", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(52, "TUNGSTEN2$ ", 183.85, 74., 19.3, .35, 10.3);
  //
  //     Ni-W-Cu 
  AliMixture(21, "Ni-W-Cu0$", aniwcu, zniwcu, 18.78, 3, wniwcu);
  AliMixture(41, "Ni-W-Cu1$", aniwcu, zniwcu, 18.78, 3, wniwcu);
  AliMixture(61, "Ni-W-Cu2$", aniwcu, zniwcu, 18.78, 3, wniwcu);
  //
  //     Lead 
  AliMaterial(13, "LEAD0$     ", 207.19, 82., 11.35, .56, 18.5);
  AliMaterial(33, "LEAD1$     ", 207.19, 82., 11.35, .56, 18.5);
  AliMaterial(53, "LEAD2$     ", 207.19, 82., 11.35, .56, 18.5);
  //
  //     Air 
  AliMixture(15, "AIR0$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(35, "AIR1$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(55, "AIR2$      ", aAir, zAir, dAir, 4, wAir);
  //
  //     Vacuum 
  AliMixture(16, "VACUUM0$ ", aAir, zAir, dAir1, 4, wAir);
  AliMixture(36, "VACUUM1$ ", aAir, zAir, dAir1, 4, wAir);
  AliMixture(56, "VACUUM2$ ", aAir, zAir, dAir1, 4, wAir);
  //
  //     Concrete 
  AliMixture(17, "CONCRETE0$", aconc, zconc, 2.35, 10, wconc);
  AliMixture(37, "CONCRETE1$", aconc, zconc, 2.35, 10, wconc);
  AliMixture(57, "CONCRETE2$", aconc, zconc, 2.35, 10, wconc);
  //
  //     Poly CH2 
  AliMixture(18, "POLYETHYLEN0$", apoly, zpoly, .95, -2, wpoly);
  //
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  //
  AliMixture(38, "POLYETHYLEN1$", apoly, zpoly, .95, 2, wpoly);
  AliMixture(58, "POLYETHYLEN2$", apoly, zpoly, .95, 2, wpoly);
  //
  //     stainless Steel 
  AliMixture(19, "STAINLESS STEEL0$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(39, "STAINLESS STEEL1$", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(59, "STAINLESS STEEL2$", asteel, zsteel, 7.88, 4, wsteel);
  //
  //     Insulation powder 
  AliMixture(14, "INSULATION0$", ains, zins, 0.41, 4, wins);
  AliMixture(34, "INSULATION1$", ains, zins, 0.41, 4, wins);
  AliMixture(54, "INSULATION2$", ains, zins, 0.41, 4, wins);
  // Polymere Concrete 
  AliMixture(20, "Poly Concrete0$", aPolyCc, zPolyCc, 3.53, -9, wPolyCc);
  AliMixture(40, "Poly Concrete1$", aPolyCc, zPolyCc, 3.53,  9, wPolyCc);
  AliMixture(60, "Poly Concrete2$", aPolyCc, zPolyCc, 3.53,  9, wPolyCc);

  //
  // **************** 
  //     Defines tracking media parameters. 
  //
  epsil  = .001;    // Tracking precision, 
  stemax = -0.01;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;    // Maximum angle due to field deflection 
  deemax = -.3;     // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // *************** 
  //
  //    Carbon 
  AliMedium(6,  "C_C0             ",  6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(26, "C_C1             ", 26, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(46, "C_C2             ", 46, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Aluminum 
  AliMedium(9,  "ALU_C0          ",  9, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(29, "ALU_C1          ", 29, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(49, "ALU_C2          ", 49, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Magnesium
  AliMedium(7,  "MG_C0           ",  7, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Iron 
  AliMedium(10, "FE_C0           ", 10, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(30, "FE_C1           ", 30, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(50, "FE_C2           ", 50, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Copper
  AliMedium(11, "Cu_C0            ", 11, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(31, "Cu_C1            ", 31, 0, isxfld, sxmgmx, tmaxfd, -stemax, deemax, epsil, stmin);
  AliMedium(51, "Cu_C2            ", 51, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Tungsten 
  AliMedium(12, "W_C0            ", 12, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(32, "W_C1            ", 32, 0, isxfld, sxmgmx, tmaxfd, -stemax, deemax, epsil, stmin);
  AliMedium(52, "W_C2            ", 52, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //  
  //    Ni/Tungsten 
  AliMedium(21, "Ni/W0           ", 21, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(41, "Ni/W1           ", 41, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(61, "Ni/W3           ", 61, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Lead 
  AliMedium(13, "PB_C0           ", 13, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(33, "PB_C1           ", 33, 0, isxfld, sxmgmx, tmaxfd, -stemax, deemax, epsil, stmin);
  AliMedium(53, "PB_C2           ", 53, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Insulation Powder 
  AliMedium(14, "INS_C0          ", 14, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(34, "INS_C1          ", 34, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(54, "INS_C2          ", 54, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Air 
  AliMedium(15, "AIR_C0          ", 15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_C1          ", 35, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(55, "AIR_C2          ", 55, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Vacuum 
  AliMedium(16, "VA_C0           ", 16, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(36, "VA_C1           ", 36, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(56, "VA_C2           ", 56, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Concrete 
  AliMedium(17, "CC_C0           ", 17, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(37, "CC_C1           ", 37, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(57, "CC_C2           ", 57, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Polyethilene 
  AliMedium(18, "CH2_C0          ", 18, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(38, "CH2_C1          ", 38, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(58, "CH2_C2          ", 58, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Steel 
  AliMedium(19, "ST_C0           ", 19, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(39, "ST_C1           ", 39, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(59, "ST_C3           ", 59, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  // Polymer Concrete 
  AliMedium(20, "PCc_C0           ", 20, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(40, "PCc_C1           ", 40, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(60, "PCc_C3           ", 60, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
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
 
Int_t  AliABSO::GetMatId(Int_t imat) const 
{
// Get geant material number
    Int_t kmat=(*fIdmate)[imat]; 
    return kmat;
}
