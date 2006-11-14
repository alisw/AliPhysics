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
//  Muon Shield Class                                                        //
//  This class contains a description of the muon shield                     //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliSHILClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliSHIL.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h"
#include "AliLog.h"

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
  : AliModule(name,title)
{
  //
  // Standard constructor for muon shield
  //
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliSHIL::CreateGeometry()
{
  //
  // Build muon shield geometry
  //
}

//_____________________________________________________________________________
void AliSHIL::CreateMaterials()
{
  //
  // Defines materials for the muon shield
  //
  Int_t   isxfld1 = gAlice->Field()->Integ();
  Int_t   isxfld2 = gAlice->Field()->PrecInteg();

  Float_t sxmgmx = gAlice->Field()->Max();
// Steel  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
// PbW
  Float_t apbw[2]   = { 207.2,183.85 };
  Float_t zpbw[2]   = { 82.,74. };
  Float_t wpbw[2]   = { .5,.5 };
// Concrete
  Float_t aconc[10] = { 1.,12.01,15.994,22.99,24.305,26.98,
			28.086,39.1,40.08,55.85 };
  Float_t zconc[10] = { 1.,6.,8.,11.,12.,13.,14.,19.,20.,26. };
  Float_t wconc[10] = { .01,.001,.529107,.016,.002,.033872,
			.337021,.013,.044,.014 };
// Ni-Cu-W alloy
  Float_t aniwcu[3] ={58.6934, 183.84, 63.546};
  Float_t zniwcu[3] ={28.,      74.,   29.};
  Float_t wniwcu[3] ={ 0.015,    0.95,  0.035};
//
// Insulation powder
//                    Si         O       Ti     Al
  Float_t ains[4] ={28.0855, 15.9994, 47.867,  26.982};
  Float_t zins[4] ={14.,      8.    , 22.   ,  13.   };
  Float_t wins[4] ={ 0.3019,  0.4887,  0.1914,  0.018};
//
// Air
//
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  Float_t dAir1 = 1.20479E-10;

  
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  
  //     STEEL 
  
  
  //     LEAD/TUNGSTEN MIXTURE 
  
  
  // --- Define the various materials for GEANT --- 
  //     Aluminum 
  AliMaterial(9,  "ALU1      ", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(29, "ALU2      ", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(49, "ALU3      ", 26.98, 13., 2.7, 8.9, 37.2);
  
  //     Iron 
  AliMaterial(10, "IRON1     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(30, "IRON2     ", 55.85, 26., 7.87, 1.76, 17.1);
  AliMaterial(50, "IRON3     ", 55.85, 26., 7.87, 1.76, 17.1);

  //
  //     Copper
  AliMaterial(11, "COPPER1   ", 63.55, 29., 8.96, 1.43, 15.1);
  AliMaterial(31, "COPPER2   ", 63.55, 29., 8.96, 1.43, 15.1);
  AliMaterial(51, "COPPER3   ", 63.55, 29., 8.96, 1.43, 15.1);
  
  //     Tungsten 
  AliMaterial(12, "TUNGSTEN1 ", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(32, "TUNGSTEN2 ", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(52, "TUNGSTEN3 ", 183.85, 74., 19.3, .35, 10.3);
  
  //     Lead 
  AliMaterial(13, "LEAD1     ", 207.19, 82., 11.35, .56, 18.5);
  AliMaterial(33, "LEAD2     ", 207.19, 82., 11.35, .56, 18.5);
  AliMaterial(53, "LEAD3     ", 207.19, 82., 11.35, .56, 18.5);
  
  //     Air 
  AliMixture(15, "AIR1      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(35, "AIR2      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(55, "AIR3      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(75, "AIR_MUON  ", aAir, zAir, dAir, 4, wAir);

  //     Vacuum 
  AliMixture(16, "VACUUM1 ", aAir, zAir, dAir1, 4, wAir);
  AliMixture(36, "VACUUM2 ", aAir, zAir, dAir1, 4, wAir);
  AliMixture(56, "VACUUM3 ", aAir, zAir, dAir1, 4, wAir);
  
  //     Stainless Steel 
  AliMixture(19, "STAINLESS STEEL1", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(39, "STAINLESS STEEL2", asteel, zsteel, 7.88, 4, wsteel);
  AliMixture(59, "STAINLESS STEEL3", asteel, zsteel, 7.88, 4, wsteel);
  
  //     Lead/Tungsten 
  AliMixture(20, "LEAD/TUNGSTEN1", apbw, zpbw, 15.325, 2, wpbw);
  AliMixture(40, "LEAD/TUNGSTEN2", apbw, zpbw, 15.325, 2, wpbw);
  AliMixture(60, "LEAD/TUNGSTEN3", apbw, zpbw, 15.325, 2, wpbw);

  //     Ni-W-Cu 
  AliMixture(21, "Ni-W-Cu1", aniwcu, zniwcu, 18.78, 3, wniwcu);
  AliMixture(41, "Ni-W-Cu2", aniwcu, zniwcu, 18.78, 3, wniwcu);
  AliMixture(61, "Ni-W-Cu3", aniwcu, zniwcu, 18.78, 3, wniwcu);

  //     Concrete 
  AliMixture(17, "CONCRETE1", aconc, zconc, 2.35, 10, wconc);
  AliMixture(37, "CONCRETE2", aconc, zconc, 2.35, 10, wconc);
  AliMixture(57, "CONCRETE3", aconc, zconc, 2.35, 10, wconc);

  //     Insulation powder 
  AliMixture(14, "INSULATION1", ains, zins, 0.41, 4, wins);
  AliMixture(34, "INSULATION2", ains, zins, 0.41, 4, wins);
  AliMixture(54, "INSULATION3", ains, zins, 0.41, 4, wins);
  
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
  AliMedium(9,  "ALU_C0          ", 9, 0,  isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(29, "ALU_C1          ", 29, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(49, "ALU_C2          ", 49, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Iron 
  AliMedium(10, "FE_C0           ", 10, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(30, "FE_C1           ", 30, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(50, "FE_C2           ", 50, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //    Copper 
  AliMedium(11, "Cu_C0           ", 11, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(31, "Cu_C1           ", 31, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(51, "Cu_C2           ", 51, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Tungsten 
  AliMedium(12, "W_C0            ", 12, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(32, "W_C1            ", 32, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(52, "W_C2            ", 52, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Lead 
  AliMedium(13, "PB_C0           ", 13, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(33, "PB_C1           ", 33, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(53, "PB_C2           ", 53, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //    Insulation Powder 
  AliMedium(14, "INS_C0          ", 14, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(34, "INS_C1          ", 34, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(54, "INS_C2          ", 54, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Air 
  AliMedium(15, "AIR_C0          ", 15, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_C1          ", 35, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(55, "AIR_C2          ", 55, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(75, "AIR_MUON        ", 75, 0, isxfld2, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Vacuum 
  AliMedium(16, "VA_C0           ", 16, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(36, "VA_C1           ", 36, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(56, "VA_C2           ", 56, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Steel 
  AliMedium(19, "ST_C0           ", 19, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(39, "ST_C1           ", 39, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(59, "ST_C3           ", 59, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Lead/Tungsten 
  AliMedium(20, "PB/W0           ", 20, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(40, "PB/W1           ", 40, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(60, "PB/W3           ", 60, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //    Ni/Tungsten 
  AliMedium(21, "Ni/W0           ", 21, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(41, "Ni/W1           ", 41, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(61, "Ni/W3           ", 61, 0, isxfld1, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

//    Concrete 
  AliMedium(17, "CC_C0           ", 17, 0, 0,      sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(37, "CC_C1           ", 37, 0, 0,      sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(57, "CC_C2           ", 57, 0, 0,      sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

}

//_____________________________________________________________________________
void AliSHIL::DrawModule () const
{
    // Drawing options
}

//_____________________________________________________________________________
void AliSHIL::Init()
{
  //
  // Initialise the muon shield after it has been built
  //
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" SHIL_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    // Here the SHIL initialisation code (if any!)
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}

 
