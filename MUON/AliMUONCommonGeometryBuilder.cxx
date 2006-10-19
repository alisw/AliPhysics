/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *      SigmaEffect_thetadegrees                                                                  *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
//
// Class AliMUONCommonGeometryBuilder
// ----------------------------------
// Geometry construction common to all stations
// (material definition).
// separated from AliMUONGeometryBuilder


#include <TVirtualMC.h>

#include "AliMUONCommonGeometryBuilder.h"
#include "AliMUON.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONCommonGeometryBuilder)
/// \endcond
 
//______________________________________________________________________________//___________________________________________
AliMUONCommonGeometryBuilder::AliMUONCommonGeometryBuilder(AliMUON* muon)
  : AliMUONVGeometryBuilder(-1, 0),
    fMUON(muon)
{
/// Standard constructor
}

//______________________________________________________________________________//___________________________________________
AliMUONCommonGeometryBuilder::AliMUONCommonGeometryBuilder() 
  : AliMUONVGeometryBuilder(),
    fMUON(0)
{
/// Default constructor
} 

//______________________________________________________________________________
AliMUONCommonGeometryBuilder::~AliMUONCommonGeometryBuilder()
{
/// Destructor
}

//
// public functions
//

//_____________________________________________________________________________
void AliMUONCommonGeometryBuilder::CreateMaterials()
{
/// Definition of common materials

  //
  //     Ar-CO2 gas (80%+20%)
  Float_t ag1[3]   = { 39.95,12.01,16. };
  Float_t zg1[3]   = { 18.,6.,8. };
  Float_t wg1[3]   = { .8,.0667,.13333 };
  Float_t dg1      = .001821;
  //
  //     Ar-buthane-freon gas -- trigger chambers 
  Float_t atr1[4]  = { 39.95,12.01,1.01,19. };
  Float_t ztr1[4]  = { 18.,6.,1.,9. };
  Float_t wtr1[4]  = { .56,.1262857,.2857143,.028 };
  Float_t dtr1     = .002599;
  //
  //     Ar-CO2 gas 
  Float_t agas[3]  = { 39.95,12.01,16. };
  Float_t zgas[3]  = { 18.,6.,8. };
  Float_t wgas[3]  = { .74,.086684,.173316 };
  Float_t dgas     = .0018327;
  //
  //     Ar-Isobutane gas (80%+20%) -- tracking 
  Float_t ag[3]    = { 39.95,12.01,1.01 };
  Float_t zg[3]    = { 18.,6.,1. };
  Float_t wg[3]    = { .8,.057,.143 };
  Float_t dg       = .0019596;
  //
  //     Ar-Isobutane-Forane-SF6 gas (49%+7%+40%+4%) -- trigger 
  Float_t atrig[5] = { 39.95,12.01,1.01,19.,32.066 };
  Float_t ztrig[5] = { 18.,6.,1.,9.,16. };
  Float_t wtrig[5] = { .49,1.08,1.5,1.84,0.04 };
  Float_t dtrig    = .0031463;
  //
  //     bakelite: C6 H6 O
  Float_t abak[3] = {12.01 , 1.01 , 16.};
  Float_t zbak[3] = {6.     , 1.   , 8.};
  Float_t wbak[3] = {6.     , 6.   , 1.}; 
  Float_t dbak = 1.4;

  Int_t iSXFLD   = gAlice->Field()->PrecInteg();
  Float_t sXMGMX = gAlice->Field()->Max();
  //
  // --- Define the various materials for GEANT --- 
  fMUON->AliMaterial(9, "ALUMINIUM0$", 26.98, 13., 2.7, 8.9, 37.2);
  fMUON->AliMaterial(10, "ALUMINIUM1$", 26.98, 13., 2.7, 8.9, 37.2);
  fMUON->AliMaterial(49, "Kapton$", 12.01,6,1.42,-28.6,999);          // from DPG
  fMUON->AliMaterial(42, "Copper$", 63.546,29.,8.96,-1.43,9.6);

  //fMUON->AliMaterial(43, "FR4$", 17.749, 8.875, 1.7, -19.4, 999.);    // from DPG
  Float_t aFR[4] = {16.0, 28.09, 12.011, 1.00794} ;
  Float_t zFR[4] = {8.0, 14.0, 6.0, 1.0} ;
  Float_t wFR[4] = {292.0, 68.0, 462.0, 736.0} ;
  Float_t dFR = 1.8 ; 
  fMUON->AliMixture(43, "FR4$", aFR, zFR, dFR, -4, wFR);

  fMUON->AliMaterial(44, "FrameEpoxy",12.24,6.0,1.85,-19.14,999);// use 16.75cm


  // Air
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  fMUON->AliMixture(15, "AIR$      ", aAir,  zAir, dAir,4, wAir);
  //    fMUON->AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
  fMUON->AliMixture(19, "Bakelite$", abak, zbak, dbak, -3, wbak);
  fMUON->AliMixture(20, "ArC4H10 GAS$", ag, zg, dg, 3, wg);
  fMUON->AliMixture(21, "TRIG GAS$", atrig, ztrig, dtrig, -5, wtrig);
  fMUON->AliMixture(22, "ArCO2 80%$", ag1, zg1, dg1, 3, wg1);
  fMUON->AliMixture(23, "Ar-freon $", atr1, ztr1, dtr1, 4, wtr1);
  fMUON->AliMixture(24, "ArCO2 GAS$", agas, zgas, dgas, 3, wgas);

  // materials for slat: 
  //     Sensitive area: gas (already defined) 
  //     PCB: copper 
  //     insulating material: vetronite -> replacing by G10 Ch. Finck
  //     spacer: noryl Ch. Finck
  //     panel sandwich: carbon, nomex, carbon replacing rohacell by nomex Ch. Finck

  // G10: SiO2(60%) + C8H14O4(40%)
  Float_t aglass[5] = {12.01, 28.09, 16., 1.01,  16.};
  Float_t zglass[5] = { 6.,   14.,    8., 1.,    8.};
  Float_t wglass[5] = { 0.22, 0.28, 0.32, 0.03,  0.15};
  Float_t dglass    = 1.7;

  // rohacell: C9 H13 N1 O2
  Float_t arohac[4] = {12.01,  1.01, 14.010, 16.};
  Float_t zrohac[4] = { 6.,    1.,    7.,     8.};
  Float_t wrohac[4] = { 9.,   13.,    1.,     2.};
  Float_t drohac    = 0.03;

  // Nomex: C22 H10 N2 O5
  Float_t aNomex[4] = {12.01,  1.01, 14.010, 16.};
  Float_t zNomex[4] = { 6.,    1.,    7.,     8.};
  Float_t wNomex[4] = { 22.,   10.,   2.,     5.};
  Float_t dNomex    = 0.024; //honey comb
  Float_t dNomex2   = 1.43;  //bulk material


  // Noryl: C8 H8 O polyphenylene oxyde (di-methyl not sure)
  Float_t aNoryl[3] = {12.01,  1.01, 16.};
  Float_t zNoryl[3] = { 6.,    1.,    8.};
  Float_t wNoryl[3] = { 8.,    8.,    1.};
  Float_t dNoryl    = 1.06;

  fMUON->AliMaterial(31, "COPPER$",   63.54,    29.,   8.96,   1.4, 0.);
  fMUON->AliMixture( 32, "G10$",      aglass, zglass, dglass, -5, wglass);
  fMUON->AliMaterial(33, "Carbon$",   12.01,     6.,  2.265,  18.8, 49.9);
  fMUON->AliMixture( 34, "Rohacell$", arohac, zrohac, drohac, -4, wrohac); 
  fMUON->AliMixture( 35, "Nomex$",    aNomex, zNomex, dNomex, -4, wNomex); 
  fMUON->AliMixture( 36, "Noryl$",    aNoryl, zNoryl, dNoryl, -3, wNoryl); 
  fMUON->AliMixture( 37, "Nomex_bulk$",aNomex, zNomex, dNomex2, -4, wNomex); 

  Float_t  epsil  = .001; // Tracking precision, 
  Float_t  stemax = -1.;  // Maximum displacement for multiple scat 
  Float_t  tmaxfd = -20.; // Maximum angle due to field deflection 
  Float_t  deemax = -.3;  // Maximum fractional energy loss, DLS 
  Float_t  stmin  = -.8;
  Float_t  maxDestepAlu = fMUON->GetMaxDestepAlu();
  Float_t  maxDestepGas = fMUON->GetMaxDestepGas();
  Float_t  maxStepAlu = fMUON->GetMaxStepAlu();
  Float_t  maxStepGas = fMUON->GetMaxStepGas();

  //
  //    Air 
  fMUON->AliMedium(1, "AIR_CH_US         ", 15, 1, iSXFLD, sXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
 
  //
  //    Aluminum 
  fMUON->AliMedium(4, "ALU_CH_US0         ", 9, 0, iSXFLD, sXMGMX, tmaxfd, maxStepAlu, 
		   maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(5, "ALU_CH_US1         ", 10, 0, iSXFLD, sXMGMX, tmaxfd, maxStepAlu, 
		   maxDestepAlu, epsil, stmin);
  //
  //    Ar-isoC4H10 gas 
  fMUON->AliMedium(6, "AR_CH_US          ", 20, 1, iSXFLD, sXMGMX, tmaxfd, maxStepGas, 
		   maxDestepGas, epsil, stmin);
  //
  //    Ar-Isobuthane-Forane-SF6 gas 
  fMUON->AliMedium(7, "GAS_CH_TRIGGER    ", 21, 1, iSXFLD, sXMGMX, tmaxfd, stemax, deemax, epsil, stmin);

  fMUON->AliMedium(8, "BAKE_CH_TRIGGER   ", 19, 0, iSXFLD, sXMGMX, tmaxfd, maxStepAlu, 
		   maxDestepAlu, epsil, stmin);
  //
  // slat medium
  fMUON->AliMedium(9, "ARG_CO2   ", 22, 1, iSXFLD, sXMGMX, tmaxfd, maxStepGas, 
		   maxDestepAlu, epsil, stmin);
  //
  // tracking media for slats: check the parameters!! 
  fMUON->AliMedium(11, "PCB_COPPER        ", 31, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(12, "G10               ", 32, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(13, "CARBON            ", 33, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(14, "Rohacell          ", 34, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(15, "Nomex             ", 35, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(16, "Noryl             ", 36, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(17, "Nomex bulk        ", 37, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);

  // for station 2 only
		   // was med: 4  mat: 9
  fMUON->AliMedium(22, "COPPER_II        ", 42, 0, iSXFLD, sXMGMX, 
                   tmaxfd, maxStepAlu, maxDestepAlu, epsil, stmin);
		   // was med: 10  mat: 30
  fMUON->AliMedium(23, "FR4_CH           ", 43, 0, iSXFLD, sXMGMX, 
                   10.0, 0.01, 0.1, 0.003, 0.003);
  fMUON->AliMedium(24, "FrameCH$",   44, 1, iSXFLD, sXMGMX, 
                   10.0, 0.001, 0.001, 0.001, 0.001);
  fMUON->AliMedium(29, "Kapton            ", 49, 0, iSXFLD, sXMGMX,  
                   10.0, 0.01, 1.0, 0.003, 0.003);
		   // was med: 18  mat: 34 
}


