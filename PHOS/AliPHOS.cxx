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

//_________________________________________________________________________
// Base Class of PHOS.
// Only creates the materials
//*-- Author : Laurent Aphecetche  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOS.h"
#include "AliMC.h"
//#include "TGeant3.h"
#include "AliRun.h"

ClassImp(AliPHOS)

//____________________________________________________________________________
AliPHOS::AliPHOS(const char* name, const char* title) 
  : AliDetector(name,title) 
{
}

//____________________________________________________________________________
AliPHOS::AliPHOS() : AliDetector()
{
}

//____________________________________________________________________________
AliPHOS::~AliPHOS()
{
  delete fHits ;
  delete fDigits ;
}

//____________________________________________________________________________
void AliPHOS::CreateMaterials()
{
  // DEFINITION OF PHOS MATERIALS

  // --- The PbWO4 crystals ---
  Float_t AX[3] = {207.19, 183.85, 16.0} ;
  Float_t ZX[3] = {82.0, 74.0, 8.0} ;
  Float_t WX[3] = {1.0, 1.0, 4.0} ;
  Float_t DX = 8.28 ;

  AliMixture(0, "PbWO4$", AX, ZX, DX, -3, WX) ;


  // --- The polysterene scintillator (CH) ---
  Float_t AP[2] = {12.011, 1.00794} ;
  Float_t ZP[2] = {6.0, 1.0} ;
  Float_t WP[2] = {1.0, 1.0} ;
  Float_t DP = 1.032 ;

  AliMixture(1, "Polystyrene$", AP, ZP, DP, -2, WP) ;

  // --- Aluminium ---
  AliMaterial(2, "Al$", 26.98, 13., 2.7, 8.9, 999., 0, 0) ;
  // ---         Absorption length is ignored ^

 // --- Tyvek (CnH2n) ---
  Float_t AT[2] = {12.011, 1.00794} ;
  Float_t ZT[2] = {6.0, 1.0} ;
  Float_t WT[2] = {1.0, 2.0} ;
  Float_t DT = 0.331 ;

  AliMixture(3, "Tyvek$", AT, ZT, DT, -2, WT) ;

  // --- Polystyrene foam ---
  Float_t AF[2] = {12.011, 1.00794} ;
  Float_t ZF[2] = {6.0, 1.0} ;
  Float_t WF[2] = {1.0, 1.0} ;
  Float_t DF = 0.12 ;

  AliMixture(4, "Foam$", AF, ZF, DF, -2, WF) ;

 // --- Titanium ---
  Float_t ATIT[3] = {47.88, 26.98, 54.94} ;
  Float_t ZTIT[3] = {22.0, 13.0, 25.0} ;
  Float_t WTIT[3] = {69.0, 6.0, 1.0} ;
  Float_t DTIT = 4.5 ;

  AliMixture(5, "Titanium$", ATIT, ZTIT, DTIT, -3, WTIT);

 // --- Silicon ---
  AliMaterial(6, "Si$", 28.0855, 14., 2.33, 9.36, 42.3, 0, 0) ;



  // --- Foam thermo insulation ---
  Float_t ATI[2] = {12.011, 1.00794} ;
  Float_t ZTI[2] = {6.0, 1.0} ;
  Float_t WTI[2] = {1.0, 1.0} ;
  Float_t DTI = 0.1 ;

  AliMixture(7, "Thermo Insul.$", ATI, ZTI, DTI, -2, WTI) ;

  // --- Textolitn ---
  Float_t ATX[4] = {16.0, 28.09, 12.011, 1.00794} ;
  Float_t ZTX[4] = {8.0, 14.0, 6.0, 1.0} ;
  Float_t WTX[4] = {292.0, 68.0, 462.0, 736.0} ;
  Float_t DTX    = 1.75 ;

  AliMixture(8, "Textolit$", ATX, ZTX, DTX, -4, WTX) ;

  //--- FR4  ---
  Float_t AFR[3] = {28.0855, 15.9994, 17.749} ; 
  Float_t ZFR[3] = {14., 8., 8.875} ; 
  Float_t WFR[3] = {.28, .32, .4} ;
  Float_t DFR = 1.8 ; 

  AliMixture(9, "FR4$", AFR, ZFR, DFR, -3, WFR) ;

  // --- The Composite Material for  micromegas (so far polyetylene) ---                                       
  Float_t ACM[2] = {12.01, 1.} ; 
  Float_t ZCM[2] = {6., 1.} ; 
  Float_t WCM[2] = {1., 2.} ; 
  Float_t DCM = 0.935 ; 

  AliMixture(10, "Compo Mat$", ACM, ZCM, DCM, -2, WCM) ;

  // --- Copper ---                                                                    
  AliMaterial(11, "Cu$", 63.546, 29, 8.96, 1.43, 14.8, 0, 0) ;
 
  // --- G10 : Printed Circuit material ---                                                  
  Float_t AG10[4] = { 12., 1., 16., 28.} ;
  Float_t ZG10[4] = { 6., 1., 8., 14.} ;
  Float_t WG10[4] = { .259, .288, .248, .205} ;
  Float_t DG10  = 1.7 ;
  
  AliMixture(12, "G10$", AG10, ZG10, DG10, -4, WG10);

  // --- Lead ---                                                                     
  AliMaterial(13, "Pb$", 207.2, 82, 11.35, 0.56, 0., 0, 0) ;

 // --- The gas mixture ---                                                                
 // Co2
  Float_t ACO[2] = {12.0, 16.0} ; 
  Float_t ZCO[2] = {6.0, 8.0} ; 
  Float_t WCO[2] = {1.0, 2.0} ; 
  Float_t DCO = 0.001977 ; 

  AliMixture(14, "CO2$", ACO, ZCO, DCO, -2, WCO);

 // Ar
  Float_t DAr = 0.001782 ; 
  AliMaterial(15, "Ar$", 39.948, 18.0, DAr, 14.0, 0., 0, 0) ;   
 
 // ArCo2
  Char_t namate[21];
  Float_t AGM[2] ; 
  Float_t ZGM[2] ; 
  Float_t WGM[2] ; 
  Float_t DGM ; 

  Float_t AbsL, RadL, Density ;
  Float_t buf[1] ;
  Int_t nbuf ;

  gMC->Gfmate((*fIdmate)[15], namate, AGM[0], ZGM[0], Density, RadL, AbsL, buf, nbuf) ; // Get properties of Ar 
  gMC->Gfmate((*fIdmate)[14], namate, AGM[1], ZGM[1], Density, RadL, AbsL, buf, nbuf) ; // Get properties of CO2 


  // Create gas mixture 

  Float_t ArContent    = 0.80 ;  // Ar-content of the Ar/CO2-mixture (80% / 20%) 
 
  WGM[0] = ArContent;
  WGM[1] = 1. - ArContent ;
  DGM    = WGM[0] * DAr + WGM[1] * DCO;

  AliMixture(16, "ArCO2$", AGM, ZGM, DGM,  2, WGM) ;

 
  // --- Air ---
  AliMaterial(99, "Air$", 14.61, 7.3, 0.001205, 30420., 67500., 0, 0) ;
  
 
  // DEFINITION OF THE TRACKING MEDIA

  // for PHOS: idtmed[699->798] equivalent to fIdtmed[0->100]
  Int_t * idtmed = fIdtmed->GetArray() - 699 ; 
  Int_t   ISXFLD = gAlice->Field()->Integ() ;
  Float_t SXMGMX = gAlice->Field()->Max() ;

  // The scintillator of the calorimeter made of PBW04                              -> idtmed[699]
  AliMedium(0, "PHOS Xtal    $", 0, 1,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[700]
  AliMedium(1, "CPV scint.   $", 1, 1,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Various Aluminium parts made of Al                                             -> idtmed[701]
  AliMedium(2, "Al parts     $", 2, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // The Tywek which wraps the calorimeter crystals                                 -> idtmed[702]
  AliMedium(3, "Tyvek wrapper$", 3, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // The Polystyrene foam around the calorimeter module                             -> idtmed[703]
  AliMedium(4, "Polyst. foam $", 4, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The Titanium around the calorimeter crystal                                    -> idtmed[704]
  AliMedium(5, "Titan. cover $", 5, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.0001, 0.0001, 0, 0) ;

  // The Silicon of the pin diode to read out the calorimeter crystal               -> idtmed[705] 
 AliMedium(6, "Si PIN       $", 6, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.01, 0.01, 0, 0) ;

 // The thermo insulating material of the box which contains the calorimeter module -> idtmed[706]
  AliMedium(7, "Thermo Insul.$", 7, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The Textolit which makes up the box which contains the calorimeter module      -> idtmed[707]
  AliMedium(8, "Textolit     $", 8, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // FR4: The Plastic which makes up the frame of micromegas                        -> idtmed[708]
  AliMedium(9, "FR4 $", 9, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.0001, 0, 0) ; 


  // The Composite Material for  micromegas                                         -> idtmed[709]
  AliMedium(10, "CompoMat   $", 10, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Copper                                                                         -> idtmed[710]
  AliMedium(11, "Copper     $", 11, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.0001, 0, 0) ;

  // G10: Printed Circuit material                                                  -> idtmed[711]
 
  AliMedium(12, "G10        $", 12, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.01, 0, 0) ;

  // The Lead                                                                       -> idtmed[712]
 
  AliMedium(13, "Lead      $", 13, 0,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The gas mixture: ArCo2                                                         -> idtmed[715]
 
  AliMedium(16, "ArCo2      $", 16, 1,
	    ISXFLD, SXMGMX, 10.0, 0.1, 0.1, 0.1, 0.01, 0, 0) ;
 
  // Air                                                                            -> idtmed[798] 
  AliMedium(99, "Air          $", 99, 0,
	    ISXFLD, SXMGMX, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;

  // --- Set decent energy thresholds for gamma and electron tracking

  // Tracking threshold for photons and electrons in the scintillator crystal 
  gMC->Gstpar(idtmed[699], "CUTGAM",0.5E-4) ; 
  gMC->Gstpar(idtmed[699], "CUTELE",1.0E-4) ;
 
  // --- Generate explicitly delta rays in the titan cover ---
  gMC->Gstpar(idtmed[704], "LOSS",3.) ;
  gMC->Gstpar(idtmed[704], "DRAY",1.) ;

  // --- and in aluminium parts ---
  gMC->Gstpar(idtmed[701], "LOSS",3.) ;
  gMC->Gstpar(idtmed[701], "DRAY",1.) ;

// Tracking threshold for photons and electrons in the gas ArC02 
  //  TGeant3 *geant3 = (TGeant3*)gMC;
  //geant3->SetERAN(5.e-8, 1.e1,90);

  gMC->Gstpar(idtmed[715], "CUTGAM",1.E-8) ; 
  gMC->Gstpar(idtmed[715], "CUTELE",1.E-8) ;
  gMC->Gstpar(idtmed[715], "CUTNEU",1.E-8) ;
  gMC->Gstpar(idtmed[715], "CUTHAD",1.E-8) ;
  gMC->Gstpar(idtmed[715], "CUTMUO",1.E-8) ;
  gMC->Gstpar(idtmed[715], "BCUTE",1.E-8) ;
  gMC->Gstpar(idtmed[715], "BCUTM",1.E-8) ;
  gMC->Gstpar(idtmed[715], "DCUTE",1.E-8) ;
  gMC->Gstpar(idtmed[715], "DCUTM",1.E-8) ;
  gMC->Gstpar(idtmed[715], "PPCUTM",1.E-8) ;
  gMC->Gstpar(idtmed[715], "LOSS",2.) ;
  gMC->Gstpar(idtmed[715], "DRAY",0.) ;
  gMC->Gstpar(idtmed[715], "STRA",2.) ;

  
 

}
