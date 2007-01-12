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
/* History of cvs commits:
 *
 * $Log$
 * Revision 1.104  2006/11/23 13:40:44  hristov
 * Common class for raw data reading and ALTRO mappiing for PHOS and EMCAL (Gustavo, Cvetan)
 *
 * Revision 1.103  2006/11/14 17:11:15  hristov
 * Removing inheritances from TAttLine, TAttMarker and AliRndm in AliModule. The copy constructor and assignment operators are moved to the private part of the class and not implemented. The corresponding changes are propagated to the detectors
 *
 * Revision 1.102  2006/10/27 17:14:27  kharlov
 * Introduce AliDebug and AliLog (B.Polichtchouk)
 *
 * Revision 1.101  2006/10/13 06:47:29  kharlov
 * Simulation of RAW data applies real mapping (B.Polichtchouk)
 *
 * Revision 1.100  2006/08/11 12:36:26  cvetan
 * Update of the PHOS code needed in order to read and reconstruct the beam test raw data (i.e. without an existing galice.root)
 *
 * Revision 1.99  2006/06/28 11:36:09  cvetan
 * New detector numbering scheme (common for DAQ/HLT/Offline). All the subdetectors shall use the AliDAQ class for the sim and rec of the raw data. The AliDAQ and raw reader classes now provide all the necessary interfaces to write and select the detector specific raw-data payload. Look into the AliDAQ.h and AliRawReader.h for more details.
 *
 * Revision 1.98  2006/05/11 11:30:48  cvetan
 * Major changes in AliAltroBuffer. Now it can be used only for writing of raw data. All the corresponding read method are removed. It is based now on AliFstream in order to avoid endianess problems. The altro raw data is written always with little endian
 *
 * Revision 1.97  2006/04/22 10:30:17  hristov
 * Add fEnergy to AliPHOSDigit and operate with EMC amplitude in energy units (Yu.Kharlov)
 *
 * Revision 1.96  2006/04/07 08:41:59  hristov
 * Follow AliAlignObj framework and remove AliPHOSAlignData (Yu.Kharlov)
 *
 * Revision 1.95  2006/03/14 19:40:41  kharlov
 * Remove De-digitizing of raw data and digitizing the raw data fit
 *
 * Revision 1.94  2006/03/07 18:56:25  kharlov
 * CDB is passed via environment variable
 *
 * Revision 1.93  2005/11/22 08:45:11  kharlov
 * Calibration is read from CDB if any (Boris Polichtchouk)
 *
 * Revision 1.92  2005/11/03 13:09:19  hristov
 * Removing meaningless const declarations (linuxicc)
 *
 * Revision 1.91  2005/07/27 15:08:53  kharlov
 * Mixture ArCO2 is corrected
 *
 * Revision 1.90  2005/06/17 07:39:07  hristov
 * Removing GetDebug and SetDebug from AliRun and AliModule. Using AliLog for the messages
 *
 * Revision 1.89  2005/05/28 12:10:07  schutz
 * Copy constructor is corrected (by T.P.)
 *
 */

//_________________________________________________________________________
// Base Class for PHOS description:
//   PHOS consists of a PbWO4 calorimeter (EMCA) and a gazeous charged 
//    particles detector (CPV or PPSD).
//   The only provided method here is CreateMaterials, 
//    which defines the materials common to all PHOS versions.   
// 
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
class TFile;
#include <TFolder.h> 
#include <TTree.h>
#include <TVirtualMC.h> 
#include <TH1F.h> 
#include <TF1.h> 
#include <TRandom.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliMagF.h"
#include "AliPHOS.h"
#include "AliPHOSGetter.h"
#include "AliRun.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSDigit.h"
#include "AliAltroBuffer.h"
#include "AliAltroMapping.h"
#include "AliCaloAltroMapping.h"
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliPHOSCalibData.h"
#include "AliDAQ.h"

ClassImp(AliPHOS)

Double_t AliPHOS::fgCapa        = 1.;        // 1pF 
Int_t    AliPHOS::fgOrder       = 2 ;
Double_t AliPHOS::fgTimeMax     = 2.56E-5 ;  // each sample is over 100 ns fTimeMax/fTimeBins
Double_t AliPHOS::fgTimePeak    = 4.1E-6 ;   // 4 micro seconds
Double_t AliPHOS::fgTimeTrigger = 100E-9 ;      // 100ns, just for a reference

Double_t AliPHOS::fgHighCharge  = 8.2;       // adjusted for a high gain range of 5.12 GeV (10 bits)
Double_t AliPHOS::fgHighGain    = 6.64;
Double_t AliPHOS::fgHighLowGainFactor = 16.; // adjusted for a low gain range of 82 GeV (10 bits) 

//____________________________________________________________________________
  AliPHOS:: AliPHOS() : AliDetector()
{
  // Default ctor
  fName   = "PHOS" ;

}

//____________________________________________________________________________
AliPHOS::AliPHOS(const char* name, const char* title): AliDetector(name, title)
{
  //   ctor : title is used to identify the layout
}

//____________________________________________________________________________
AliPHOS::~AliPHOS() 
{  
}

//____________________________________________________________________________
AliDigitizer* AliPHOS::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliPHOSDigitizer(manager);
}

//____________________________________________________________________________
void AliPHOS::CreateMaterials()
{
  // Definitions of materials to build PHOS and associated tracking media.
  // media number in idtmed are 699 to 798.

  // --- The PbWO4 crystals ---
  Float_t aX[3] = {207.19, 183.85, 16.0} ;
  Float_t zX[3] = {82.0, 74.0, 8.0} ;
  Float_t wX[3] = {1.0, 1.0, 4.0} ;
  Float_t dX = 8.28 ;

  AliMixture(0, "PbWO4$", aX, zX, dX, -3, wX) ;


  // --- The polysterene scintillator (CH) ---
  Float_t aP[2] = {12.011, 1.00794} ;
  Float_t zP[2] = {6.0, 1.0} ;
  Float_t wP[2] = {1.0, 1.0} ;
  Float_t dP = 1.032 ;

  AliMixture(1, "Polystyrene$", aP, zP, dP, -2, wP) ;

  // --- Aluminium ---
  AliMaterial(2, "Al$", 26.98, 13., 2.7, 8.9, 999., 0, 0) ;
  // ---         Absorption length is ignored ^

 // --- Tyvek (CnH2n) ---
  Float_t aT[2] = {12.011, 1.00794} ;
  Float_t zT[2] = {6.0, 1.0} ;
  Float_t wT[2] = {1.0, 2.0} ;
  Float_t dT = 0.331 ;

  AliMixture(3, "Tyvek$", aT, zT, dT, -2, wT) ;

  // --- Polystyrene foam ---
  Float_t aF[2] = {12.011, 1.00794} ;
  Float_t zF[2] = {6.0, 1.0} ;
  Float_t wF[2] = {1.0, 1.0} ;
  Float_t dF = 0.12 ;

  AliMixture(4, "Foam$", aF, zF, dF, -2, wF) ;

 // --- Titanium ---
  Float_t aTIT[3] = {47.88, 26.98, 54.94} ;
  Float_t zTIT[3] = {22.0, 13.0, 25.0} ;
  Float_t wTIT[3] = {69.0, 6.0, 1.0} ;
  Float_t dTIT = 4.5 ;

  AliMixture(5, "Titanium$", aTIT, zTIT, dTIT, -3, wTIT);

 // --- Silicon ---
  AliMaterial(6, "Si$", 28.0855, 14., 2.33, 9.36, 42.3, 0, 0) ;



  // --- Foam thermo insulation ---
  Float_t aTI[2] = {12.011, 1.00794} ;
  Float_t zTI[2] = {6.0, 1.0} ;
  Float_t wTI[2] = {1.0, 1.0} ;
  Float_t dTI = 0.04 ;

  AliMixture(7, "Thermo Insul.$", aTI, zTI, dTI, -2, wTI) ;

  // --- Textolith ---
  Float_t aTX[4] = {16.0, 28.09, 12.011, 1.00794} ;
  Float_t zTX[4] = {8.0, 14.0, 6.0, 1.0} ;
  Float_t wTX[4] = {292.0, 68.0, 462.0, 736.0} ;
  Float_t dTX    = 1.75 ;

  AliMixture(8, "Textolit$", aTX, zTX, dTX, -4, wTX) ;

  //--- FR4  ---
  Float_t aFR[4] = {16.0, 28.09, 12.011, 1.00794} ;
  Float_t zFR[4] = {8.0, 14.0, 6.0, 1.0} ;
  Float_t wFR[4] = {292.0, 68.0, 462.0, 736.0} ;
  Float_t dFR = 1.8 ; 

  AliMixture(9, "FR4$", aFR, zFR, dFR, -4, wFR) ;

  // --- The Composite Material for  micromegas (so far polyetylene) ---                                       
  Float_t aCM[2] = {12.01, 1.} ; 
  Float_t zCM[2] = {6., 1.} ; 
  Float_t wCM[2] = {1., 2.} ; 
  Float_t dCM = 0.935 ; 

  AliMixture(10, "Compo Mat$", aCM, zCM, dCM, -2, wCM) ;

  // --- Copper ---                                                                    
  AliMaterial(11, "Cu$", 63.546, 29, 8.96, 1.43, 14.8, 0, 0) ;
 
  // --- G10 : Printed Circuit material ---                                                  
  Float_t aG10[4] = { 12., 1., 16., 28.} ;
  Float_t zG10[4] = { 6., 1., 8., 14.} ;
  Float_t wG10[4] = { .259, .288, .248, .205} ;
  Float_t dG10  = 1.7 ;
  
  AliMixture(12, "G10$", aG10, zG10, dG10, -4, wG10);

  // --- Lead ---                                                                     
  AliMaterial(13, "Pb$", 207.2, 82, 11.35, 0.56, 0., 0, 0) ;

 // --- The gas mixture ---                                                                
 // Co2
  Float_t aCO[2] = {12.0, 16.0} ; 
  Float_t zCO[2] = {6.0, 8.0} ; 
  Float_t wCO[2] = {1.0, 2.0} ; 
  Float_t dCO = 0.001977 ; 

  AliMixture(14, "CO2$", aCO, zCO, dCO, -2, wCO);

 // Ar
  Float_t dAr = 0.001782 ; 
  AliMaterial(15, "Ar$", 39.948, 18.0, dAr, 14.0, 0., 0, 0) ;   
 
  // Ar+CO2 Mixture (80% / 20%)
  Float_t arContent = 0.80 ;  // Ar-content of the ArCO2-mixture
  Float_t aArCO[3]  = {39.948, 12.0, 16.0} ;
  Float_t zArCO[3]  = {18.0  ,  6.0,  8.0} ;
  Float_t wArCO[3];
  wArCO[0] = arContent;
  wArCO[1] = (1-arContent)*1;
  wArCO[2] = (1-arContent)*2;
  Float_t dArCO = arContent*dAr + (1-arContent)*dCO ;
  AliMixture(16, "ArCO2$", aArCO, zArCO, dArCO,  -3, wArCO) ;

  // --- Stainless steel (let it be pure iron) ---
  AliMaterial(17, "Steel$", 55.845, 26, 7.87, 1.76, 0., 0, 0) ;


  // --- Fiberglass ---
  Float_t aFG[4] = {16.0, 28.09, 12.011, 1.00794} ;
  Float_t zFG[4] = {8.0, 14.0, 6.0, 1.0} ;
  Float_t wFG[4] = {292.0, 68.0, 462.0, 736.0} ;
  Float_t dFG    = 1.9 ;

  AliMixture(18, "Fibergla$", aFG, zFG, dFG, -4, wFG) ;

  // --- Cables in Air box  ---
  // SERVICES

  Float_t aCA[4] = { 1.,12.,55.8,63.5 };
  Float_t zCA[4] = { 1.,6.,26.,29. }; 
  Float_t wCA[4] = { .014,.086,.42,.48 };
  Float_t dCA    = 0.8 ;  //this density is raw estimation, if you know better - correct

  AliMixture(19, "Cables  $", aCA, zCA, dCA, -4, wCA) ;


  // --- Air ---
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
 
  AliMixture(99, "Air$", aAir, zAir, dAir, 4, wAir) ;

  // DEFINITION OF THE TRACKING MEDIA

  // for PHOS: idtmed[699->798] equivalent to fIdtmed[0->100]
  Int_t * idtmed = fIdtmed->GetArray() - 699 ; 
  Int_t   isxfld = gAlice->Field()->Integ() ;
  Float_t sxmgmx = gAlice->Field()->Max() ;

  // The scintillator of the calorimeter made of PBW04                              -> idtmed[699]
  AliMedium(0, "PHOS Xtal    $", 0, 1,
	    isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[700]
  AliMedium(1, "CPV scint.   $", 1, 1,
	    isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Various Aluminium parts made of Al                                             -> idtmed[701]
  AliMedium(2, "Al parts     $", 2, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // The Tywek which wraps the calorimeter crystals                                 -> idtmed[702]
  AliMedium(3, "Tyvek wrapper$", 3, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // The Polystyrene foam around the calorimeter module                             -> idtmed[703]
  AliMedium(4, "Polyst. foam $", 4, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The Titanium around the calorimeter crystal                                    -> idtmed[704]
  AliMedium(5, "Titan. cover $", 5, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.0001, 0.0001, 0, 0) ;

  // The Silicon of the pin diode to read out the calorimeter crystal               -> idtmed[705] 
 AliMedium(6, "Si PIN       $", 6, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.01, 0.01, 0, 0) ;

 // The thermo insulating material of the box which contains the calorimeter module -> idtmed[706]
  AliMedium(7, "Thermo Insul.$", 7, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The Textolit which makes up the box which contains the calorimeter module      -> idtmed[707]
  AliMedium(8, "Textolit     $", 8, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // FR4: The Plastic which makes up the frame of micromegas                        -> idtmed[708]
  AliMedium(9, "FR4 $", 9, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, 0, 0) ; 


  // The Composite Material for  micromegas                                         -> idtmed[709]
  AliMedium(10, "CompoMat   $", 10, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Copper                                                                         -> idtmed[710]
  AliMedium(11, "Copper     $", 11, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, 0, 0) ;

  // G10: Printed Circuit material                                                  -> idtmed[711]
 
  AliMedium(12, "G10        $", 12, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.01, 0, 0) ;

  // The Lead                                                                       -> idtmed[712]
 
  AliMedium(13, "Lead      $", 13, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // The gas mixture: ArCo2                                                         -> idtmed[715]
 
  AliMedium(16, "ArCo2      $", 16, 1,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.01, 0, 0) ;
 
  // Stainless steel                                                                -> idtmed[716]
  AliMedium(17, "Steel     $", 17, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.0001, 0, 0) ;

  // Fibergalss                                                                     -> idtmed[717]
  AliMedium(18, "Fiberglass$", 18, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Cables in air                                                                  -> idtmed[718]
  AliMedium(19, "Cables    $", 19, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Air                                                                            -> idtmed[798] 
  AliMedium(99, "Air          $", 99, 0,
	     isxfld, sxmgmx, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;

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
  // --- and in PIN diode
  gMC->Gstpar(idtmed[705], "LOSS",3) ;
  gMC->Gstpar(idtmed[705], "DRAY",1) ;
  // --- and in the passive convertor
  gMC->Gstpar(idtmed[712], "LOSS",3) ;
  gMC->Gstpar(idtmed[712], "DRAY",1) ;
  // Tracking threshold for photons and electrons in the gas ArC02 
  gMC->Gstpar(idtmed[715], "CUTGAM",1.E-5) ; 
  gMC->Gstpar(idtmed[715], "CUTELE",1.E-5) ;
  gMC->Gstpar(idtmed[715], "CUTNEU",1.E-5) ;
  gMC->Gstpar(idtmed[715], "CUTHAD",1.E-5) ;
  gMC->Gstpar(idtmed[715], "CUTMUO",1.E-5) ;
  gMC->Gstpar(idtmed[715], "BCUTE",1.E-5) ;
  gMC->Gstpar(idtmed[715], "BCUTM",1.E-5) ;
  gMC->Gstpar(idtmed[715], "DCUTE",1.E-5) ;
  gMC->Gstpar(idtmed[715], "DCUTM",1.E-5) ;
  gMC->Gstpar(idtmed[715], "PPCUTM",1.E-5) ;
  gMC->Gstpar(idtmed[715], "LOSS",2.) ;
  gMC->Gstpar(idtmed[715], "DRAY",0.) ;
  gMC->Gstpar(idtmed[715], "STRA",2.) ;

}

//____________________________________________________________________________
void AliPHOS::Digits2Raw()
{
// convert digits of the current event to raw data
  
  AliPHOSLoader * loader = dynamic_cast<AliPHOSLoader*>(fLoader) ; 

  // get the digits
  loader->LoadDigits();
  TClonesArray* digits = loader->Digits() ;

  if (!digits) {
    AliError(Form("No digits found !"));
    return;
  }

  // get the geometry
  AliPHOSGeometry* geom = GetGeometry();
  if (!geom) {
    AliError(Form("No geometry found !"));
    return;
  }

  // some digitization constants
//   const Int_t    kThreshold = 1; // skip digits below this threshold // YVK
  const Float_t    kThreshold = 0.001; // skip digits below 1 MeV
  const Int_t      kAdcThreshold = 1;  // Lower ADC threshold to write to raw data

  AliAltroBuffer* buffer = NULL;
  Int_t prevDDL = -1;
  Int_t adcValuesLow[fkTimeBins];
  Int_t adcValuesHigh[fkTimeBins];


  //!!!!for debug!!!
  Int_t modMax=-111;
  Int_t colMax=-111;
  Int_t rowMax=-111;
  Float_t eMax=-333;
  //!!!for debug!!!

  // loop over digits (assume ordered digits)
  for (Int_t iDigit = 0; iDigit < digits->GetEntries(); iDigit++) {
    AliPHOSDigit* digit = dynamic_cast<AliPHOSDigit *>(digits->At(iDigit)) ;
    if (digit->GetEnergy() < kThreshold) 
      continue;
    Int_t relId[4];
    geom->AbsToRelNumbering(digit->GetId(), relId);
    Int_t module = relId[0];
 
   // Begin FIXME 
    if (relId[1] != 0) 
      continue;    // ignore digits from CPV
   // End FIXME 

    Int_t row = relId[2]-1;
    Int_t col = relId[3]-1;

    Int_t iRCU = -111;

    //RCU0
    if(0<=row&&row<32 && 0<=col&&col<28) iRCU=0;

    //RCU1
    if(0<=row&&row<32 && 28<=col&&col<56) iRCU=1;

    //RCU2
    if(32<=row&&row<64 && 0<=col&&col<28) iRCU=2;

    //RCU3
    if(32<=row&&row<64 && 28<=col&&col<56) iRCU=3;


    // PHOS EMCA has 4 DDL per module. Splitting is based on the (row,column) numbers.
    // PHOS internal convention: 1<module<5.
    Int_t iDDL = 4 * (module - 1) + iRCU;

    // new DDL
    if (iDDL != prevDDL) {
      // write real header and close previous file
      if (buffer) {
	buffer->Flush();
	buffer->WriteDataHeader(kFALSE, kFALSE);
	delete buffer;
      }

      // open new file and write dummy header
      TString fileName = AliDAQ::DdlFileName("PHOS",iDDL);

      TString path = gSystem->Getenv("ALICE_ROOT");
      path += "/PHOS/mapping/RCU";
      path += iRCU;
      path += ".data";

      AliAltroMapping* mapping = new AliCaloAltroMapping(path.Data());
      buffer = new AliAltroBuffer(fileName.Data(),mapping);
      buffer->WriteDataHeader(kTRUE, kFALSE);  //Dummy;

      prevDDL = iDDL;
    }

    // out of time range signal (?)
    if (digit->GetTimeR() > GetRawFormatTimeMax() ) {
      AliInfo("Signal is out of time range.\n");
      buffer->FillBuffer((Int_t)digit->GetEnergy());
      buffer->FillBuffer(GetRawFormatTimeBins() );  // time bin
      buffer->FillBuffer(3);          // bunch length      
      buffer->WriteTrailer(3, relId[3], relId[2], module);  // trailer
      
    // calculate the time response function
    } else {
      Double_t energy = 0 ;
      Int_t   module = relId[0];
      if ( digit->GetId() <= geom->GetNModules() *  geom->GetNCristalsInModule()) {
	energy=digit->GetEnergy();
	AliDebug(2,Form("digit energy: %f\n",digit->GetEnergy()));
	if(energy>eMax) {eMax=energy; modMax=module; colMax=col; rowMax=row;}
      }
      else {
 	energy = 0; // CPV raw data format is now know yet
      }        
      Bool_t lowgain = RawSampledResponse(digit->GetTimeR(), energy, adcValuesHigh, adcValuesLow) ; 
      
	buffer->WriteChannel(relId[3]-1, relId[2]-1, 0, 
			     GetRawFormatTimeBins(), adcValuesLow , kAdcThreshold);
 	buffer->WriteChannel(relId[3]-1, relId[2]-1, 1, 
 			     GetRawFormatTimeBins(), adcValuesHigh, kAdcThreshold);
      
    }
  }
  
  // write real header and close last file
  if (buffer) {
    buffer->Flush();
    buffer->WriteDataHeader(kFALSE, kFALSE);
    delete buffer;
  }
  
  AliDebug(1,Form("Digit with max. energy:  modMax %d colMax %d rowMax %d  eMax %f\n",
	 modMax,colMax,rowMax,eMax));

  loader->UnloadDigits();
}

//____________________________________________________________________________
void AliPHOS::Hits2SDigits()  
{ 
// create summable digits

  AliPHOSSDigitizer phosDigitizer(fLoader->GetRunLoader()->GetFileName().Data()) ;
  phosDigitizer.SetEventRange(0, -1) ; // do all the events
  phosDigitizer.ExecuteTask("all") ; 
}

//____________________________________________________________________________
AliLoader* AliPHOS::MakeLoader(const char* topfoldername)
{
//different behaviour than standard (singleton getter)
// --> to be discussed and made eventually coherent
 fLoader = new AliPHOSLoader(GetName(),topfoldername);
 return fLoader;
}

//__________________________________________________________________
Double_t AliPHOS::RawResponseFunction(Double_t *x, Double_t *par) 
{
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function of the general form
  // v(t) = n**n * Q * A**n / C *(t/tp)**n * exp(-n * t/tp) with 
  // tp : peaking time par[0]
  // n  : order of the function
  // C  : integrating capacitor in the preamplifier
  // A  : open loop gain of the preamplifier
  // Q  : the total APD charge to be measured Q = C * energy
  
  Double_t signal ;
  Double_t xx = x[0] - ( fgTimeTrigger + par[3] ) ; 

  if (xx < 0 || xx > fgTimeMax) 
    signal = 0. ;  
  else { 
    Double_t fac = par[0] * TMath::Power(fgOrder, fgOrder) * TMath::Power(par[1], fgOrder) / fgCapa ; 
    signal = fac * par[2] * TMath::Power(xx / fgTimePeak, fgOrder) * TMath::Exp(-fgOrder * (xx / fgTimePeak)) ; 
  }
  return signal ;  
}

//__________________________________________________________________
Double_t AliPHOS::RawResponseFunctionMax(Double_t charge, Double_t gain) 
{
  return ( charge * TMath::Power(fgOrder, fgOrder) * TMath::Power(gain, fgOrder) 
     / ( fgCapa * TMath::Exp(fgOrder) ) );  

}

//__________________________________________________________________
Bool_t AliPHOS::RawSampledResponse(Double_t dtime, Double_t damp, Int_t * adcH, Int_t * adcL) const 
{
  // for a start time dtime and an amplitude damp given by digit, 
  // calculates the raw sampled response AliPHOS::RawResponseFunction
  // Input: dtime - signal start time
  //        damp  - signal amplitude (energy)
  // Output: adcH - array[fkTimeBins] of 10-bit samples for high-gain channel
  //         adcL - array[fkTimeBins] of 10-bit samples for low-gain channel

  const Int_t kRawSignalOverflow = 0x3FF ; 
  Bool_t lowGain = kFALSE ; 

  TF1 signalF("signal", RawResponseFunction, 0, GetRawFormatTimeMax(), 4);

  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    signalF.SetParameter(0, GetRawFormatHighCharge() ) ; 
    signalF.SetParameter(1, GetRawFormatHighGain() ) ; 
    signalF.SetParameter(2, damp) ; 
    signalF.SetParameter(3, dtime) ; 
    Double_t time = iTime * GetRawFormatTimeMax() / GetRawFormatTimeBins() ;
    Double_t signal = signalF.Eval(time) ;     
    if ( static_cast<Int_t>(signal+0.5) > kRawSignalOverflow ){  // larger than 10 bits 
      signal = kRawSignalOverflow ;
      lowGain = kTRUE ; 
    }
    adcH[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    AliDebug(4,Form("iTime: %d Energy: %f HG signal: %f adcH: %d ",iTime,damp,signal,adcH[iTime]));

    signalF.SetParameter(0, GetRawFormatLowCharge() ) ;     
    signalF.SetParameter(1, GetRawFormatLowGain() ) ; 
    signal = signalF.Eval(time) ;  
    if ( static_cast<Int_t>(signal+0.5) > kRawSignalOverflow)  // larger than 10 bits 
      signal = kRawSignalOverflow ;
    adcL[iTime] = static_cast<Int_t>(0.5 + signal ) ; 
    AliDebug(4,Form("..LG: %f adcL: %d\n",signal,adcL[iTime]));

  }
  return lowGain ; 
}

//____________________________________________________________________________
void AliPHOS::SetTreeAddress()
{ 
  // Links Hits in the Tree to Hits array
  TBranch *branch;
  char branchname[20];
  sprintf(branchname,"%s",GetName());
  // Branch address for hit tree
    TTree *treeH = TreeH();
  if (treeH) {
    branch = treeH->GetBranch(branchname);
    if (branch) 
     { 
       if (fHits == 0x0) fHits= new TClonesArray("AliPHOSHit",1000);
       //AliInfo(Form("<%s> Setting Hits Address",GetName()));
       branch->SetAddress(&fHits);
     }
  }
}

