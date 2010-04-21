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
//_________________________________________________________________________
// Base Class for EMCAL description:
// This class contains material definitions    
// for the EMCAL - It does not place the detector in Alice
//*-- Author: Yves Schutz (SUBATECH) 
//
//*-- Additional Contributions: Sahal Yacoob (LBNL/UCT)
//
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
class TFile;
#include <TFolder.h> 
#include <TGeoGlobalMagField.h>
#include <TGraph.h> 
#include <TH1F.h> 
#include <TRandom.h> 
#include <TTree.h>
#include <TVirtualMC.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliMagF.h"
#include "AliEMCAL.h"
#include "AliRun.h"
#include "AliEMCALLoader.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRawUtils.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

ClassImp(AliEMCAL)

//____________________________________________________________________________
AliEMCAL::AliEMCAL()
  : AliDetector(),
    fBirkC0(0),
    fBirkC1(0.),
    fBirkC2(0.),
    fGeometry(0)
{
  // Default ctor 
  fName = "EMCAL" ;
  InitConstants();

  // Should call  AliEMCALGeometry::GetInstance(EMCAL->GetTitle(),"") for getting EMCAL geometry
}

//____________________________________________________________________________
AliEMCAL::AliEMCAL(const char* name, const char* title)
  : AliDetector(name,title),
    fBirkC0(0),
    fBirkC1(0.),
    fBirkC2(0.),
    fGeometry(0)
{
  //   ctor : title is used to identify the layout
  InitConstants();
}

//____________________________________________________________________________
AliEMCAL::~AliEMCAL()
{
  //dtor
}

//____________________________________________________________________________
void AliEMCAL::InitConstants()
{
  //initialize EMCAL values
  fBirkC0 = 1;
  fBirkC1 = 0.013/1.032;
  fBirkC2 = 9.6e-6/(1.032 * 1.032);
  }

//Not needed, modify $ALICE_ROOT/data/galice.cuts instead.
//Load the modified one in the configuration file with SetTransPar
// //____________________________________________________________________________
// void AliEMCAL::DefineMediumParameters()
// {
//   //
//   // EMCAL cuts (Geant3) 
//   // 
//   Int_t * idtmed = fIdtmed->GetArray() - 1599 ; 
// // --- Set decent energy thresholds for gamma and electron tracking

//   // Tracking threshold for photons and electrons in Lead 
//   Float_t cutgam=10.e-5; // 100 kev;
//   Float_t cutele=10.e-5; // 100 kev;
//   TString ntmp(GetTitle()); 
//   ntmp.ToUpper();
//   if(ntmp.Contains("10KEV")) {
//     cutele = cutgam = 1.e-5;
//   } else if(ntmp.Contains("50KEV")) {
//     cutele = cutgam = 5.e-5;
//   } else if(ntmp.Contains("100KEV")) {
//     cutele = cutgam = 1.e-4;
//   } else if(ntmp.Contains("200KEV")) {
//     cutele = cutgam = 2.e-4;
//   } else if(ntmp.Contains("500KEV")) {
//     cutele = cutgam = 5.e-4;
//   }

//   gMC->Gstpar(idtmed[1600],"CUTGAM", cutgam);
//   gMC->Gstpar(idtmed[1600],"CUTELE", cutele); // 1MEV -> 0.1MEV; 15-aug-05
//   gMC->Gstpar(idtmed[1600],"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
//   gMC->Gstpar(idtmed[1600],"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
//   // --- Generate explicitly delta rays in Lead ---
//   gMC->Gstpar(idtmed[1600], "LOSS", 3) ;
//   gMC->Gstpar(idtmed[1600], "DRAY", 1) ;
//   gMC->Gstpar(idtmed[1600], "DCUTE", cutele) ;
//   gMC->Gstpar(idtmed[1600], "DCUTM", cutele) ;

// // --- in aluminium parts ---
//   gMC->Gstpar(idtmed[1602],"CUTGAM", cutgam) ;
//   gMC->Gstpar(idtmed[1602],"CUTELE", cutele) ;
//   gMC->Gstpar(idtmed[1602],"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
//   gMC->Gstpar(idtmed[1602],"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
//   gMC->Gstpar(idtmed[1602], "LOSS",3.) ;
//   gMC->Gstpar(idtmed[1602], "DRAY",1.) ;
//   gMC->Gstpar(idtmed[1602], "DCUTE", cutele) ;
//   gMC->Gstpar(idtmed[1602], "DCUTM", cutele) ;

// // --- and finally thresholds for photons and electrons in the scintillator ---
//   gMC->Gstpar(idtmed[1601],"CUTGAM", cutgam) ;
//   gMC->Gstpar(idtmed[1601],"CUTELE", cutele) ;// 1MEV -> 0.1MEV; 15-aug-05
//   gMC->Gstpar(idtmed[1601],"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
//   gMC->Gstpar(idtmed[1601],"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
//   gMC->Gstpar(idtmed[1601], "LOSS",3) ; // generate delta rays 
//   gMC->Gstpar(idtmed[1601], "DRAY",1) ;
//   gMC->Gstpar(idtmed[1601], "DCUTE", cutele) ;
//   gMC->Gstpar(idtmed[1601], "DCUTM", cutele) ;

//   // S steel - 
//   gMC->Gstpar(idtmed[1603],"CUTGAM", cutgam);
//   gMC->Gstpar(idtmed[1603],"CUTELE", cutele);
//   gMC->Gstpar(idtmed[1603],"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
//   gMC->Gstpar(idtmed[1603],"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
//   // --- Generate explicitly delta rays 
//   gMC->Gstpar(idtmed[1603], "LOSS",3);
//   gMC->Gstpar(idtmed[1603], "DRAY",1);
//   gMC->Gstpar(idtmed[1603], "DCUTE", cutele) ;
//   gMC->Gstpar(idtmed[1603], "DCUTM", cutele) ;

//   AliEMCALGeometry* geom = GetGeometry();
//   if(geom->GetILOSS()>=0) {
//     for(int i=1600; i<=1603; i++) gMC->Gstpar(idtmed[i], "LOSS", geom->GetILOSS()) ; 
//   } 
//   if(geom->GetIHADR()>=0) {
//     for(int i=1600; i<=1603; i++) gMC->Gstpar(idtmed[i], "HADR", geom->GetIHADR()) ; 
//   }
// }

//____________________________________________________________________________
AliDigitizer* AliEMCAL::CreateDigitizer(AliRunDigitizer* manager) const
{
  //create and return the digitizer
  return new AliEMCALDigitizer(manager);
}

//____________________________________________________________________________
void AliEMCAL::CreateMaterials()
{
  // Definitions of materials to build EMCAL and associated tracking media.
  // media number in idtmed are 1599 to 1698.
  // --- Air ---               
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  AliMixture(0, "Air$", aAir, zAir, dAir, 4, wAir) ;

  // --- Lead ---                                                                     
  AliMaterial(1, "Pb$", 207.2, 82, 11.35, 0.56, 0., 0, 0) ;


  // --- The polysterene scintillator (CH) ---
  Float_t aP[2] = {12.011, 1.00794} ;
  Float_t zP[2] = {6.0, 1.0} ;
  Float_t wP[2] = {1.0, 1.0} ;
  Float_t dP = 1.032 ;

  AliMixture(2, "Polystyrene$", aP, zP, dP, -2, wP) ;

  // --- Aluminium ---
  AliMaterial(3, "Al$", 26.98, 13., 2.7, 8.9, 999., 0, 0) ;
  // ---         Absorption length is ignored ^

  // 25-aug-04 by PAI - see  PMD/AliPMDv0.cxx for STEEL definition
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  AliMixture(4, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);

  // DEFINITION OF THE TRACKING MEDIA

  // for EMCAL: idtmed[1599->1698] equivalent to fIdtmed[0->100]
  Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ() ;
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max() ;

  // Air                                                                         -> idtmed[1599]
 AliMedium(0, "Air$", 0, 0,
	     isxfld, sxmgmx, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;

  // The Lead                                                                      -> idtmed[1600]
 
  AliMedium(1, "Lead$", 1, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

 // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[1601]
  float deemax = 0.1; // maximum fractional energy loss in one step (0 < DEEMAX < deemax )
  AliMedium(2, "Scintillator$", 2, 1,
            isxfld, sxmgmx, 10.0, 0.001, deemax, 0.001, 0.001, 0, 0) ;

  // Various Aluminium parts made of Al                                            -> idtmed[1602]
  AliMedium(3, "Al$", 3, 0,
             isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // 25-aug-04 by PAI : see  PMD/AliPMDv0.cxx for STEEL definition                 -> idtmed[1603]
  AliMedium(4, "S steel$", 4, 0, 
             isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;


  //set constants for Birk's Law implentation
  fBirkC0 =  1;
  fBirkC1 =  0.013/dP;
  fBirkC2 =  9.6e-6/(dP * dP);

}

//____________________________________________________________________________
void  AliEMCAL::Init()
{ 
  // Init
  //Not needed, modify $ALICE_ROOT/data/galice.cuts instead.
  //Load the modified one in the configuration file with SetTransPar
  //DefineMediumParameters(); 
}     

//____________________________________________________________________________
void AliEMCAL::Digits2Raw() {

  static AliEMCALRawUtils rawUtils;
  rawUtils.Digits2Raw();

}
//____________________________________________________________________________
void AliEMCAL::Hits2SDigits()  
{ 
// create summable digits

  GetGeometry();
  AliEMCALSDigitizer emcalDigitizer(fLoader->GetRunLoader()->GetFileName().Data()) ;
  emcalDigitizer.SetEventRange(0, -1) ; // do all the events
  emcalDigitizer.ExecuteTask() ;
}

//____________________________________________________________________________

AliLoader* AliEMCAL::MakeLoader(const char* topfoldername)
{
//different behaviour than standard (singleton getter)
// --> to be discussed and made eventually coherent
 fLoader = new AliEMCALLoader(GetName(),topfoldername);
 return fLoader;
}
