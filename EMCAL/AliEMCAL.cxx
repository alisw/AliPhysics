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
//    
// 
//*-- Author: Yves Schutz (SUBATECH) 
//
//*-- Additional Contributions: Sahal Yacoob (LBNL/UCT)
//
//////////////////////////////////////////////////////////////////////////////

// --- Standard library ---
#include <strstream.h>

// --- ROOT system ---
class TFile;
#include "TBranch.h" 
#include "TClonesArray.h" 
#include "TTree.h" 

// --- AliRoot header files ---

#include "AliEMCAL.h"
#include "AliEMCALGeometry.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliMagF.h"

ClassImp(AliEMCAL)

//____________________________________________________________________________
AliEMCAL::AliEMCAL():AliDetector()
{
  // ctor 
  //We do not create objects, because these pointers will be overwritten during reading from file.
  fGeom     = 0;
}
//____________________________________________________________________________
AliEMCAL::AliEMCAL(const char* name, const char* title): AliDetector(name,title) {
//   ctor : title is used to identify the layout
  fGeom     = 0;
//   gets an instance of the geometry parameters class  
}
//____________________________________________________________________________
AliEMCAL::~AliEMCAL(){
  // dtor
    delete fHits;
}

//____________________________________________________________________________
void AliEMCAL::CreateMaterials(){
  // Definitions of materials to build EMCAL and associated tracking media.
  // media number in idtmed are 1599 to 1698.

  // --- Air ---
  AliMaterial(0, "Air$", 14.61, 7.3, 0.001205, 30420., 67500., 0, 0) ;


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



  // DEFINITION OF THE TRACKING MEDIA

  // for EMCAL: idtmed[1599->1698] equivalent to fIdtmed[0->100]
  Int_t * idtmed = fIdtmed->GetArray() - 1599 ; 
  Int_t   isxfld = gAlice->Field()->Integ() ;
  Float_t sxmgmx = gAlice->Field()->Max() ;
 


   // Air                                                                           -> idtmed[1599] 
  AliMedium(0, "Air          $", 0, 0,
	     isxfld, sxmgmx, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;

  // The Lead                                                                       -> idtmed[1600]
 
  AliMedium(1, "Lead      $", 1, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

 
 // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[1601]
  AliMedium(2, "CPV scint.   $", 2, 1,
            isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

  // Various Aluminium parts made of Al                                             -> idtmed[1602]
  AliMedium(3, "Al parts     $", 3, 0,
             isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;




// --- Set decent energy thresholds for gamma and electron tracking

  // Tracking threshold for photons and electrons in Lead 
  gMC->Gstpar(idtmed[1600], "CUTGAM",0.5E-4) ; 
  gMC->Gstpar(idtmed[1600], "CUTELE",1.0E-4) ;

  // --- Generate explicitly delta rays in Lead ---
  gMC->Gstpar(idtmed[1600], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1600], "DRAY",1.) ;
 
// --- and in aluminium parts ---
  gMC->Gstpar(idtmed[1602], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1602], "DRAY",1.) ;

}
//____________________________________________________________________________
void AliEMCAL::SetTreeAddress()
{ 

  TBranch *branch;
  char branchname[20];
  sprintf(branchname,"%s",GetName());
  
  // Branch address for hit tree
  TTree *treeH = gAlice->TreeH();
  if (treeH && fHits) {
    branch = treeH->GetBranch(branchname);
    if (branch) branch->SetAddress(&fHits);
  }
 
}


