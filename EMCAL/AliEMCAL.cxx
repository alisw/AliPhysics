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
//////////////////////////////////////////////////////////////////////////////

// --- Standard library ---
#include <strstream.h>

// --- ROOT system ---
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
  fSDigits       = 0 ;
}
//____________________________________________________________________________
AliEMCAL::AliEMCAL(const char* name, const char* title): AliDetector(name,title) 
{
  // ctor : title is used to identify the layout
  
  // gets an instance of the geometry parameters class  
  
  if (strcmp(GetTitle(),"") != 0 ) 
    fGeom =  AliEMCALGeometry::GetInstance(GetTitle(), "") ; 
}
//____________________________________________________________________________
AliEMCAL::~AliEMCAL()
{
  // dtor

}

//____________________________________________________________________________
void AliEMCAL::CreateMaterials()
{
  // Definitions of materials to build EMCAL and associated tracking media.
  // media number in idtmed are 1599 to 1698.

  // --- Air ---
  AliMaterial(0, "Air$", 14.61, 7.3, 0.001205, 30420., 67500., 0, 0) ;

  // --- Lead ---                                                                     
  AliMaterial(1, "Pb$", 207.2, 82, 11.35, 0.56, 0., 0, 0) ;

  // --- Average properties of the active material ---                                                                     
  AliMaterial(2, "EmcalMat$", fGeom->GetAmat(), 
	                      fGeom->GetZmat(),
	                      fGeom->GetDmat(),
	                      fGeom->GetRmat(), 
	                      0) ;

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

  // The Average properties of the active material                                  -> idtmed[1601]

  AliMedium(2, "EmcalMat  $", 2, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

 // --- Set decent energy thresholds for gamma and electron tracking

  // Tracking threshold for photons and electrons in Lead 
  gMC->Gstpar(idtmed[1600], "CUTGAM",0.5E-4) ; 
  gMC->Gstpar(idtmed[1600], "CUTELE",1.0E-4) ;
  gMC->Gstpar(idtmed[1601], "CUTGAM",0.5E-4) ; 
  gMC->Gstpar(idtmed[1601], "CUTELE",1.0E-4) ;

  // --- Generate explicitly delta rays in Lead ---
  gMC->Gstpar(idtmed[1600], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1600], "DRAY",1.) ;
  gMC->Gstpar(idtmed[1601], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1601], "DRAY",1.) ;
 
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
 
  // Branch address for digit tree
  TTree *treeD = gAlice->TreeD();
  
  if(fDigits)
    fDigits->Clear();

  if (treeD && fDigits) {
    branch = treeD->GetBranch(branchname);
    if (branch) branch->SetAddress(&fDigits);
  }

  if(fSDigits)
    fSDigits->Clear();

  if (gAlice->TreeS()  && fSDigits ) {
    branch = gAlice->TreeS()->GetBranch("EMCAL");
    if (branch) branch->SetAddress(&fSDigits) ;
  } 
}


