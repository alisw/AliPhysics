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
#include <TTree.h>
#include <TVirtualMC.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliMagF.h"
#include "AliEMCAL.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALLoader.h"

ClassImp(AliEMCAL)
//____________________________________________________________________________
AliEMCAL::AliEMCAL():AliDetector()
{
  // Default ctor 
  fName="EMCAL";
  //fQATask = 0;
  fTreeQA = 0;
  fGeom = 0 ; 
}

//____________________________________________________________________________
AliEMCAL::AliEMCAL(const char* name, const char* title): AliDetector(name,title)
{
  //   ctor : title is used to identify the layout
  
  //fQATask = 0;
  fTreeQA = 0;
  fGeom = 0;
}

//____________________________________________________________________________
AliEMCAL::~AliEMCAL()
{

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


  // --- The polysterene scintillator (CH) ---
  Float_t aP[2] = {12.011, 1.00794} ;
  Float_t zP[2] = {6.0, 1.0} ;
  Float_t wP[2] = {1.0, 1.0} ;
  Float_t dP = 1.032 ;

  AliMixture(2, "Polystyrene$", aP, zP, dP, -2, wP) ;

  // --- Aluminium ---
  AliMaterial(3, "Al$", 26.98, 13., 2.7, 8.9, 999., 0, 0) ;
  // ---         Absorption length is ignored ^

  // --- Copper ---
  AliMaterial(4, "Cu$", 63.546, 29, 8.96, 1.43, 14.8, 0, 0) ; 
  // ---         Absorption length is ignored ^


  // DEFINITION OF THE TRACKING MEDIA

  // for EMCAL: idtmed[1599->1698] equivalent to fIdtmed[0->100]
  Int_t * idtmed = fIdtmed->GetArray() - 1599 ; 
  Int_t   isxfld = gAlice->Field()->Integ() ;
  Float_t sxmgmx = gAlice->Field()->Max() ;
 


   // Air                                                                           -> idtmed[1599] 
  AliMedium(0, "Air          $", 0, 0,
	     isxfld, sxmgmx, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;

  // The Lead                                                                      -> idtmed[1600]
 
  AliMedium(1, "Lead      $", 1, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

 
 // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[1601]
  AliMedium(2, "CPV scint.   $", 2, 1,
            isxfld, sxmgmx, 10.0, 0.001, 0.1, 0.001, 0.001, 0, 0) ;

  // Various Aluminium parts made of Al                                            -> idtmed[1602]
  AliMedium(3, "Al parts     $", 3, 0,
             isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // Copper for HCal (post shower)                                                 -> idtmed[1603]
  AliMedium(4, "Copper       $", 4, 0,
             isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;



// --- Set decent energy thresholds for gamma and electron tracking

  // Tracking threshold for photons and electrons in Lead 
  gMC->Gstpar(idtmed[1600],"CUTGAM",0.00008) ;
  gMC->Gstpar(idtmed[1600],"CUTELE",0.001) ;
  gMC->Gstpar(idtmed[1600],"BCUTE",0.0001) ;

  // --- Generate explicitly delta rays in Lead ---
  gMC->Gstpar(idtmed[1600], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1600], "DRAY",1.) ;
  gMC->Gstpar(idtmed[1600], "DCUTE",0.00001) ;
  gMC->Gstpar(idtmed[1600], "DCUTM",0.00001) ;

// --- in aluminium parts ---
  gMC->Gstpar(idtmed[1602], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1602], "DRAY",1.) ;
  gMC->Gstpar(idtmed[1602], "DCUTE",0.00001) ;
  gMC->Gstpar(idtmed[1602], "DCUTM",0.00001) ;

// --- in copper parts ---
  gMC->Gstpar(idtmed[1603], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1603], "DRAY",1.) ;
  gMC->Gstpar(idtmed[1603], "DCUTE",0.00001) ;
  gMC->Gstpar(idtmed[1603], "DCUTM",0.00001) ;



// --- and finally thresholds for photons and electrons in the scintillator ---
  gMC->Gstpar(idtmed[1601],"CUTGAM",0.00008) ;
  gMC->Gstpar(idtmed[1601],"CUTELE",0.001) ;
  gMC->Gstpar(idtmed[1601],"BCUTE",0.0001) ;


}

//____________________________________________________________________________
AliEMCALGeometry * AliEMCAL::GetGeometry() const 
{  
  // gets the pointer to the AliEMCALGeometry unique instance 

  if (fGeom) 
    return fGeom ; 
  else 
    return AliEMCALGeometry::GetInstance(GetTitle(),"") ;  
}

//____________________________________________________________________________
void AliEMCAL::SetTreeAddress()
{ 
  // Linking Hits in Tree to Hits array
  TBranch *branch;
  char branchname[20];
  sprintf(branchname,"%s",GetName());
  
  // Branch address for hit tree
  TTree *treeH = TreeH();
  if (treeH) {
    branch = treeH->GetBranch(branchname);
    if (branch) 
      { 
	if (fHits == 0x0) 
	  fHits= new TClonesArray("AliEMCALHit",1000);
	//Info("SetTreeAddress","<%s> Setting Hits Address",GetName());
	branch->SetAddress(&fHits);
      }
    else
      {
	Warning("SetTreeAddress","<%s> Failed",GetName());
      }
  }
}
//____________________________________________________________________________
void AliEMCAL::WriteQA()
{
  
  // Make TreeQA in the output file. 
  
  if(fTreeQA == 0)
    fTreeQA = new TTree("TreeQA", "QA Alarms") ;    
  // Create Alarms branches
  Int_t bufferSize = 32000 ;    
  Int_t splitlevel = 0 ; 
  
  TFolder* topfold = GetLoader()->GetTopFolder(); //get top aliroot folder; skowron
  TString emcalqafn(AliConfig::Instance()->GetQAFolderName()+"/"); //get name of QAaut folder relative to top event; skowron
  emcalqafn+=GetName(); //hard wired string!!! add the detector name to the pathname; skowron 
  TFolder * alarmsF = (TFolder*)topfold->FindObjectAny(emcalqafn); //get the folder
  
  if (alarmsF == 0x0)
    {
      Error("WriteQA","Can not find folder with qa alarms");
      return;
    }
  TString branchName(alarmsF->GetName());
  TBranch * alarmsBranch = fTreeQA->Branch(branchName,"TFolder", &alarmsF, bufferSize, splitlevel);
  TString branchTitle = branchName + " QA alarms" ; 
  alarmsBranch->SetTitle(branchTitle);
  alarmsBranch->Fill() ; 
  
  //fTreeQA
}

//____________________________________________________________________________
AliLoader* AliEMCAL::MakeLoader(const char* topfoldername)
{
//different behaviour than standard (singleton getter)
// --> to be discussed and made eventually coherent
 fLoader = new AliEMCALLoader(GetName(),topfoldername);
 return fLoader;
}
