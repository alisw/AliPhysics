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
// Implementation version v0 of EMCAL Manager class 
// An object of this class does not produce hits nor digits
// It is the one to use if you do not want to produce outputs in TREEH or TREED
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

#include "TTUBS.h"
#include "TNode.h"
#include "TRandom.h"
#include "TGeometry.h"


// --- Standard library ---

#include <iostream.h> 

// --- AliRoot header files ---

#include "AliEMCALv0.h"
#include "AliEMCALGeometry.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliEMCALv0)

//____________________________________________________________________________
AliEMCALv0::AliEMCALv0(const char *name, const char *title):
  AliEMCAL(name,title)
{
  // ctor

}

//____________________________________________________________________________
void AliEMCALv0::BuildGeometry()
{

  const Int_t kColorArm1   = kBlue ;
  const Int_t kColorArm2   = kBlue ;
  const Int_t kColorArm1Active   = kRed ;
  const Int_t kColorArm2Active   = kRed ;

  // make the container of  Arm1
  
  new TTUBS("Envelop1", "Tubs that contains arm 1", "void", 
	    fGeom->GetEnvelop(0),     // rmin 
	    fGeom->GetEnvelop(1),     // rmax
	    fGeom->GetEnvelop(2)/2.0, // half length in Z
	    fGeom->GetArm1PhiMin(),   // minimun phi angle
	    fGeom->GetArm1PhiMax()   // maximun phi angle
	    ) ; 
   // Active material of  Arm1
 
  new TTUBS("Arm1", "Active material of  arm 1", "void", 
	    fGeom->GetEnvelop(0),     // rmin 
	    fGeom->GetEnvelop(1),     // rmax
	    fGeom->GetLmat()/2.0,     // half length in Z
	    fGeom->GetArm1PhiMin(),   // minimun phi angle
	    fGeom->GetArm1PhiMax()   // maximun phi angle
	    ) ; 
  // make the container of  Arm2
 
  new TTUBS("Envelop2", "Tubs that contains arm 2", "void", 
	    fGeom->GetEnvelop(0),     // rmin 
	    fGeom->GetEnvelop(1),     // rmax
	    fGeom->GetEnvelop(2)/2.0, // half length in Z
	    fGeom->GetArm2PhiMin(),   // minimun phi angle
	    fGeom->GetArm2PhiMax()   // maximun phi angle
	    ) ;
	   
  // Active material of  Arm2
 
  new TTUBS("Arm2", "Active material of  arm 2", "void", 
	    fGeom->GetEnvelop(0),     // rmin 
	    fGeom->GetEnvelop(1),     // rmax
	    fGeom->GetLmat()/2.0,     // half length in Z	    ) ; 
	    fGeom->GetArm2PhiMin(),   // minimun phi angle
	    fGeom->GetArm2PhiMax()   // maximun phi angle
	    ) ;

  TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
  top->cd();
  
  // Arm 1 inside alice
  TNode * envelop1node = new TNode("Arm1 Envelop", "Arm1 Envelop", "Envelop1") ;
  envelop1node->SetLineColor(kColorArm1) ;
  fNodes->Add(envelop1node) ;

  // Arm 2 inside alice
  TNode * envelop2node = new TNode("Arm2 Envelop", "Arm2 Envelop", "Envelop2") ;
  envelop2node->SetLineColor(kColorArm2) ;
  fNodes->Add(envelop2node) ;

  // active material inside Arm 1
  envelop1node->cd() ; 
  TNode * arm1node = new TNode("Arm1 Mat", "Arm1 Mat", "Arm1Mat") ;
  arm1node->SetLineColor(kColorArm1Active) ;
  fNodes->Add(arm1node) ; 

  // active material inside Arm 2
  envelop2node->cd() ; 
  TNode * arm2node = new TNode("Arm2 Mat", "Arm2 Mat", "Arm2Mat") ;
  arm2node->SetLineColor(kColorArm2Active) ;
  fNodes->Add(arm2node) ; 
  
}


//____________________________________________________________________________
void AliEMCALv0::CreateGeometry()
{
  // Create the EMCAL geometry for Geant

  AliEMCALv0 *emcaltmp = (AliEMCALv0*)gAlice->GetModule("EMCAL") ;

  if ( emcaltmp == NULL ) {
    
    fprintf(stderr, "EMCAL detector not found!\n") ;
    return;
    
  }
  // Get pointer to the array containing media indices
  Int_t *idtmed = fIdtmed->GetArray() - 1599 ;


  // Create a tube sector that contains Arm 1 
 
  Float_t envelopA[5] ; 
  envelopA[0] = fGeom->GetEnvelop(0) ;         // rmin
  envelopA[1] = fGeom->GetEnvelop(1) ;         // rmax
  envelopA[2] = fGeom->GetEnvelop(2) / 2.0 ;   // dz
  envelopA[3] = fGeom->GetArm1PhiMin() ;       // minimun phi angle
  envelopA[4] = fGeom->GetArm1PhiMax() ;       // maximun phi angle

  gMC->Gsvolu("XEN1", "TUBS ", idtmed[1599], envelopA, 5) ; // filled with air

  // Create a tube sector that contains active material Arm 1 
 
  envelopA[2] = fGeom->GetLmat() / 2.0 ;       // dz

  gMC->Gsvolu("XAR1", "TUBS ", idtmed[1601], envelopA, 5) ; // filled with active material (average)


  // Create a tube sector that contains Arm 2 
 
  envelopA[3] = fGeom->GetArm2PhiMin() ;       // minimun phi angle
  envelopA[4] = fGeom->GetArm2PhiMax() ;       // maximun phi angle

  gMC->Gsvolu("XEN2", "TUBS ", idtmed[1599], envelopA, 5) ; // filled with air

  // Create a tube sector that contains active material Arm 1 
 
  envelopA[2] = fGeom->GetLmat() / 2.0 ;       // dz

  gMC->Gsvolu("XAR2", "TUBS ", idtmed[1601], envelopA, 5) ; // filled with active material (average)
  
  Int_t idrotm =1;
  AliMatrix(idrotm, 90.0, 0., 90.0, 90.0, 0.0, 0.0) ;

  // Position  ENV1 container in ALIC
  gMC->Gspos("XEN1", 1, "ALIC", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
  
  // Position  ARM1  into ENV1
  gMC->Gspos("XAR1", 1, "XEN1", 0.0, 0.0, 0.0, idrotm, "ONLY") ;

  // Position  ENV2 container in ALIC  
  gMC->Gspos("XEN2", 1, "ALIC", 0.0, 0.0, 0.0, idrotm, "ONLY") ;

  // Position  ARM2  into ENV2
  gMC->Gspos("XAR2", 1, "XEN2", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
  
}

//____________________________________________________________________________
void AliEMCALv0::Init(void)
{
  // Just prints an information message
  
  Int_t i;

  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" EMCAL_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");

  // Here the EMCAL initialisation code (if any!)

  if (fGeom!=0)  
    cout << "AliEMCAL" << Version() << " : EMCAL geometry intialized for " << fGeom->GetName() << endl ;
  else
    cout << "AliEMCAL" << Version() << " : EMCAL geometry initialization failed !" << endl ;   
  
  for(i=0;i<80;i++) printf("*");
  printf("\n");
  
}

