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
// This class places a Geometry of the EMCAL in the ALICE Detector as defined in AliEMCALGeometry.cxx                 
//*-- Author: Yves Schutz (SUBATECH)
//*-- and   : Sahal Yacoob (LBL / UCT)

// This Version of AliEMCALv0 reduces the number of volumes placed in XEN1 (the envelope) to less than five hundred
// The Envelope is Placed in Alice, And the Aluminium layer. Mini envelopes (XU) are then placed in XEN1.
// Each mini envelope contains 2 scintillator, and 2 lead layers, except the last one which contains just one scintillator layer.
// At the moment I cannot place the 36 and above layers in the mini envelopes so all layers are still placed in XEN1


// --- ROOT system ---
#include "TPGON.h"
#include "TTUBS.h"
#include "TNode.h"
#include "TRandom.h"
#include "TGeometry.h"
//#include "Tstring.h"

// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <strstream.h>
#include <iostream.h>

// --- AliRoot header files ---

#include "AliEMCALv0.h"
#include "AliEMCALGeometry.h"
#include "AliConst.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliEMCALv0)

//______________________________________________________________________
AliEMCALv0::AliEMCALv0(const char *name, const char *title):
    AliEMCAL(name,title){
    // Standard Constructor

    if (strcmp(GetTitle(),"") != 0 )
	fGeom =  AliEMCALGeometry::GetInstance(GetTitle(), "") ;

}
//______________________________________________________________________
void AliEMCALv0::BuildGeometry(){
    // Display Geometry for display.C

    const Int_t kColorArm1   = kBlue ;

    // Difine the shape of the Calorimeter 
  
    new TTUBS("Envelop1", "Tubs that contains arm 1", "void", 
	      fGeom->GetEnvelop(0),     // rmin 
	      fGeom->GetEnvelop(1) +30 ,     // rmax
	      fGeom->GetEnvelop(2)/2.0, // half length in Z
	      fGeom->GetArm1PhiMin(),   // minimun phi angle
	      fGeom->GetArm1PhiMax()    // maximun phi angle
	);
    // Place the Node
    TNode * envelop1node = new TNode("Envelop1", "Arm1 Envelop", "Envelop1") ;
    envelop1node->SetLineColor(kColorArm1) ;
    fNodes->Add(envelop1node) ;
}
//______________________________________________________________________
void AliEMCALv0::CreateGeometry(){
    // Create the EMCAL geometry for Geant

    AliEMCALv0 *emcaltmp = (AliEMCALv0*)gAlice->GetModule("EMCAL") ;

    if ( emcaltmp == NULL ) {
	Warning("CreateGeometry","detector not found!");
	return;
    } // end if
    // Get pointer to the array containing media indices
    Int_t *idtmed = fIdtmed->GetArray() - 1599 ;

    // Create an Envelope within which to place the Detector 
 
    Float_t envelopA[5]  ; 
      envelopA[0] = fGeom->GetEnvelop(0) ;         // rmin
      envelopA[1] = fGeom->GetEnvelop(1) + 30 ;    // rmax
      envelopA[2] = fGeom->GetEnvelop(2) / 2.0 ;   // dz
      envelopA[3] = fGeom->GetArm1PhiMin() ;       // minimun phi angle
      envelopA[4] = fGeom->GetArm1PhiMax() ;       // maximun phi angle
    
    
     // create XEN1
    gMC->Gsvolu("XEN1", "TUBS ", idtmed[1599], envelopA, 5) ; //filled with air
    
    Int_t idrotm = 1;
    AliMatrix(idrotm, 90.0, 0., 90.0, 90.0, 0.0, 0.0) ;

    // Position the Envelope in Alice  
    gMC->Gspos("XEN1", 1, "ALIC", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
     envelopA[1]  =envelopA[0] + fGeom->GetGap2Active(); 
    TString label = "XU0";

         envelopA[0] = envelopA[1]+ 3.18 ; //rmin Start mini envelopes after the aluminium layer
	envelopA[1] = envelopA[0] + 2.2    ;  //rmax larger for first two layers (preshower)

    gMC->Gsvolu(label.Data(), "TUBS ", idtmed[1599], envelopA, 5) ; //filled with air
    gMC->Gspos(label.Data(), 1, "XEN1", 0.0, 0.0, 0.0, idrotm, "ONLY") ; // Place XU0 in to XEN1

 for (int i = 1; i < ((fGeom->GetNLayers()-1)/2) + 1 ; i++ ){
       label = "XU" ;
       label += i ;	
       
       envelopA[0] = envelopA[1] ; //rmin
       envelopA[1] = envelopA[0] + 2.0  ;  //rmax 

    gMC->Gsvolu(label.Data(), "TUBS ", idtmed[1599], envelopA, 5) ; //filled with air
    gMC->Gspos(label.Data(), 1, "XEN1", 0.0, 0.0, 0.0, idrotm, "ONLY") ;

	} // end  i

    

    // Create the shapes of active material (LEAD/Aluminium/Scintillator) to be placed 

    Float_t envelopB[10]; // First Layer of Aluminium
    Float_t envelopC[10]; // Scintillator Layers
    Float_t envelopD[10]; //  Lead Layers
    
    envelopC[0] = envelopD[0] =  envelopB[0] = fGeom->GetArm1PhiMin(); //starting position in Phi
    envelopC[1] = envelopD[1] =  envelopB[1] = fGeom->GetArm1PhiMax() -  // Angular size of the Detector in Phi
	                                       fGeom->GetArm1PhiMin();
    envelopC[2] = envelopD[2] =  envelopB[2] = fGeom->GetNPhi() ; // Number of Section in Phi
    envelopD[3] = envelopC[3] = envelopB[3] = 2; // each section will be passed 2 z coordinates    

    envelopB[4] = (fGeom->GetEnvelop(0) + fGeom->GetGap2Active()) /
                  (tan(2*atan(exp(0.7)))) ;                          // z co-ordinate 1
    envelopB[5] = fGeom->GetEnvelop(0) + fGeom->GetGap2Active(); //rmin at z1
    envelopD[6] = envelopB[6] = envelopB[5] + 3.18;  //rmax at z1
    envelopB[7] = (fGeom->GetEnvelop(0) + fGeom->GetGap2Active()) /
	          (tan(2*atan(exp(-0.7)))) ;              // z co-ordinate 2
    envelopB[8] = envelopB[5] ; //
    envelopB[9] = envelopB[6] ; // radii are the same.

    // filled shapes wit hactive material 
    gMC->Gsvolu("XALU", "PGON", idtmed[1602], envelopB, 10); // Define Aluminium volume completely
    gMC->Gsvolu("XPST", "PGON", idtmed[1601], 0, 0) ;  // The polystyrene layers will be defined when placed 
    gMC->Gsvolu("XPBX", "PGON", idtmed[1600], 0, 0) ; //  as will the lead layers 
    gMC->Gsdvn("XPHI", "XPST", fGeom->GetNPhi(), 2) ; //  Dividind eta polystyrene divisions into phi segments.
    
    // Position Aluminium Layer in the Envelope 
    gMC->Gspos("XALU", 1, "XEN1", 0.0, 0.0, 0.0 , idrotm, "ONLY") ;

// The loop below places the scintillator in Lead Layers alternately.
 
    for (int i = 0; i < fGeom->GetNLayers() ; i++ ){

       label = "XU" ;
       
        label += (int) i/2  ; // we will place two layers (i = one layer) in each mini envelope)	
        envelopC[5] = envelopD[6] ; //rmin
	envelopC[6] = envelopD[6] + ((i > 1)  ? 0.5 : 0.6)  ;  //rmax larger for first two layers (preshower)
	envelopC[8] = envelopD[6] ; //rmin
	envelopC[9] = envelopD[6] + ((i > 1 ) ? 0.5 : 0.6)  ;  //rmax larger for first two layers (preshower)
	for (int j =0; j < (fGeom->GetNZ()) ; j++){
	    envelopC[4] = envelopD[6]/tan(2*atan(exp(0.7-(j*1.4/
                                                      (fGeom->GetNZ()))))); //z begin  
	    envelopC[7] = envelopD[6]/tan(2*atan(exp(0.7-((j+1)*1.4/
                                                      (fGeom->GetNZ()))))); // z end 
	    gMC->Gsposp("XPST",1+j+i*(fGeom->GetNZ()), label.Data(), // should be used but there's a weird crash above i = 18, 
			0.0, 0.0, 0.0 , idrotm, "ONLY", envelopC, 10); // Position and define layer
	} // end for j
	if (i < (fGeom->GetNLayers()-1)){
	   
             envelopD[5] = envelopC[6] ; //rmin
	    envelopD[6] = envelopC[6] + 0.5;  //rmax
	    envelopD[8] = envelopC[6] ; //rmin
	    envelopD[9] = envelopC[6] + 0.5;  //rmax
	    for (int j =0; j < (fGeom->GetNZ()) ; j++){
		envelopD[4] = envelopC[6]/tan(2*atan(exp(0.7-(j*1.4/
                                                       (fGeom->GetNZ()))))); // z begin 
		envelopD[7] = envelopC[6]/tan(2*atan(exp(0.7-((j+1)*1.4/
                                                        (fGeom->GetNZ()))))); // z end
		gMC->Gsposp("XPBX",1+ j+i*(fGeom->GetNZ()), label.Data(), 
			    0.0, 0.0, 0.0 , idrotm, "ONLY", envelopD, 10) ; // Position and Define Layer
	    } // end for j
	} // end if i
    }  // for i
}
//______________________________________________________________________
void AliEMCALv0::Init(void){
    // Just prints an information message
    Int_t i;

    cout << endl;
    for(i=0;i<35;i++) cout <<"*";
    cout << " EMCAL_INIT ";
    for(i=0;i<35;i++) cout << "*";
    cout << endl;

    // Here the EMCAL initialisation code (if any!)

    if (fGeom!=0)  
	cout << "AliEMCAL" << Version() << " : EMCAL geometry intialized for "
	     << fGeom->GetName() << endl ;
    else
	cout << "AliEMCAL" << Version() << 
	    " : EMCAL geometry initialization failed !" << endl ;
    for(i=0;i<80;i++) printf("*");
    cout << endl;
}
