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
  AliEMCAL(name,title)
{
  // ctor : title is used to identify the layout
  GetGeometry() ; 

}

//______________________________________________________________________
void AliEMCALv0::BuildGeometry()
{
    // Display Geometry for display.C

    const Int_t kColorArm1   = kBlue ;

    AliEMCALGeometry * geom = GetGeometry() ; 

    // Define the shape of the Calorimeter 
    TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
    new TTUBS("Envelop1", "Tubs that contains arm 1", "void", 
	      geom->GetEnvelop(0),     // rmin 
	      geom->GetEnvelop(1) +30 ,     // rmax
	      geom->GetEnvelop(2)/2.0, // half length in Z
	      geom->GetArm1PhiMin(),   // minimun phi angle
	      geom->GetArm1PhiMax()    // maximun phi angle
	);

    // Place the Node
    top->cd();
    TNode * envelop1node = new TNode("Envelop1", "Arm1 Envelop", "Envelop1"
				     ,0., 0., 0., "") ;
    envelop1node->SetLineColor(kColorArm1) ;
    fNodes->Add(envelop1node) ;
}

//______________________________________________________________________
void AliEMCALv0::CreateGeometry()
{
    // Create the EMCAL geometry for Geant

    Float_t etamin,etamax;
    Float_t *dum=0;

    AliEMCALGeometry * geom = GetGeometry() ; 

    if(!(geom->IsInitialized())){
	Error("CreateGeometry","EMCAL Geometry class has not been set up.");
    } // end if

    // Get pointer to the array containing media indices
    Int_t *idtmed = fIdtmed->GetArray() - 1599 ;

    // Create an Envelope within which to place the Detector 
    Float_t envelopA[5];
    envelopA[0] = geom->GetEnvelop(0);     // rmin
    envelopA[1] = geom->GetEnvelop(1);     // rmax
    envelopA[2] = geom->GetEnvelop(2)/2.0; // dz
    envelopA[3] = geom->GetArm1PhiMin();   // minimun phi angle
    envelopA[4] = geom->GetArm1PhiMax();   // maximun phi angle

    // create XEN1
    gMC->Gsvolu("XEN1", "TUBS ", idtmed[1599], envelopA, 5) ; //filled with air
    Int_t idrotm = 1;
    AliMatrix(idrotm, 90.0, 0., 90.0, 90.0, 0.0, 0.0) ;

    // Position the EMCAL Mother Volume in Alice  
    gMC->Gspos("XEN1", 1, "ALIC", 0.0, 0.0, 0.0, idrotm, "ONLY") ;

    //  
    TString label = "XU0";

    //rmin Start mini envelopes after the aluminium layer
    envelopA[0] = geom->GetEnvelop(0) + geom->GetGap2Active() + 
	          geom->GetAlFrontThickness();

    //rmax larger for first two layers (preshower);
    Float_t tseg = geom->GetPreSintThick()+geom->GetPbRadThick();
    envelopA[1] = envelopA[0] + 2.0*tseg;
    envelopA[2] = geom->GetEnvelop(2)/2.0; // dz
    envelopA[3] = geom->GetArm1PhiMin();   // minimun phi angle
    envelopA[4] = geom->GetArm1PhiMax();   // maximun phi angle

    //filled with air
    gMC->Gsvolu(label.Data(), "TUBS ", idtmed[1599], envelopA, 5);

    // Place XU0 in to XEN1
    gMC->Gspos(label.Data(), 1, "XEN1", 0.0, 0.0, 0.0, idrotm, "ONLY");

    tseg = geom->GetFullSintThick()+geom->GetPbRadThick();
    for (int i = 1; i < ((geom->GetNLayers()-1)/2) + 1 ; i++ ){
	label = "XU" ;
	label += i ;
	envelopA[0] = envelopA[1]; //rmin
	envelopA[1] = envelopA[0] + 2.0*tseg;  //rmax

	//filled with air
	gMC->Gsvolu(label.Data(), "TUBS ", idtmed[1599], envelopA, 5);
	gMC->Gspos(label.Data(), 1, "XEN1", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
    } // end  i

    // Create the shapes of active material (LEAD/Aluminium/Scintillator)
    // to be placed
    Float_t envelopB[10]; // First Layer of Aluminium
    Float_t envelopC[10]; // Scintillator Layers
    Float_t envelopD[10]; // Lead Layers

    //starting position in Phi
    envelopC[0] = envelopD[0] =  envelopB[0] = geom->GetArm1PhiMin();

    // Angular size of the Detector in Phi
    envelopB[1] = geom->GetArm1PhiMax() - geom->GetArm1PhiMin();
    envelopC[1] = envelopD[1] = envelopB[1];

    // Number of Section in Phi
    envelopC[2] = envelopD[2] = envelopB[2] = geom->GetNPhi();

    // each section will be passed 2 z coordinates    
    envelopD[3] = envelopC[3] = envelopB[3] = 2;
    envelopB[4] = geom->ZFromEtaR(geom->GetEnvelop(0)+geom->GetGap2Active(),
				   geom->GetArm1EtaMin());// z co-ordinate 1
    envelopB[5] = geom->GetEnvelop(0) + geom->GetGap2Active(); //rmin at z1
    envelopB[6] = envelopB[5] + geom->GetAlFrontThickness();//rmax at z1
    envelopD[6] = envelopB[6];
    envelopB[7] = geom->ZFromEtaR(geom->GetEnvelop(0)+geom->GetGap2Active(),
				  geom->GetArm1EtaMax()); // z co-ordinate 2
    envelopB[8] = envelopB[5] ; //
    envelopB[9] = envelopB[6] ; // radii are the same.

    // filled shapes wit hactive material 

    // Define Aluminium volume completely
    gMC->Gsvolu("XALU", "PGON", idtmed[1602], envelopB, 10);

    // The polystyrene layers will be defined when placed 
    gMC->Gsvolu("XPST", "PGON", idtmed[1601], dum, 0);
    gMC->Gsvolu("XPBX", "PGON", idtmed[1600], dum, 0);//  as will the lead layers

    //  Dividind eta polystyrene divisions into phi segments.
    gMC->Gsdvn("XPHI", "XPST", geom->GetNPhi(), 2);

    // Position Aluminium Layer in the Envelope 
    gMC->Gspos("XALU", 1, "XEN1", 0.0, 0.0, 0.0 , idrotm, "ONLY") ;

    // The loop below places the scintillator in Lead Layers alternately.
    for (int i = 0; i < geom->GetNLayers() ; i++ ){
	label = "XU" ;
        label += (int) i/2  ; // we will place two layers (i = one layer) in each mini envelope)	
        envelopC[5] = envelopD[6] ; //rmin
	envelopC[6] = envelopD[6] + ((i > 1)  ? geom->GetFullSintThick() : 
				     geom->GetPreSintThick());//rmax larger for first two layers (preshower)
	envelopC[8] = envelopD[6] ; //rmin
	envelopC[9] = envelopD[6] + ((i > 1 ) ? geom->GetFullSintThick() :
				     geom->GetPreSintThick());//rmax larger for first two layers (preshower)
	for (int j =0; j < (geom->GetNEta()) ; j++){
	    etamin = geom->GetArm1EtaMin()+
		(j*geom->GetDeltaEta());
	    etamax = geom->GetArm1EtaMin()+
		((j+1)*geom->GetDeltaEta());
	    envelopC[4] = geom->ZFromEtaR(envelopD[6],etamin); //z begin  
	    envelopC[7] = geom->ZFromEtaR(envelopD[6],etamax);// z end 
	    gMC->Gsposp("XPST",1+j+i*(geom->GetNEta()), label.Data(), // should be used but there's a weird crash above i = 18, 
			0.0, 0.0, 0.0 , idrotm, "ONLY", envelopC, 10); // Position and define layer
	} // end for j

	if (i < (geom->GetNLayers()-1)){
	    envelopD[5] = envelopC[6] ; //rmin
	    envelopD[6] = envelopC[6] + geom->GetPbRadThick();  //rmax
	    envelopD[8] = envelopC[6] ; //rmin
	    envelopD[9] = envelopC[6] + geom->GetPbRadThick();  //rmax
	    for (int j =0; j < (geom->GetNEta()) ; j++){
		etamin = geom->GetArm1EtaMin()+
		    (j*geom->GetDeltaEta());
		etamax = geom->GetArm1EtaMin()+
		    ((j+1)*geom->GetDeltaEta());
		envelopD[4] = geom->ZFromEtaR(envelopC[6],etamin);//z begin  
		envelopD[7] = geom->ZFromEtaR(envelopC[6],etamax);// z end

		// Position and Define Layer
		gMC->Gsposp("XPBX",1+ j+i*(geom->GetNEta()), label.Data(), 
			    0.0, 0.0, 0.0 , idrotm, "ONLY", envelopD, 10);
	    } // end for j
	} // end if i
    }  // for i
}

//______________________________________________________________________
void AliEMCALv0::Init(void)
{
    // Just prints an information message

    Int_t i;
    
    if(fDebug) { 
      cout << endl;
      for(i=0;i<35;i++) 
	cout <<"*";
    cout << "INFO: " << ClassName() << "::Init ";
    for(i=0;i<35;i++) 
      cout << "*";
    cout << endl;

    // Here the EMCAL initialisation code (if any!)

    AliEMCALGeometry * geom = GetGeometry() ; 
 
    if (geom!=0)  
	cout << "AliEMCAL" << Version() << " : EMCAL geometry intialized for "
	     << geom->GetName() << endl ;
    else
	cout << "AliEMCAL" << Version() << 
	    " : EMCAL geometry initialization failed !" << endl ;

    for(i=0;i<80;i++) 
      cout << "*" ;
    cout << endl;
    }
}
