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

/* $Id$*/

//_________________________________________________________________________
// Geometry class  for EMCAL : singleton  
// EMCAL consists of layers of scintillator and lead
// Places the the Barrel Geometry of The EMCAL at Midrapidity
// between 0 and 120 degrees of Phi and
// -0.7 to 0.7 in eta 
// Number of Modules and Layers may be controlled by 
// the name of the instance defined               
// EMCALArch2x has more modules along both phi and eta
// EMCALArchxa has less Layers in the Radial Direction
//*-- Author: Sahal Yacoob (LBL / UCT)
//     and  : Yves Schutz (SUBATECH)
//     and  : Jennifer Klay (LBL)

// --- ROOT system ---

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliEMCALGeometry.h"
#include "AliConst.h"

ClassImp(AliEMCALGeometry);

AliEMCALGeometry *AliEMCALGeometry::fgGeom = 0;
Bool_t            AliEMCALGeometry::fgInit = kFALSE;

//______________________________________________________________________
AliEMCALGeometry::~AliEMCALGeometry(void){
    // dtor
}
//______________________________________________________________________
void AliEMCALGeometry::Init(void){
    // Initializes the EMCAL parameters

    if(!(  (strcmp( fName, "EMCALArch1a" ) == 0) |
	   (strcmp( fName, "EMCALArch1b" ) == 0) | 
	   (strcmp( fName, "EMCALArch2a" ) == 0) | 
	   (strcmp( fName, "EMCALArch2b" ) == 0) )){
	fgInit = kFALSE;
	cout <<"Instance " << fName << " undefined" << endl;
    } // end if
    fgInit = kTRUE; 

    // geometry 
    fAirGap     = 5.0; 
    fArm1PhiMin = 0.0; 
    fArm1PhiMax = 120.0; 

    fIPDistance = 454.0; 
    fZLength = 817.0; 
    fEnvelop[0] = fIPDistance; 
    fEnvelop[2] = fZLength; 
    fGap2Active = 1.0; 
    fShellThickness = 3.18 + 1.2 + (double)((2*fNLayers -3)/2);   
    fEnvelop[1] = fIPDistance + fShellThickness;

    if (((strcmp( fName, "EMCALArch1a" ))    == 0) |
	((strcmp( fName, "EMCALArch1b" ))    == 0)){
	fNZ         = 96;
	fNPhi       = 144;
    } // end if
    if (((strcmp( fName, "EMCALArch2a" ))    == 0) |
	((strcmp( fName, "EMCALArch2b" ))    == 0)){
	fNZ         = 112;
	fNPhi       = 168;
    } // end if
    if (((strcmp( fName, "EMCALArch1a" ))    == 0) |
	((strcmp( fName, "EMCALArch2a" ))    == 0)){
	fNLayers    = 21;
    } // end if
    if (((strcmp( fName, "EMCALArch1b" ))    == 0) |
	((strcmp( fName, "EMCALArch2b" ))    == 0)){
	fNLayers    = 25;
    } // end if
}
//______________________________________________________________________
AliEMCALGeometry *  AliEMCALGeometry::GetInstance(){ 
    // Returns the pointer of the unique instance

    return (AliEMCALGeometry *) fgGeom; 
}
//______________________________________________________________________
AliEMCALGeometry* AliEMCALGeometry::GetInstance(const Text_t* name,
						const Text_t* title){
    // Returns the pointer of the unique instance

    AliEMCALGeometry * rv = 0; 
    if ( fgGeom == 0 ) {
	if ( strcmp(name,"") == 0 ) rv = 0;
	else {    
	    fgGeom = new AliEMCALGeometry(name, title);
	    if ( fgInit ) rv = (AliEMCALGeometry * ) fgGeom;
	    else {
		rv = 0; 
		delete fgGeom; 
		fgGeom = 0; 
	    } // end if fgInit
	} // end if strcmp(name,"")
    }else{
	if ( strcmp(fgGeom->GetName(), name) != 0 ) {
	    cout << "AliEMCALGeometry <E> : current geometry is " 
		 << fgGeom->GetName() << endl
		 << "                      you cannot call     " << name 
		 << endl; 
	}else{
	    rv = (AliEMCALGeometry *) fgGeom; 
	} // end if
    }  // end if fgGeom
    return rv; 
}
